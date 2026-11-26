"""Suite de pruebas unitarias para el simulador y análisis estadístico del Tablero de Galton."""

from __future__ import annotations

import numpy as np
import pytest

from src.simulator import (
    GaltonConfig,
    compute_exact_distribution,
    compute_exact_random_defect_distribution,
    run_monte_carlo,
    simulate_single_run,
)
from src.statistics import (
    chi_square_discretized_normal_test,
    discretized_normal_expected,
    evaluate_row_count_normality,
    pool_tails_cochran,
)


def test_galton_config_validation() -> None:
    """Verifica validación de parámetros físicos del tablero y ubicación del defecto."""
    with pytest.raises(ValueError):
        GaltonConfig(num_rows=0)
    with pytest.raises(ValueError):
        GaltonConfig(num_rows=5, num_balls=0)
    with pytest.raises(ValueError):
        GaltonConfig(num_rows=5, defect_row=4, defect_col=None)
    with pytest.raises(ValueError):
        GaltonConfig(num_rows=5, defect_row=4, defect_col=4)
    with pytest.raises(ValueError):
        GaltonConfig(num_rows=3, defect_row=4, defect_col=0)
    with pytest.raises(ValueError):
        GaltonConfig(num_rows=10, defect_row=4, defect_col=1, defect_p_left=1.2)

    cfg = GaltonConfig(num_rows=10, num_balls=100, defect_row=4, defect_col="random")
    assert cfg.num_rows == 10


def test_simulate_single_run_conservation_and_bounds() -> None:
    """La suma de bolas por receptáculo debe ser exactamente N y caer en 0..X."""
    cfg = GaltonConfig(num_rows=10, num_balls=100)
    rng = np.random.default_rng(123)
    res = simulate_single_run(cfg, rng)

    assert len(res["counts"]) == 11
    assert int(np.sum(res["counts"])) == 100
    assert np.all(res["positions"] >= 0) and np.all(res["positions"] <= 10)
    assert res["row4_visits"] is not None
    assert int(np.sum(res["row4_visits"])) == 100
    assert res["chosen_defect_col"] is None


def test_exact_distribution_standard_and_row4_traffic() -> None:
    """En tablero estándar las probabilidades siguen Binomial(X, 0.5) y tráfico renglón 4 = [1/8, 3/8, 3/8, 1/8]."""
    probs, row4_probs = compute_exact_distribution(num_rows=10)
    assert np.isclose(np.sum(probs), 1.0)
    assert row4_probs is not None
    assert np.allclose(row4_probs, np.array([0.125, 0.375, 0.375, 0.125]))

    # Media exacta = X / 2 = 5.0, varianza exacta = X / 4 = 2.5
    k = np.arange(11, dtype=np.float64)
    mean_exact = float(np.sum(k * probs))
    var_exact = float(np.sum(((k - mean_exact) ** 2) * probs))
    assert np.isclose(mean_exact, 5.0)
    assert np.isclose(var_exact, 2.5)


def test_defect_shift_matches_analytical_expectation() -> None:
    """Un desvío 70/30 en el nodo j* desplaza la media exactamente en -0.20 * P(visita j*)."""
    probs_base, row4_probs = compute_exact_distribution(num_rows=10)
    assert row4_probs is not None
    k = np.arange(11, dtype=np.float64)
    mu_base = float(np.sum(k * probs_base))

    for j_star in range(4):
        probs_def, _ = compute_exact_distribution(
            num_rows=10, defect_row=4, defect_col=j_star, defect_p_left=0.70
        )
        assert np.isclose(np.sum(probs_def), 1.0)
        mu_def = float(np.sum(k * probs_def))
        expected_shift = -0.20 * float(row4_probs[j_star])
        assert np.isclose(mu_def - mu_base, expected_shift, atol=1e-12)


def test_random_defect_exact_mixture_and_mc_convergence() -> None:
    """El promedio Monte Carlo converge a la distribución exacta para defecto fijo o aleatorio."""
    mix_probs, _ = compute_exact_random_defect_distribution(num_rows=10, defect_row=4, defect_p_left=0.70)
    assert np.isclose(np.sum(mix_probs), 1.0)

    mc = run_monte_carlo(
        GaltonConfig(num_rows=10, num_balls=100, defect_row=4, defect_col="random", defect_p_left=0.70),
        num_replications=500,
        seed=99,
    )
    emp_probs = mc["mean_counts"] / 100.0
    assert np.allclose(emp_probs, mix_probs, atol=0.015)


def test_discretized_normal_and_cochran_pooling() -> None:
    """La normal discretizada suma N y el agrupamiento garantiza todas las clases E >= 5."""
    exp = discretized_normal_expected(num_rows=10, num_balls=100)
    assert len(exp) == 11
    assert np.isclose(np.sum(exp), 100.0)

    obs = np.round(exp).astype(np.int64)
    obs[5] += 100 - int(np.sum(obs))
    obs_p, exp_p, m = pool_tails_cochran(obs, exp, min_expected=5.0)
    assert np.all(exp_p >= 5.0)
    assert np.isclose(np.sum(obs_p), 100.0)
    assert np.isclose(np.sum(exp_p), 100.0)
    assert m >= 5

    chi_res = chi_square_discretized_normal_test(obs, exp)
    assert chi_res["valid"] is True
    assert chi_res["p_value"] > 0.90


def test_stage1_normality_monotone_convergence() -> None:
    """Con X=3 no cumple normalidad aceptable; con X=10 sí cumple y Shapiro-Wilk rechaza por empates."""
    mc_low = run_monte_carlo(GaltonConfig(num_rows=3, num_balls=100), num_replications=300, seed=10)
    res_low = evaluate_row_count_normality(mc_low)
    assert res_low["is_acceptable"] is False

    mc_ten = run_monte_carlo(GaltonConfig(num_rows=10, num_balls=100), num_replications=300, seed=10)
    res_ten = evaluate_row_count_normality(mc_ten)
    assert res_ten["chi2_accept_rate"] >= 0.90
    assert res_ten["dagostino_accept_rate"] >= 0.90
    assert res_ten["shapiro_accept_rate"] < 0.10
    assert res_ten["is_acceptable"] is True
