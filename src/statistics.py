"""Métricas estadísticas y pruebas de bondad de ajuste para el Tablero de Galton."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import stats


def discretized_normal_expected(
    num_rows: int,
    num_balls: int = 100,
    mu: float | None = None,
    sigma: float | None = None,
) -> np.ndarray:
    """Calcula las frecuencias esperadas por receptáculo bajo una Normal discretizada.

    Utiliza corrección por continuidad asignando a cada receptáculo `k` la probabilidad
    normal acumulada en `[k - 0.5, k + 0.5]`, extendiendo las colas extremas hacia
    `(-inf, 0.5]` y `[X - 0.5, +inf)` para conservar suma exacta igual a `num_balls`.

    Args:
        num_rows: Cantidad de renglones `X`.
        num_balls: Cantidad total de bolas `N`.
        mu: Media teórica (por defecto `X / 2`).
        sigma: Desvío estándar teórico (por defecto `sqrt(X) / 2`).

    Returns:
        Vector de longitud `X + 1` con las frecuencias esperadas `E_k`.
    """
    if mu is None:
        mu = num_rows / 2.0
    if sigma is None:
        sigma = float(np.sqrt(num_rows) / 2.0)
    if sigma <= 0.0:
        raise ValueError(f"sigma debe ser > 0, recibido: {sigma}")

    if num_rows == 0:
        return np.array([float(num_balls)], dtype=np.float64)

    edges = np.arange(num_rows, dtype=np.float64) + 0.5
    cdfs = stats.norm.cdf(edges, loc=mu, scale=sigma)

    probs = np.empty(num_rows + 1, dtype=np.float64)
    probs[0] = cdfs[0]
    if num_rows > 1:
        probs[1:-1] = np.diff(cdfs)
    probs[-1] = 1.0 - cdfs[-1]

    return probs * float(num_balls)


def pool_tails_cochran(
    observed: np.ndarray,
    expected: np.ndarray,
    min_expected: float = 5.0,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Agrupa colas izquierda y derecha para satisfacer la condición de Cochran (E_i >= 5).

    Args:
        observed: Frecuencias observadas por receptáculo.
        expected: Frecuencias esperadas por receptáculo.
        min_expected: Frecuencia esperada mínima requerida por clase agrupada.

    Returns:
        Tupla `(observed_pooled, expected_pooled, effective_bins)`.
    """
    obs = observed.astype(np.float64).tolist()
    exp = expected.astype(np.float64).tolist()

    # Agrupar colas extremas hacia el centro mientras alguna tenga E_i < min_expected.
    while len(exp) > 1 and (exp[0] < min_expected or exp[-1] < min_expected):
        if exp[0] < min_expected and (exp[0] <= exp[-1] or exp[-1] >= min_expected):
            obs[1] += obs[0]
            exp[1] += exp[0]
            obs.pop(0)
            exp.pop(0)
        else:
            obs[-2] += obs[-1]
            exp[-2] += exp[-1]
            obs.pop()
            exp.pop()

    # Resguardo adicional si quedara alguna clase interior por debajo del umbral.
    idx = 0
    while len(exp) > 1 and idx < len(exp):
        if exp[idx] < min_expected:
            if idx == 0:
                target = 1
            elif idx == len(exp) - 1:
                target = len(exp) - 2
            else:
                target = idx - 1 if exp[idx - 1] <= exp[idx + 1] else idx + 1
            obs[target] += obs[idx]
            exp[target] += exp[idx]
            obs.pop(idx)
            exp.pop(idx)
            idx = max(0, idx - 1)
        else:
            idx += 1

    return np.asarray(obs, dtype=np.float64), np.asarray(exp, dtype=np.float64), len(exp)


def chi_square_discretized_normal_test(
    counts: np.ndarray,
    expected: np.ndarray,
    ddof: int = 0,
    min_expected: float = 5.0,
) -> dict[str, float | int | bool]:
    """Ejecuta el test Chi-cuadrado de bondad de ajuste con agrupamiento de colas de Cochran.

    Args:
        counts: Frecuencias observadas en los receptáculos.
        expected: Frecuencias esperadas bajo el modelo normal discretizado.
        ddof: Grados de libertad perdidos por parámetros estimados desde la muestra
            (0 para H0 con parámetros teóricos conocidos `mu=X/2, sigma=sqrt(X)/2`).
        min_expected: Mínimo esperado por contenedor agrupado (5.0 según regla de Cochran).

    Returns:
        Diccionario con `statistic`, `p_value`, `df`, `effective_bins` y `valid`.
    """
    obs_p, exp_p, effective_bins = pool_tails_cochran(counts, expected, min_expected=min_expected)
    df = effective_bins - 1 - ddof
    if df <= 0:
        return {
            "statistic": float("nan"),
            "p_value": float("nan"),
            "df": int(df),
            "effective_bins": int(effective_bins),
            "valid": False,
        }

    # Normalizar suma exacta para evitar errores de redondeo de punto flotante en scipy.chisquare
    exp_p = exp_p * (np.sum(obs_p) / np.sum(exp_p))
    stat, p_val = stats.chisquare(f_obs=obs_p, f_exp=exp_p, ddof=ddof)
    return {
        "statistic": float(stat),
        "p_value": float(p_val),
        "df": int(df),
        "effective_bins": int(effective_bins),
        "valid": True,
    }


def evaluate_single_sample_metrics(
    positions: np.ndarray,
    counts: np.ndarray,
    num_rows: int,
    mu_ref: float | None = None,
    sigma_ref: float | None = None,
    ddof: int = 0,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Calcula estadísticos descriptivos y pruebas de normalidad para una corrida individual."""
    n = len(positions)
    expected = discretized_normal_expected(
        num_rows=num_rows, num_balls=n, mu=mu_ref, sigma=sigma_ref
    )
    chi2_res = chi_square_discretized_normal_test(counts, expected, ddof=ddof)

    mean_val = float(np.mean(positions))
    std_val = float(np.std(positions, ddof=1)) if n > 1 else 0.0

    if std_val > 0.0:
        skew_val = float(stats.skew(positions, bias=False))
        kurt_val = float(stats.kurtosis(positions, fisher=True, bias=False))
        shapiro_stat, shapiro_p = stats.shapiro(positions)
        if n >= 8:
            dag_stat, dag_p = stats.normaltest(positions)
        else:
            dag_stat, dag_p = float("nan"), float("nan")
    else:
        skew_val = 0.0
        kurt_val = 0.0
        shapiro_stat, shapiro_p = float("nan"), float("nan")
        dag_stat, dag_p = float("nan"), float("nan")

    obs_probs = counts / float(n)
    exp_probs = expected / float(np.sum(expected))
    tvd = 0.5 * float(np.sum(np.abs(obs_probs - exp_probs)))

    return {
        "mean": mean_val,
        "std": std_val,
        "skewness": skew_val,
        "excess_kurtosis": kurt_val,
        "chi2_stat": chi2_res["statistic"],
        "chi2_p": chi2_res["p_value"],
        "chi2_df": chi2_res["df"],
        "effective_bins": chi2_res["effective_bins"],
        "chi2_accepted": bool(chi2_res["valid"] and chi2_res["p_value"] > alpha),
        "shapiro_stat": float(shapiro_stat),
        "shapiro_p": float(shapiro_p),
        "shapiro_accepted": bool(not np.isnan(shapiro_p) and shapiro_p > alpha),
        "dagostino_stat": float(dag_stat),
        "dagostino_p": float(dag_p),
        "dagostino_accepted": bool(not np.isnan(dag_p) and dag_p > alpha),
        "tvd": tvd,
    }


def evaluate_row_count_normality(
    mc_result: dict[str, Any],
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Evalúa la normalidad Monte Carlo para un tablero estándar de X renglones (Etapa 1).

    Criterio compuesto de normalidad aceptable con N = 100 bolas:
      1. Criterio estructural de De Moivre-Laplace / Cochran: `X * p >= 5` y `X * (1 - p) >= 5`
         (equivale a `X >= 10`, `sigma >= 1.581` y exceso de curtosis teórico `|gamma_2| <= 0.20`),
         junto con al menos 6 contenedores efectivos (`E_i >= 5`).
      2. Bondad de ajuste Chi-cuadrado contra Normal Discretizada: tasa de aceptación >= 90%.
      3. Test omnibus de forma (D'Agostino-Pearson K²): tasa de no-rechazo >= 90%.

    Nota: Shapiro-Wilk crudo sobre enteros discretos se reporta como referencia comparativa para
    evidenciar la penalización por empates en datos agrupados en contenedores.
    """
    config = mc_result["config"]
    num_rows = config.num_rows
    num_balls = config.num_balls
    counts_matrix = mc_result["counts_matrix"]
    positions_matrix = mc_result["positions_matrix"]
    num_reps = mc_result["num_replications"]

    per_rep = [
        evaluate_single_sample_metrics(
            positions=positions_matrix[i],
            counts=counts_matrix[i],
            num_rows=num_rows,
            alpha=alpha,
        )
        for i in range(num_reps)
    ]

    chi2_p_arr = np.array([m["chi2_p"] for m in per_rep], dtype=np.float64)
    chi2_stat_arr = np.array([m["chi2_stat"] for m in per_rep], dtype=np.float64)
    shapiro_p_arr = np.array([m["shapiro_p"] for m in per_rep], dtype=np.float64)
    shapiro_stat_arr = np.array([m["shapiro_stat"] for m in per_rep], dtype=np.float64)
    dag_p_arr = np.array([m["dagostino_p"] for m in per_rep], dtype=np.float64)
    tvd_arr = np.array([m["tvd"] for m in per_rep], dtype=np.float64)
    means_arr = np.array([m["mean"] for m in per_rep], dtype=np.float64)
    stds_arr = np.array([m["std"] for m in per_rep], dtype=np.float64)
    skews_arr = np.array([m["skewness"] for m in per_rep], dtype=np.float64)
    kurts_arr = np.array([m["excess_kurtosis"] for m in per_rep], dtype=np.float64)

    chi2_accept_rate = float(np.mean([m["chi2_accepted"] for m in per_rep]))
    shapiro_accept_rate = float(np.mean([m["shapiro_accepted"] for m in per_rep]))
    dagostino_accept_rate = float(np.mean([m["dagostino_accepted"] for m in per_rep]))

    mu_theo = num_rows / 2.0
    sigma_theo = float(np.sqrt(num_rows) / 2.0)
    kurt_theo = -2.0 / float(num_rows)
    np_criterion = (num_rows * 0.5) >= 5.0
    effective_bins = int(per_rep[0]["effective_bins"])

    is_acceptable = bool(
        np_criterion
        and effective_bins >= 6
        and chi2_accept_rate >= 0.90
        and dagostino_accept_rate >= 0.90
    )
    if is_acceptable:
        verdict = "NORMAL ACEPTABLE"
    elif num_rows >= 7 and chi2_accept_rate >= 0.90 and dagostino_accept_rate >= 0.90:
        verdict = "TRANSICIÓN"
    else:
        verdict = "NO ACEPTABLE"

    return {
        "num_rows": num_rows,
        "num_bins": num_rows + 1,
        "num_balls": num_balls,
        "mu_theo": mu_theo,
        "sigma_theo": sigma_theo,
        "kurt_theo": kurt_theo,
        "effective_bins": per_rep[0]["effective_bins"],
        "chi2_df": per_rep[0]["chi2_df"],
        "mean_obs": float(np.mean(means_arr)),
        "std_obs": float(np.mean(stds_arr)),
        "skew_obs": float(np.mean(skews_arr)),
        "kurt_obs": float(np.mean(kurts_arr)),
        "chi2_stat_mean": float(np.nanmean(chi2_stat_arr)),
        "chi2_p_median": float(np.nanmedian(chi2_p_arr)),
        "chi2_accept_rate": chi2_accept_rate,
        "shapiro_w_mean": float(np.nanmean(shapiro_stat_arr)),
        "shapiro_p_median": float(np.nanmedian(shapiro_p_arr)),
        "shapiro_accept_rate": shapiro_accept_rate,
        "dagostino_p_median": float(np.nanmedian(dag_p_arr)),
        "dagostino_accept_rate": dagostino_accept_rate,
        "tvd_mean": float(np.mean(tvd_arr)),
        "np_criterion": np_criterion,
        "is_acceptable": is_acceptable,
        "verdict": verdict,
    }
