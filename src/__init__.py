"""Paquete de simulación y análisis estadístico del Tablero de Galton (Las copas de champagne)."""

from src.simulator import (
    GaltonConfig,
    compute_exact_distribution,
    run_monte_carlo,
    simulate_single_run,
)
from src.statistics import (
    chi_square_discretized_normal_test,
    discretized_normal_expected,
    evaluate_row_count_normality,
    evaluate_single_sample_metrics,
)

__all__ = [
    "GaltonConfig",
    "compute_exact_distribution",
    "simulate_single_run",
    "run_monte_carlo",
    "discretized_normal_expected",
    "chi_square_discretized_normal_test",
    "evaluate_single_sample_metrics",
    "evaluate_row_count_normality",
]
