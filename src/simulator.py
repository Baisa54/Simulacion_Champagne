"""Motor de simulación estocástica y cálculo exacto del Tablero de Galton ("Las copas de champagne")."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class GaltonConfig:
    """Configuración física y experimental de una corrida en el Tablero de Galton.

    Attributes:
        num_rows: Cantidad de renglones (filas de triángulos), X >= 1.
        num_balls: Cantidad de bolas lanzadas por corrida (por defecto N = 100).
        defect_row: Índice del renglón (base 1, 1..X) que contiene un triángulo defectuoso, o None.
        defect_col: Índice del triángulo dentro de `defect_row` (base 0, 0..defect_row-1),
            `"random"` para sortear uniformemente uno en cada corrida, o None si no hay defecto.
        defect_p_left: Probabilidad de desvío a la izquierda en el triángulo defectuoso (por defecto 0.70).
    """

    num_rows: int
    num_balls: int = 100
    defect_row: int | None = None
    defect_col: int | str | None = None
    defect_p_left: float = 0.70

    def __post_init__(self) -> None:
        if self.num_rows < 1:
            raise ValueError(f"num_rows debe ser >= 1, recibido: {self.num_rows}")
        if self.num_balls < 1:
            raise ValueError(f"num_balls debe ser >= 1, recibido: {self.num_balls}")
        if not (0.0 <= self.defect_p_left <= 1.0):
            raise ValueError(
                f"defect_p_left debe estar en [0, 1], recibido: {self.defect_p_left}"
            )

        has_row = self.defect_row is not None
        has_col = self.defect_col is not None
        if has_row != has_col:
            raise ValueError(
                "defect_row y defect_col deben especificarse juntos o ser ambos None."
            )

        if has_row:
            assert self.defect_row is not None
            if not (1 <= self.defect_row <= self.num_rows):
                raise ValueError(
                    f"defect_row debe estar entre 1 y {self.num_rows}, recibido: {self.defect_row}"
                )
            if self.defect_col != "random":
                if (
                    not isinstance(self.defect_col, int)
                    or isinstance(self.defect_col, bool)
                    or not (0 <= self.defect_col < self.defect_row)
                ):
                    raise ValueError(
                        f"defect_col debe ser 'random' o un entero en [0, {self.defect_row - 1}], "
                        f"recibido: {self.defect_col!r}"
                    )


def simulate_single_run(
    config: GaltonConfig, rng: np.random.Generator | None = None
) -> dict[str, Any]:
    """Simula una corrida completa de N bolas cayendo por X renglones del Tablero de Galton.

    Args:
        config: Parámetros del tablero y del experimento.
        rng: Generador de números pseudoaleatorios de NumPy.

    Returns:
        Diccionario con:
            - `positions`: vector `(N,)` con el receptáculo final `k in {0, ..., X}` de cada bola.
            - `counts`: vector `(X + 1,)` con la cantidad de bolas por receptáculo.
            - `row4_visits`: vector `(4,)` con el paso por cada triángulo del 4.º renglón (o None si X < 4).
            - `chosen_defect_col`: índice del triángulo defectuoso usado en la corrida (o None).
    """
    if rng is None:
        rng = np.random.default_rng()

    chosen_defect_col: int | None = None
    if config.defect_row is not None:
        if config.defect_col == "random":
            chosen_defect_col = int(rng.integers(0, config.defect_row))
        else:
            chosen_defect_col = int(config.defect_col)  # type: ignore[arg-type]

    col_idx = np.zeros(config.num_balls, dtype=np.int64)
    row4_visits: np.ndarray | None = None

    for row in range(1, config.num_rows + 1):
        if row == 4:
            row4_visits = np.bincount(col_idx, minlength=4)

        p_left = np.full(config.num_balls, 0.5, dtype=np.float64)
        if config.defect_row == row and chosen_defect_col is not None:
            p_left[col_idx == chosen_defect_col] = config.defect_p_left

        # Desvío a la izquierda mantiene col_idx; desvío a la derecha incrementa en 1.
        go_right = rng.random(config.num_balls) >= p_left
        col_idx += go_right.astype(np.int64)

    counts = np.bincount(col_idx, minlength=config.num_rows + 1)

    return {
        "positions": col_idx,
        "counts": counts,
        "row4_visits": row4_visits,
        "chosen_defect_col": chosen_defect_col,
    }


def compute_exact_distribution(
    num_rows: int,
    defect_row: int | None = None,
    defect_col: int | None = None,
    defect_p_left: float = 0.70,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Calcula por programación dinámica las probabilidades exactas en los X + 1 receptáculos.

    Args:
        num_rows: Cantidad de renglones X >= 1.
        defect_row: Renglón con defecto (base 1), o None.
        defect_col: Índice de triángulo defectuoso (base 0), o None.
        defect_p_left: Probabilidad de desvío a la izquierda del triángulo defectuoso.

    Returns:
        Tupla `(final_probs, row4_probs)` donde `final_probs` tiene longitud `X + 1` y
        `row4_probs` tiene longitud 4 con la probabilidad de paso por cada triángulo del 4.º renglón
        (o None si `num_rows < 4`).
    """
    if defect_row is not None and defect_col is None:
        raise ValueError("defect_col no puede ser None cuando defect_row está definido.")
    if defect_col is not None and defect_row is None:
        raise ValueError("defect_row no puede ser None cuando defect_col está definido.")
    GaltonConfig(
        num_rows=num_rows,
        num_balls=1,
        defect_row=defect_row,
        defect_col=defect_col,
        defect_p_left=defect_p_left,
    )

    probs = np.array([1.0], dtype=np.float64)
    row4_probs: np.ndarray | None = None

    for row in range(1, num_rows + 1):
        if row == 4:
            row4_probs = probs.copy()

        next_probs = np.zeros(row + 1, dtype=np.float64)
        p_left_vec = np.full(row, 0.5, dtype=np.float64)
        if defect_row == row and defect_col is not None:
            p_left_vec[defect_col] = defect_p_left

        next_probs[:-1] += probs * p_left_vec
        next_probs[1:] += probs * (1.0 - p_left_vec)
        probs = next_probs

    return probs, row4_probs


def compute_exact_random_defect_distribution(
    num_rows: int,
    defect_row: int = 4,
    defect_p_left: float = 0.70,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Calcula la distribución exacta marginal cuando el triángulo defectuoso se elige al azar."""
    if not (1 <= defect_row <= num_rows):
        raise ValueError(f"defect_row debe estar entre 1 y {num_rows}, recibido: {defect_row}")

    mixture = np.zeros(num_rows + 1, dtype=np.float64)
    row4_probs_ref: np.ndarray | None = None
    for col in range(defect_row):
        probs, row4_probs = compute_exact_distribution(
            num_rows=num_rows,
            defect_row=defect_row,
            defect_col=col,
            defect_p_left=defect_p_left,
        )
        mixture += probs / defect_row
        if row4_probs_ref is None:
            row4_probs_ref = row4_probs

    return mixture, row4_probs_ref


def run_monte_carlo(
    config: GaltonConfig,
    num_replications: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Ejecuta múltiples réplicas Monte Carlo independientes de la configuración dada.

    Args:
        config: Configuración del tablero de Galton.
        num_replications: Número de corridas independientes de `config.num_balls` bolas.
        seed: Semilla para reproducibilidad.

    Returns:
        Diccionario con matrices de conteos, posiciones, visitas al 4.º renglón y resúmenes.
    """
    if num_replications < 1:
        raise ValueError(f"num_replications debe ser >= 1, recibido: {num_replications}")

    rng = np.random.default_rng(seed)
    counts_matrix = np.zeros((num_replications, config.num_rows + 1), dtype=np.int64)
    positions_matrix = np.zeros((num_replications, config.num_balls), dtype=np.int64)
    row4_visits_matrix = (
        np.zeros((num_replications, 4), dtype=np.int64) if config.num_rows >= 4 else None
    )
    chosen_defect_cols: list[int | None] = []

    for rep in range(num_replications):
        res = simulate_single_run(config, rng)
        counts_matrix[rep] = res["counts"]
        positions_matrix[rep] = res["positions"]
        if row4_visits_matrix is not None and res["row4_visits"] is not None:
            row4_visits_matrix[rep] = res["row4_visits"]
        chosen_defect_cols.append(res["chosen_defect_col"])

    return {
        "config": config,
        "num_replications": num_replications,
        "seed": seed,
        "counts_matrix": counts_matrix,
        "positions_matrix": positions_matrix,
        "mean_counts": np.mean(counts_matrix, axis=0),
        "std_counts": np.std(counts_matrix, axis=0, ddof=1) if num_replications > 1 else np.zeros(config.num_rows + 1),
        "row4_visits_matrix": row4_visits_matrix,
        "chosen_defect_cols": chosen_defect_cols,
    }
