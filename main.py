"""Script principal para ejecutar las Etapas 1 y 2 del TP 'Las copas de champagne' (MSTD)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.simulator import (
    GaltonConfig,
    compute_exact_distribution,
    compute_exact_random_defect_distribution,
    run_monte_carlo,
)
from src.statistics import evaluate_row_count_normality, evaluate_single_sample_metrics
from src.visualization import (
    plot_stage1_column_comparison,
    plot_stage1_normality_evolution,
    plot_stage1_xmin_detail,
    plot_stage2_defect_comparison,
    plot_stage2_traffic_and_impact,
)

NUM_BALLS = 100
NUM_REPLICATIONS = 1000
ALPHA = 0.05
ROW_RANGE = range(3, 16)
BASE_SEED = 42
OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def _pick_representative_index(
    mc: dict[str, Any],
    mu_ref: float | None = None,
    sigma_ref: float | None = None,
) -> tuple[int, dict[str, Any]]:
    """Selecciona la corrida cuya estadística χ² está más cerca de la mediana Monte Carlo."""
    num_rows = mc["config"].num_rows
    num_reps = mc["num_replications"]
    metrics_list = [
        evaluate_single_sample_metrics(
            positions=mc["positions_matrix"][i],
            counts=mc["counts_matrix"][i],
            num_rows=num_rows,
            mu_ref=mu_ref,
            sigma_ref=sigma_ref,
            alpha=ALPHA,
        )
        for i in range(num_reps)
    ]
    p_vals = np.array([m["chi2_p"] for m in metrics_list], dtype=np.float64)
    med_p = float(np.nanmedian(p_vals))
    rep_idx = int(np.nanargmin(np.abs(p_vals - med_p)))
    return rep_idx, metrics_list[rep_idx]


def _exact_moments(probs: np.ndarray) -> tuple[float, float, float]:
    """Calcula media, desvío estándar y asimetría exactos de una distribución discreta."""
    k = np.arange(len(probs), dtype=np.float64)
    mu = float(np.sum(k * probs))
    var = float(np.sum(((k - mu) ** 2) * probs))
    sigma = float(np.sqrt(var))
    skew = float(np.sum(((k - mu) ** 3) * probs) / (sigma**3)) if sigma > 0 else 0.0
    return mu, sigma, skew


def run_stage1() -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]], int]:
    """Ejecuta Etapa 1: búsqueda del mínimo de renglones X para normalidad aceptable."""
    stage1_data: dict[int, dict[str, Any]] = {}
    summaries: list[dict[str, Any]] = []

    for x_val in ROW_RANGE:
        cfg = GaltonConfig(num_rows=x_val, num_balls=NUM_BALLS)
        mc = run_monte_carlo(cfg, num_replications=NUM_REPLICATIONS, seed=BASE_SEED + x_val)
        summary = evaluate_row_count_normality(mc, alpha=ALPHA)
        rep_idx, rep_metrics = _pick_representative_index(mc)

        stage1_data[x_val] = {
            "mc": mc,
            "summary": summary,
            "rep_idx": rep_idx,
            "rep_counts": mc["counts_matrix"][rep_idx],
            "rep_metrics": rep_metrics,
        }
        summaries.append(summary)

    acceptable = [s["num_rows"] for s in summaries if s["is_acceptable"]]
    if not acceptable:
        raise RuntimeError("Ningún valor de X cumplió el criterio de normalidad aceptable.")
    x_min = int(min(acceptable))

    # Imprimir tabla resumen en consola
    print("\n" + "=" * 112)
    print(f"ETAPA 1: DETERMINACIÓN DEL MÍNIMO DE RENGLONES (X) PARA NORMALIDAD ACEPTABLE (N={NUM_BALLS}, R={NUM_REPLICATIONS})")
    print("=" * 112)
    df1 = pd.DataFrame(
        [
            {
                "X": s["num_rows"],
                "Recept.": s["num_bins"],
                "Bins_Eff": s["effective_bins"],
                "gl": s["chi2_df"],
                "Chi2_med": round(s["chi2_stat_mean"], 2),
                "p_med_Chi2": round(s["chi2_p_median"], 3),
                "Acep_Chi2(%)": round(s["chi2_accept_rate"] * 100, 1),
                "Acep_K2(%)": round(s["dagostino_accept_rate"] * 100, 1),
                "Acep_SW_crudo(%)": round(s["shapiro_accept_rate"] * 100, 1),
                "Kurt_Teor": round(s["kurt_theo"], 3),
                "TVD": round(s["tvd_mean"], 4),
                "Veredicto": s["verdict"],
            }
            for s in summaries
        ]
    )
    print(df1.to_string(index=False))
    print("-" * 112)
    print(
        f"--> MÍNIMO DETERMINADO: X_min = {x_min} renglones ({x_min + 1} receptáculos). "
        f"Cumple De Moivre-Laplace (X*p={x_min * 0.5:.1f} >= 5), |curtosis| <= 0.20, "
        f"Bins efectivos >= 6, χ² >= 90% y D'Agostino K² >= 90%."
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df1.to_csv(OUTPUT_DIR / "resumen_etapa1.csv", index=False, encoding="utf-8-sig")

    return stage1_data, summaries, x_min


def run_stage2(x_min: int) -> tuple[list[dict[str, Any]], np.ndarray]:
    """Ejecuta Etapa 2: impacto de un triángulo defectuoso (70% izq / 30% der) en el 4.º renglón."""
    mu_base = x_min / 2.0
    sigma_base = float(np.sqrt(x_min) / 2.0)

    scenario_defs = [
        {
            "key": "baseline",
            "label": "Sin defecto (Tablero estándar 50/50)",
            "short_label": "Sin defecto",
            "defect_row": None,
            "defect_col": None,
            "traffic_prob": 1.0,
        },
        {
            "key": "pos_0",
            "label": "Defecto Renglón 4, Pos. 0 (Ext. Izq.)",
            "short_label": "Pos 0 (12.5%)",
            "defect_row": 4,
            "defect_col": 0,
            "traffic_prob": 0.125,
        },
        {
            "key": "pos_1",
            "label": "Defecto Renglón 4, Pos. 1 (Centro-Izq.)",
            "short_label": "Pos 1 (37.5%)",
            "defect_row": 4,
            "defect_col": 1,
            "traffic_prob": 0.375,
        },
        {
            "key": "pos_2",
            "label": "Defecto Renglón 4, Pos. 2 (Centro-Der.)",
            "short_label": "Pos 2 (37.5%)",
            "defect_row": 4,
            "defect_col": 2,
            "traffic_prob": 0.375,
        },
        {
            "key": "pos_3",
            "label": "Defecto Renglón 4, Pos. 3 (Ext. Der.)",
            "short_label": "Pos 3 (12.5%)",
            "defect_row": 4,
            "defect_col": 3,
            "traffic_prob": 0.125,
        },
        {
            "key": "random",
            "label": "Defecto Renglón 4, Posición Aleatoria",
            "short_label": "Aleatoria (Prom.)",
            "defect_row": 4,
            "defect_col": "random",
            "traffic_prob": 0.25,
        },
    ]

    scenarios: list[dict[str, Any]] = []
    row4_traffic_sim = np.zeros(4, dtype=np.float64)

    for idx, sdef in enumerate(scenario_defs):
        cfg = GaltonConfig(
            num_rows=x_min,
            num_balls=NUM_BALLS,
            defect_row=sdef["defect_row"],
            defect_col=sdef["defect_col"],
            defect_p_left=0.70,
        )
        mc = run_monte_carlo(cfg, num_replications=NUM_REPLICATIONS, seed=BASE_SEED + 100 + idx)
        if sdef["key"] == "baseline" and mc["row4_visits_matrix"] is not None:
            row4_traffic_sim = np.mean(mc["row4_visits_matrix"], axis=0) / float(NUM_BALLS)

        if sdef["defect_col"] is None:
            exact_probs, _ = compute_exact_distribution(x_min)
        elif sdef["defect_col"] == "random":
            exact_probs, _ = compute_exact_random_defect_distribution(
                x_min, defect_row=4, defect_p_left=0.70
            )
        else:
            exact_probs, _ = compute_exact_distribution(
                x_min, defect_row=4, defect_col=int(sdef["defect_col"]), defect_p_left=0.70
            )

        exact_mean, exact_std, exact_skew = _exact_moments(exact_probs)
        delta_mu = exact_mean - mu_base

        # Evaluar todas las réplicas contra la Normal estándar base N(X/2, X/4)
        rep_metrics_all = [
            evaluate_single_sample_metrics(
                positions=mc["positions_matrix"][r],
                counts=mc["counts_matrix"][r],
                num_rows=x_min,
                mu_ref=mu_base,
                sigma_ref=sigma_base,
                alpha=ALPHA,
            )
            for r in range(NUM_REPLICATIONS)
        ]

        chi2_p_vals = np.array([m["chi2_p"] for m in rep_metrics_all], dtype=np.float64)
        med_p = float(np.nanmedian(chi2_p_vals))
        rep_idx = int(np.nanargmin(np.abs(chi2_p_vals - med_p)))
        rep_m = rep_metrics_all[rep_idx]

        scenarios.append(
            {
                **sdef,
                "mc": mc,
                "exact_probs": exact_probs,
                "exact_mean": exact_mean,
                "exact_std": exact_std,
                "exact_skew": exact_skew,
                "delta_mu": delta_mu,
                "mc_mean": float(np.mean([m["mean"] for m in rep_metrics_all])),
                "mc_std": float(np.mean([m["std"] for m in rep_metrics_all])),
                "mc_skew": float(np.mean([m["skewness"] for m in rep_metrics_all])),
                "chi2_p_median": med_p,
                "chi2_reject_rate": float(1.0 - np.mean([m["chi2_accepted"] for m in rep_metrics_all])),
                "dagostino_reject_rate": float(1.0 - np.mean([m["dagostino_accepted"] for m in rep_metrics_all])),
                "rep_counts": mc["counts_matrix"][rep_idx],
                "rep_chi2_stat": rep_m["chi2_stat"],
                "rep_chi2_p": rep_m["chi2_p"],
            }
        )

    print("\n" + "=" * 116)
    print(f"ETAPA 2: IMPACTO DE TRIÁNGULO DEFECTUOSO (70% IZQ / 30% DER) EN EL 4.º RENGLÓN (X={x_min}, N={NUM_BALLS}, R={NUM_REPLICATIONS})")
    print("=" * 116)
    df2 = pd.DataFrame(
        [
            {
                "Escenario": sc["short_label"],
                "Trafico_Nodo(%)": round(sc["traffic_prob"] * 100, 1),
                "Bolas_Afect(/100)": round(sc["traffic_prob"] * 100, 1) if sc["key"] != "baseline" else 0.0,
                "Media_Exacta": round(sc["exact_mean"], 4),
                "Delta_Mu": round(sc["delta_mu"], 4),
                "Media_MC": round(sc["mc_mean"], 3),
                "Desvio_MC": round(sc["mc_std"], 3),
                "Asim_Exacta": round(sc["exact_skew"], 4),
                "Asim_MC": round(sc["mc_skew"], 3),
                "p_med_Chi2_vs_Base": round(sc["chi2_p_median"], 3),
                "Rechazo_Chi2(%)": round(sc["chi2_reject_rate"] * 100, 1),
            }
            for sc in scenarios
        ]
    )
    print(df2.to_string(index=False))
    print("-" * 116)
    print(
        "Tráfico empírico simulado por el 4.º renglón: "
        + ", ".join(f"Pos {i}={row4_traffic_sim[i]*100:.2f}%" for i in range(4))
        + " (Teórico: 12.50%, 37.50%, 37.50%, 12.50%)."
    )

    df2.to_csv(OUTPUT_DIR / "resumen_etapa2.csv", index=False, encoding="utf-8-sig")
    return scenarios, row4_traffic_sim


def main() -> None:
    """Punto de entrada principal del proyecto MSTD."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stage1_data, summaries, x_min = run_stage1()
    scenarios, row4_traffic_sim = run_stage2(x_min)

    x_panel = [3, 4, 6, 8, x_min, 14]
    if x_min in {3, 4, 6, 8, 14}:
        x_panel = [3, 4, 6, 8, 10, 14]

    saved_paths = [
        plot_stage1_column_comparison(
            stage1_data=stage1_data,
            x_values_to_plot=x_panel,
            x_min=x_min,
            output_path=OUTPUT_DIR / "etapa1_comparacion_columnas.png",
        ),
        plot_stage1_normality_evolution(
            summaries=summaries,
            x_min=x_min,
            output_path=OUTPUT_DIR / "etapa1_evolucion_normalidad.png",
        ),
        plot_stage1_xmin_detail(
            entry=stage1_data[x_min],
            output_path=OUTPUT_DIR / "etapa1_detalle_x_min.png",
        ),
        plot_stage2_defect_comparison(
            scenarios=scenarios,
            x_min=x_min,
            output_path=OUTPUT_DIR / "etapa2_comparacion_columnas_defecto.png",
        ),
        plot_stage2_traffic_and_impact(
            scenarios=scenarios,
            row4_traffic_sim=row4_traffic_sim,
            x_min=x_min,
            output_path=OUTPUT_DIR / "etapa2_analisis_trafico_e_impacto.png",
        ),
    ]

    print("\nArchivos generados en /output:")
    for p in saved_paths:
        print(f" - {p}")
    print(f" - {OUTPUT_DIR / 'resumen_etapa1.csv'}")
    print(f" - {OUTPUT_DIR / 'resumen_etapa2.csv'}")


if __name__ == "__main__":
    main()
