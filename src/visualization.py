"""Generación de gráficos de columnas y análisis visual para ambas etapas del TP."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

from src.statistics import discretized_normal_expected


def _configure_style() -> None:
    """Configura un estilo visual limpio y académico para todos los gráficos."""
    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "#FAFAFA",
            "axes.edgecolor": "#CCCCCC",
            "axes.titleweight": "bold",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
        }
    )


def plot_stage1_column_comparison(
    stage1_data: dict[int, dict[str, Any]],
    x_values_to_plot: list[int],
    x_min: int,
    output_path: Path,
) -> Path:
    """Genera panel 2x3 de gráficos de columnas comparando X insuficiente vs. X mínimo y superior."""
    _configure_style()
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.8), sharey=False)
    axes_flat = axes.flatten()

    verdict_colors = {
        "NO ACEPTABLE": ("#D62728", "#FDEDEC"),
        "TRANSICIÓN": ("#D97706", "#FEF3C7"),
        "NORMAL ACEPTABLE": ("#15803D", "#DCFCE7"),
    }

    for ax, x_val in zip(axes_flat, x_values_to_plot):
        entry = stage1_data[x_val]
        mc = entry["mc"]
        summary = entry["summary"]
        rep_metrics = entry["rep_metrics"]
        rep_counts = entry["rep_counts"]

        bins = np.arange(x_val + 1)
        mean_counts = mc["mean_counts"]
        std_counts = mc["std_counts"]

        verdict = summary["verdict"]
        accent_color, bg_badge = verdict_colors[verdict]

        # Columnas de la corrida representativa (N = 100 bolas)
        ax.bar(
            bins,
            rep_counts,
            width=0.72,
            color="#4C78A8",
            alpha=0.78,
            edgecolor="#1F3A5F",
            linewidth=0.9,
            label="Corrida N=100",
            zorder=2,
        )

        # Promedio Monte Carlo +- 1 desvío estándar
        ax.errorbar(
            bins,
            mean_counts,
            yerr=std_counts,
            fmt="o",
            color="#111827",
            markersize=3.5,
            elinewidth=1.1,
            capsize=2.5,
            alpha=0.85,
            label="Media MC ±1σ",
            zorder=4,
        )

        # Curva normal continua escalada a N = 100 y puntos de normal discretizada
        mu_t = summary["mu_theo"]
        sigma_t = summary["sigma_theo"]
        x_dense = np.linspace(-0.5, x_val + 0.5, 300)
        y_norm_cont = stats.norm.pdf(x_dense, loc=mu_t, scale=sigma_t) * mc["config"].num_balls
        y_norm_disc = discretized_normal_expected(x_val, mc["config"].num_balls)

        ax.plot(
            x_dense,
            y_norm_cont,
            color="#E45756",
            linewidth=2.0,
            label=f"N({mu_t:.1f}, {sigma_t**2:.2f})",
            zorder=3,
        )
        ax.plot(
            bins,
            y_norm_disc,
            "s--",
            color="#B91C1C",
            markersize=3.8,
            linewidth=0.9,
            alpha=0.8,
            label="Normal discretizada",
            zorder=3,
        )

        tag_min = " [X mínimo]" if x_val == x_min else ""
        ax.set_title(
            f"X = {x_val} renglones ({x_val + 1} receptáculos){tag_min}",
            color=accent_color,
            pad=8,
        )
        ax.set_xlabel("Receptáculo (k)")
        ax.set_ylabel("Cantidad de bolas")
        ax.set_xticks(bins if x_val <= 10 else bins[::2])
        ax.set_xlim(-0.7, x_val + 0.7)
        y_max = max(float(np.max(rep_counts)), float(np.max(y_norm_cont)), float(np.max(mean_counts + std_counts)))
        ax.set_ylim(0, y_max * 1.38)

        info_text = (
            f"{verdict}\n"
            f"χ²={rep_metrics['chi2_stat']:.2f} (p={rep_metrics['chi2_p']:.3f})\n"
            f"D'Agostino K²={rep_metrics['dagostino_stat']:.2f} (p={rep_metrics['dagostino_p']:.3f})\n"
            f"Acep. MC: χ²={summary['chi2_accept_rate']*100:.0f}% | K²={summary['dagostino_accept_rate']*100:.0f}% | SW={summary['shapiro_accept_rate']*100:.0f}%\n"
            f"Curtosis teór.: {summary['kurt_theo']:+.3f}"
        )
        ax.text(
            0.03,
            0.96,
            info_text,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=7.7,
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": bg_badge,
                "edgecolor": accent_color,
                "linewidth": 1.1,
                "alpha": 0.94,
            },
            zorder=5,
        )
        ax.legend(loc="upper right", frameon=True, fontsize=7.5)

    fig.suptitle(
        "Etapa 1: Comparación de Distribuciones por Receptáculo (N = 100 bolas por corrida)\n"
        "Transición desde distribuciones discretas/platicúrticas hacia normalidad aceptable",
        fontsize=13.5,
        fontweight="bold",
        y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.935))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def plot_stage1_normality_evolution(
    summaries: list[dict[str, Any]],
    x_min: int,
    output_path: Path,
) -> Path:
    """Genera gráfico de evolución de tasas de aceptación y métricas de forma para X = 3..15."""
    _configure_style()
    x_vals = [s["num_rows"] for s in summaries]
    chi2_rates = [s["chi2_accept_rate"] * 100.0 for s in summaries]
    sw_rates = [s["shapiro_accept_rate"] * 100.0 for s in summaries]
    dag_rates = [s["dagostino_accept_rate"] * 100.0 for s in summaries]
    tvds = [s["tvd_mean"] for s in summaries]
    abs_kurts = [abs(s["kurt_theo"]) for s in summaries]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.5, 8.2), sharex=True)

    ax1.plot(x_vals, chi2_rates, "o-", linewidth=2.2, color="#1D4ED8", label="χ² Normal Discretizada (p > 0.05)")
    ax1.plot(x_vals, dag_rates, "^-", linewidth=2.0, color="#7C3AED", label="D'Agostino-Pearson K² (p > 0.05)")
    ax1.plot(x_vals, sw_rates, "s--", linewidth=1.9, color="#059669", label="Shapiro-Wilk crudo (penalizado por empates discretos)")

    ax1.axhline(90.0, color="#1D4ED8", linestyle=":", linewidth=1.3, alpha=0.8, label="Umbral aceptación χ² y K² (90%)")
    ax1.axvline(x_min, color="#DC2626", linestyle="--", linewidth=1.8, label=f"Mínimo óptimo X = {x_min}")
    ax1.axvspan(x_min - 0.3, max(x_vals) + 0.3, color="#DCFCE7", alpha=0.45, zorder=0)

    ax1.set_ylabel("Tasa de aceptación Monte Carlo (%)")
    ax1.set_ylim(-2, 106)
    ax1.set_title("Evolución de la Bondad de Ajuste Normal con R = 1000 réplicas (N = 100 bolas)")
    ax1.legend(loc="center right", frameon=True)

    ax2_twin = ax2.twinx()
    l1 = ax2.plot(x_vals, tvds, "D-", linewidth=2.1, color="#B45309", label="Distancia de Variación Total (TVD media)")
    l2 = ax2_twin.plot(x_vals, abs_kurts, "o--", linewidth=2.0, color="#BE123C", label="|Exceso de curtosis teórico| = 2/X")
    ax2.axvline(x_min, color="#DC2626", linestyle="--", linewidth=1.8)
    ax2.axvspan(x_min - 0.3, max(x_vals) + 0.3, color="#DCFCE7", alpha=0.45, zorder=0)

    ax2.set_xlabel("Cantidad de renglones (X)")
    ax2.set_ylabel("TVD media respecto a Normal Discretizada", color="#B45309")
    ax2_twin.set_ylabel("|Exceso de curtosis teórico| (2/X)", color="#BE123C")
    ax2.set_xticks(x_vals)

    lines = l1 + l2
    labels = [line.get_label() for line in lines]
    ax2.set_ylim(min(tvds) * 0.92, max(tvds) * 1.12)
    ax2.legend(lines, labels, loc="center right", frameon=True)
    ax2.set_title("Discrepancia muestral (TVD con N=100) y reducción de platicurtosis al aumentar X")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def plot_stage1_xmin_detail(
    entry: dict[str, Any],
    output_path: Path,
) -> Path:
    """Genera gráfico de columnas detallado para el valor óptimo X_min."""
    _configure_style()
    summary = entry["summary"]
    mc = entry["mc"]
    rep_counts = entry["rep_counts"]
    rep_metrics = entry["rep_metrics"]
    x_val = summary["num_rows"]
    num_balls = summary["num_balls"]

    bins = np.arange(x_val + 1)
    expected = discretized_normal_expected(x_val, num_balls)
    mu_t = summary["mu_theo"]
    sigma_t = summary["sigma_theo"]

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    width = 0.38

    bars_obs = ax.bar(
        bins - width / 2,
        rep_counts,
        width=width,
        color="#2563EB",
        alpha=0.88,
        edgecolor="#1E3A8A",
        label="Frecuencia observada O_k (corrida N=100)",
        zorder=3,
    )
    bars_exp = ax.bar(
        bins + width / 2,
        expected,
        width=width,
        color="#93C5FD",
        alpha=0.90,
        edgecolor="#1D4ED8",
        hatch="//",
        label="Esperada Normal Discretizada E_k",
        zorder=2,
    )

    x_dense = np.linspace(-0.5, x_val + 0.5, 350)
    y_norm = stats.norm.pdf(x_dense, loc=mu_t, scale=sigma_t) * num_balls
    ax.plot(
        x_dense,
        y_norm,
        color="#DC2626",
        linewidth=2.4,
        label=f"Densidad Normal continua N({mu_t:.1f}, {sigma_t**2:.2f}) × 100",
        zorder=4,
    )

    for bar in bars_obs:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.4, f"{int(h)}", ha="center", va="bottom", fontsize=8, color="#1E3A8A", fontweight="bold")
    for bar in bars_exp:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.4, f"{h:.1f}", ha="center", va="bottom", fontsize=7.5, color="#1D4ED8")

    ax.set_title(
        f"Detalle en X_mín = {x_val} renglones ({x_val + 1} receptáculos) — N = {num_balls} bolas\n"
        f"χ² = {rep_metrics['chi2_stat']:.2f} (gl={rep_metrics['chi2_df']}, p={rep_metrics['chi2_p']:.3f}) | "
        f"D'Agostino K² = {rep_metrics['dagostino_stat']:.2f} (p={rep_metrics['dagostino_p']:.3f}) | "
        f"Aceptación MC χ² = {summary['chi2_accept_rate']*100:.1f}%",
        fontsize=11.5,
    )
    ax.set_xlabel("Receptáculo (k)")
    ax.set_ylabel("Frecuencia (cantidad de bolas)")
    ax.set_xticks(bins)
    ax.set_ylim(0, max(float(np.max(rep_counts)), float(np.max(y_norm))) * 1.22)
    ax.legend(loc="upper right", frameon=True)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def plot_stage2_defect_comparison(
    scenarios: list[dict[str, Any]],
    x_min: int,
    output_path: Path,
) -> Path:
    """Genera panel 2x3 de gráficos de columnas comparando sin defecto vs. cada ubicación en renglón 4."""
    _configure_style()
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.8), sharey=True)
    axes_flat = axes.flatten()

    bins = np.arange(x_min + 1)
    mu_base = x_min / 2.0
    sigma_base = np.sqrt(x_min) / 2.0
    x_dense = np.linspace(-0.5, x_min + 0.5, 300)
    y_norm_base = stats.norm.pdf(x_dense, loc=mu_base, scale=sigma_base) * 100.0

    for ax, sc in zip(axes_flat, scenarios):
        rep_counts = sc["rep_counts"]
        mean_counts = sc["mc"]["mean_counts"]
        exact_counts = sc["exact_probs"] * 100.0
        is_baseline = sc["key"] == "baseline"

        bar_color = "#2563EB" if is_baseline else "#EA580C"
        edge_color = "#1E3A8A" if is_baseline else "#9A3412"

        ax.bar(
            bins,
            rep_counts,
            width=0.72,
            color=bar_color,
            alpha=0.74,
            edgecolor=edge_color,
            linewidth=0.9,
            label="Corrida N=100",
            zorder=2,
        )
        ax.plot(
            x_dense,
            y_norm_base,
            color="#374151",
            linestyle="--",
            linewidth=1.7,
            label=f"Normal base N({mu_base:.1f}, {sigma_base**2:.2f})",
            zorder=3,
        )
        ax.plot(
            bins,
            exact_counts,
            "o-",
            color="#B91C1C" if not is_baseline else "#1D4ED8",
            linewidth=1.8,
            markersize=4.2,
            label="Esperada exacta (con modelo)",
            zorder=4,
        )
        ax.axvline(sc["mc_mean"], color="#7C2D12", linestyle=":", linewidth=1.6, label=f"Media MC={sc['mc_mean']:.2f}")

        ax.set_title(sc["label"], fontsize=10.8, color="#111827" if is_baseline else "#9A3412")
        ax.set_xlabel("Receptáculo (k)")
        ax.set_ylabel("Cantidad de bolas")
        ax.set_xticks(bins)
        ax.set_xlim(-0.7, x_min + 0.7)

        badge_bg = "#EFF6FF" if is_baseline else "#FFF7ED"
        badge_edge = "#2563EB" if is_baseline else "#EA580C"
        info_txt = (
            f"Tráfico nodo: {sc['traffic_prob']*100:.1f}%\n"
            f"μ_exacta={sc['exact_mean']:.3f} (Δμ={sc['delta_mu']:+.3f})\n"
            f"σ_MC={sc['mc_std']:.3f} | Asim.={sc['mc_skew']:+.3f}\n"
            f"χ² vs N_base (corrida): {sc['rep_chi2_stat']:.2f} (p={sc['rep_chi2_p']:.3f})\n"
            f"Rechazo χ² MC: {sc['chi2_reject_rate']*100:.1f}%"
        )
        ax.text(
            0.03,
            0.96,
            info_txt,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=7.6,
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": badge_bg,
                "edgecolor": badge_edge,
                "linewidth": 1.0,
                "alpha": 0.94,
            },
            zorder=5,
        )
        ax.legend(loc="upper right", fontsize=7.2, frameon=True)

    y_max = max(float(np.max(sc["rep_counts"])) for sc in scenarios)
    for ax in axes_flat:
        ax.set_ylim(0, max(y_max, float(np.max(y_norm_base))) * 1.48)

    fig.suptitle(
        f"Etapa 2: Impacto de un Triángulo Defectuoso (70% Izq. / 30% Der.) en el 4.º Renglón (X = {x_min}, N = 100 bolas)",
        fontsize=13.5,
        fontweight="bold",
        y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def plot_stage2_traffic_and_impact(
    scenarios: list[dict[str, Any]],
    row4_traffic_sim: np.ndarray,
    x_min: int,
    output_path: Path,
) -> Path:
    """Genera panel tripartito explicando tráfico en renglón 4, Δμ/asimetría y curvas exactas."""
    _configure_style()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(17.0, 5.6))

    # Panel 1: Tráfico por cada uno de los 4 triángulos del 4.º renglón
    positions = np.arange(4)
    theo_traffic = np.array([12.5, 37.5, 37.5, 12.5])
    sim_traffic = row4_traffic_sim * 100.0
    w = 0.36

    b1 = ax1.bar(positions - w / 2, theo_traffic, width=w, color="#0EA5E9", edgecolor="#0369A1", label="Teórico Bin(3, 0.5)")
    b2 = ax1.bar(positions + w / 2, sim_traffic, width=w, color="#F97316", edgecolor="#C2410C", alpha=0.88, label="Simulado MC")
    for bar in list(b1) + list(b2):
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2, h + 0.8, f"{h:.1f}%", ha="center", va="bottom", fontsize=8)

    ax1.set_title("Tráfico de bolas por el 4.º renglón")
    ax1.set_xticks(positions)
    ax1.set_xticklabels(["Pos 0\n(Ext. Izq)", "Pos 1\n(Centro-Izq)", "Pos 2\n(Centro-Der)", "Pos 3\n(Ext. Der)"])
    ax1.set_ylabel("Porcentaje de bolas que impactan el nodo (%)")
    ax1.set_ylim(0, 48)
    ax1.legend(loc="upper right", frameon=True)

    # Panel 2: Corrimiento de media (Δμ) y asimetría por posición fija
    pos_scenarios = [sc for sc in scenarios if sc["key"].startswith("pos_")]
    pos_labels = ["Pos 0\n(12.5%)", "Pos 1\n(37.5%)", "Pos 2\n(37.5%)", "Pos 3\n(12.5%)"]
    delta_mus = [sc["delta_mu"] for sc in pos_scenarios]
    skews = [sc["exact_skew"] for sc in pos_scenarios]

    x_idx = np.arange(len(pos_scenarios))
    ax2.bar(x_idx - 0.18, delta_mus, width=0.35, color="#DC2626", alpha=0.85, label="Corrimiento media Δμ (exacto)")
    ax2.bar(x_idx + 0.18, skews, width=0.35, color="#7C3AED", alpha=0.85, label="Asimetría exacta (Skewness)")
    ax2.axhline(0.0, color="#374151", linewidth=1.0)
    for idx_v, (dm, sk) in enumerate(zip(delta_mus, skews)):
        ax2.text(idx_v - 0.18, dm - 0.004, f"{dm:+.3f}", ha="center", va="top", fontsize=7.8, color="#991B1B", fontweight="bold")
        va_sk = "bottom" if sk >= 0 else "top"
        off_sk = 0.003 if sk >= 0 else -0.004
        ax2.text(idx_v + 0.18, sk + off_sk, f"{sk:+.3f}", ha="center", va=va_sk, fontsize=7.8, color="#5B21B6", fontweight="bold")

    ax2.set_title("Impacto en Δμ y Asimetría según ubicación")
    ax2.set_xticks(x_idx)
    ax2.set_xticklabels(pos_labels)
    ax2.set_ylabel("Magnitud del efecto")
    y_low = min(min(delta_mus), min(skews)) - 0.035
    y_high = max(0.03, max(skews) + 0.03)
    ax2.set_ylim(y_low, y_high)
    ax2.legend(loc="lower right", frameon=True, fontsize=8)

    # Panel 3: Curvas de probabilidad exacta P(K = k) comparadas
    bins = np.arange(x_min + 1)
    baseline_sc = next(sc for sc in scenarios if sc["key"] == "baseline")
    ax3.plot(bins, baseline_sc["exact_probs"] * 100.0, "k--", linewidth=2.2, label="Sin defecto (Base)")
    palette = ["#0284C7", "#DC2626", "#EA580C", "#16A34A"]
    for sc, col in zip(pos_scenarios, palette):
        ax3.plot(bins, sc["exact_probs"] * 100.0, "o-", linewidth=1.7, markersize=4, color=col, label=sc["short_label"])

    ax3.set_title(f"Distribuciones exactas esperadas (X = {x_min})")
    ax3.set_xlabel("Receptáculo (k)")
    ax3.set_ylabel("Frecuencia esperada (%)")
    ax3.set_xticks(bins)
    ax3.legend(loc="upper right", frameon=True, fontsize=8)

    fig.suptitle(
        "Análisis Estructural del Defecto en el 4.º Renglón: Peso del Tráfico Nodal sobre la Distribución Final",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path
