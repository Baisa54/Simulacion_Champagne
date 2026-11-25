"""Lanzador de la vista gráfica interactiva del TP 'Las copas de champagne' (MSTD)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import webbrowser

ROOT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT_DIR / "output"
HTML_PATH = ROOT_DIR / "demo_interactiva.html"

REQUIRED_ASSETS = [
    OUTPUT_DIR / "etapa1_comparacion_columnas.png",
    OUTPUT_DIR / "etapa1_evolucion_normalidad.png",
    OUTPUT_DIR / "etapa1_detalle_x_min.png",
    OUTPUT_DIR / "etapa2_comparacion_columnas_defecto.png",
    OUTPUT_DIR / "etapa2_analisis_trafico_e_impacto.png",
    OUTPUT_DIR / "resumen_etapa1.csv",
    OUTPUT_DIR / "resumen_etapa2.csv",
]


def ensure_outputs_ready() -> None:
    """Genera los gráficos y tablas en /output si alguno todavía no existe."""
    missing = [p for p in REQUIRED_ASSETS if not p.exists()]
    if missing:
        print("Faltan archivos en /output; ejecutando main.py para generarlos...")
        import main

        main.main()


def verify_demo_bundle() -> bool:
    """Verifica que la demo HTML y todos los gráficos existan."""
    required = [HTML_PATH, *REQUIRED_ASSETS]
    missing = [p for p in required if not p.exists()]
    if missing:
        print("Archivos faltantes:")
        for item in missing:
            print(f" - {item}")
        return False

    print("Verificación OK: todos los archivos de la demo interactiva están listos:")
    for item in required:
        print(f" - {item.name} ({item.stat().st_size:,} bytes)")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Abre la demo gráfica interactiva del Tablero de Galton en el navegador."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verifica que todos los archivos de la demo existan sin abrir el navegador.",
    )
    args = parser.parse_args()

    ensure_outputs_ready()
    if not verify_demo_bundle():
        return 1

    if args.check:
        return 0

    url = HTML_PATH.resolve().as_uri()
    print("\nAbriendo demo interactiva en tu navegador predeterminado:")
    print(f" -> {url}")
    print("\nVistas disponibles en la aplicación:")
    print(" 1) Pestaña 1 ('Simulador en Vivo'): probá X=4 (enunciado) vs X=10 (óptimo) y clic en los triángulos del 4.º renglón.")
    print(" 2) Pestaña 2 ('Resultados Etapa 1'): tabla Monte Carlo (R=1000) y gráficos comparativos.")
    print(" 3) Pestaña 3 ('Resultados Etapa 2'): árbol de tráfico nodal, tabla comparativa y gráficos en alta resolución.")
    webbrowser.open(url)
    return 0


if __name__ == "__main__":
    sys.exit(main())
