# TP: Las Copas de Champagne (Tablero de Galton) - MSTD

Trabajo práctico para la materia **Modelado, Simulación y Teoría de la Decisión (UNLu)** sobre la simulación del Tablero de Galton con $N = 100$ bolas por tirada y $1.000$ corridas de Monte Carlo por escenario.

## Qué se hizo

1. **Etapa 1 (Cantidad mínima de renglones):**  
   Se simularon tableros de 3 a 15 renglones tirando 100 bolas por corrida (1.000 corridas de Monte Carlo por cada tamaño) para ver a partir de cuántos renglones la distribución en las copas se puede considerar Normal. Usando el test Chi-cuadrado con agrupamiento de extremos (regla de Cochran, mínimo 5 bolas esperadas), la condición `X * p >= 5` y el test de D'Agostino-Pearson, llegamos a que el mínimo es **`X = 10` renglones (11 copas)**. También se dejó registrado por qué Shapiro-Wilk no sirve acá (al caer las 100 bolas en pocas copas enteras hay demasiados valores repetidos y rechaza siempre).

2. **Etapa 2 (Triángulo defectuoso en el 4.º renglón):**  
   Sobre el tablero de 10 renglones se puso un triángulo fallado en el 4.º renglón que desvía 70% a la izquierda y 30% a la derecha. Se comparó qué pasa en las 4 posiciones posibles y al azar: por los triángulos del centro pasa el triple de bolas (37,5%) que por los de las puntas (12,5%), corriendo el triple la media (`-0,075` vs `-0,025`), aunque en una sola tirada de 100 bolas el desvío es tan chico que casi no se detecta.

## Cómo probar la simulación

1. Instalar las librerías necesarias:
   ```powershell
   pip install -r requirements.txt
   ```

2. Correr la simulación por consola (imprime las tablas de las dos etapas y guarda los gráficos y CSV en `output/`):
   ```powershell
   python main.py
   ```

3. Abrir el simulador visual en el navegador:
   ```powershell
   python demo.py
   ```
   *(También se puede abrir directamente `demo_interactiva.html` haciendo doble clic en el archivo).*

4. Correr las pruebas unitarias:
   ```powershell
   python -m pytest -v
   ```

## Archivos del proyecto

- `main.py`: ejecuta las Etapas 1 y 2, muestra las tablas en consola y guarda gráficos y CSV en `output/`.
- `demo.py` y `demo_interactiva.html`: simulador visual interactivo en el navegador.
- `src/simulator.py`: simulación de las corridas Monte Carlo y cálculo exacto de probabilidades.
- `src/statistics.py`: test Chi-cuadrado con regla de Cochran, D'Agostino-Pearson y Shapiro-Wilk.
- `src/visualization.py`: generación de los gráficos de ambas etapas.
- `tests/test_galton.py`: pruebas unitarias con `pytest`.
- `INFORME_TECNICO.md` y `Informe_TP_Copas_de_Champagne.docx`: informe completo del TP con tablas, gráficos y conclusiones.
