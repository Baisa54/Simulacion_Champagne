# TP: Las Copas de Champagne (Tablero de Galton) - MSTD

Trabajo práctico para la materia **Modelado, Simulación y Teoría de la Decisión (UNLu)** sobre la simulación del Tablero de Galton con $N = 100$ bolas por tirada y $1.000$ corridas de Monte Carlo por escenario.

## Qué se hizo

1. **Etapa 1 (Cantidad mínima de renglones):**  
   Se simularon tableros de 3 a 15 renglones tirando 100 bolas por corrida (1.000 corridas de Monte Carlo por cada tamaño) para ver a partir de cuántos renglones la distribución en las copas se puede considerar Normal. Usando el test Chi-cuadrado con agrupamiento de extremos (regla de Cochran, mínimo 5 bolas esperadas), la condición `X * p >= 5` y el test de D'Agostino-Pearson, llegamos a que el mínimo es **`X = 10` renglones (11 copas)**. También se dejó registrado por qué Shapiro-Wilk no sirve acá (al caer las 100 bolas en pocas copas enteras hay demasiados valores repetidos y rechaza siempre).

2. **Etapa 2 (Triángulo defectuoso en el 4.º renglón):**  
   Sobre el tablero de 10 renglones se puso un triángulo fallado en el 4.º renglón que desvía 70% a la izquierda y 30% a la derecha. Se comparó qué pasa en las 4 posiciones posibles y al azar: por los triángulos del centro pasa el triple de bolas (37,5%) que por los de las puntas (12,5%), corriendo el triple la media (`-0,075` vs `-0,025`), aunque en una sola tirada de 100 bolas el desvío es tan chico que casi no se detecta.

## Cómo probar la simulación (paso a paso)

Podés probar el simulador de dos formas. Si solo querés verlo funcionando en la pantalla sin instalar nada, usá la **Opción 1**.

### Opción 1: La forma más fácil (sin instalar nada, solo el simulador visual)

1. En esta misma página de GitHub, buscá el botón verde que dice **`<> Code`** (está arriba a la derecha de la lista de archivos), hacé clic ahí y elegí **`Download ZIP`**.
2. Andá a la carpeta **Descargas** de tu computadora. Vas a ver un archivo comprimido (con ícono de cierre). Hacé **clic derecho** sobre ese archivo y tocá **"Extraer todo..."** (o "Descomprimir aquí") para pasarlo a una carpeta normal.
3. Entrá a la carpeta que se acaba de crear, buscá el archivo **`demo_interactiva.html`** y hacele **doble clic**.
4. Se va a abrir solo en tu navegador (Chrome, Edge, Firefox o el que uses). Adentro ya podés tirar bolas, cambiar la cantidad de renglones, activar el triángulo defectuoso en el 4.º renglón y mirar los gráficos y tablas de las dos etapas, sin tocar ninguna consola ni instalar ningún programa.

### Opción 2: Si querés correr el código en Python (para volver a generar los gráficos y tablas)

1. **Instalar Python:** si no lo tenés en tu computadora, descargalo desde https://www.python.org/downloads/ y abrilo para instalarlo. Muy importante: en la primera pantalla del instalador, tildá la casilla que dice **"Add Python to PATH"** antes de tocar *Install Now*.
2. **Bajar el proyecto:** descargá y descomprimí la carpeta siguiendo los pasos 1 y 2 de la Opción 1 (o usá `git clone` si ya sabés usar Git).
3. **Abrir la ventana de comandos en la carpeta:** abrí la carpeta del proyecto en el Explorador de Windows, hacé clic arriba de todo en la barra donde figura la dirección de la carpeta, escribí **`cmd`** y apretá **Enter**. Se va a abrir una ventana negra ya parada dentro del proyecto.
4. **Instalar las librerías (se hace una sola vez):** en esa ventana negra escribí el siguiente comando y apretá **Enter** (esperá a que termine de descargar):
   ```powershell
   pip install -r requirements.txt
   ```
5. **Ejecutar lo que quieras probar:**
   - Para correr todas las simulaciones, ver las tablas en consola y volver a generar los gráficos y CSV dentro de `output/`:
     ```powershell
     python main.py
     ```
   - Para abrir el simulador visual en el navegador desde Python:
     ```powershell
     python demo.py
     ```
   - Para correr las pruebas automáticas y verificar que todo el código funciona bien:
     ```powershell
     python -m pytest -v
     ```

## Archivos del proyecto

- `demo_interactiva.html` y `demo.py`: simulador visual interactivo en el navegador.
- `main.py`: ejecuta las Etapas 1 y 2, muestra las tablas en consola y guarda gráficos y CSV en `output/`.
- `src/simulator.py`: simulación de las corridas Monte Carlo y cálculo exacto de probabilidades.
- `src/statistics.py`: test Chi-cuadrado con regla de Cochran, D'Agostino-Pearson y Shapiro-Wilk.
- `src/visualization.py`: generación de los gráficos de ambas etapas.
- `tests/test_galton.py`: pruebas unitarias con `pytest`.
- `INFORME_TECNICO.md` y `Informe_TP_Copas_de_Champagne.docx`: informe completo del TP con tablas, gráficos y conclusiones.
