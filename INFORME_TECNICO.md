# Informe Técnico: Simulación y Análisis Estadístico del Tablero de Galton ("Las Copas de Champagne")

**Materia:** Modelado, Simulación y Teoría de la Decisión (MSTD)  
**Tamaño de muestra por experimento:** $N = 100$ bolas  
**Réplicas Monte Carlo por escenario:** $R = 1\,000$ corridas independientes ($\alpha = 0.05$)

---

## 1. Introducción y Modelado del Sistema

El experimento "Las copas de champagne" corresponde a un **Tablero de Galton** compuesto por $X$ renglones (filas horizontales de triángulos clavados, numeradas $r = 1, 2, \dots, X$) donde el renglón $r$ posee exactamente $r$ triángulos indexados $j \in \{0, 1, \dots, r-1\}$ de izquierda a derecha. En la base se encuentran $X + 1$ receptáculos numerados $k \in \{0, 1, \dots, X\}$.

Cada bola ingresa centrada sobre el triángulo superior $(r=1, j=0)$. Al impactar un triángulo estándar, la bola cae hacia la izquierda ($j \to j$) con probabilidad $p_L = 0.50$ o hacia la derecha ($j \to j + 1$) con probabilidad $p_R = 0.50$. En ausencia de defectos, el receptáculo final $K$ equivale al número total de desvíos a la derecha en $X$ etapas independientes:

$$K \sim \text{Binomial}(X, p = 0.5)$$

Con parámetros teóricos exactos:
$$\mu = \frac{X}{2}, \qquad \sigma = \frac{\sqrt{X}}{2}, \qquad \gamma_1 = 0 \text{ (asimetría)}, \qquad \gamma_2 = -\frac{2}{X} \text{ (exceso de curtosis)}$$

---

## 2. Etapa 1: Determinación del Mínimo de Renglones ($X_{\min}$) para Normalidad Aceptable

### 2.1 Justificación de la Herramienta de Bondad de Ajuste

Para evaluar si la distribución observada de $N = 100$ bolas en $X + 1$ receptáculos constituye una **distribución normal aceptable**, es crucial distinguir entre pruebas diseñadas para variables continuas y pruebas adecuadas para conteos discretos agrupados en contenedores:

1. **Por qué NO debe usarse Shapiro-Wilk crudo como único criterio decisor:**  
   El test de Shapiro-Wilk presupone una variable continua sin valores repetidos (*ties*). En el Tablero de Galton, las $N = 100$ bolas caen únicamente en $X + 1$ valores enteros ($4$ a $16$ receptáculos), generando entre $84\%$ y $96\%$ de empates exactos. Como verificaron las simulaciones Monte Carlo, Shapiro-Wilk aplicado directamente sobre los enteros discreto-agrupados rechaza casi el $100\%$ de las muestras ($0.0\%$ de aceptación para $X \le 11$ y apenas $4.9\%$ en $X = 15$), lo cual no refleja falta de convergencia binomial sino un desajuste metodológico del test frente a datos discretizados.

2. **Herramienta elegida — Test $\chi^2$ de Bondad de Ajuste contra la Normal Discretizada (con corrección por continuidad y regla de Cochran):**  
   Para datos discretos y agrupados en contenedores $k \in \{0, \dots, X\}$, la referencia probabilística correcta es la **Normal Discretizada** $\mathcal{N}(\mu = X/2, \sigma^2 = X/4)$ integrada mediante corrección por continuidad sobre cada receptáculo $[k - 0.5, k + 0.5]$ (extendiendo las colas extremas a $(-\infty, 0.5]$ y $[X - 0.5, +\infty)$ para conservar $\sum E_k = 100$).  
   Para garantizar la validez asintótica del estadístico de Pearson $\chi^2 = \sum_{i=1}^{m} \frac{(O_i - E_i)^2}{E_i}$, se aplica el **agrupamiento de colas según la regla de Cochran** ($E_i \ge 5.0$), reduciendo los $X+1$ receptáculos originales a $m$ clases efectivas con $gl = m - 1$ grados de libertad.

3. **Por qué el test $\chi^2$ solo tampoco basta en tableros muy pequeños (Criterio Compuesto):**  
   Con $N = 100$ bolas, la diferencia máxima entre $\text{Binomial}(X, 0.5)$ y la Normal Discretizada es diminuta incluso para $X = 3$, por lo que $\chi^2$ ya acepta el $95.1\%$ de las corridas con $X = 3$ renglones (solo 4 columnas y curtosis severamente platicúrtica $\gamma_2 = -0.667$, lejos de una campana gaussiana). Para evitar este falso positivo estructural, definimos que $X$ alcanza una **normalidad aceptable** cuando cumple simultáneamente cuatro condiciones:
   - **Condición estructural de De Moivre-Laplace:** $X \cdot p \ge 5$ y $X \cdot (1 - p) \ge 5 \implies X \ge 10$ ($\sigma \ge 1.581$ y exceso de curtosis acotado $|\gamma_2| = 2/X \le 0.20$).
   - **Resolución de clases (Cochran):** al menos $m \ge 6$ contenedores efectivos con $E_i \ge 5$.
   - **Bondad de ajuste $\chi^2$ (Normal Discretizada):** tasa de aceptación ($p > 0.05$) $\ge 90\%$ en $R = 1\,000$ réplicas.
   - **Test omnibus de forma D'Agostino-Pearson $K^2$:** tasa de aceptación ($p > 0.05$) $\ge 90\%$ evaluando conjuntamente asimetría y curtosis.

### 2.2 Resultados Numéricos de la Etapa 1 ($X = 3 \dots 15$)

| $X$ | Receptáculos | Clases Eff. ($m$) | $gl$ | $\chi^2$ medio | $p_{\text{med}}$ ($\chi^2$) | Acep. $\chi^2$ (%) | Acep. $K^2$ (%) | Acep. SW crudo (%) | Curtosis teór. ($\gamma_2$) | TVD media | Veredicto |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 3 | 4 | 4 | 3 | 3.05 | 0.483 | 95.1% | 81.3% | 0.0% | -0.667 | 0.0660 | NO ACEPTABLE |
| 4 | 5 | 5 | 4 | 4.13 | 0.497 | 94.3% | 93.3% | 0.0% | -0.500 | 0.0749 | NO ACEPTABLE |
| 5 | 6 | 4 | 3 | 3.05 | 0.482 | 94.8% | 95.3% | 0.0% | -0.400 | 0.0813 | NO ACEPTABLE |
| 6 | 7 | 5 | 4 | 4.02 | 0.489 | 95.4% | 95.3% | 0.0% | -0.333 | 0.0867 | NO ACEPTABLE |
| 7 | 8 | 6 | 5 | 5.09 | 0.481 | 94.4% | 96.9% | 0.0% | -0.286 | 0.0917 | TRANSICIÓN |
| 8 | 9 | 5 | 4 | 3.97 | 0.480 | 95.1% | 96.1% | 0.0% | -0.250 | 0.0967 | TRANSICIÓN |
| 9 | 10 | 6 | 5 | 4.94 | 0.502 | 95.8% | 96.9% | 0.0% | -0.222 | 0.0979 | TRANSICIÓN |
| **10** | **11** | **7** | **6** | **6.02** | **0.515** | **95.0%** | **96.4%** | **0.0%** | **-0.200** | **0.1025** | **NORMAL ACEPTABLE** |
| 11 | 12 | 6 | 5 | 4.90 | 0.494 | 95.7% | 94.7% | 0.0% | -0.182 | 0.1060 | NORMAL ACEPTABLE |
| 12 | 13 | 7 | 6 | 6.36 | 0.451 | 93.6% | 95.6% | 0.4% | -0.167 | 0.1109 | NORMAL ACEPTABLE |
| 13 | 14 | 6 | 5 | 4.90 | 0.504 | 94.9% | 97.0% | 0.5% | -0.154 | 0.1083 | NORMAL ACEPTABLE |
| 14 | 15 | 7 | 6 | 5.94 | 0.500 | 95.5% | 96.5% | 1.9% | -0.143 | 0.1124 | NORMAL ACEPTABLE |
| 15 | 16 | 8 | 7 | 6.92 | 0.501 | 94.8% | 95.8% | 4.9% | -0.133 | 0.1138 | NORMAL ACEPTABLE |

*Nota sobre TVD:* La Distancia de Variación Total muestral (`TVD`) aumenta ligeramente con $X$ porque el tamaño muestral está fijo en $N = 100$ bolas repartidas en cada vez más receptáculos (de 4 a 16), incrementando el ruido multinomial por celda mientras el sesgo estructural platicúrtico $|\gamma_2| = 2/X$ disminuye monótonamente.

### 2.3 Conclusión de la Etapa 1

El mínimo necesario para obtener una distribución normal aceptable con $N = 100$ bolas es **$X_{\min} = 10$ renglones** ($11$ receptáculos, $\mathcal{N}(5.0, 2.50)$):
- A partir de $X = 10$ se satisface la regla clásica de De Moivre-Laplace ($X p = 5.0 \ge 5$), con una curtosis platicúrtica ya acotada ($|\gamma_2| = 0.200$), **7 clases activas de Cochran** ($gl = 6$), **$95.0\%$ de aceptación en $\chi^2$** ($p_{\text{med}} = 0.515$) y **$96.4\%$ de aceptación en D'Agostino-Pearson $K^2$**.
- Gráficos de respaldo:
  - [`output/etapa1_comparacion_columnas.png`](output/etapa1_comparacion_columnas.png)
  - [`output/etapa1_evolucion_normalidad.png`](output/etapa1_evolucion_normalidad.png)
  - [`output/etapa1_detalle_x_min.png`](output/etapa1_detalle_x_min.png)

---

## 3. Etapa 2: Impacto de un Triángulo Defectuoso (70% Izq. / 30% Der.) en el 4.º Renglón

Tomando el tablero óptimo $X_{\min} = 10$ renglones, introducimos un único triángulo defectuoso en el **4.º renglón** (donde existen $4$ triángulos posibles: posiciones $j^* \in \{0, 1, 2, 3\}$ de izquierda a derecha), el cual desvía las bolas con probabilidad $p_L^* = 0.70$ hacia la izquierda y $p_R^* = 0.30$ hacia la derecha.

### 3.1 Fundamento Analítico: Probabilidad de Paso y Desplazamiento Exacto de la Media

Para alcanzar el 4.º renglón, toda bola debe atravesar primero los renglones estándar $r = 1, 2, 3$ ($3$ bifurcaciones simétricas $50/50$). Por lo tanto, el índice de llegada al 4.º renglón sigue una ley $\text{Binomial}(3, 0.5)$:

$$P(\text{visitar nodo } j^* \text{ en renglón 4}) = \binom{3}{j^*} \left(\frac{1}{2}\right)^3 = \begin{cases} 1/8 = 12.5\% & \text{si } j^* \in \{0, 3\} \text{ (extremos)} \\ 3/8 = 37.5\% & \text{si } j^* \in \{1, 2\} \text{ (centrales)} \end{cases}$$

La simulación Monte Carlo confirmó estos flujos nodales con alta precisión: `Pos 0 = 12.40%`, `Pos 1 = 37.20%`, `Pos 2 = 37.71%` y `Pos 3 = 12.68%`.

Cada bola estándar aporta $+0.50$ al índice esperado final, mientras que una bola que impacta el triángulo defectuoso aporta solo $+0.30$ en esa etapa (una diferencia de $\delta = 0.30 - 0.50 = -0.20$ receptáculos). Dado que los renglones posteriores ($r = 5 \dots 10$) son nuevamente simétricos, por linealidad de la esperanza el corrimiento exacto de la media poblacional es:

$$\Delta \mu(j^*) = \mu_{\text{def}}(j^*) - \mu_{\text{base}} = -0.20 \cdot P(\text{visitar } j^*)$$

- **Si el defecto está en los extremos ($j^* = 0$ ó $j^* = 3$):** solo $12.5$ de cada $100$ bolas pasan por él, produciendo un desplazamiento de $\Delta \mu = -0.20 \times 0.125 = -0.025$ ($\mu = 4.975$).
- **Si el defecto está en el centro ($j^* = 1$ ó $j^* = 2$):** pasan $37.5$ de cada $100$ bolas (**el triple de tráfico**), triplicando el desplazamiento medio a $\Delta \mu = -0.20 \times 0.375 = -0.075$ ($\mu = 4.925$).
- **Si la posición es aleatoria uniforme en $\{0, 1, 2, 3\}$:** el desplazamiento esperado marginal es el promedio $\Delta \mu = -0.050$ ($\mu = 4.950$).

### 3.2 Resultados Comparativos de la Etapa 2 ($X = 10$, $N = 100$, $R = 1\,000$)

| Escenario | Tráfico nodo (%) | Bolas afectadas (/100) | $\mu$ exacta | $\Delta\mu$ exacto | $\bar{x}$ MC | $s$ MC | Asimetría exacta ($\gamma_1$) | Asimetría MC | $p_{\text{med}}$ ($\chi^2$ vs. base) | Rechazo $\chi^2$ (%) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Sin defecto (Base) | 100.0% | 0.0 | 5.000 | 0.000 | 5.007 | 1.576 | 0.0000 | -0.002 | 0.504 | 6.2% |
| Posición 0 (Ext. Izq.) | 12.5% | 12.5 | 4.975 | -0.025 | 4.974 | 1.600 | -0.0229 | -0.035 | 0.486 | 5.5% |
| Posición 1 (Centro-Izq.) | 37.5% | 37.5 | 4.925 | -0.075 | 4.929 | 1.599 | +0.0403 | +0.032 | 0.457 | 6.2% |
| Posición 2 (Centro-Der.) | 37.5% | 37.5 | 4.925 | -0.075 | 4.923 | 1.553 | +0.0352 | +0.040 | 0.457 | 6.6% |
| Posición 3 (Ext. Der.) | 12.5% | 12.5 | 4.975 | -0.025 | 4.973 | 1.554 | -0.0280 | -0.023 | 0.490 | 4.1% |
| Posición Aleatoria | 25.0% | 25.0 | 4.950 | -0.050 | 4.951 | 1.577 | +0.0063 | -0.003 | 0.474 | 6.1% |

### 3.3 Análisis Morfológico e Interpretación por Ubicación

1. **Magnitud del impacto (Centrales vs. Extremos):**  
   Los nodos centrales (`Pos 1` y `Pos 2`) interceptan exactamente $3$ veces más bolas que los extremos (`Pos 0` y `Pos 3`), por lo que triplican el corrimiento de la media ($\Delta \mu = -0.075$ vs. $-0.025$) y reducen más el p-valor mediano frente a la distribución normal original ($p_{\text{med}}$ cae de $0.504$ a $0.457$ en centrales vs. $0.486\text{--}0.490$ en extremos).

2. **Efecto opuesto sobre la varianza según el hemisferio del defecto:**  
   Aunque `Pos 0` y `Pos 3` comparten el mismo $\Delta \mu = -0.025$ (y `Pos 1` y `Pos 2` comparten $\Delta \mu = -0.075$), su efecto sobre la dispersión y la forma difiere cualitativamente:
   - **Defecto en la mitad izquierda (`Pos 0` y `Pos 1`):** empuja bolas que ya venían por la izquierda aún más hacia la cola izquierda, **aumentando el desvío estándar** ($s_{\text{MC}} = 1.600$ y $1.599$ frente a $1.576$ en la base). En `Pos 0`, esto engrosa la cola izquierda ($\gamma_1 = -0.0229$); en `Pos 1`, concentra masa entre los receptáculos $k = 3$ y $k = 4$ dejando una cola derecha relativamente más larga ($\gamma_1 = +0.0403$).
   - **Defecto en la mitad derecha (`Pos 2` y `Pos 3`):** desvía hacia la izquierda a bolas que venían hacia la derecha, devolviéndolas hacia el centro del tablero y **reduciendo el desvío estándar** ($s_{\text{MC}} = 1.553$ y $1.554$). En `Pos 2`, concentra masa justo sobre el pico $k = 5$ generando ligera asimetría positiva ($\gamma_1 = +0.0352$); en `Pos 3`, comprime la cola derecha dejando mayor peso relativo a la izquierda ($\gamma_1 = -0.0280$).

3. **Detectabilidad estadística con $N = 100$ bolas:**  
   Dado que el mayor corrimiento medio ($\Delta \mu = -0.075$) equivale apenas al $4.7\%$ de un desvío estándar ($\approx 0.047\sigma$), en una corrida individual de $N = 100$ bolas el ruido muestral propio del experimento domina sobre el sesgo local: la tasa de rechazo $\chi^2$ contra $\mathcal{N}(5.0, 2.50)$ se mantiene entre $4.1\%$ y $6.6\%$ (cercana al nivel nominal $\alpha = 5\%$). Esto explica que un único triángulo defectuoso en el 4.º renglón sea **difícil de detectar en un ensayo único de 100 bolas**, pero su huella sistemática sobre $\Delta \mu$, $s$ y $\gamma_1$ queda plenamente demostrada al analizar la distribución exacta y múltiples réplicas Monte Carlo.

- Gráficos de respaldo:
  - [`output/etapa2_comparacion_columnas_defecto.png`](output/etapa2_comparacion_columnas_defecto.png)
  - [`output/etapa2_analisis_trafico_e_impacto.png`](output/etapa2_analisis_trafico_e_impacto.png)

---

## 4. Estructura del Proyecto y Reproducción

```text
Simulacion_Champagne/
├── main.py                  # Script principal: ejecuta Etapas 1 y 2, imprime tablas y exporta gráficos/CSV
├── requirements.txt         # Dependencias del proyecto
├── INFORME_TECNICO.md       # Informe técnico estructurado listo para presentar
├── src/
│   ├── __init__.py
│   ├── simulator.py         # Motor Monte Carlo vectorizado y programación dinámica exacta
│   ├── statistics.py        # Normal discretizada, test Chi2 (Cochran), D'Agostino K2 y Shapiro-Wilk
│   └── visualization.py     # Gráficos de columnas comparativos y de evolución (300 DPI)
├── tests/
│   └── test_galton.py       # Suite de pruebas unitarias (pytest)
└── output/                  # Gráficos PNG (300 DPI) y tablas resumen CSV
```

Para ejecutar las pruebas unitarias y reproducir íntegramente los experimentos:

```powershell
python -m pytest -v
python main.py
```
