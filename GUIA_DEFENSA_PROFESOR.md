# Guía Completa para Defender el TP "Las Copas de Champagne" (Desde Cero)

> **Cómo usar esta guía:** Está pensada para que entiendas **absolutamente todo desde cero** (sin dar por sabido nada de estadística ni de programación), sepas qué mostrarle en pantalla al profesor con la demo interactiva (`python demo.py`) y tengas respuestas listas para cualquier pregunta que te hagan en la defensa oral.

---

## 0. El Resumen de 30 Segundos (Si solo te acordás de 4 ideas, que sean estas)

1. **¿Qué es la máquina?** Es un **Tablero de Galton**: caen $100$ bolas de a una desde arriba, chocan contra filas ("renglones") de triángulos que las desvían $50\%$ a la izquierda o $50\%$ a la derecha, y caen en $X + 1$ copas/receptáculos abajo.
2. **¿Cuál es la respuesta de la Etapa 1?** Se necesitan como mínimo **$X = 10$ renglones** ($11$ receptáculos) para tener una distribución normal aceptable con $100$ bolas.
3. **¿Cómo lo demostraste estadísticamente?** Usamos el **Test Chi-cuadrado ($\chi^2$) de bondad de ajuste contra la Normal Discretizada** (agrupando las colas con regla de Cochran $E_i \ge 5$), complementado con el test **D'Agostino-Pearson $K^2$** y la regla clásica de **De Moivre-Laplace ($X \cdot p \ge 5 \implies X \ge 10$)**, repitiendo cada experimento **$1\,000$ veces** (Monte Carlo) para no dejarnos engañar por la suerte de una sola corrida.
4. **¿Cuál es la respuesta de la Etapa 2?** Un defecto $70\%$ izq. / $30\%$ der. en el 4.º renglón **impacta el triple si cae en el centro** (Posiciones `1` y `2`, donde pasa el **$37.5\%$** de las bolas y corre la media en **$\Delta \mu = -0.075$**) que si cae en los extremos (Posiciones `0` y `3`, donde pasa solo el **$12.5\%$** y la media se corre **$\Delta \mu = -0.025$**).

---

## 1. El Experimento Físico Explicado Como si Nunca Hubieras Visto Estadística

Imaginá una caja vertical transparente (entre dos vidrios) llena de clavitos en forma de triángulo:

- En el **renglón 1** (arriba de todo) hay **1 triángulo**.
- En el **renglón 2** hay **2 triángulos**.
- En el **renglón 3** hay **3 triángulos**.
- ...
- En el **renglón $X$** hay **$X$ triángulos**.
- Debajo del último renglón hay **$X + 1$ copas o receptáculos** numerados $0, 1, 2, \dots, X$ de izquierda a derecha.

*(Por ejemplo, en el dibujo de la hoja del TP hay $X = 4$ renglones y $5$ receptáculos).*

### ¿Por qué una bola termina en una copa y no en otra?
Cada vez que una bola golpea un triángulo sano, tiene la misma chance de ir hacia la izquierda o hacia la derecha: **$50\%$ y $50\%$**, exactamente igual que tirar una moneda al aire.

Pensá qué significa caer en cada copa después de pasar por $X$ renglones:
- Para caer en la **copa 0** (extremo izquierdo), la bola tuvo que doblar **siempre a la izquierda** ($0$ veces a la derecha). Hay **un solo camino** posible para hacer eso, por eso caen muy pocas bolas ahí.
- Para caer en la **copa del medio** (por ejemplo, la copa `5` cuando hay $X = 10$ renglones), la bola tuvo que doblar $5$ veces a la derecha y $5$ veces a la izquierda, en **cualquier orden**. ¡Hay **252 caminos distintos** para hacer eso! Por eso las copas del medio se llenan muchísimo más que los extremos.

### La conexión mágica con la Estadística
Contar cuántas veces dobló a la derecha una bola en $X$ intentos con probabilidad $0.5$ se llama **Distribución Binomial** $\text{Binomial}(X, p=0.5)$.

Cuando $X$ es grande, esa montaña de barras apiladas en el medio y bajita a los costados forma una **Campana de Gauss (Distribución Normal)**:
- **Centro de la campana (Media $\mu$):** es exactamente la mitad de los renglones: $\mu = \frac{X}{2}$ (para $X=10$, el centro está en la copa $5$).
- **Ancho de la campana (Desvío estándar $\sigma$):** mide qué tan desparramadas caen las bolas: $\sigma = \frac{\sqrt{X}}{2}$ (para $X=10$, $\sigma \approx 1.58$).

---

## 2. Traductor de Conceptos Estadísticos (Sin Miedo a las Palabras Difíciles)

Si el profesor nombra alguno de estos términos, acá tenés qué significa en español simple:

| Término que suena difícil | Qué significa realmente | Para qué lo usamos en el TP |
|---|---|---|
| **Muestra de $N = 100$ bolas** | Cada vez que hacés una corrida tirando 100 bolitas. | Es el tamaño fijo que pide el enunciado. |
| **Simulación Monte Carlo ($R = 1\,000$ réplicas)** | En vez de tirar 100 bolas una sola vez (donde podés tener mala o buena suerte), la computadora repite el experimento $1\,000$ veces. | Evita sacar conclusiones falsas por azar de una sola tirada. |
| **Normal Discretizada (con corrección por continuidad)** | La campana normal es una curva continua (con decimales), pero nuestras copas son números enteros ($0, 1, 2\dots$). Entonces calculamos cuánta área de la campana cae dentro de cada copa (por ej. copa `5` toma de `4.5` a `5.5`). | Permite comparar correctamente la campana continua teórica contra las barras enteras de las copas. |
| **Test Chi-cuadrado ($\chi^2$)** | Mide la diferencia total entre **lo que observamos ($O_i$)** en las copas y **lo que esperábamos ($E_i$)** según la campana normal. | Es la herramienta principal de bondad de ajuste para datos en categorías/contenedores. |
| **Regla de Cochran ($E_i \ge 5$)** | El test $\chi^2$ se rompe si una copa de la esquina espera muy poquitas bolas (por ej. $0.2$ bolas). Cochran dice: *"juntá las copas de los extremos hasta que cada grupo espere al menos $5$ bolas"*. | Hace que el test $\chi^2$ sea matemáticamente válido en las colas. |
| **p-valor ($p > 0.05$)** | Es la probabilidad de ver una diferencia igual o más grande que la observada si la campana fuera realmente normal. **Ojo:** acá **NO queremos rechazar** la normalidad; o sea, **si $p > 0.05$, el ajuste es ACEPTABLE**. | Si en $\ge 90\%$ de las 1.000 réplicas da $p > 0.05$, la distribución ajusta bien. |
| **Test de Shapiro-Wilk** | Es un test de normalidad muy famoso, pero está hecho para **números continuos con decimales** (como alturas o pesos), donde casi nunca se repite un número exacto. | En nuestras copas solo hay entre $4$ y $16$ números posibles para $100$ bolas: ¡se repiten muchísimo! Por eso Shapiro-Wilk crudo rechaza casi siempre ($0\%$). |
| **Asimetría (*Skewness* $\gamma_1$)** | Mide si la campana está torcida hacia un lado ($0$ = perfectamente simétrica). | En Etapa 1 es $0$. En Etapa 2 cambia porque el defecto desvía más hacia un lado. |
| **Exceso de Curtosis ($\gamma_2$)** | Mide si la curva es más "puntiaguda" o más "achatada" que una normal ($\gamma_2 = 0$). En la binomial es negativa: $\gamma_2 = -2/X$ (más chata). | Con $X = 3$, $\gamma_2 = -0.667$ (súper chata). Recién en $X = 10$ baja a $|\gamma_2| \le 0.20$. |
| **Regla de De Moivre-Laplace ($X \cdot p \ge 5$)** | Teorema clásico de probabilidad que dice cuándo una Binomial puede aproximarse bien por una Normal: exige que $X \cdot p \ge 5$ y $X \cdot (1-p) \ge 5$. | Como $p = 0.5$, $X \cdot 0.5 \ge 5$ da directamente **$X \ge 10$**. |

---

## 3. Cómo Explicar la Etapa 1 (Determinación de $X_{\min} = 10$)

### El argumento ganador frente al profesor
1. **¿Qué probamos?** Simulamos tableros desde $X = 3$ hasta $X = 15$ renglones, haciendo $1\,000$ corridas de $N = 100$ bolas para cada $X$.
2. **¿Qué herramienta elegimos y por qué?**
   - Elegimos el **Test $\chi^2$ de bondad de ajuste contra la Normal Discretizada** (con agrupamiento de colas de Cochran $E_i \ge 5$) porque las bolas caen en **contenedores discretos**, no en una recta continua.
   - Descartamos **Shapiro-Wilk crudo** como criterio principal porque al haber $100$ bolas en pocos receptáculos enteros se generan más del $85\%$ de empates (*ties*), haciendo que Shapiro-Wilk dé $0.0\%$ de aceptación hasta $X = 11$ y apenas $4.9\%$ en $X = 15$. *(Decir esto demuestra que realmente entendiste qué test va para datos discretos vs. continuos).*
3. **¿Por qué no alcanza con $X = 3$ o $X = 4$ aunque $\chi^2$ dé $p > 0.05$?**
   - Este es el punto más inteligente del análisis: si mirás solo $\chi^2$ contra la normal discretizada en $X = 3$, ya acepta el $95.1\%$ de las veces porque con apenas $100$ bolas y $4$ barras, $\chi^2$ tiene poca potencia para distinguir un bloque platicúrtico de una campana real.
   - Pero visual y estructuralmente, en $X = 3$ y $X = 4$ (el dibujo del TP) la distribución tiene muy pocos escalones ($4$ o $5$ copas) y un exceso de curtosis muy negativo ($\gamma_2 = -0.667$ y $-0.500$, muy achatada).
4. **El criterio compuesto que resuelve el problema:**
   Exigimos cuatro condiciones simultáneas para hablar de **"normalidad aceptable"**:
   - Bondad de ajuste **$\chi^2$** $\ge 90\%$ de aceptación ($p > 0.05$).
   - Test de forma **D'Agostino-Pearson $K^2$** $\ge 90\%$ de aceptación (mide asimetría y curtosis juntas).
   - Al menos **$6$ clases efectivas** tras agrupar colas con Cochran ($E_i \ge 5$).
   - Regla estructural de **De Moivre-Laplace**: $X \cdot p \ge 5$ (y $|\gamma_2| = 2/X \le 0.20$).
5. **Conclusión exacta:**
   - De $X = 3$ a $6$: **NO ACEPTABLE**.
   - De $X = 7$ a $9$: **ZONA DE TRANSICIÓN** (ya pasan $\chi^2$ y $K^2$, pero la curtosis sigue por encima de $0.20$ y $X \cdot p < 5$).
   - Desde **$X_{\min} = 10$ renglones ($11$ receptáculos)**: **NORMAL ACEPTABLE** ($\chi^2$ acepta el **$95.0\%$**, $p_{\text{med}} = 0.515$; D'Agostino $K^2$ acepta el **$96.4\%$**; $7$ clases activas de Cochran; $|\gamma_2| = 0.200$ y $X \cdot p = 5.0 \ge 5$).

---

## 4. Cómo Explicar la Etapa 2 (Triángulo Defectuoso 70/30 en el 4.º Renglón)

Una vez fijado $X = 10$, nos piden poner **un triángulo defectuoso en el 4.º renglón** que desvía **$70\%$ a la izquierda** y **$30\%$ a la derecha** (en lugar de $50/50$). En el 4.º renglón hay **4 triángulos posibles** (llamémoslos `Pos 0`, `Pos 1`, `Pos 2` y `Pos 3` de izquierda a derecha).

### La clave física de toda la Etapa 2: ¿Cuántas bolas tocan cada triángulo?
Antes de llegar al 4.º renglón, las bolas atraviesan $3$ renglones sanos ($r = 1, 2, 3$). Es como tirar $3$ monedas:
- Hay $2^3 = 8$ caminos posibles:
  - `Pos 0` (extremo izquierdo, 3 izquierdas seguidas): **$1/8 = 12.5\%$** de las bolas.
  - `Pos 1` (centro-izquierda, 2 izq. y 1 der.): **$3/8 = 37.5\%$** de las bolas.
  - `Pos 2` (centro-derecha, 1 izq. y 2 der.): **$3/8 = 37.5\%$** de las bolas.
  - `Pos 3` (extremo derecho, 3 derechas seguidas): **$1/8 = 12.5\%$** de las bolas.

*(Nuestra simulación lo verificó casi al decimal: `12.40%`, `37.20%`, `37.71%` y `12.68%`).*

### ¿Cuánto corre la media según la posición real del triángulo?
- En un triángulo sano, ir a la derecha ($+1$) ocurre $50\%$ de las veces $\implies$ aporte promedio $+0.50$.
- En el triángulo defectuoso, ir a la derecha ($+1$) ocurre solo $30\%$ de las veces $\implies$ aporte promedio $+0.30$.
- O sea, cada bola que toca el defecto pierde en promedio $0.30 - 0.50 = -0.20$ posiciones (se corre $0.20$ hacia la izquierda).

Multiplicando esa pérdida por la fracción de bolas que pasan por el nodo:

$$\Delta \mu = -0.20 \times P(\text{visitar el triángulo})$$

- **En los extremos (`Pos 0` o `Pos 3`):** $\Delta \mu = -0.20 \times 0.125 = \mathbf{-0.025}$ (la media baja apenas de $5.000$ a $4.975$).
- **En el centro (`Pos 1` o `Pos 2`):** $\Delta \mu = -0.20 \times 0.375 = \mathbf{-0.075}$ (la media baja de $5.000$ a $4.925$ — **¡impacta exactamente el TRIPLE!**).
- **Si el defecto cae en un lugar al azar del 4.º renglón:** el promedio de los cuatro casos es $\Delta \mu = \mathbf{-0.050}$ ($\mu = 4.950$).

### ¡Detalle de 10 puntos! Efecto en la dispersión (varianza) y en la forma
Aunque `Pos 0` y `Pos 3` corren igual la media (`-0.025`), y `Pos 1` y `Pos 2` también empatan en media (`-0.075`), **no hacen lo mismo con el ancho de la campana**:
- **Si el defecto está en la mitad izquierda (`Pos 0` o `Pos 1`):** manda hacia la izquierda a bolas que *ya venían por la izquierda*, alejándolas todavía más del centro. Resultado: **aumenta el desvío estándar** ($s_{\text{MC}} = 1.600$ y $1.599$ frente a $1.576$ sin defecto).
- **Si el defecto está en la mitad derecha (`Pos 2` o `Pos 3`):** manda hacia la izquierda a bolas que venían yéndose hacia la derecha, es decir, **las devuelve hacia el centro**. Resultado: **achica el desvío estándar** ($s_{\text{MC}} = 1.553$ y $1.554$).

### ¿Por qué el test $\chi^2$ con $N = 100$ bolas casi nunca "salta" frente a este defecto?
Porque el desplazamiento máximo ($\Delta \mu = -0.075$) es apenas el **$4.7\%$ de un desvío estándar** ($0.075 / 1.581 \approx 0.047\sigma$), y debajo del 4.º renglón todavía quedan **6 renglones sanos** que vuelven a mezclar las bolas. Con solo $N = 100$ bolas, el ruido natural del azar tapa ese pequeño corrimiento (la tasa de rechazo de $\chi^2$ contra la normal original queda entre $4.1\%$ y $6.6\%$). Recién al analizar la distribución exacta y promediar $1\,000$ réplicas aparece con claridad el efecto real.

---

## 5. Paso a Paso para Hacer la Demo en Vivo Frente al Profesor

1. Abrí una terminal en la carpeta del proyecto y ejecutá:
   ```powershell
   python demo.py
   ```
   Se va a abrir `demo_interactiva.html` en tu navegador (también podés abrir directamente la pestaña **4. Guion para Mostrar al Profesor** dentro de la misma página, que tiene botones que cargan cada caso automáticamente).
2. **En la Pestaña 1 (`Simulador en Vivo`):**
   - Seleccioná **`X = 4 (TP)`** y apretá **`🔍 Lanzar 1 bola en cámara lenta`** para mostrarle visualmente cómo cae y elige izquierda/derecha en cada triángulo. Luego apretá **`▶ Lanzar 100 bolas`** para ver que con 4 renglones quedan barras muy gruesas y achatadas.
   - Seleccioná **`X = 10 ⭐`** y apretá **`⚡ Simular 1.000 corridas (Monte Carlo)`**: mostrale cómo la barra verde dice **`NORMAL ACEPTABLE`** ($95\%$ de aceptación $\chi^2$ y alineamiento perfecto con la curva roja).
   - Sin salir de $X = 10$, hacé clic en **`Pos 0 (Ext. Izq 12.5%)`** (o directamente sobre el triángulo izquierdo del 4.º renglón en el dibujo) y corré **`⚡ Simular 1.000 corridas`**: mostrale en la consola de abajo que pasa ~12.5% de las bolas y $\Delta \mu = -0.025$.
   - Ahora hacé clic en **`Pos 1 (Centro 37.5%)`** y corré **`⚡ Simular 1.000 corridas`**: mostrale que el tráfico sube a ~37.5% y $\Delta \mu = -0.075$ (el triple).
3. **Mostrá las Pestañas 2 y 3:** tienen las tablas completas y los gráficos de 300 DPI listos para hacer clic y ampliar.

---

## 6. Las 10 Preguntas Clásicas del Profesor (y Qué Responder)

1. **¿Por qué hiciste $1\,000$ réplicas por cada $X$ si el enunciado dice $100$ bolas?**  
   *Respuesta:* Porque $N = 100$ bolas es una muestra relativamente chica. En una sola corrida, por pura suerte puede pasarte que $X = 4$ parezca normal o que $X = 12$ rechace por azar (un $5\%$ de las veces). Con $1\,000$ réplicas medimos la **tasa consistente** de aceptación de cada $X$.

2. **¿Por qué elegiste el test $\chi^2$ y no Kolmogorov-Smirnov o Shapiro-Wilk?**  
   *Respuesta:* Porque la salida del tablero son **conteos discretos en $X+1$ recipientes**, no una variable continua. Shapiro-Wilk y KS suponen continuidad sin empates; al tener $100$ bolas agrupadas en pocos enteros hay $>85\%$ de valores repetidos y Shapiro-Wilk crudo da $0\%$ de aceptación hasta $X = 11$. En cambio $\chi^2$ contra la **Normal Discretizada** (integrada en $[k-0.5, k+0.5]$) y agrupando colas con Cochran ($E_i \ge 5$) está hecho exactamente para frecuencias en contenedores.

3. **¿Qué es la corrección por continuidad en la Normal Discretizada?**  
   *Respuesta:* Como el receptáculo $k$ es discreto y la Normal es continua, a la copa $k$ le asignamos la probabilidad de la campana entre $k - 0.5$ y $k + 0.5$ (y en las dos puntas extendemos a $-\infty$ y $+\infty$) para que la suma de probabilidades esperadas dé exactamente $100\%$.

4. **Si $\chi^2$ ya da $95\%$ de aceptación en $X = 3$, ¿por qué no dijiste que $X_{\min} = 3$?**  
   *Respuesta:* Porque con $N = 100$ bolas y solo $4$ contenedores, $\chi^2$ no tiene potencia para detectar que la forma real tiene una curtosis de $-0.667$ (muy achatada). Por eso exigimos además la condición clásica de **De Moivre-Laplace** ($X \cdot p \ge 5 \implies X \ge 10$), que acota $|\gamma_2| \le 0.20$, junto con al menos $6$ clases efectivas de Cochran y $\ge 90\%$ en D'Agostino-Pearson $K^2$.

5. **¿Qué pasa si $N$ fuera mucho más grande, por ejemplo $N = 10\,000$ bolas?**  
   *Respuesta:* Con $N = 10\,000$, el ruido muestral se achica muchísimo y el test $\chi^2$ detectaría incluso la mínima diferencia entre la Binomial($X, 0.5$) y la Normal Discretizada, aumentando el poder estadístico del test (tanto en Etapa 1 como para detectar el defecto en Etapa 2).

6. **En la Etapa 2, ¿por qué los triángulos del medio del 4.º renglón afectan más que los de las puntas?**  
   *Respuesta:* Porque para llegar al 4.º renglón las bolas pasan antes por $3$ filas sanas. La distribución de llegadas al 4.º renglón es $\text{Binomial}(3, 0.5)$: al centro (`Pos 1` y `Pos 2`) llega el **$37.5\%$** de las bolas en cada nodo, mientras que a los extremos (`Pos 0` y `Pos 3`) solo llega el **$12.5\%$**. Como reciben el triple de tráfico, triplican el corrimiento de la media ($\Delta \mu = -0.075$ vs. $-0.025$).

7. **¿Por qué $\Delta \mu$ da negativo en todas las posiciones del defecto?**  
   *Respuesta:* Porque el defecto desvía el $70\%$ hacia la izquierda (índice menor) y solo $30\%$ hacia la derecha (índice mayor). Siempre empuja el promedio hacia la izquierda respecto de $\mu = 5.0$.

8. **¿Qué diferencia hay entre poner el defecto en `Pos 1` o en `Pos 2` si ambos tienen el mismo tráfico ($37.5\%$) y el mismo $\Delta \mu$ ($-0.075$)?**  
   *Respuesta:* Cambia la **varianza**. `Pos 1` está a la izquierda del eje central y desvía bolas aún más hacia la izquierda, aumentando la dispersión ($s = 1.599$). `Pos 2` está a la derecha del eje central y desvía bolas hacia la izquierda (es decir, de vuelta hacia el centro), reduciendo la dispersión ($s = 1.553$).

9. **¿Qué pasa si la ubicación del defecto dentro del 4.º renglón es aleatoria?**  
   *Respuesta:* Cada una de las 4 posiciones tiene $25\%$ de probabilidad de ser elegida. El efecto promedio esperado es la media de los cuatro casos: $\Delta \mu = -0.050$ ($\mu = 4.950$), afectando en promedio al $25\%$ de las bolas por corrida.

10. **¿Cómo verificaste que el simulador no tiene errores de programación?**  
    *Respuesta:* Además de las $1\,000$ réplicas Monte Carlo, implementamos el **cálculo exacto por programación dinámica** (propagando las probabilidades renglón por renglón) y una suite de **7 tests unitarios con `pytest`** que comprueban la conservación exacta de las $100$ bolas, la coincidencia analítica $\Delta \mu = -0.20 \cdot P(\text{visita } j^*)$ y la convergencia del promedio Monte Carlo a la distribución exacta.
