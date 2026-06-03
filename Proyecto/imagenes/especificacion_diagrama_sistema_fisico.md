# Especificación del diagrama del sistema físico — PMSM con reductora y carga péndulo

Este documento describe completamente el contenido y el estilo del diagrama del sistema físico del Proyecto Global Integrador (Automática y Máquinas Eléctricas — UNCuyo). Está pensado para servir como prompt o guía para modelos de generación de imágenes, herramientas de diagramación (Mermaid, Draw.io, TikZ, Inkscape, etc.) o reescritura manual del SVG.

---

## 1. Descripción general

Diagrama técnico, estilo libro de ingeniería eléctrica/mecatrónica, formato horizontal (16:9 o similar), fondo blanco, **sin título**.

Muestra el sistema físico completo de un accionamiento de CA con motor sincrónico de imanes permanentes (PMSM) trifásico, acoplado a una carga péndulo de un grado de libertad a través de una reductora planetaria, **sin ningún bloque de control** — solo la planta física más los sensores.

El diagrama está dividido visualmente en **tres subsistemas**, agrupados con marcos de fondo de colores suaves distintos y bordes punteados, con sus respectivos títulos arriba:

1. **Subsistema Eléctrico (Potencia)** — fondo gris muy claro
2. **Subsistema Electromagnético + Térmico** — fondo celeste muy claro
3. **Subsistema Mecánico (Transmisión + Carga)** — fondo amarillo muy claro

---

## 2. Componentes (de izquierda a derecha)

### 2.1 Subsistema Eléctrico

**(a) Bus DC / batería**
- Símbolo estándar de batería con 3 celdas verticales
- Terminales marcados con "+" arriba y "−" abajo
- Etiqueta a la derecha: $V_{DC}$
- Texto debajo: "Bus DC"

**(b) Cableado DC**
- Dos rieles horizontales: uno superior (+) y uno inferior (−) que conectan la batería al inversor

**(c) Inversor trifásico**
- Caja rectangular blanca con borde negro
- Etiqueta superior: "Inversor 3φ"
- Etiqueta inferior: "6 llaves IGBT/MOSFET · PWM"
- Internamente muestra: 3 piernas verticales, cada una con 2 llaves (rectángulos grises) en serie, formando la topología clásica de inversor de 6 llaves
- Entradas: 2 terminales en el lado izquierdo conectados al bus DC
- Salidas: 3 terminales (puntos negros) en el lado derecho, a tres alturas distintas, correspondientes a las fases a, b, c

### 2.2 Subsistema Electromagnético + Térmico

**(d) Cables trifásicos de fase**
- Tres líneas horizontales paralelas que van del inversor al PMSM, a tres alturas distintas
- Etiqueta arriba de las líneas: $v_{as}(t), v_{bs}(t), v_{cs}(t)$

**(e) Sensores de corriente (3 CTs)**
- Tres pequeños círculos color ámbar/dorado con borde anaranjado, atravesando cada una de las tres líneas de fase
- Cada uno marcado con "CT" en su interior
- Simbolizan sensores de efecto Hall

**(f) Motor PMSM (representación interna como 3 bobinados en Y)**
- Círculo grande blanco con borde negro (representa el estator)
- En su interior, en la mitad izquierda: **tres bobinados de fase en estrella (Y)**, cada uno compuesto por:
  - Una resistencia $R_s(T_s)$ (rectángulo) en serie con una inductancia $L_s$ (símbolo de bobina con 4 lazos consecutivos)
- Los tres bobinados están dispuestos paralelos horizontalmente (uno por fase)
- A la derecha de los bobinados, una **línea vertical interna** une los tres extremos formando el punto neutro
- El **punto neutro** se indica con un pequeño círculo abierto (no rellenado) etiquetado "n (flotante)"
- **Tres terminales** marcados con puntos negros en el borde izquierdo del estator, identificados como "a", "b", "c"
- Etiqueta encima del motor: "PMSM"
- Etiqueta debajo del motor: "estator en Y - neutro NO accesible"
- **NO hay imanes ni rotor dibujado**, solo el estator con sus bobinados

**(g) Sensor de temperatura (RTD)**
- Símbolo de termómetro pequeño (cuerpo rectangular con bulbo circular abajo) apoyado en la cara superior del estator
- Etiqueta "RTD" al lado

### 2.3 Subsistema Mecánico

**(h) Eje del motor**
- Rectángulo horizontal gris que sale por el lado derecho del estator
- Pequeño acoplamiento al final (rectángulo más oscuro)
- Etiquetas debajo del eje: $T_m, \omega_m, \theta_m$ y debajo $J_m, b_m$ (eje motor)

**(i) Encoder / Resolver**
- Caja rectangular color ámbar claro con borde anaranjado, ubicada arriba del eje del motor
- Etiqueta: "Encoder / Resolver"
- Conexión al eje mediante una línea punteada vertical corta (representa que mide la rotación del eje)

**(j) Reductora planetaria**
- Forma trapezoidal gris
- Pequeñas marcas en la cara superior representando dientes de engranaje
- Etiqueta interior: "Reductora planetaria"
- Etiqueta debajo en rojo: "r = 120 : 1"

**(k) Eje de salida (hacia la carga)**
- Rectángulo horizontal gris similar al del motor, del lado de salida de la reductora
- Termina en un círculo gris oscuro (acoplamiento al pivote del péndulo)
- Etiquetas debajo: $T_q, \omega_l, \theta_l$ y "(eje carga)"

**(l) Carga péndulo (1 g.d.l.)**
- Pivote (círculo negro pequeño) en la base superior, con una pequeña barra horizontal de apoyo
- **Brazo rígido** rotado aproximadamente 28° respecto a la vertical, hacia abajo y a la derecha
  - Rectángulo alargado gris claro con borde negro
  - Pequeño círculo negro en el medio del brazo (centro de masa)
- **Masa de carga útil** $m_l$ al extremo del brazo: círculo grande amarillo con etiqueta "$m_l$" en su centro
- **Referencia vertical**: línea punteada gris vertical que baja desde el pivote
- **Ángulo $\theta_l$**: arco rojo desde la vertical hasta el brazo, con etiqueta "$\theta_l$" en rojo
- **Etiqueta del brazo**: callout que apunta a la mitad del brazo con texto "$m, l_{cm}, J_{cm}$"
- **Dimensión $l_l$**: línea con marcas en los extremos paralela al brazo (a la izquierda), midiendo desde el pivote hasta el centro de la masa, etiqueta "$l_l$"
- **Vector gravedad**: flecha azul vertical hacia abajo a la derecha del péndulo, etiqueta "$g$"
- **Fricción $b_l$**: etiqueta gris al lado del pivote
- **Torque de perturbación $T_{ld}$**: flecha curvada roja en la masa, etiqueta "$T_{ld}$"

---

## 3. Sensores físicos (3 totales) y sus salidas

Cada uno conectado por **líneas punteadas color naranja/ámbar** desde el sensor físico hacia un pequeño cuadro de salida en la parte superior del diagrama.

| Sensor | Mide | Posición |
|---|---|---|
| 3× CT (Hall) | $i_{as}(t), i_{bs}(t), i_{cs}(t)$ | Sobre las 3 líneas de fase entre inversor y motor |
| 1× RTD (sonda térmica) | $T_s°(t)$ | Apoyado sobre el estator del motor |
| 1× Encoder/Resolver | $\theta_m(t)$ **(solo posición, NO velocidad)** | Acoplado al eje del motor |

Los tres cuadros de salida en la parte superior (fondo amarillo claro, borde anaranjado) muestran:

- "SENSORES DE CORRIENTE" → $i_{as}(t), i_{bs}(t), i_{cs}(t)$
- "SONDA TÉRMICA" → $T_s°(t)$
- "ENCODER / RESOLVER" → $\theta_m(t)$

**Importante**: la velocidad $\omega_m$ NO se mide; se reconstruye en el controlador mediante un observador reducido de Luenberger.

---

## 4. Conexiones (flujo del sistema)

```
Bus DC (V_DC)
   ↓ (rieles + y −)
Inversor 3φ (PWM)
   ↓ (3 fases: v_as, v_bs, v_cs)
3× CT (medición de corriente)
   ↓
PMSM (3 bobinados R+L en Y con neutro flotante)
   ↓ (eje motor: T_m, ω_m, θ_m)
Encoder mide θ_m
   ↓
Reductora planetaria (r = 120:1)
   ↓ (eje carga: T_q, ω_l, θ_l)
Carga péndulo 1 g.d.l. (brazo + masa, sometido a g y a T_ld)
```

Adicionalmente:

- La sonda RTD mide $T_s°$ en el estator del PMSM
- El controlador (NO mostrado) recibiría las señales de los 3 CTs, la RTD y el encoder

---

## 5. Estilo visual

### Paleta de colores (limitada y plana, sin gradientes ni sombras)

- **Negro / gris oscuro**: trazos y componentes principales
- **Gris claro**: rellenos de componentes mecánicos (eje, reductora, brazo del péndulo)
- **Naranja / ámbar**: sensores y sus líneas de salida (señales medidas)
- **Azul**: vector gravedad
- **Rojo**: ángulo $\theta_l$ y torque de perturbación $T_{ld}$
- **Amarillo**: masa de carga del péndulo

### Tipografía

- Familia serif (estilo Times New Roman)
- Variables en itálica
- Etiquetas y títulos en regular

### Líneas

- Trazos negros continuos: conexiones eléctricas y mecánicas
- Trazos punteados naranjas: señales medidas por sensores
- Trazos punteados grises: referencias geométricas (vertical de referencia del péndulo, dimensiones)

### Reglas importantes

- **Sin imanes** dibujados (representación simplificada del estator)
- **Sin bloques de control**, **sin lazos cerrados**, **sin diagramas de control**: es solo la planta física + sensores
- **Sin título principal** del diagrama
- **Sin caja de "parámetros del sistema"** (los parámetros se incluyen por separado en una tabla del informe)
- **Sin caja de leyenda de sensores** (los sensores están directamente en sus posiciones físicas con etiquetas claras)

---

## 6. Resumen de variables y símbolos

| Símbolo | Significado |
|---|---|
| $V_{DC}$ | Tensión del bus de continua del inversor |
| $v_{as}, v_{bs}, v_{cs}$ | Tensiones de fase entregadas por el inversor |
| $i_{as}, i_{bs}, i_{cs}$ | Corrientes de fase (medidas por sensores CT) |
| $R_s(T_s)$ | Resistencia del estator dependiente de la temperatura |
| $L_s$ | Inductancia del bobinado del estator |
| $n$ | Punto neutro (no accesible, flotante) |
| $T_s°$ | Temperatura del estator (medida por sonda RTD) |
| $T_m$ | Torque electromagnético del motor |
| $\omega_m$ | Velocidad angular del eje del motor (reconstruida por observador) |
| $\theta_m$ | Posición angular del eje del motor (medida por encoder) |
| $J_m, b_m$ | Inercia y fricción viscosa del rotor + caja |
| $r$ | Relación de transmisión de la reductora (120 : 1) |
| $T_q, \omega_l, \theta_l$ | Torque, velocidad y posición en el eje de la carga |
| $m, l_{cm}, J_{cm}$ | Masa, longitud al CM, e inercia del brazo |
| $m_l, l_l$ | Masa de carga útil y longitud al extremo |
| $b_l$ | Fricción viscosa de la articulación |
| $g$ | Aceleración de la gravedad |
| $T_{ld}$ | Torque de perturbación externa en la carga |

---

## 7. Notas finales para herramientas de generación

- Para **modelos de imagen generativos** (DALL-E 3, Midjourney, Stable Diffusion): este tipo de diagrama técnico exacto suele dar resultados imprecisos. Conviene usar herramientas de diagramación específicas (Mermaid, Draw.io, TikZ) o generar el SVG por código.
- Para **GPT-4 con SVG/HTML**: pegarle esta especificación completa y pedirle que produzca un SVG limpio en una sola pasada suele funcionar bien.
- Para **TikZ / LaTeX**: convendría adaptar componentes como el inversor (CircuitTikZ tiene macros de IGBT) y el motor (símbolos estándar de máquina eléctrica).
- Para **Inkscape / Draw.io manual**: usar esta especificación como checklist al momento de armar el diagrama componente por componente.
