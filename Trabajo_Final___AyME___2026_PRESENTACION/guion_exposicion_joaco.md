# Guion cronometrado — Exposición en mesa
### Joaquín Calderón · Proyecto Global Integrador AyME 2026

**Tus páginas: 1–31, 66–117, 153–165 (96 de 182) · Tu tiempo: ~21 min de los 40**

> ✅ **Numeración verificada** contra `main_presentacion.tex` compilado el 30/07: **182 páginas**.
> ⚠️ El `main_presentacion.pdf` que está en la carpeta es del **29/07 18:47 y quedó viejo** (184 páginas: todavía tiene las dos láminas de métricas de transitorio). **Recompilalo y aseguráte de que Francisco presente desde el mismo archivo**, o los números de lámina no van a coincidir entre ustedes.

---

## Regla de oro

13 segundos por lámina es el promedio bruto. **No hables de cada diapositiva.** El guion está armado en **26 beats**: cada beat es una idea, y las láminas del beat se pasan mientras hablás. Ensayá los beats, no las láminas.

Tres decisiones que conviene respetar:

- **Primero el concepto, después el número.** Nunca leas una ecuación en voz alta. Señalá qué representa y pasá.
- **Las láminas de "Tabla de contenidos" (4, 7, 15, 22, 26, 66, 90) se pasan en silencio.** Son 7 páginas gratis.
- Los números que **sí** decís están en **negrita**. Son ~15 en total. El resto están en pantalla y no hace falta pronunciarlos.

**Total: 20:40**, con 20 s de margen dentro de tus 21 minutos.

---

# BLOQUE A — Modelado del sistema físico
### Páginas 1–31 · **7:15**

---

### A1 · pág. 1–3 · `0:40` — Apertura

> Buenos días. Vamos a presentar el proyecto integrador: el control de un accionamiento de corriente alterna con motor sincrónico de imanes permanentes.
>
> El sistema es un accionamiento de cuatro cuadrantes que mueve la articulación de un brazo manipulador. La cadena es: bus de continua, inversor trifásico, la PMSM, una reductora planetaria **120 a 1**, y el brazo.
>
> El trabajo recorre cinco etapas: modelamos, linealizamos, analizamos la estructura del sistema, diseñamos el control y verificamos, hasta llegar al código embebido. Yo cubro el modelado y la primera parte del control; Francisco continúa después.

`pág. 4 → pasar en silencio`

---

### A2 · pág. 5–6 · `0:55` — ⚓ El sistema y la carga

> Este es el sistema completo. Tres cosas importan acá.
>
> Primero: se divide en **tres subsistemas acoplados** — electromagnético, mecánico y térmico.
>
> Segundo: medimos tres cosas — las corrientes de fase, la posición del eje del motor con un encoder, y la temperatura del bobinado con un RTD.
>
> Y tercero, esto es clave: **lo que queremos controlar, la posición de la articulación, no se mide.** Se mide del lado del motor y se divide por 120.

*(pág. 6)*

> La carga es un péndulo rígido de un grado de libertad en plano vertical. La gravedad entra como perturbación permanente que depende de la posición, y la masa útil es variable **entre 0 y 1,5 kg y no se conoce a priori**. Esa es la incertidumbre principal de todo el proyecto.

`pág. 7 → pasar en silencio`

---

### A3 · pág. 8–11 · `0:35` — Los tres balances *(pasar rápido las 4 láminas)*

> El modelado mecánico son tres balances de torque: el del péndulo, el del tren de transmisión y el del rotor. Como asumimos **transmisión rígida** —sin backlash y sin elasticidad torsional— el torque de reacción del reductor se elimina algebraicamente y los tres se combinan en uno solo.

---

### A4 · pág. 12–14 · `1:10` — ⚓ El modelo mecánico equivalente

> Y queda esta única ecuación referida al eje del motor. Conceptualmente: motor, caja y carga se comportan como **una sola inercia con una sola fricción**, más un término no lineal de gravedad.
>
> El detalle que vale la pena mirar es el factor **1 sobre r cuadrado**. La inercia y la fricción de la carga se reflejan al eje del motor divididas por 120 al cuadrado, o sea por **14.400**. Eso tiene dos consecuencias.
>
> La carga queda enormemente diluida: la fricción del brazo aporta 7 millonésimas frente a las 15 del motor. Despreciable.
>
> Pero la inercia **sí** importa: la inercia equivalente pasa de **1,98 a 4,58 por diez a la menos cinco** cuando cargamos el brazo al máximo. Se multiplica por 2,3. Y este es el parámetro que al final de mi segundo bloque nos va a degradar el amortiguamiento del PID, así que téngalo presente.

*(pág. 14)*

> En espacio de estados, el seno de theta sobre r impide escribir A por x más B por u. El subsistema mecánico, por sí solo, ya es no lineal.

`pág. 15 → pasar en silencio`

---

### A5 · pág. 16–19 · `0:55` — Por qué Park

> Ahora el subsistema electromagnético. La idea de la transformación de Park es proyectar las tres variables de fase sobre un marco que gira con el rotor.
>
> ¿Para qué? Porque en régimen las variables trifásicas son senoidales, y en el marco rotórico se vuelven **cuasi-continuas**. En vez de controlar tres senoidales acopladas entre sí, controlamos dos magnitudes casi constantes: una asociada al flujo y otra al torque.

*(pág. 17–18, pasar mientras hablás)*

> El inversor sintetiza tres tensiones desfasadas 120 grados. Sus variables manipulables son el módulo de tensión de línea y la frecuencia eléctrica, cuyo **signo define el sentido de giro** del campo rotante.

*(pág. 19)*

> Aplicando Park a las ecuaciones del estator quedan tres ecuaciones con la misma estructura: caída resistiva, caída inductiva, y un término de acoplamiento proporcional a la velocidad.

---

### A6 · pág. 20–21 · `1:10` — ⚓⚓ Los acoplamientos y el torque

> **Esta lámina define todo el diseño de control que viene después**, así que me detengo un momento.
>
> En el eje q: disipación resistiva, término inductivo, la **fuerza contraelectromotriz de los imanes**, y el acoplamiento cruzado que trae el eje d. En el eje d, lo mismo más el acoplamiento cruzado que trae el eje q.
>
> Dos observaciones. Una: la resistencia depende de la temperatura, así que el subsistema eléctrico queda acoplado con el térmico. Dos: todos los términos proporcionales a omega acoplan la dinámica eléctrica con la mecánica.
>
> **Estos son exactamente los términos que el modulador de torque va a cancelar por realimentación.**

*(pág. 21)*

> El torque tiene dos contribuciones: el de imanes, proporcional a la corriente de cuadratura, y el de **reluctancia**, que existe porque L de d es distinta de L de q y es proporcional al producto de las dos corrientes.
>
> La estrategia base del trabajo es forzar la corriente de eje directo a cero: ahí el torque queda lineal y proporcional a una sola corriente. Pero dejarla libre es lo que después habilita reforzar o debilitar el campo.

`pág. 22 → pasar en silencio`

---

### A7 · pág. 23–25 · `0:40` — El subsistema térmico

> El modelo térmico es concentrado y de primer orden: consideramos solo pérdidas Joule en el bobinado, despreciamos las magnéticas del núcleo, y disipamos al ambiente por conducción y convección natural, sin ventilación forzada.
>
> Lo interesante es que hay un **lazo cerrado propio**: la temperatura sube la resistencia, y la resistencia mayor genera más pérdidas. Es una de las no linealidades del modelo.
>
> Pero es un lazo **muy lento**: la constante térmica es de **120 segundos**, tres órdenes de magnitud más lenta que la dinámica eléctrica. Eso es lo que después justifica desacoplarlo.

`pág. 26 → pasar en silencio`

---

### A8 · pág. 27–31 · `1:10` — ⚓ El modelo no lineal global y cierre

> Integrando los tres subsistemas queda el modelo no lineal global: **seis estados** —posición, velocidad, las tres corrientes en qd0 y la temperatura— y seis entradas contando las perturbaciones.

*(pág. 28)*

> No voy a leer las seis ecuaciones. Lo que quiero que se vea es **dónde están las no linealidades**: el seno de la gravedad, el producto de las dos corrientes en el torque, los términos cuadráticos de las pérdidas, y la dependencia de la resistencia con la temperatura. Cuatro fuentes.

*(pág. 29–30)*

> Y una simplificación importante. Como el estator está en estrella con **neutro flotante**, Kirchhoff en el centro de estrella impone que las tres corrientes de fase sumen cero. Y esto es **topológico**: no depende de que el sistema esté equilibrado. Por lo tanto la componente homopolar es idénticamente nula, la ecuación de secuencia cero se satisface de forma trivial, y el modelo efectivo queda en **cinco estados**.

*(pág. 31)*

> Y así se implementó en Simulink: la planta recibe las tensiones de fase y el torque de carga, y los tres subsistemas interactúan simultáneamente. Sobre este modelo se hace todo lo que sigue. Le paso a Francisco.

---
---

# BLOQUE B — Simulación a lazo abierto y control en cascada
### Páginas 66–117 · **10:25**

> ⚠️ Tu bloque crítico: 50 láminas de contenido. La disciplina de tiempo se gana o se pierde acá.

`pág. 66 → pasar en silencio`

---

### B1 · pág. 67–71 · `0:55` — El ensayo

> Retomo con la simulación dinámica a lazo abierto, que compara el modelo no lineal completo contra el LTI aumentado que acaba de presentar Francisco.
>
> El ensayo excita con un pulso de tensión en el eje q de **±19,6 volt** —que equivale a alimentar con 24 volt de línea, un 80 % de la tensión nominal— y con un doble pulso de torque de carga. Entre las dos señales quedan **diez transitorios** consecutivos, y la combinación recorre los cuatro cuadrantes.
>
> Los dos modelos corren en paralelo con las mismas entradas. En el LTI usamos resistencia constante de **1,058 ohm**, evaluada a la temperatura media de operación y no a la de referencia de 20 grados: es lo consistente con haber desacoplado lo térmico.

---

### B2 · pág. 72–73 · `0:40` — ⚓ Posición

> Posición: crece como rampa mientras dura el pulso. Es lo esperable, porque theta es la integral de omega — hay un **integrador puro** en la planta.
>
> El error entre modelos se mantiene por debajo de 0,1 radianes sobre una excursión de **242 radianes**: menos del 0,05 %. Pero no se anula, deriva lentamente.
>
> Y eso es justamente **por** el integrador: acumula las pequeñas discrepancias instantáneas de velocidad. No es un error de modelado, es la naturaleza del integrador.

---

### B3 · pág. 74–75 · `1:00` — ⚓⚓ Velocidad: el resultado más interesante

> Velocidad. Acá está el resultado más interesante del ensayo.
>
> Uno **esperaría** ver dos escalas de tiempo separadas: una eléctrica rápida y una mecánica lenta, con la constante mecánica J sobre b, que da **0,9 segundos**.
>
> No es lo que pasa. La realimentación de fuerza contraelectromotriz **fusiona los dos modos** en una única resonancia electromecánica subamortiguada, con zeta **0,527** y frecuencia natural **174**.
>
> El motor llega a su velocidad de vacío —**405 radianes por segundo**, que es donde la contraelectromotriz equilibra la tensión aplicada— con **14 % de sobrepico**, y se establece en **50 milisegundos**. Cincuenta milisegundos, no novecientos. Y esos números —14,3 % de sobrepico, 21 ms de pico, 50 ms de establecimiento— son consistentes en los diez transitorios.
>
> La coincidencia entre modelos es excelente: menos de medio radián por segundo en régimen.

---

### B4 · pág. 76–80 · `0:50` — Corriente, temperatura, tensión d, fases

> Corriente de eje q. El error clásico sería decir que se establece en tensión sobre resistencia — eso serían 18 amperes. **No.** Lo que se ve es un **pico de torque transitorio de 10,4 amperes** que decae a 0,1, porque a medida que el motor acelera la contraelectromotriz absorbe la tensión, y en vacío la corriente solo tiene que vencer la fricción.
>
> Ese pico de 10,4 es el número relevante frente a los **2 amperes** máximos del motor: este ensayo los supera ampliamente, porque es un ensayo de caracterización a lazo abierto, no de operación.

*(pág. 78–80, pasar mientras hablás)*

> La temperatura casi no se mueve: la constante térmica de 120 segundos es cincuenta veces más larga que el ensayo. Eso valida la hipótesis de resistencia constante en el LTI.
>
> Esta es la tensión que impone internamente la ley de control mínima; es no nula solo cuando hay corriente y velocidad simultáneamente. Y por Park inversa recuperamos las corrientes de fase: senoidales moduladas en amplitud, con frecuencia proporcional a la velocidad. Cuando el rotor está quieto, la frecuencia eléctrica es cero y las corrientes de fase son constantes.

---

### B5 · pág. 81 · `0:25` — ⚓ Los cuatro cuadrantes

> Esta es la curva paramétrica torque-velocidad. Cada color es uno de los diez transitorios y los rombos son los puntos de equilibrio.
>
> Se recorren **los cuatro cuadrantes**: motorización en los dos sentidos y frenado regenerativo en los dos sentidos. Esto confirma que el accionamiento es efectivamente de cuatro cuadrantes, como pedía la especificación.

---

### B6 · pág. 82–84 · `0:20` — Ángulo de carga

> El ángulo de carga es el desfasaje entre el ángulo eléctrico de la tensión y la posición eléctrica del rotor. Para calcularlo tomamos la tensión de fase a como referencia, proyectamos con **Clarke** al plano alfa-beta, y recuperamos el ángulo con **atan2**.

---

### B7 · pág. 85–89 · `0:50` — Condición inicial de i_ds y field weakening

> Efecto de la condición inicial de la corriente de eje directo. Arrancamos con cero, más 0,5 y menos 0,5 amperes. En los tres casos converge a cero con constante de tiempo L de d sobre resistencia, **6,2 milisegundos**, y se extingue en unos 30 milisegundos sin afectar nada más.
>
> Esto es lo que **justifica** que la hipótesis de corriente de eje d inicial nula, usada en la linealización, no sea restrictiva en la práctica.

*(pág. 87–88)*

> Después inyectamos una consigna de tensión en el eje d, de un décimo de la nominal, y aparecen los dos regímenes.
>
> **Reforzamiento de campo**, con corriente d positiva: suma torque de reluctancia, pero aumenta la contraelectromotriz y baja la velocidad máxima alcanzable. **Debilitamiento**, con corriente d negativa: resigna torque por amper, pero permite girar más rápido dentro del límite de tensión del inversor.

*(pág. 89)*

> Y algo importante: el modelo LTI **no captura** esto. Solo lo ve el no lineal, porque son términos de producto entre estados.

`pág. 90 → pasar en silencio`

---

### B8 · pág. 91–96 · `1:15` — ⚓⚓ El modulador de torque: la idea central

> Pasamos al diseño del controlador. La primera pieza es el **modulador de torque**, que es el controlador interno vectorial de corriente y torque.
>
> La idea central es esta. La planta tiene **realimentaciones físicas naturales**: los acoplamientos cruzados, la contraelectromotriz de los imanes y las caídas resistivas. Lo que hacemos es calcular esos términos y sumárselos por prealimentación a la tensión de comando, para que se **cancelen exactamente** contra los de la planta.

*(pág. 94 — la lámina clave)*

> Cuando uno sustituye y simplifica, queda esto: la tensión de comando de cada eje es igual a la inductancia por la derivada de la corriente. O sea, **cada eje se convierte en un integrador puro, perfectamente desacoplado del otro.**
>
> Ese es el resultado clave de todo el bloque: convertimos una planta de tres ecuaciones acopladas y no lineales en tres integradores independientes. Todo el control que viene después es fácil gracias a esto.

*(pág. 95–96, pasar rápido)*

---

### B9 · pág. 97–100 · `0:40` — Diferencias y el efecto de R_s

> ¿En qué se diferencia esto de la linealización por realimentación no lineal que ya habíamos hecho? En tres cosas. Primero, se **libera** la restricción de corriente d nula, lo que habilita el debilitamiento de campo. Segundo, antes solo cancelábamos los acoplamientos cruzados; ahora cancelamos **también** la contraelectromotriz y las caídas resistivas. Y tercero, lo más importante: ahora podemos **asignar los polos** del lazo de corriente donde queramos.

*(pág. 98–100)*

> Sobre la resistencia hicimos una comparación: estimarla en línea con el sensor de temperatura, contra usar el valor nominal constante. Con valor constante, al cabo de 4 segundos aparece un **déficit de torque del 15 %**, que integrado se traduce en un **4 % de error de posición**. Conclusión: en régimen prolongado hay que estimarla.

---

### B10 · pág. 101–104 · `0:50` — ⚓ Lazos de corriente y la acción integral

> Con la planta reducida a integradores puros, los lazos de corriente son **simplemente proporcionales**. La transferencia queda de primer orden con constante de tiempo L sobre K p, y la consigna pide los polos en **menos 5000**, así que las ganancias son directamente la inductancia por 5000.

*(pág. 102 — ⚠️ la lámina con la pregunta abierta. Esto lo decís vos, no está escrito)*

> Acá hay algo que la consigna pide investigar explícitamente: **¿conviene agregar acción integral?**
>
> La respuesta es **no**, y la razón es el principio del modelo interno: la planta de cada eje **ya es de tipo 1** —tiene su propio integrador—, así que el proporcional solo ya alcanza error nulo ante escalón.
>
> Agregar un integrador la volvería tipo 2 sin mejorar el seguimiento, y traería solo desventajas: sobrepico, más sensibilidad al ruido de medición, y riesgo de **windup** cuando el inversor satura. Por eso el lazo queda puramente proporcional.

---

### B11 · pág. 105–107 · `0:35` — Torque precalculado

> Falta traducir una consigna de torque a una consigna de corriente. Como la única forma de variar el torque es por la corriente, invertimos la ecuación de torque.
>
> Y aprovechamos para aplicar **torque precalculado**: el modulador le suma a la consigna la compensación de fricción viscosa y la de gravedad. Así el lazo externo no tiene que pelear contra la gravedad — ya viene compensada por prealimentación. Eso simplifica muchísimo el diseño del PID.
>
> Y noten que la corriente de eje d queda como **parámetro libre** en el denominador: ahí es donde se enchufaría el debilitamiento de campo. Para este proyecto la dejamos en cero.

---

### B12 · pág. 108–114 · `1:10` — ⚓⚓ El PID externo

> El lazo externo es un PID que manipula directamente el torque, sintonizado por el método de **sintonía serie** con zeta **0,75** y frecuencia natural **800**.
>
> La clave del diseño es **qué planta ve el PID**. Como los lazos de corriente son más de seis veces más rápidos, a la escala del lazo de movimiento el torque sigue instantáneamente a su consigna. Y como el torque precalculado ya compensó fricción y gravedad, lo que queda es simplemente **un doble integrador** con la perturbación de carga. Una planta trivial.

*(pág. 110–112, pasar mientras hablás)*

> El método es: se saca el polinomio característico del lazo cerrado en función de las ganancias —el **estructural**—, se construye el **requerido** a partir de los polos deseados, y se igualan coeficiente a coeficiente.
>
> La elección de polos es un polo real y un par complejo que **comparten la misma frecuencia natural**: los tres quedan sobre la circunferencia de radio 800.
>
> Y lo elegante del método es que **cada ganancia aparece en un único coeficiente**, así que el despeje es directo, sin resolver ningún sistema.

*(pág. 113–114)*

> Quedan estos tres números, y así se implementó en Simulink. Un detalle que vale la pena: el factor s en el numerador de la transferencia de perturbación muestra que la acción integral **rechaza exactamente** perturbaciones de carga constantes.

---

### B13 · pág. 115–116 · `0:30` — ⚓ Separación de escalas

> Este mapa de polos resume la arquitectura. En verde los lazos de corriente en **menos 5000**; en rojo el lazo cerrado del PID alrededor de **800**; en azul la planta original con su modo en **174**.

*(pág. 116)*

> Y acá está escrito el argumento que valida todo el diseño en cascada: **174, mucho menor que 800, mucho menor que 5000**. Un orden de magnitud entre niveles. Eso es lo que permite diseñar cada lazo suponiendo que el interior ya está resuelto y que el exterior no lo perturba.

---

### B14 · pág. 117 · `0:25` — ⚓ Cierre: qué pasa con la carga

> Y cierro con la contracara. Las ganancias se sintonizaron con la inercia nominal, sin carga útil, y **no se reajustan en línea**. Al barrer la masa de carga de 0 a 1,5 kg, los polos de lazo cerrado migran hacia el eje imaginario: el amortiguamiento **cae de 0,75 a 0,30**.
>
> El sistema sigue estable en todo el rango, pero se vuelve notablemente más oscilatorio en el peor caso. Si hubiera que conservar el amortiguamiento especificado, haría falta agendado de ganancias o un esquema de control robusto. Sigue Francisco.

---
---

# BLOQUE C — Verificación de desempeño y no idealidades
### Páginas 153–165 · **3:00**

---

### C1 · pág. 153–155 · `0:50` — ⚓⚓ El ensayo térmico

> Retomo con la verificación del desempeño térmico. El ensayo consiste en repetir cíclicamente la consigna de trabajo: un trapecio de **15 segundos** —sube 5, mantiene 5, baja 5— con 25 grados de ambiente y el torque de perturbación máximo de 5 newton-metro.

*(pág. 154)*

> El resultado es que la temperatura crece de forma **acumulativa ciclo a ciclo y nunca alcanza el equilibrio térmico**. Cruza los 115 grados de la clase de aislación a los **225,9 segundos**, o sea a los 15 ciclos.

*(pág. 155)*

> Y la conclusión es negativa, y conviene decirla como tal: **el accionamiento no satisface el requerimiento térmico de operación continua repetitiva** en esas condiciones. O se agrega disipación activa, o hay que reformular el ciclo de trabajo.

---

### C2 · pág. 156–159 · `0:40` — Sensores como filtros

> Sensores no ideales. Hasta acá los habíamos supuesto ideales, con ganancia unitaria. Ahora los reemplazamos por **filtros pasa bajos en espacio de estados**, según lo que especifica la consigna: los de corriente y el de posición, de segundo orden críticamente amortiguados; el térmico, de primer orden con 20 segundos de constante.
>
> La derivación a espacio de estados es directa: se antitransforma la transferencia y se eligen la salida y su derivada como variables de estado.

---

### C3 · pág. 160–162 · `0:50` — ⚓⚓ El hallazgo de esta sección

> Y acá aparece el hallazgo más interesante de esta parte. **Con los polos que especifica la consigna, el sistema directamente no logra seguir la consigna de posición.**
>
> La razón es de separación de dinámicas. Nuestro lazo de corriente tiene el polo en menos 5000, y los sensores de corriente en menos 6000: son apenas un **20 % más rápidos**. Eso no es separación de escalas — el retardo del sensor entra de lleno en el lazo. Y con el observador pasa lo mismo.
>
> La solución fue escalar los polos de todos los sensores en conjunto por un factor N. Con **N igual a 3** el desempeño vuelve a ser aceptable.
>
> Y la conclusión de ingeniería es concreta: los sensores que efectivamente habría que comprar tienen que ser **por lo menos tres veces más rápidos** que los que sugiere la especificación de la guía.

---

### C4 · pág. 163–165 · `0:40` — El inversor y cierre del bloque

> Con el inversor pasa exactamente lo mismo. Un inversor real tiene dos restricciones que el bloque ideal no captura: **ancho de banda finito** por los semiconductores, y **saturación de tensión** por la fuente de continua. Lo modelamos igual, como un pasa bajos de segundo orden.

*(pág. 164)*

> Con los 6000 radianes por segundo iniciales no alcanza — otra vez, es solo un 20 % más rápido que el lazo de corriente.

*(pág. 165)*

> Duplicando a **12000** se consigue un desempeño satisfactorio cumpliendo el resto de las condiciones de operación.
>
> Y el cierre del bloque es que los tres elementos —sensores, observador e inversor— comparten el mismo diagnóstico: **la especificación de la guía no da suficiente separación de escalas frente a un lazo de corriente ubicado en 5000**. Hay que exigirles entre dos y tres veces más rapidez. Francisco sigue con la discretización.

---

## Resumen de tiempos

| Bloque | Páginas | Beats | Tiempo |
|---|---|---|---|
| A — Modelado físico | 1–31 | 8 | 7:15 |
| B — Simulación y control | 66–117 | 14 | 10:25 |
| C — Verificación | 153–165 | 4 | 3:00 |
| **Total** | **96** | **26** | **20:40** |

---

## Los 4 anclajes que no podés perder

Si te quedás sin tiempo, sacrificá todo menos esto:

1. **pág. 20–21** — Los acoplamientos físicos son los que después se cancelan. Es la bisagra entre modelado y control.
2. **pág. 74–75** — Los modos eléctrico y mecánico se **fusionan** por la contraelectromotriz. Es el resultado no obvio del trabajo.
3. **pág. 94** — El desacople convierte cada eje en un **integrador puro**. Es lo que hace posible todo el control.
4. **pág. 160–162** — Los sensores especificados **no alcanzan**; hacen falta 3× más rápidos. Es la conclusión de ingeniería más concreta.

---

## Antes de la mesa

- [ ] **Recompilar `main_presentacion.pdf`** — el de la carpeta es del 29/07 y tiene 184 páginas. El actual tiene **182**.
- [ ] Confirmar con Francisco que ambos presentan desde el **mismo PDF**.
- [ ] Avisarle que su segundo bloque ahora termina en **152** (la tabla de contenidos de verificación) y que vos arrancás en el ensayo térmico.
- [ ] Ensayar la pág. **102** (la lámina "¿Acción integral?" está vacía a propósito: el argumento lo decís vos).
