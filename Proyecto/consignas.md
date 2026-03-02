# UNCUYO -- Ingeniería Mecatrónica
## Espacio Curricular 311: Automática y Máquinas Eléctricas
**Prof. Titular:** Ing. Gabriel L. Julián
**Proyecto Global Integrador (2025): Guía de Trabajo y Especificaciones**
**Rev.0:** 18/02/2026
**Lugar:** Mendoza, Argentina

---

# Proyecto Global Integrador: Guía de Trabajo
## Control de Accionamiento de CA con Motor Sincrónico de Imanes Permanentes (PMSM)

---

## 1. Objetivo y Alcances

Proyecto didáctico, con el objetivo de integrar los conocimientos y competencias fundamentales del Espacio Curricular "Automática y Máquinas Eléctricas" en una aplicación mecatrónica concreta, a partir de especificaciones de requisitos simplificadas.

Se trabajará sobre el **modelado, simulación, diseño y análisis de desempeño** de un **Sistema de Control Automático de Posición y Movimientos** para un accionamiento eléctrico de **4 cuadrantes**, compuesto por:

- 1 máquina eléctrica de CA trifásica sincrónica con excitación por imanes permanentes (**PMSM**)
- 1 inversor electrónico trifásico desde fuente de tensión continua (**CC**)
- 1 reductor de velocidad con engranajes planetarios conectado a la carga mecánica
- Sensores:
  - 1 sensor de posición en eje del motor
  - 3 sensores de corriente instantánea de fases (salida del inversor hacia el estator)
  - 1 sensor de temperatura en el bobinado del estator

---

## 2. Lineamientos generales

- Trabajo colaborativo en equipo de **dos (2) alumnos**; ambos deben dominar el proyecto completo.
- Modelado, simulación, análisis y diseño en **Matlab/Simulink**. Respetar nomenclatura indicada.
- El contenido debe ser **producción propia y original** (no copia ni adaptación total/parcial) y cumplir la guía.
- Horario de consulta semanal para dudas y avances.
- Presentación de **Informe Técnico escrito**, completo y breve, incluyendo:
  - Resumen
  - Introducción
  - Desarrollo: modelado y esquemas conceptuales; análisis; diseño e implementación; simulación; resultados
  - Conclusiones
  - Referencias
  - Anexos
  (Ver documentos complementarios: *Guía para preparar Informe Técnico y Plantilla*.)
- Exposición oral presencial y demostración:
  - Coloquio basado en el Proyecto (al menos 1 semana antes del Examen)
  - Presentación en Mesa de Examen Final (Ordinaria o Especial)
- Plazo de presentación final: **Hasta Junio de 2026 (fecha a coordinar)** si desean cursar 317 AyCD en 2026.
  Nota: para fechas posteriores, usar la última versión vigente de esta guía.

---

## 3. Especificaciones del Sistema Físico (Accionamiento) y Modelado Dinámico

### Convención de signos (torque y velocidad)
- **Iguales** (`++` o `--`) => motorización (**cuadrantes I y III**)
- **Opuestos** (`+-` o `-+`) => frenado regenerativo (**cuadrantes II y IV**)

---

### 3.1 Carga mecánica

Aplicación mecatrónica simplificada (adaptada de [1]): control de movimiento de un eje (control descentralizado) para la articulación de un brazo manipulador robótico elemental (péndulo rígido actuado), con un solo grado de libertad (1 g.d.l.), eje de rotación horizontal fijo a la base en un sistema de referencia inercial; sometido a la acción externa de la aceleración de gravedad \(g\) (perturbación externa vertical constante).

Modelo simplificado equivalente (No lineal con parámetros variables), referido al eje de salida del tren de transmisión:

- Coordenada articular:
  \[
  q(t)\equiv \theta_l(t)
  \]
  referida al eje vertical hacia abajo (equilibrio estable), positiva antihoraria.
- Torque impulsor:
  \[
  \tau(t)\equiv T_q(t)
  \]
- Perturbación externa: gravedad constante hacia abajo
  \[
  g = 9.80665~\text{m/s}^2
  \]

**Ecuación dinámica (Ec. 1.1):**
\[
J_l \frac{d\omega_l(t)}{dt} = T_q(t) - b_l \omega_l(t) - T_l(t)
\]
donde:
\[
T_l(t)= g\,k_l \sin(\theta_l(t)) + T_{ld}(t)
\]

**Relación posición-velocidad (Ec. 1.2):**
\[
\frac{d\theta_l(t)}{dt}\equiv \omega_l(t)\;\;\Longleftrightarrow\;\;
\theta_l(t) = \int_{t_0}^{t}\omega_l(\xi)\,d\xi + \theta_l(t_0)
\]

#### Parámetros equivalentes variables (valor nominal +/- variación máx.)

- Coeficiente de fricción viscosa en la articulación:
  \[
  b_l \approx (0.1\pm 0.03)\; \text{N·m/(rad/s)} \quad \text{(incertidumbre)}
  \]
- Masa del brazo manipulador: \(m=1.0\) kg
- Longitud e inercia equivalente (centro de masa):
  \(l_{cm}=0.25\) m; \(J_{cm}=0.0208\) kg·m²
- Longitud total (extremo): \(l_l=0.50\) m
- Masa de carga útil en extremo (variable): \(m_l\in[0,1.5]\) kg (no conocida a priori)
- Momento de inercia total (a eje de rotación):
  \[
  J_l = (m\,l_{cm}^2 + J_{cm}) + m_l\,l_l^2 = 0.0833 + [0\ldots 0.375]\;\text{kg·m}^2
  \]
- Coeficiente \(k_l\) en Torque de carga \(T_l(t)\), Ec. 1.1:
  \[
  k_l = m\,l_{cm} + m_l\,l_l = 0.25 + [0\ldots 0.75]\;\text{kg·m}
  \]
- Torque recuperador gravitacional:
  \[
  T_g(t) = g\,k_l\,\sin(\theta_l(t))\quad [\text{N·m}]
  \]

#### Especificaciones de operación (perturbación adicional)

- Torque de perturbación por contacto:
  \[
  T_{ld}(t)\approx (0\pm 5.0)\;\text{N·m} \quad \text{(asumir escalón)}
  \]

---

### 3.2 Tren de transmisión

Caja reductora reversible con sistema de engranajes planetarios, asumiendo acoplamiento rígido (sin elasticidad torsional ni juego/holgura/backlash). Momento de inercia equivalente y pérdidas por fricción interna reflejados al eje de entrada y considerados junto al motor (parte de \(J_m\) y \(b_m\), Ec. 3.1).

**Modelo equivalente (rígido) (Ec. 2.1 / 2.2):**
\[
\omega_l(t)=\frac{1}{r}\,\omega_m(t)
\]
\[
T_q(t)= r\,T_d(t)
\]

- Relación de reducción total: \(r=120.0:1\)

**Límites (no sobrepasar):**
- Velocidad nominal salida:
  \(n_{l,nom}=60\) rpm (\(\omega_{l,nom}=6.28\) rad/s)
- Torque nominal salida: \(T_{q,nom}=17.0\) N·m (régimen continuo o rms)
- Torque pico salida: \(T_{q,max}=45.0\) N·m (corta duración, aceleración)

---

### 3.3 Máquina eléctrica PMSM

Máquina eléctrica de CA trifásica sincrónica con excitación por imanes permanentes (PMSM), con estator simétrico y equilibrado, bobinas abc conectadas en estrella (Y), centro de estrella \(O\) (neutro) flotante (no accesible) y bornes de línea ABC accesibles desde el inversor.

#### Subsistema mecánico (rotor, referido a estator estacionario)

(Ec. 3.1)
\[
J_m \frac{d\omega_m(t)}{dt} = T_m(t) - b_m\,\omega_m(t) - T_d(t)
\]

(Ec. 3.2)
\[
\frac{d\theta_m(t)}{dt}\equiv \omega_m(t)
\;\;\Longleftrightarrow\;\;
\theta_m(t)=\int_{t_0}^{t}\omega_m(\xi)\,d\xi + \theta_m(t_0)
\]

#### Subsistema electromagnético (qd0 rotor, Park)

(Ec. 3.3)
\[
\frac{d\theta_r(t)}{dt}\equiv \omega_r(t)
\;\;\Longleftrightarrow\;\;
\theta_r(t)=\int_{t_0}^{t}\omega_r(\xi)\,d\xi + \theta_r(t_0)
\]

(Ec. 3.4)
\[
\theta_r(t)=P_p\,\theta_m(t)\quad\Rightarrow\quad \omega_r(t)=P_p\,\omega_m(t)
\]

**Torque electromagnético (Ec. 3.5):**
\[
T_m(t)=\frac{3}{2}P_p\lambda'_{mr}\, i^r_{qs}(t)
+\frac{3}{2}P_p(L_d-L_q)\, i^r_{ds}(t)\, i^r_{qs}(t)
\]

#### Balance de tensiones (qd0 rotor)

(Ec. 3.6)
\[
v^r_{qs}(t)=R_s(T^\circ_s(t))\,i^r_{qs}(t)+L_q\frac{di^r_{qs}(t)}{dt}
+\left[\lambda'_{mr}+L_d i^r_{ds}(t)\right]\omega_r(t)
\]

(Ec. 3.7)
\[
v^r_{ds}(t)=R_s(T^\circ_s(t))\,i^r_{ds}(t)+L_d\frac{di^r_{ds}(t)}{dt}
- L_q i^r_{qs}(t)\,\omega_r(t)
\]

(Ec. 3.8)
\[
v_{0s}(t)=R_s(T^\circ_s(t))\,i_{0s}(t)+L_{ls}\frac{di_{0s}(t)}{dt}
\]

Resistencia con temperatura (Ec. 3.9):
\[
R_s(T^\circ_s(t))=R_{sREF}\left(1+\alpha_{Cu}\left(T^\circ_s(t)-T^\circ_{sREF}\right)\right)
\]

**Nota sobre Ec. 3.8:** Para un sistema trifásico equilibrado y circuito simétrico, en el caso general (retorno de neutro) la Ec. 3.8 no influye: \(i_{0s}(t)\equiv 0\); \(v_{0s}(t)\equiv 0\). Además, al tener el neutro flotante (no accesible), la Ec. 3.8 tampoco influye (demostrar mediante Transformación de Park y Leyes de Kirchoff).

---

### 3.4 Subsistema térmico (primer orden)

Pérdidas Joule reales abc (Ec. 3.10):
\[
P_{s,perd}(t)=R_s(T^\circ_s(t))\left(i_{as}^2(t)+i_{bs}^2(t)+i_{cs}^2(t)\right)
\]

O equivalente virtual qd0 (para análisis):
\[
P_{s,perd}(t) = R_s(T^\circ_s(t)) \frac{3}{2}\left(i^{r\,2}_{qs}(t)+i^{r\,2}_{ds}(t)+2\,i_{0s}^2(t)\right)
\]

Balance térmico (Ec. 3.11):
\[
P_{s,perd}(t)=C_{ts}\frac{dT^\circ_s(t)}{dt} + \frac{1}{R_{ts-amb}}\left(T^\circ_s(t)-T^\circ_{amb}(t)\right)
\]

#### Parámetros nominales (tolerancia +/-1%, salvo aclaración)

- \(J_m \approx 14.0\times 10^{-6}\) kg·m²
- \(b_m \approx 15.0\times 10^{-6}\) N·m/(rad/s)
- Pares de polos: \(P_p=3\) (i.e. 6 polos)
- Flujo concatenado: \(\lambda'_{mr}\approx 0.016\) Wb-turn (V/(rad/s))
- \(L_q\approx 5.8\) mH
- \(L_d\approx 6.6\) mH
- \(L_{ls}\approx 0.8\) mH
- \(R_s\approx 1.02~\Omega = R_{sREF}\) (@ \(T^\circ_{sREF}=20^\circ C\))
- \(\alpha_{Cu}=3.9\times 10^{-3}~1/^\circ C\)
- \(C_{ts}\approx 0.818\) W/(°C/s) (almacenamiento interno)
- \(R_{ts-amb}\approx 146.7\) °C/W (disipación al ambiente)
- \(\tau_{ts-amb}=R_{ts-amb}C_{ts}\approx 120\) s

#### Límites operativos (estator abc)

- Velocidad nominal rotor: \(n_{m,nom}=6600\) rpm (\(\omega_{m,nom}=691.15\) rad/s)
- Tensión nominal línea: \(V_{sl,nom}=30\) Vca rms
  (fase: \(V_{sf,nom}=V_{sl,nom}/\sqrt{3}\))
- Corriente nominal: \(I_{s,nom}=0.4\) Aca rms (continuo)
- Corriente máxima: \(I_{s,max}=2.0\) Aca rms (corta duración)
- Temperatura máxima bobinado: \(T^\circ_{s,max}=115^\circ C\)
- Temperatura ambiente: \(-15^\circ C \le T^\circ_{amb}(t)\le 40^\circ C\) (perturbación)

**Nota:** \(T^\circ_{amb}(t)\) es perturbación externa (puede asumirse aprox. "constante" para intervalos cortos, pero podría variar). Elegir \(T^\circ_{amb}(t)\) adecuadamente frente a \(T^\circ_s(t_0)\), según escenario.

---

### 3.5 Inversor electrónico trifásico (modulador de tensión)

No se analiza detalle PWM. Se considera inversor + fuente CC como **modulador idealizado**.

Inicialmente:
- Sin saturación de tensión
- Respuesta ideal: \(G(s)\equiv 1\) (filtro "pasatodo", ganancia unitaria y ancho de banda "infinito")

Finalmente: limitar ancho de banda y agregar saturación por tensión finita (ver 5.2.5.e).

#### Modelo promediado de tensiones (componente fundamental)

(Ec. 4.1-4.3)
\[
v_{as}(t)\approx \sqrt{2}\frac{V_{sl}(t)}{\sqrt{3}}\cos(\theta_{ev}(t))
\]
\[
v_{bs}(t)\approx \sqrt{2}\frac{V_{sl}(t)}{\sqrt{3}}\cos\left(\theta_{ev}(t)-\frac{2\pi}{3}\right)
\]
\[
v_{cs}(t)\approx \sqrt{2}\frac{V_{sl}(t)}{\sqrt{3}}\cos\left(\theta_{ev}(t)+\frac{2\pi}{3}\right)
\]

(Ec. 4.4)
\[
\omega_e(t)\equiv 2\pi f_e(t)\equiv \frac{d\theta_{ev}(t)}{dt}
\;\;\Longleftrightarrow\;\;
\theta_{ev}(t)=\int_{t_0}^{t}\omega_e(\xi)\,d\xi + \theta_{ev}(t_0)
\]

- Variables manipulables: \(V_{sl}(t)\) y \(f_e(t)\)
- Límites:
  - \(V_{sl}\in[0,48]\) Vca rms
  - \(f_e\in[-330,0,+330]\) Hz
    El signo de \(f_e\) define secuencia abc o acb y sentido de giro.

---

### 3.6 Ángulo de torque del rotor \(\delta(t)\)

(Ec. 4.5)
\[
\delta(t)\equiv \theta_{ev}(t)-\theta_r(t)
=\int_{t_0}^{t}\left(\omega_e(\xi)-\omega_r(\xi)\right)\,d\xi
+\left[\theta_{ev}(t_0)-\theta_r(t_0)\right]
\]

Interpretación: rotor "atrasado" en modos de motorización y "adelantado" en modos de frenado regenerativo.

---

### 3.7 Sensores de retroalimentación

- 1 sensor de posición angular (encoder incremental) en eje motor (homing y decodificación idealizados)
  Variable medida: \(\theta_m(t)\), posición angular absoluta "rectificada" (al girar más de una revolución), tal que:
  \[
  \theta_l(t)=\frac{1}{r}\theta_m(t)
  \]
  y ajuste \(\theta_m\equiv 0\) (Home) para \(\theta_l\equiv 0\).
- 3 sensores de corriente instantánea de fase: \(i_{as}(t), i_{bs}(t), i_{cs}(t)\)
- 1 sensor de temperatura (ej. RTD) en bobinado: \(T^\circ_s(t)\)
  (para estimación de \(R_s(T^\circ_s(t))\), monitoreo y protección)

Nota: asumir inicialmente todos los sensores con respuesta ideal \(G(s)\equiv 1\). Luego limitar BW (ver 5.2.5.d).

---

### 3.8 Variables principales en el Modelo Dinámico completo

**a) Excitaciones (entradas) externas:**
- Variable manipulada (vectorial): \(\mathbf{v}_{abcs}(t)\) => (Park directa) => \(\mathbf{v}^r_{qd0s}(t)\) virtual equivalente interna.
- Variables de perturbación: \(T_l(t)\) (torque de carga) y \(T^\circ_{amb}(t)\).

**b) Estado interno:**
- Posición \(\theta_m(t)\) y velocidad \(\omega_m(t)\) en eje del motor.
- Corrientes virtuales equivalentes \(\mathbf{i}^r_{qd0s}(t)\); temperatura \(T^\circ_s(t)\).

**c) Respuestas (salidas) externas:**
- Variable controlada (no medida directamente): \(q(t)\equiv\theta_l(t) = \frac{1}{r}\theta_m(t)\)
- Variables medidas (retroalimentación): \(\theta_m(t)\), \(\mathbf{i}_{abcs}(t)\), \(T^\circ_s(t)\).

---

## 4. Requisitos generales de Análisis y Diseño

a) Diseño en tiempo continuo \(t[s]\in\mathbb{R}\). Inicialmente asumir disponibilidad de todas las variables de estado requeridas por sensores; luego reemplazar donde corresponda por estimación mediante Observador reducido (5.2.3).

b) Valores nominales de carga:
- \(J_{l,nom}\approx 0.0833\) kg·m²
- \(b_{l,nom}\approx 0.1\) N·m/(rad/s)
Inicialmente sin carga útil: \(m_l\equiv 0\).
Luego evaluar migración variando \(m_l\in[0,1.5]\) kg y \(b_l\approx(0.1\pm 0.03)\), manteniendo invariante el control diseñado con valores nominales (robustez).

c) Estado inicial consistente (ejemplo):
\(\theta_l(t_0)=0\), \(\omega_m(t_0)=0\) rad/s, \(\mathbf{i}^r_{qd0s}(t_0)=\mathbf{0}\) A, \(T^\circ_s(t_0)=T^\circ_{amb}(t_0)=T^\circ_{amb,max}=40^\circ C\).

d) Estrategia de control vectorial con campo orientado ("desacoplamiento" flujo/torque) forzando:
\[
i^r_{ds}(t)\equiv 0
\]
Imponer restricción/ley de control NL sobre \(\mathbf{v}^r_{qd0s}(t)\) (o equivalente en \(\mathbf{v}_{abcs}(t)\) vía Park).

**Nota para diagramas de bloques:** Separar claramente sistema físico del sistema de control. Indicar correctamente todas las Transformaciones de Park necesarias (directas/inversas). Indicar dónde se implementan los "desacoplamientos". No repetir diagramas innecesariamente; usar colores y leyendas inteligentemente, compactar y sintetizar.

---

# 5. Tareas a desarrollar

## 5.1 Modelado, Análisis y Simulación dinámica del SISTEMA FÍSICO a "Lazo Abierto" (Sin Controlador externo de Movimiento)

**1)** Modelo matemático equivalente (1 g.d.l.) del subsistema mecánico completo: subsistema mecánico del motor (Ec. 3.1/3.2) + transmisión rígida (Ec. 2.1/2.2) + carga (Ec. 1.1/1.2), referido al eje del motor. ¿Por qué se puede realizar esta simplificación o compactación?

**2)** Modelo dinámico del sistema físico completo, incorporando los subsistemas electromagnético y térmico acoplados de la máquina eléctrica al subsistema mecánico completo obtenido en ítem 1:

> **a)** Modelo global no lineal (NL), para \(i^r_{ds}(t)\neq 0\) (caso gral.):
> - I) Ecuaciones vectoriales NL de estado y de salida (con estado inicial genérico, en coordenadas virtuales qd0).
> - II) Diagrama de bloques de estado (forma desagregada o escalar); incorporar las Transformaciones de Park virtuales según corresponda, para acceso físico en bornes a las tensiones y corrientes de fase reales de estator (coordenadas abcs).
>
> **b)** Linealización Jacobiana: Modelo global linealizado con parámetros variables (LPV), para \(i^r_{ds}(t)\neq 0\) (caso gral.), a partir de modelo NL mediante aproximación con serie de Taylor truncada de 1° orden en punto genérico de operación. Ecuaciones:
> - I) Espacio de operación global NL (cuasi-estacionario).
> - II) Modelo dinámico LPV (pequeñas variaciones locales), función de parámetros variables según el punto de operación. Indicar estado inicial genérico.
>
> **c)** Linealización por Retroalimentación NL: Modelo simplificado lineal invariante (LTI) equivalente (sin tener en cuenta acoplamiento NL con el subsistema térmico, pero sí considerar su dinámica lineal), imponiendo directamente el requisito \(i^r_{ds}(t)\equiv 0\) (estrategia de Control Vectorial con campo orientado), a partir del modelo NL original:
> - **I)** Ecuaciones vectoriales/matriciales LTI de estado y de salida (con estado inicial genérico) => matrices del modelo LTI equivalente.
> - **II)** Diagrama de bloques de estado (forma desagregada o escalar).
> - **III)** Determinación de la **Restricción o Ley de Control mínima** que es necesario aplicar sobre la variable manipulada virtual \(\mathbf{v}^r_{qd0s}(t)\), a través de las Transformaciones de Park necesarias con \(\mathbf{v}_{abcs}(t)\), para cumplir \(i^r_{ds}(t)\equiv 0\) ("desacoplamiento" de canales de flujo magnético y torque); ¿qué hipótesis se asume para el estado inicial de \(i^r_{ds}(t)\)?
> - **IV)** Implementación, en el modelo global NL completo (ítem 2.a), de esta Ley de control ("desacoplamiento" o compensación, y linealización por retroalimentación directa NL de estado parcial) mediante un controlador parcial; incorporando el inversor (modulador de tensión trifásico equivalente), las Transformaciones de Park reales según corresponda, sensores de retroalimentación ideal de variables de estado, etc. (separar claramente el controlador de la planta, e indicar la/s entrada/s de manipulación y variables medidas para retroalimentación resultantes).
> - **V)** Modelo de la **dinámica residual equivalente** para \(i^r_{ds}(t)\) (eje d) al aplicar esta ley de control mínima, para el caso general en que no se cumple la hipótesis asumida para el estado inicial de \(i^r_{ds}(t)\) => incorporar al modelo LTI equivalente (ítems 2.c.I y 2.c.II), despreciando el acoplamiento residual NL con el eje q. ¿Cuál es este acoplamiento residual NL y por qué se puede despreciar sin error significativo en régimen forzado?
> - **VI)** ¿Se puede implementar alguna **Restricción o Ley de Control complementaria mínima** en el eje q para eliminar completamente este acoplamiento residual NL aún en régimen natural y obtener un modelo equivalente completamente lineal, independiente del estado inicial de \(i^r_{ds}(t)\)?; en tal caso, agregar al controlador parcial en ítem 2.c.IV. Mostrar claramente ambos modelos resultantes: **NL desacoplado con Ley de control NL** y **LTI equivalente aumentado**.
>
> **d)** Comparación del modelo dinámico LTI equivalente aumentado vs. el modelo dinámico global LPV forzando \(I^r_{ds,o}\equiv 0\). Evaluación del modelo LPV para otros puntos de operación con \(I^r_{ds,o}< 0\) (debilitamiento de campo), \(I^r_{ds,o}> 0\) (reforzamiento de campo): migración de propiedades ante cambios de punto de operación variando \(I^r_{ds,o}\).
>
> **e)** Funciones de Transferencia para el modelo LTI equivalente aumentado (ítem 2.c.VI), desde ambas entradas \(v^r_{qs}(t)\) y \(T_l(t)\) hacia la salida \(\theta_m(t)\). Indicar el estado inicial considerado, y si hay estados internos que no aportan ni se ven reflejados en las funciones de transferencia.

**3)** Análisis de Estabilidad a lazo abierto para el modelo LTI equivalente aumentado (ítem 2.c.VI):
> a) Determinar autovalores = polos y ceros (valores numéricos y mapa en plano s); identificar correspondencia con modos de oscilación, subsistemas, entradas, etc.; calcular sus frecuencias naturales y amortiguamientos, o constantes de tiempo correspondientes.
> b) Evaluar estabilidad parcial y completa, y dinámica o efecto de los ceros.
> Nota: considerar migración de propiedades ante variación de parámetros de carga, ítem 4.b.

**4)** Análisis de **Observabilidad** completa de estado para el modelo LTI equivalente aumentado (ítem 2.c.VI) desde salida medida \(\theta_m(t)\). ¿Existe algún estado no observable desde \(\theta_m(t)\)?
> Alternativa: medir velocidad \(\omega_m(t)\) con tacogenerador, en vez de medir posición con encoder y reevaluar Observabilidad. ¿Existe algún estado no observable desde \(\omega_m(t)\)?

**5)** Análisis de **Controlabilidad** completa de estado para el modelo LTI equivalente aumentado (ítem 2.c.VI) desde entrada manipulada \(v^r_{qs}(t)\), sin considerar la perturbación de la carga mecánica. ¿Existe algún estado no controlable desde \(v^r_{qs}(t)\)?; en tal caso, ¿podría controlarse agregando alguna otra entrada de control razonable?

**6)** Simulación dinámica en DT, comparando el modelo NL completo desacoplado con Ley de control NL vs LTI equivalente aumentado (ítem 2.c.VI), comparando \(i^r_{ds}(0)=\pm 0.5\) A vs \(i^r_{ds}(0)= 0\) A:

> **a)** Respuesta del estado interno \(\{\theta_m(t); \omega_m(t); \mathbf{i}^r_{qd0s}(t); T^\circ_s(t)\}\) (y \(v^r_{ds}(t)\) forzada) a **pulso de consigna de tensión** de estator en eje q:
> \[
> v^{r*}_{qs}(t) = 0\text{V} \to (+19.596\text{Vcc en }t_{step1}=0.1\text{s}) \to (0\text{V en }t_{step4}=0.7\text{s}) \to (-19.596\text{Vcc en }t_{step6}=1.1\text{s}) \to (0\text{V en }t_{step9}=1.7\text{s})
> \]
> superpuesto con **doble pulso de torque de carga**:
> \[
> T_l(t) = 0 \to (+6.28\text{Nm en }t_{step2}=0.3\text{s}) \to (-6.28\text{Nm en }t_{step3}=0.5\text{s}) \to (0\text{Nm en }t_{step5}=0.9\text{s}) \to (+6.28\text{Nm en }t_{step7}=1.3\text{s}) \to (-6.28\text{Nm en }t_{step8}=1.5\text{s}) \to (0\text{Nm en }t_{step10}=1.9\text{s})
> \]
> Graficar evolución de variables vs. tiempo (tensiones y corrientes en ambas coordenadas qd0 <=> abcs) y curvas paramétricas torque vs. velocidad vs. curvas características cuasi-estacionarias, evaluar cuadrantes de operación y evolución del ángulo de torque del rotor. También \(i^r_{ds}\) vs \(i^r_{qs}\).
>
> **b)** Determinar velocidad y corriente final de establecimiento luego de cada transitorio, tiempos de crecimiento (10% al 90%), tiempo de establecimiento (±1%), sobrepico, etc. ¿Qué influencia relativa tienen cada una de las dos acciones externas? ¿A qué se debe?
>
> **c)** Comparar comportamiento de \(i^r_{ds}(t)\) para \(i^r_{ds}(0)=\pm 0.5\) A vs \(i^r_{ds}(0)= 0\) A. ¿Qué efecto tiene \(i^r_{ds}(0)\neq 0\) A sobre el sistema, en ambos modelos?
>
> **d)** Agregar luego una consigna de tensión en eje d, \(v^{r*}_{ds}(t) = 0 \to V^r_{qs,nom}/10 = \pm 1.9596\) Vcc en \(t_{step1}=0.5\) s (field forcing/weakening a lazo abierto), sumada a restricción o ley de control NL. ¿Qué efecto tiene sobre el sistema, en ambos modelos?

---

## 5.2 Diseño, Análisis, Simulación e "Implementación" de CONTROLADOR de Movimiento en Cascada con Modulador de Torque (Control Vectorial)

**1)** Modulador de Torque equivalente (Controlador interno vectorial de corriente/torque), con su diagrama de bloques completo a partir del modelo NL completo y valores de parámetros correspondientes, basado en los siguientes lineamientos:

> **a)** "Desacoplamiento" o compensación de todas las retroalimentaciones físicas naturales de estado hacia la entrada. Comparar con la Linealización por Retroalimentación NL completa realizada en ítem 5.1.2.c.VI. Evaluar efecto de estimar y utilizar, en el Controlador, la dependencia \(R_s(T^\circ_s(t))\) vs. usar valor "nominal" constante de \(R_s\) (la planta física siempre tiene dependencia \(R_s(T^\circ_s(t))\)). Comparar desempeño en ambos casos.
>
> **b)** Diseño de lazos de control de corrientes \(\mathbf{i}^r_{qd0s}(t)\) desacoplados entre sí y de la velocidad, con **control proporcional solamente**, con polos en \(p_i = -5000\) rad/s (\(BW \cong 796\) Hz) para todos los ejes. ¿En qué cambia la dinámica, comparada con la obtenida en ítem 5.1.2.c.VI? Investigar el agregado o no de acción integral en los lazos de control de corrientes.
>
> **c)** Incorporación adecuada y completa de **consigna de torque** (nueva variable manipulada) y "desacoplamiento" o compensación de fricción viscosa equivalente.
>
> **d)** "Desacoplamiento" o **Compensación del Torque de carga por gravedad** ("Precomputed Torque", o "torque precalculado" = Linealización por retroalim. NL del Torque de carga).
>
> **Nota:** Simular con Modulador de Torque en DT y comparar con el comportamiento dinámico obtenido en ítem 5.1.6, asumiendo torque de carga y señal consigna de torque necesaria. Comparar requisitos e implementación simple para estrategia base con \(i^{r*}_{ds}(t)\equiv 0\), vs estrategia extendida para \(i^{r*}_{ds}(t)\neq 0\) en general (debilitamiento/reforzamiento de campo magnético).

**2)** Controlador externo de movimientos: posición/velocidad (con "acceso directo" a manipular el torque motor a través de la consigna de torque al modulador interno (para estrategia base con \(i^{r*}_{ds}(t)\equiv 0\)), diseñado utilizando el **método de sintonía serie con acción integral ("PID")**, con:
\[
\zeta = 0.75 \;;\; \omega_n = 800\;\text{rad/s}
\]
considerando valores nominales de \(J_l, b_l\); con su diagrama de bloques completo y valores de parámetros correspondientes. Indicar en el plano s dónde quedan ubicados los polos correspondientes, en comparación con los polos de reguladores de corriente y los polos de la planta original (evaluar la influencia de variación extrema de parámetros de carga \(J_l, b_l\) => migración de polos, etc.).

Incorporar entrada de referencia o setpoint de posición \(q^*_1(t)\equiv\frac{1}{r}\theta^*_m(t)\) al diagrama de bloques del sistema.

**3)** Incorporación y diseño de **Observador de Estado de orden reducido** sólo para la parte mecánica de este controlador, que estime la posición y velocidad a partir de sensor de posición \(\theta_m(t)\) (no es necesario estimar las corrientes, ya que se dispone de sensores de corriente para el control vectorial, modulador de torque). Ubicar los dos polos del observador reales iguales en:
\[
p_{obs\,1,2} = -3200\;\text{rad/s}
\]
para no interferir demasiado con el controlador de estado. Adecuar todas las retroalimentaciones de velocidad a los valores estimados.

**4)** Simulación en tiempo continuo con modelo completo NL, mostrando:

> **a)** Seguimiento de consignas de movimiento \(q^*_1(t)\equiv\frac{1}{r}\theta^*_m(t)\) con **perfil trapezoidal de posición**:
> \[
> q^*_1(t) = 0 \to (\Delta t_{ramp}=5\text{s}) \to 2\pi\;\text{[rad]} \to (\Delta t_{ramp}=5\text{s}) \to 0
> \]
>
> **b)** Rechazo a perturbaciones (cambios en escalón): considerando valores nominales y variación máx. de los parámetros de carga mecánica física.

**5)** Verificación de desempeño y/o mejoras:

> **a)** Verificar si se supera/n alguna/s de las Especificaciones de operación (valores límite) de velocidad, torque, corriente y tensión de los componentes del sistema físico (caja reductora, motor, inversor), o si existe margen para aumentar el desempeño del sistema.
> **Nota:** En caso de superar los valores límites, determinar el origen o causa primera; evaluar qué restricciones es necesario imponer al controlador o consignas para respetar dichos límites => realizar los ajustes y verificar.
>
> **b)** Observador: Verificar si existe error de estimación de régimen permanente distinto de cero ante perturbaciones de carga, o si en este caso también converge asintóticamente la estimación al valor real no medido.
> **Nota:** En caso de tener error de estimación estacionario no nulo, proponer esquema alternativo o agregado para compensar este error y llevarlo a cero (para perturbación constante) => realizar los ajustes y verificar.
>
> **c)** Comportamiento térmico del motor: Verificar si la temperatura del bobinado se mantiene dentro de los valores admisibles para operación continua repetitiva con ciclo de operación especificado, al llegar al "equilibrio térmico".
>
> **d)** Evaluar si aparece alguna degradación de desempeño del sistema cuando se considera la **respuesta no ideal (ancho de banda limitado) de los sensores y acondicionadores de señal**, reemplazando el modelo ideal \(G(s)\equiv 1\) por modelos aproximados equivalentes con características de **filtro Pasa Bajos (LP) con ganancia unitaria**, implementados en el Espacio de Estados (SS):
> - Corrientes \(i_{as}(t), i_{bs}(t), i_{cs}(t)\): modelo LP en SS 2° orden, \(\omega_n = 6000\) rad/s, \(\zeta = 1\).
> - Posición angular \(\theta_m(t)\): modelo LP en SS 2° orden, \(\omega_n = 2000\) rad/s, \(\zeta = 1\).
> - Temperatura \(T^\circ_s(t)\): modelo LP en SS 1° orden, \(\tau = 20\) s.
>
> **Nota:** En la implementación en SS de cada filtro, calibrar sus condiciones iniciales propias en forma consistente con el valor inicial de la señal respectiva a medir, a fin de evitar transitorios iniciales de medición.
> Evaluar el efecto de aumentar \(\omega_n \to \omega_n\times 2 \to \omega_n\times 3\), etc., de los filtros posición y corriente, y su interacción con \(\omega_n\) del observador, transformación de Park y lazos de corriente, respectivamente.
>
> **e)** Evaluar adicionalmente si aparece alguna degradación de desempeño del sistema cuando se considera la **respuesta no ideal (Saturación y ancho de banda limitado) en el modulador trifásico de tensión** que representa como modelo promediado el Inversor, reemplazando el modelo ideal \(G(s)\equiv 1\) por modelo aproximado equivalente con **saturación** y características de **filtro Pasa Bajos (LP) con ganancia unitaria**, implementados en el Espacio de Estados (SS):
> - **Saturación:** \(|v_{as}(t)|, |v_{bs}(t)|, |v_{cs}(t)| \leq \sqrt{2}\cdot\frac{V_{sl,max}}{\sqrt{3}}\), con \(V_{sl,max} = 48\) Vca rms
> - **Tensiones** \(v_{as}(t), v_{bs}(t), v_{cs}(t)\): modelo LP en SS 2° orden, \(\omega_n = 6000\) rad/s, \(\zeta = 1\).

**6)** Diseño Final: Finalmente, mostrar Controlador Completo con todas las etapas de control integradas correctamente en el DT continuo \(t\in\mathbb{R}\). **Discretizar** el Controlador completo considerando muestreo en instantes \(t_k\equiv k\cdot T_s\), \(k\in\mathbb{Z}\), utilizando el **método de Tustin** (integración numérica por Trapecios) con periodo de muestreo único y constante \(T_s\) [s] \(\in\mathbb{R}\) y retenedor de orden cero (ZOH) para las señales de control; determinar \(T_s\) para desempeño adecuado del controlador en tiempo discreto.

**7)** Codificación de Software de Control: reemplazar el modelo completo de todo el controlador discreto diseñado por **código programado en lenguaje de texto estructurado** (lenguaje Matlab) dentro de `Matlab function` y validar su comportamiento, comparativamente con el modelo en bloques obtenido en el punto 6. Utilizar funciones escalares desagregadas, precálculo de todas las variables que no dependen del instante actual, etc. y documentar código (en forma similar a programar un microcontrolador).

---

## 6. Referencias

1. R. Kelly et al, *Control of Robot Manipulators in Joint Space*, Springer, 2005. (Example & Figure 2.2).
2. P. Krause et al, *Analysis of Electric Machinery and Drive Systems*, 3rd Ed., IEEE-Wiley, 2013.
3. G. Franklin et al, *Feedback Control of Dynamic Systems*, 7th Ed., Pearson, 2015.
