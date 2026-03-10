# Instrucciones para agente — Correcciones del informe hasta 5.1

## Alcance
Trabajar **únicamente hasta la sección 5.1 inclusive**.  
**No modificar ni desarrollar 5.2 en adelante.**

El objetivo es corregir el informe en todo lo relativo a 5.1, con foco en:

- inconsistencias matemáticas y de notación,
- parámetros y tablas mal calculadas,
- revisión del código y scripts asociados,
- coherencia entre texto, ecuaciones y resultados,
- cumplimiento exacto de la consigna hasta 5.1.

No modificar imágenes salvo que una corrección del informe obligue a actualizar su pie, referencia o interpretación. La revisión visual final de imágenes la haré manualmente.

---

## Entregables esperados
Al finalizar, entregar:

1. archivo `.tex` corregido hasta 5.1,
2. tablas recalculadas,
3. scripts/código corregidos y consistentes,
4. changelog breve con todas las modificaciones,
5. checklist final de cumplimiento de 5.1 contra la consigna.

---

# Fase 0 — Regla de trabajo

## Tarea 0.1
No tocar contenido de 5.2 en adelante.

Limitar todos los cambios a:

- portada, resumen e introducción si hace falta por consistencia,
- sección 5.1 completa,
- tablas, ecuaciones, referencias cruzadas y bibliografía usada en 5.1,
- scripts y código que generan resultados de 5.1.

**Entregable:** confirmación de alcance congelado en 5.2+.

---

# Fase 1 — Auditoría global de consistencia

## Tarea 1.1 — Unificar notación
Revisar todo el documento y unificar símbolos repetidos o ambiguos. Corregir especialmente:

- `T_l(t)` = torque total de carga,
- `T_ld(t)` = perturbación adicional/disturbio mecánico, solo si realmente se quiere separar de la gravedad,
- `R_s(T_s^\circ)`, `R_s0`, `R_sREF`: elegir una sola convención clara,
- `R_th` vs `R_ts-amb`: dejar una sola,
- `T_s^\circ` vs `T_s`,
- `λ'_m`, `λ'^r_m`, `λ'_{mr}`: dejar una sola,
- `i_qs`, `i_qs^r`, `i_qs^r(t)`: normalizar,
- `θ_r`, `θ_e`, `θ_ev`: aclarar cuál es cuál y mantenerlo igual en todo el texto.

**Entregable:** tabla interna de notación unificada y reemplazo en el `.tex`.

---

## Tarea 1.2 — Corregir errores tipográficos matemáticos
Buscar y corregir errores evidentes de escritura matemática, por ejemplo:

- apariciones de `ω_n` donde debía decir `ω_m`,
- cadenas incorrectas como `\dot{θ}_r = ω_r = θ_r = ∫...`,
- signos incorrectos,
- términos con `g` faltante o duplicado,
- referencias a ecuaciones inexistentes o mal numeradas.

**Entregable:** listado de correcciones puntuales aplicadas.

---

## Tarea 1.3 — Coherencia entre texto, ecuaciones y tablas
Verificar que cada parámetro usado en ecuaciones:

- exista en la tabla de parámetros,
- tenga el mismo símbolo en todo el informe,
- tenga las mismas unidades,
- tenga el mismo valor en todas sus apariciones.

**Entregable:** tabla de chequeo `símbolo / valor / unidad / apariciones`.

---

# Fase 2 — Corrección de parámetros físicos y tablas

## Tarea 2.1 — Recalcular parámetros derivados
Recalcular desde cero, con script reproducible, todos los parámetros derivados usados hasta 5.1, especialmente:

- `J_l` nominal,
- `k_l`,
- `J_eq = J_m + J_l / r^2`,
- `b_eq = b_m + b_l / r^2`,
- constantes de tiempo eléctricas y mecánicas,
- constante de torque `k_t = (3/2) P_p λ'_m`.

Usar como nominales los valores definidos por la consigna.

**Entregable:** script de parámetros + tabla final de valores corregidos.

---

## Tarea 2.2 — Corregir la tabla general de parámetros
Rehacer la tabla de parámetros para que:

- no tenga valores inconsistentes,
- separe claramente:
  - parámetros dados por consigna,
  - parámetros calculados,
  - rangos de variación,
- use unidades homogéneas,
- marque explícitamente cuáles son nominales y cuáles variables.

**Entregable:** nueva tabla principal de parámetros.

---

## Tarea 2.3 — Corregir tablas de migración
Recalcular y rehacer las tablas de:

- migración con `R_s`,
- migración con `m_l`,
- migración con `i_ds0^r`,

asegurando que usen el mismo conjunto de parámetros nominales corregidos.

**Entregable:** tablas nuevas + script que las genera.

---

# Fase 3 — Revisión algebraica del modelo hasta 5.1

## Tarea 3.1 — Modelo mecánico equivalente
Auditar y corregir la derivación del modelo mecánico equivalente referido al eje motor:

- transformación por relación `r`,
- torque y velocidad,
- inercia y fricción equivalentes,
- inclusión de gravedad,
- definición de salida.

Asegurarse de que esta subsección responda directamente al ítem 5.1.1 de la consigna.

**Entregable:** subsección corregida y compacta.

---

## Tarea 3.2 — Modelo NL completo
Revisar el modelo no lineal completo y dejarlo formalmente presentado como:

- vector de estado,
- vector de entradas,
- vector de perturbaciones,
- vector de salidas medidas,
- ecuaciones de estado,
- ecuaciones de salida.

Alinear esto con la consigna y con las variables físicas realmente usadas en simulación.

**Entregable:** bloque formal del modelo NL limpio y consistente.

---

## Tarea 3.3 — Secuencia cero
Revisar la demostración de irrelevancia de la componente homopolar para:

- neutro flotante,
- sistema equilibrado.

Mantener la idea, pero:

- compactarla,
- eliminar redundancias,
- asegurar que no haya afirmaciones demasiado fuertes sobre `v_0` si no están justificadas.

**Entregable:** demostración más limpia y rigurosa.

---

## Tarea 3.4 — Modelo LPV por Jacobiano
Revisar línea por línea la linealización Jacobiana:

- punto de operación,
- variables incrementales,
- matrices `A_0`, `B_0^c`, `B_0^d`,
- términos térmicos,
- signos,
- términos cruzados,
- consistencia dimensional.

Corregir cualquier error en:

- derivadas parciales,
- dependencia de `R_s(T_s)`,
- notación de `R_s,ref`,
- constante térmica `Γ_th`.

**Entregable:** versión auditada del modelo LPV, con script simbólico o numérico de verificación si es posible.

---

## Tarea 3.5 — Modelo LTI equivalente y aumentado
Revisar y corregir la sección de linealización por realimentación no lineal:

- ley mínima en eje `d`,
- dinámica residual de `i_ds`,
- ley complementaria en eje `q`,
- modelo LTI nominal,
- modelo aumentado.

Aclarar expresamente en el texto que:

- el subsistema dinámico es lineal,
- la gravedad quedó separada como perturbación no lineal exógena si así se mantiene,
- no llamar “LTI completo” a algo que aún depende de `sin(θ_m/r)` en un canal externo.

**Entregable:** subsección corregida, con terminología rigurosa.

---

## Tarea 3.6 — Funciones de transferencia
Revisar la derivación de las funciones de transferencia desde:

- `v_qs^r(t)` hacia `θ_m(t)`,
- perturbación de carga hacia `θ_m(t)`,

y confirmar:

- orden correcto,
- cancelación correcta del polo asociado a `i_ds`,
- denominador común,
- ceros,
- interpretación física.

Corregir la nomenclatura para no mezclar `T_l` con `T_ld`.

**Entregable:** derivación depurada + chequeo con herramienta automática (`ss2tf` o equivalente).

---

# Fase 4 — Revisión de análisis estructural y dinámico

## Tarea 4.1 — Estabilidad
Recalcular autovalores del modelo aumentado con parámetros corregidos y rehacer:

- tabla de autovalores,
- interpretación física,
- `ω_n`, `ζ`, `τ`,
- estabilidad interna,
- estabilidad BIBO,
- estabilidad parcial.

**Entregable:** tabla y texto corregidos.

---

## Tarea 4.2 — Migración de polos y ceros
Recalcular migración de polos y ceros al variar:

- `R_s`,
- `m_l`,
- `i_ds0^r`,

usando scripts reproducibles.

Verificar que los valores sean físicamente plausibles y que correspondan exactamente al modelo usado en el informe.

**Entregable:** datos corregidos + script + texto interpretativo actualizado.

---

## Tarea 4.3 — Observabilidad
Verificar formalmente el análisis de observabilidad desde:

- `θ_m`,
- `ω_m` como alternativa.

Corroborar rango con cálculo automático y corregir el texto si hubiera alguna imprecisión.

**Entregable:** subsección validada.

---

## Tarea 4.4 — Controlabilidad
Verificar formalmente el análisis de controlabilidad desde:

- `v_qs^r`,
- alternativa incorporando `v_ds^r`.

Corroborar rango con cálculo automático y revisar coherencia con la interpretación física.

**Entregable:** subsección validada.

---

# Fase 5 — Revisión completa del código y resultados numéricos de 5.1

## Tarea 5.1 — Auditoría del código base
Localizar todos los scripts y modelos que generan resultados de 5.1 y verificar:

- que usen el mismo set de parámetros,
- que las unidades sean consistentes,
- que no mezclen rpm con rad/s,
- que no mezclen ángulo mecánico con eléctrico,
- que no mezclen torque total con disturbio.

**Entregable:** informe corto de auditoría de código.

---

## Tarea 5.2 — Corregir resultados físicamente imposibles
Detectar y corregir el origen de resultados no plausibles, especialmente:

- velocidades absurdamente altas,
- tiempos de establecimiento incoherentes,
- sobrepicos inconsistentes,
- tablas con magnitudes fuera de rango del motor/inversor.

**Entregable:** lista de errores encontrados en código + corrección aplicada.

---

## Tarea 5.3 — Unificar scripts de generación de tablas y gráficos
Crear o consolidar scripts para que cada tabla/resultado de 5.1 salga automáticamente desde código, evitando valores tipeados a mano.

**Entregable:** carpeta de scripts organizada con nombres claros.

---

## Tarea 5.4 — Verificación cruzada analítico vs simulación
Comparar, donde corresponda:

- autovalores analíticos vs `eig(A)`,
- funciones de transferencia analíticas vs `tf(ss(...))`,
- respuesta libre de `i_ds` vs exponencial teórica,
- constantes de tiempo eléctricas vs respuesta simulada.

**Entregable:** cuadro de validación cruzada.

---

# Fase 6 — Corrección de la sección 5.1.6 de simulación dinámica

## Tarea 6.1 — Verificar señales de excitación
Asegurar que las señales usadas en simulación coincidan exactamente con la consigna para 5.1.6:

- pulso de `v_qs^{r*}(t)` en 0.1, 0.7, 1.1, 1.7 s,
- doble pulso de torque en 0.3, 0.5, 0.9, 1.3, 1.5, 1.9 s.

**Entregable:** bloque de generación de entradas validado.

---

## Tarea 6.2 — Completar variables requeridas por la consigna
Verificar que la sección de simulación incluya, al menos en resultados numéricos y descripción, todas las variables que pide 5.1.6.a:

- `θ_m(t)`,
- `ω_m(t)`,
- `i_qd0^r(t)`,
- `T_s^\circ(t)`,
- `v_ds^r(t)` forzada,
- tensiones y corrientes en `qd0` y `abc`,
- curvas paramétricas torque vs velocidad,
- `i_d` vs `i_q`,
- evolución del ángulo de torque `δ(t)`.

**Entregable:** checklist de variables cubiertas; agregar texto o referencias a figuras si faltan.

---

## Tarea 6.3 — Corregir el análisis de métricas transitorias
Rehacer la tabla de métricas de transitorio para que los valores:

- sean físicamente plausibles,
- usen unidades correctas,
- correspondan al estado correcto,
- distingan claramente cuándo hay y cuándo no hay nuevo régimen.

**Entregable:** nueva tabla de métricas.

---

## Tarea 6.4 — Revisar el caso `i_ds^r(0) ≠ 0`
Verificar por simulación y por fórmula la respuesta de `i_ds^r(t)` para:

- `0`,
- `+0.5 A`,
- `-0.5 A`,

y chequear que el texto describa exactamente el comportamiento observado.

**Entregable:** subsección validada.

---

## Tarea 6.5 — Corregir y separar bien field forcing / field weakening
Reescribir el caso de 5.1.6.d para que quede claro que son **dos simulaciones distintas** o dos casos distintos:

- `v_ds^* > 0`: field forcing,
- `v_ds^* < 0`: field weakening.

Corregir la definición por tramos de `v_ds^*(t)`, que hoy es ambigua.

**Entregable:** subsección corregida + ecuación válida.

---

## Tarea 6.6 — Cerrar el ángulo de carga `δ(t)`
Completar la parte del ángulo de carga:

- verificar definición,
- obtenerlo desde simulación,
- agregar interpretación física,
- conectarlo con la operación en cuatro cuadrantes si corresponde.

**Entregable:** subsección terminada o, si no aporta valor real, eliminarla del texto y dejar solo lo exigido por la consigna.

---

# Fase 7 — Limpieza editorial de 5.1

## Tarea 7.1 — Compactar repeticiones
Reducir repeticiones innecesarias en:

- LPV,
- comparación LPV vs LTI,
- dinámica residual,
- explicaciones sobre desacoplamiento.

**Entregable:** texto más limpio, sin pérdida de contenido técnico.

---

## Tarea 7.2 — Corregir referencias internas
Revisar todas las referencias:

- `\label`,
- `\ref`,
- `\eqref`,
- figuras,
- tablas,
- subsecciones,

para evitar referencias rotas o equivocadas.

**Entregable:** compilación sin warnings de referencias rotas.

---

## Tarea 7.3 — Revisar bibliografía usada en 5.1
Verificar que todas las citas usadas en 5.1:

- existan en `references.bib`,
- tengan formato consistente,
- correspondan al contenido citado.

**Entregable:** bibliografía mínima correcta para 5.1.

---

# Fase 8 — Validación final de la sección 5.1

## Tarea 8.1 — Checklist contra la consigna
Construir un checklist final `ítem de consigna vs evidencia en el informe` para:

- 5.1.1,
- 5.1.2.a,
- 5.1.2.b,
- 5.1.2.c,
- 5.1.2.d,
- 5.1.2.e,
- 5.1.3,
- 5.1.4,
- 5.1.5,
- 5.1.6.a,
- 5.1.6.b,
- 5.1.6.c,
- 5.1.6.d.

**Entregable:** tabla de cumplimiento de 5.1.

---

## Tarea 8.2 — Entrega técnica final
Entregar:

- archivo `.tex` corregido hasta 5.1,
- tablas recalculadas,
- scripts corregidos,
- changelog de todo lo modificado,
- lista de figuras que podrían requerir revisión manual, sin modificarlas salvo que una corrección del informe obligue a actualizar su pie o referencia.

**Entregable:** paquete final listo para revisión manual.

---

# Orden recomendado de ejecución

1. Unificar notación y parámetros.  
2. Recalcular tabla de parámetros y derivados.  
3. Auditar modelo mecánico, NL, LPV y LTI.  
4. Revisar código y corregir resultados imposibles.  
5. Rehacer tablas de polos, migraciones y métricas.  
6. Cerrar 5.1.6 completo.  
7. Limpiar redacción, referencias y bibliografía.  
8. Hacer checklist final contra la consigna.

---

# Prompt breve alternativo para agente

Revisá y corregí el informe únicamente hasta la sección 5.1 inclusive. No modifiques ni desarrolles 5.2 en adelante. Tu trabajo debe enfocarse en:  
1) unificar notación y símbolos;  
2) recalcular todos los parámetros derivados y corregir tablas inconsistentes;  
3) auditar algebraicamente el modelo mecánico equivalente, el modelo NL completo, el LPV Jacobiano, el modelo LTI equivalente y aumentado, y las funciones de transferencia;  
4) revisar el código y scripts asociados a 5.1 para detectar errores de unidades, parámetros, señales o resultados físicamente imposibles;  
5) rehacer tablas de autovalores, migración de polos/ceros y métricas transitorias con valores consistentes;  
6) verificar que la simulación 5.1.6 cumpla la consigna exacta y completar lo faltante en variables requeridas, especialmente tensiones/corrientes en qd0 y abc, `i_d` vs `i_q`, `v_ds` forzada y ángulo `δ(t)`;  
7) corregir redacción técnica, referencias internas y bibliografía de 5.1.  
No cambies imágenes salvo que una corrección del informe obligue a actualizar su pie o su referencia. Entregá el `.tex` corregido, scripts corregidos, tablas recalculadas y un changelog de modificaciones.

## Regla adicional de edición de archivos

No modificar directamente `informe.tex`.

Crear una copia nueva llamada `informe_revisado.tex` y realizar allí todos los cambios correspondientes hasta la sección 5.1.

Mantener `informe.tex` como respaldo intacto.

Si se corrigen o rehacen scripts auxiliares, no sobrescribir los originales sin necesidad. Crear versiones nuevas con nombres claros, por ejemplo:

- `parametros_revisados.m`
- `tablas_5_1_revisadas.m`
- `simulacion_5_1_revisada.m`

El objetivo es preservar todos los archivos originales y concentrar las correcciones en una versión revisada, fácil de comparar contra la original.