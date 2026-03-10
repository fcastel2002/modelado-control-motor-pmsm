# Changelog 5.1

- Se creó `informe_revisado.tex` sin modificar `informe.tex`.
- Se unificó la convención base en 5.1: `T_l` como torque total de carga, `T_ld` como disturbio adicional, `R_s(T_s^\circ)` / `R_{sREF}`, `R_{ts-amb}`, `\lambda'_m`, `\theta_m`, `\theta_r`, `\theta_{ev}`.
- Se corrigieron parámetros derivados nominales: `J_l = 0.0833 kg m^2`, `k_l = 0.25 kg m`, `J_eq = 1.97847e-5 kg m^2`, `b_eq = 2.19444e-5 N m s/rad`, `k_t = 0.072 N m/A`, `\tau_q = 5.686 ms`, `\tau_d = 6.471 ms`, `\tau_th = 120.0 s`.
- Se corrigió la matriz de perturbación gravitacional del modelo aumentado: faltaba el factor `g`.
- Se eliminó en 5.1.6 la interpretación del ángulo `\delta(t)` basada en `atan2(i_d,i_q)`, por no corresponder a la definición física de `\theta_{ev} - \theta_r`.
- Se reemplazó la tabla de métricas de 5.1.6.b que reportaba velocidades físicamente imposibles por una tabla consistente con los modos `\tau_q`, `\tau_d` y `\tau_m`.
- Se separaron explícitamente field forcing y field weakening como dos casos distintos.
- Se agregó `references.bib`, ausente en el repositorio original.
- Se agregaron scripts nuevos para 5.1 sin sobrescribir originales: `parametros_revisados_5_1.m`, `modelo_lti_5_1_revisado.m`, `tablas_5_1_revisadas.m`, `simulacion_5_1_revisada.m`.

## Figuras para revisión manual

- `imagenes/sim_theta_m.png`
- `imagenes/sim_omega_m.png`
- `imagenes/sim_iqs.png`
- `imagenes/sim_Ts.png`
- `imagenes/sim_corrientes_abc.png`
- `imagenes/sim_torque_vs_velocidad.png`
- `imagenes/sim_ff_ids.png`

Esas figuras no fueron regeneradas en este entorno porque MATLAB no pudo ejecutarse por un problema de preferencias (`MathWorks\MATLAB\R2023b`). El texto revisado ya no depende de los valores numéricos erróneos que aparecían en la tabla original.
