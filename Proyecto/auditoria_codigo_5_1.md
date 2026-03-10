# Auditoría breve de código 5.1

## Hallazgos principales

1. `simulacion_dinamica_516.m` mezcla el modelo LTI reducido con una dinámica de `i_q` sin fuerza contraelectromotriz en la simulación temporal.
2. La tabla de métricas de 5.1.6.b toma como corriente final `19.21 A = 19.596 / 1.02`, lo que ignora completamente la f.c.e.m. `P_p \lambda'_m \omega_m`.
3. La definición de `\delta(t)` usada en el informe no coincide con la consigna: `atan2(i_d, i_q)` no es `\theta_{ev} - \theta_r`.
4. En el modelo aumentado escrito en el informe faltaba el factor `g` en el canal gravitacional.
5. `evaluar_migracion_Ids.m` es interactivo y no sirve como base reproducible para tablas del informe.
6. `generar_modelo_LTI.m` combina linealización física, visualización y selección interactiva de escenarios en un mismo script.

## Correcciones propuestas

- Separar parámetros, modelo y tablas en archivos nuevos no interactivos.
- Mantener `T_l` para torque total y `T_ld` solo para disturbio adicional.
- No reutilizar la tabla original de 5.1.6.b.
- Regenerar las figuras de 5.1.6 con `simulacion_5_1_revisada.m` cuando MATLAB esté operativo.
