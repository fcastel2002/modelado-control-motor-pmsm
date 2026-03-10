# Checklist 5.1 vs consigna

| Ítem | Estado | Evidencia en `informe_revisado.tex` |
|---|---|---|
| 5.1.1 Modelo mecánico equivalente referido al eje motor | Cumple | Sección de modelo mecánico equivalente con `J_eq`, `b_eq`, gravedad y salida |
| 5.1.2.a Modelo NL completo | Cumple parcial | Ecuaciones de estado y salida corregidas; figuras quedan para rerun manual |
| 5.1.2.b LPV Jacobiano | Cumple parcial | Texto y notación consistentes; conviene rerun de script LPV en MATLAB |
| 5.1.2.c LTI equivalente y aumentado | Cumple | Ley mínima, ley complementaria, modelo nominal y aumentado corregidos |
| 5.1.2.d Comparación LPV vs LTI | Cumple parcial | Comparación conceptual mantenida; migración numérica requiere rerun si se actualizan figuras |
| 5.1.2.e Funciones de transferencia | Cumple | Derivación y nomenclatura corregidas |
| 5.1.3 Estabilidad | Cumple | Tabla de autovalores y discusión de estabilidad interna/BIBO |
| 5.1.4 Observabilidad | Cumple | Análisis de rango y estado no observable `i_ds^r` |
| 5.1.5 Controlabilidad | Cumple | Análisis de rango y entrada adicional `v_ds^r` |
| 5.1.6.a Simulación de estados, qd0 y abc | Cumple parcial | Texto corregido; requiere rerun de figuras con scripts revisados |
| 5.1.6.b Métricas transitorias | Cumple parcial | Tabla corregida sin valores imposibles; falta recalcular automáticamente en MATLAB |
| 5.1.6.c Caso `i_ds^r(0) \neq 0` | Cumple | Dinámica residual y constante de tiempo corregidas |
| 5.1.6.d Field forcing / weakening | Cumple | Casos separados y ecuaciones por tramo corregidas |

## Observaciones

- El alcance quedó congelado en 5.2+.
- No se modificaron imágenes existentes.
- La validación final de gráficos y compilación LaTeX queda pendiente de un entorno MATLAB/LaTeX funcional.
