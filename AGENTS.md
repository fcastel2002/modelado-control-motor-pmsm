# Instrucciones del repositorio

## Alcance

- El trabajo activo para futuras sesiones es editar documentos LaTeX y generar sus compilados; no realizar tareas de código MATLAB/Simulink salvo solicitud explícita.
- Modificar únicamente los archivos y el contenido necesarios para la solicitud; no hacer refactors, cambios de formato ni limpiezas no solicitadas.
- No abrir, leer ni analizar archivos `.pdf`: consumen tokens y son salidas o material de referencia. Usar `.tex`, `.md` y `references.bib` como fuentes editables.
- No limpiar ni revertir artefactos MATLAB/Simulink modificados por el usuario o por una simulación si no forman parte de la tarea solicitada.

## Documentos principales

- `Trabajo_Final___AyME___2026/main.tex`: informe técnico; sus imágenes se resuelven desde `Trabajo_Final___AyME___2026`, por lo que debe compilarse allí.
- `Trabajo_Final___AyME___2026/main_presentacion.tex`: presentación Beamer; también espera `imagenes/` relativa a ese directorio.
- `Trabajo_Final___AyME___2026/registro_cambios_consigna.tex`: registro auxiliar de cambios.
- `Proyecto/consignas.md` y `Proyecto/Guía para preparar INFORME TÉCNICO.md`: especificaciones y criterios de formato; consultarlos antes de cambiar estructura, nomenclatura o referencias.

## Compilación

Ejecutar los comandos desde `Trabajo_Final___AyME___2026`:

```text
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
pdflatex main_presentacion.tex
pdflatex main_presentacion.tex
pdflatex registro_cambios_consigna.tex
```

- El informe usa `references.bib` mediante `\bibliographystyle{ieeetr}` y `\bibliography{references}`; por eso requiere la secuencia `pdflatex -> bibtex -> pdflatex -> pdflatex`.
- Si se editan referencias cruzadas, imágenes o numeración, repetir las pasadas de `pdflatex` hasta que desaparezcan las advertencias relevantes.
- Los auxiliares de LaTeX (`.aux`, `.log`, `.nav`, `.out`, `.snm`, `.toc`, `.fls`, `.fdb_latexmk`, `.synctex.gz`) están ignorados por Git. No convertirlos en fuente ni editar compilados para corregir errores.
