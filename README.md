# Metagenomics-Pipeline

El presente esquema de trabajo incluye un pipeline que permite el análisis de datos metagenómicos obtenidos por NGS. El workflow permite realizar análisis a partir de reads y análisis a partir de ensamblajes. Los análisis incluyen: 

- Preprocesamiento (Calidad y Filtrado)
- Mapeo y decontaminación
- Asignación Taxonómica
- Análisis Funcional
- Ensamblaje 
- Evaluación de Calidad de Ensamblajes
- Binning
- Refinamiento y Calidad de Bins
- Asignación Taxonómica Bins
- Evaluación y Construcción de MAGs

El presente pipeline está implementado dentro de:

- Proyecto: Interacción Parásito-Microbiota_Hospedero en modelo murino en infección por Trypanosoma cruzi (Cepa Tulahuen)
- Datos obtenidos por shotgun metagenomics (4Gb por muestra)
- Ratones BALBc y BL6
- Para cada modelo hay tomas a los 0-3-7-10-13-16 días 
- Los días 0 - 10 - 16 tienen un NI
- NI = No Infectado
. DPI: Días Post Infección
