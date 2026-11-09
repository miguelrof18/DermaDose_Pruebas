# Procesamiento TC–RDSR por lotes

Exporte cada TC desde DermaDose Surface 1.8.0. Ponga el STL, el JSON homónimo y un único RDSR DICOM en una carpeta por caso. La versión 1.8.0 propone un punto diana de tórax si identifica dos regiones pulmonares de baja densidad, cerradas y con disposición bilateral plausible. Si el TC está truncado o el criterio falla, se conserva la superficie y se documenta que falta el anclaje. Los JSON antiguos requieren una nueva exportación desde el volumen.

En el entorno Python de DermaDose:

```powershell
python -m skindose_verifier.batch_runner "C:\casos" --config "C:\config\equipo_azurion.json" --prepare-only
python -m skindose_verifier.batch_runner "C:\casos" --config "C:\config\equipo_azurion.json"
```

La preparación comprueba las entradas y escribe la colocación estimada antes del cálculo. El JSON de equipo es común y las rutas relativas a tablas se resuelven respecto a ese archivo. Cada subcarpeta deja su log y sus resultados. El resumen `dermadose_lote.json` recoge estados y el lote continúa después de un error de cálculo de un caso.

El punto pulmonar aproxima la región torácica y se sitúa en el eje del primer haz con la superficie apoyada en el plano de mesa configurado. Para una simulación sin referencia absoluta de mesa, puede usarse `reference_table_from_first_event` en el perfil Azurion de investigación como hipótesis explícita. La colocación es exploratoria. El método no estima registro a partir de proyecciones XA, movimientos del paciente respecto a la mesa ni incertidumbre dosimétrica de colocación.
