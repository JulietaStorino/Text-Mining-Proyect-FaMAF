# Trabajo de presentación para el curso: "Minería de Datos para Textos" (Text Mining)
## Impacto de la Dependencia de Datos en Modelos de Fine-Tuning para la anonimización de textos de ciencias de la vida
**Profesores**: Laura Alonso Alemany, Matias Eduardo Bordone Carranza, Milagro Teruel  
**Alumnos**: Florencia Brunello, Julieta Paola Storino, Santiago Troya.

### Índice
1. [Introducción](#introducción)
2. [Motivación](#motivación)
3. [Trabajos anteriores](#trabajos-anteriores)
4. [Modelos a testear](#modelos-a-testear)
5. [DataSet](#dataset)
6. [Planificación](#planificación)

### Introducción

El texto (lenguaje natural) es el camino más natural para codificar el conocimiento humano. Como resultado, la mayor parte del conocimiento humano se codifica en forma de datos de texto. Por ejemplo, el conocimiento científico existe casi exclusivamente en la literatura científica, mientras que los manuales técnicos contienen explicaciones detalladas sobre cómo operar los dispositivos.
El texto es, con diferencia, el tipo de información más común con el que se encuentran las personas. De hecho, la mayor parte de la información que una persona produce y consume a diario está en forma de texto. [ChengXiang Zhai & Sean Masung 2016] 

### Motivación

En el ámbito del cuidado de la salud se genera una gran cantidad de datos valiosos para la detección y caracterización de enfermedades. El acceso a estos datos se ha facilitado gracias a la digitalización del sistema de salud, permitiendo el procesamiento de grandes cantidades de datos, mediante técnicas como el procesamiento del lenguaje natural (PLN) y el aprendizaje automático, que permite analizar información sobre la salud de las poblaciones y descubrir factores sociales que influyen en los resultados de salud.

Sin embargo, estos textos médicos están compuestos, en gran medida, por datos relativos al paciente y requieren un manejo cuidadoso para evitar la identificación de las personas, ya que pueden aludir a información que históricamente ha dado lugar a diversos tipos de discriminación.

Esto plantea retos en cuanto a la protección de la privacidad de los datos de los pacientes, en la cual resalta el delicado equilibrio entre la utilidad de los datos y la privacidad, ya que, se pueden anonimizar datos hasta un punto en que no brinden información beneficiosa, volviéndose inútil, o por lo contrario, se puede tener datos que brinden información sustancial para ser vinculados a las personas a las que refieren.

Por estas razones, investigadores de todo el mundo han desarrollado técnicas y algoritmos avanzados para anonimizar datos, permitiendo su uso para los fines solicitados mientras se mantiene el anonimato del paciente. Esto nos da paso a analizar las fortalezas y debilidades de los algoritmos y software que anonimizan datos de salud estructurados, con el objetivo de identificar áreas de mejora y diseñar estrategias para optimizar su rendimiento.

### Trabajos anteriores

Este trabajo se centra en la campaña “Medical Document Anonymization (MEDDOCAN)” desarrollado por el Plan de Impulso de las tecnologías del Lenguaje del Gobierno de España. En el marco de esta campaña, se creó un corpus sintético de 1000 registros clínicos con información de salud protegida (PHI, por sus siglas en inglés), denominado “Corpus MEDDOCAN”. Este corpus se distribuyó en texto simple con codificación UTF8, donde cada caso clínico se almacena como un único archivo de texto y las anotaciones PHI se publican en formato BRAT.

<p align="center">
  <img width="auto" height="500" src="https://temu.bsc.es/meddocan/wp-content/uploads/2019/03/image-1-768x692.png" alt="Figura 1: Un ejemplo de anotación MEDDOCAN visualizada utilizando la interfaz de anotación BRAT.">
</p>

Estos casos clínicos se anotaron manualmente utilizando una versión personalizada de AnnotateIt. Las pautas de anotación de MEDDOCAN definieron un total de 29 tipos de entidades. La Tabla 1 resume la lista de tipos de entidades sensibles definidos para el seguimiento de MEDDOCAN y la cantidad de ocurrencias entre los conjuntos de entrenamiento, desarrollo y prueba.

<p align="center">
  <img width="auto" height="500" src="https://i.ibb.co/L52R0sd/Screenshot-from-2024-09-27-09-09-19.png" alt="Tabla 1: Distribución del tipo de entidad entre los conjuntos de datos.">
</p>

Durante el proyecto, se desarrollaron diversos modelos de fine-tuning, entrenados con 500 casos clínicos del corpus, seleccionados aleatoriamente. Los casos restantes se utilizaron para desarrollo y prueba. La evaluación de predicciones automáticas para este track tuvo dos escenarios o categorías diferentes:
1. *NER offset and entity type classification*: Enfocado en la identificación y clasificación de información sensible, como nombres de pacientes, teléfonos, direcciones, etc.
2. *Sensitive span detection*: Enfocado en la detección de texto sensible específico para la liberación de documentos clínicos desidentificados, donde el objetivo es identificar y enmascarar datos confidenciales, independientemente del tipo real de entidad o la correcta identificación del tipo de PHI.

Todo esto nos permite enfocarnos en evaluar la mejora de los modelos desarrollados en aspectos generales, fuera del ámbito en el que fueron entrenados. Es decir, utilizamos casos clínicos generados que no forman parte del corpus MEDDOCAN para comprobar qué tan dependientes son de los datos con los que fueron entrenados y cómo esto afecta su desempeño en otros conjuntos de datos.

### Modelos a testear

De entre todos los proyectos presentados por los 18 equipos participantes en ambas categorías, seleccionamos a los 5 proyectos con mejor desempeño basado en las métricas de evaluación primarias utilizadas para la evaluación de los modelos, entre ellos, precisión, sensibilidad y F-Score, siendo esta última la única medida de evaluación oficial de ambas categorías. Los modelos seleccionados son los presentados por las siguientes instituciones y autores:
1. Bosch Center for Artificial Intelligence
2. Universitat Rovira i Virgili, CRISES group 
3. Vicomtech
4. FSL
5. Universitat Rovira i Virgili, iTAKA Research Group

### DataSet

A los datos brindados por el corpus MEDDOCAN sumaremos historias médicas desarrolladas por el generador de pacientes sintéticos open source [Synthea](https://github.com/synthetichealth/synthea/wiki) que proporciona datos de alta calidad, realistas, pero no reales, de pacientes; en una variedad de formatos, con diferentes niveles de complejidad, cubriendo todos los aspectos de la atención médica. Los datos obtenidos no tienen coste, ni privacidad ni restricciones de seguridad.

### Planificación
* *Semana 1*: Definir los objetivos del proyecto y el alcance del testeo de los modelos. 

* *Semana 2*: Preparar el conjunto de datos. Adaptar modelos existentes para transformar los datos de entrada y salida a los formatos requeridos por los nuevos modelos a testear.

* *Semana 3*: Configurar los modelos a testear y desarrollar el modelo de evaluación para comparar los resultados obtenidos.

* *Semana 4*: Ejecutar los modelos y realizar una evaluación preliminar de los resultados obtenidos.

* *Semana 5*: Analizar los resultados obtenidos en detalle y redactar el informe final, incluyendo conclusiones y recomendaciones para futuros trabajos.

