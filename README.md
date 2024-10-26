# Trabajo de presentación para el curso: "Minería de Datos para Textos" (Text Mining)
## Impacto de la Dependencia de Datos en Modelos de Fine-Tuning para la anonimización de textos de ciencias de la vida
**Profesores**: Laura Alonso Alemany, Matias Eduardo Bordone Carranza, Milagro Teruel  
**Alumnos**: Florencia Brunello, Julieta Paola Storino, Santiago Troiano.

### Índice
1. [Introducción](#introducción)
2. [Motivación](#motivación)
3. [Objetivos](#objetivos)
4. [Trabajos anteriores](#trabajos-anteriores)
5. [Modelos a testear](#modelos-a-testear)
6. [DataSet](#dataset)
7. [Planificación](#planificación)
8. [Referencias](#referencias)

### Introducción

El texto (lenguaje natural) es el camino más natural para codificar el conocimiento humano. Como resultado, la mayor parte del conocimiento humano se codifica en forma de datos de texto. Por ejemplo, el conocimiento científico existe casi exclusivamente en la literatura científica, mientras que los manuales técnicos contienen explicaciones detalladas sobre cómo operar los dispositivos.
El texto es, con diferencia, el tipo de información más común con el que se encuentran las personas. De hecho, la mayor parte de la información que una persona produce y consume a diario está en forma de texto. [ChengXiang Zhai & Sean Masung 2016] 

### Motivación

En el ámbito del cuidado de la salud se genera una gran cantidad de datos valiosos para la detección y caracterización de enfermedades. El acceso a estos datos se ha facilitado gracias a la digitalización del sistema de salud, permitiendo el procesamiento de grandes cantidades de datos (big data), mediante técnicas como el procesamiento del lenguaje natural (NLP) y el aprendizaje automático (ML), que permite analizar información sobre la salud de las poblaciones y descubrir factores sociales que influyen en los resultados de salud. Sin embargo, estos textos médicos están compuestos, en gran medida, por datos relativos al paciente y requieren un manejo cuidadoso para evitar la identificación de las personas, ya que pueden aludir a información que históricamente ha dado lugar a diversos tipos de discriminación.

Esto plantea retos en cuanto a la protección de la privacidad de los datos de los pacientes, en la cual resalta el delicado equilibrio entre la utilidad de los datos y la privacidad, ya que, se pueden anonimizar datos hasta un punto en que no brinden información beneficiosa, volviéndose inútil o, por lo contrario, se puede tener datos que brinden información sustancial para ser vinculados a las personas a las que refieren. Por estas razones, investigadores de todo el mundo han desarrollado técnicas y algoritmos avanzados para anonimizar datos, permitiendo su uso para los fines solicitados mientras se mantiene el anonimato del paciente. Lo que nos lleva a la necesidad de evaluar la eficacia de estos algoritmos y técnicas de anonimización de datos, con el objetivo de identificar áreas de mejora y diseñar estrategias para optimizar su rendimiento.

### Objetivos

El presente informe tiene como objetivo principal comparar la dependencia de datos presenre en diferentes aproximaciones para la anonimización de textos médicos en español. Para ello, se establecen los siguientes objetivos específicos:
- Evaluar el rendimiento de modelos ya existentes: analizaremos la eficacia de los modelos de fine-tuning seleccionados en la tarea de anonimización de textos, enfocándonos en métricas clave como precisión, sensibilidad y F-Score.
- Identificar la dependencia de datos: buscaremos determinar cómo la calidad y la variedad de los conjuntos de datos utilizados en el entrenamiento de los modelos afectan su desempeño en contextos diferentes a los de entrenamiento, incluyendo el uso de datos sintéticos.
- Optimizar estrategias de anonimización: propondremos mejoras y estrategias basadas en los hallazgos del análisis, con el fin de optimizar el rendimiento de los modelos en la anonimización de datos sensibles.
- Fomentar la colaboración interdisciplinaria: el trabajo será en parte guiado por las críticas de nuestros compañeros, por lo que habrá un intercambio de conocimientos dado por la colaboración entre los diferentes equipos.

### Trabajos anteriores

Este trabajo se centra en la campaña “Medical Document Anonymization (MEDDOCAN)” desarrollado por el Plan de Impulso de las tecnologías del Lenguaje del Gobierno de España. En el marco de esta campaña, se creó un corpus sintético de 1000 registros clínicos con información de salud protegida (PHI, por sus siglas en inglés), denominado “Corpus MEDDOCAN”. Este corpus se distribuyó en texto plano con codificación UTF-8, donde cada caso clínico se almacena como un único archivo de texto y las anotaciones PHI se publican en formato BRAT para su visualización y evaluación.

<p align="center">
  <img width="auto" height="500" src="https://temu.bsc.es/meddocan/wp-content/uploads/2019/03/image-1-768x692.png">
  <figcaption align="center">Figura 1: Un ejemplo de anotación MEDDOCAN visualizada utilizando la interfaz de anotación BRAT.</figcaption>

</p>

Para llevar a cabo la anotación manual, se construyeron las primeras pautas públicas para PHI en español, a manera de proporcionar un conjunto claro y conciso de reglas y directrices para la identificación y clasificación de información sensible dentro de los historiales clínicos en español, asegurando la coherencia y la precisión en la anotación de PHI en todo el corpus MEDDOCAN. La guía se basó en la especificaciones del Reglamento General de Protección de Datos (GDPR) de la Unión Europea para garantizar que el corpus anotado cumpliera con las normas de privacidad de datos. Se consideraron las directrices de anotación de los proyectos de desidentificación i2b2, basados en la Ley de Portabilidad y Responsabilidad del Seguro Médico (HIPAA) de los Estados Unidos. Esto permitió la adaptación de sistemas y enfoques utilizados para textos en inglés al contexto del español. Las pautas de anotación de MEDDOCAN definieron un total de 29 tipos de entidades. La Tabla 1 resume la lista de tipos de entidades sensibles definidos para el seguimiento de MEDDOCAN y la cantidad de ocurrencias entre los conjuntos de entrenamiento, desarrollo y prueba.

<div align="center">
  <table>
    <thead>
      <tr>
        <th>Type</th>
        <th>Train</th>
        <th>Dev</th>
        <th>Test</th>
        <th>Total</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>TERRITORIO</td>
        <td>1875</td>
        <td>987</td>
        <td>956</td>
        <td>3818</td>
      </tr>
      <tr>
        <td>FECHAS</td>
        <td>1231</td>
        <td>724</td>
        <td>611</td>
        <td>2566</td>
      </tr>
      <tr>
        <td>EDAD_SUJETO_ASISTENCIA</td>
        <td>1035</td>
        <td>521</td>
        <td>518</td>
        <td>2074</td>
      </tr>
      <tr>
        <td>NOMBRE_SUJETO_ASISTENCIA</td>
        <td>1009</td>
        <td>503</td>
        <td>502</td>
        <td>2014</td>
      </tr>
      <tr>
        <td>NOMBRE_PERSONAL_SANITARIO</td>
        <td>1000</td>
        <td>497</td>
        <td>501</td>
        <td>1998</td>
      </tr>
      <tr>
        <td>SEXO_SUJETO_ASISTENCIA</td> 
        <td>925</td>
        <td>455</td>
        <td>461</td>
        <td>1841</td>
      </tr>
      <tr>
        <td>CALLE</td> 
        <td>862</td>
        <td>434</td>
        <td>413</td>
        <td>1709</td>
      </tr>
      <tr>
        <td>PAIS</td> 
        <td>713</td>
        <td>347</td>
        <td>363</td>
        <td>1423</td>
      </tr>
      <tr>
        <td>ID_SUJETO_ASISTENCIA</td> 
        <td>567</td>
        <td>292</td>
        <td>283</td>
        <td>1142</td>
      </tr>
      <tr>
        <td>CORREO_ELECTRONICO</td> 
        <td>469</td>
        <td>241</td>
        <td>249</td>
        <td>959</td>
      </tr>
      <tr>
        <td>ID_TITULACION_PERSONAL_SANITARIO</td> 
        <td>471</td>
        <td>226</td>
        <td>234</td>
        <td>931</td>
      </tr>
      <tr>
        <td>ID_ASEGURAMIENTO</td> 
        <td>391</td>
        <td>194</td>
        <td>198</td>
        <td>783</td>
      </tr>
      <tr>
        <td>HOSPITAL</td> 
        <td>255</td>
        <td>140</td>
        <td>130</td>
        <td>525</td>
      </tr>
      <tr>
        <td>FAMILIARES_SUJETO_ASISTENCIA</td> 
        <td>243</td>
        <td>92</td>
        <td>81</td>
        <td>416</td>
      </tr>
      <tr>
        <td>INSTITUCION</td> 
        <td>98</td>
        <td>72</td>
        <td>67</td>
        <td>237</td>
      </tr>
      <tr>
        <td>ID_CONTACTO_ASISTENCIAL</td> 
        <td>77</td>
        <td>32</td>
        <td>39</td>
        <td>148</td>
      </tr>
      <tr>
        <td>NUMERO_TELEFONO</td> 
        <td>58</td>
        <td>25</td>
        <td>26</td>
        <td>109</td>
      </tr>
      <tr>
        <td>PROFESION</td> 
        <td>24</td>
        <td>4</td>
        <td>9</td>
        <td>37</td>
      </tr>
      <tr>
        <td>NUMERO_FAX</td> 
        <td>15</td>
        <td>6</td>
        <td>7</td>
        <td>28</td>
      </tr>
      <tr>
        <td>OTROS_SUJETO_ASISTENCIA</td>
        <td>9</td>
        <td>6</td>
        <td>7</td>
        <td>22</td>
      </tr>
      <tr>
        <td>CENTRO_SALUD</td>
        <td>6</td>
        <td>2</td>
        <td>6</td>
        <td>14</td> 
      </tr>
      <tr>
        <td>ID_EMPLEO_PERSONAL_SANITARIO</td>
        <td>0</td>
        <td>1</td>
        <td>0</td>
        <td>1</td>
      </tr>
      <tr>
        <td>IDENTIF_VEHICULOS_NRSERIE_PLACAS</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
      <tr>
        <td>IDENTIF_DISPOSITIVOS_NRSERIE</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
      <tr>
        <td>NUMERO_BENEF_PLAN_SALUD</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
      <tr>
        <td>URL_WEB</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
      <tr>
        <td>DIREC_PROT_INTERNET</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
      <tr>
        <td>IDENTF_BIOMETRICOS</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
      <tr>
        <td>OTRO_NUMERO_IDENTIF</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
        <td>0</td>
      </tr>
    </tbody>
  </table>
  
  <figcaption>Tabla 1: Distribución del tipo de entidad entre los conjuntos de datos.</figcaption>
</div>



Durante el proyecto, se desarrollaron diversos modelos de fine-tuning, entrenados con 500 casos clínicos del corpus, seleccionados aleatoriamente. Los casos restantes se utilizaron para desarrollo y prueba. Para evaluar el rendimiento de los sistemas de desidentificación, el proyecto definió dos subtareas o categorías principales que constribuyeron al objetivo general del proyecto:

1. *NER offset and entity type classification*: Reconocimiento de entidades nombradas, posición y clasificación del tipo de entidad. Esta subtarea se centró en la identificación y clasificación precisa de la información sensible dentro de los textos médicos. Los sistemas participantes debían identificar las entidades de PHI, como nombres de pacientes, teléfonos, direcciones, etc., y clasificarlas correctamente según su tipo específico (NOMBRE_SUJETO_ASISTENCIA, ID ASEGURAMIENTO, ...).

2. *Sensitive span detection*: Detección de segmento sensible. Esta subtarea se enfocó en la detección general de tramos de texto sensibles, sin necesidad de clasificarlos por tipo. El objetivo principal era identificar y enmascarar los datos confidenciales para su divulgación segura. La subtarea 2 se dividió en dos escenarios, que reflejan diferentes niveles de granularidad en la identificación:
    * *Sensitive token detection with strict spans:* Los sistemas debían identificar los tramos de texto sensibles con precisión, sin fusionar tramos adyacentes.
    * *Sensitive token detection with merged spans:* Se permitía la fusión de tramos de PHI conectados por caracteres no alfanuméricos, lo que refleja un escenario más práctico de desidentificación. 

Cada equipo podía enviar hasta cinco archivos de predicción (ejecuciones) en un formato predefinido (BRAT o i2b2). La evaluación se basó en la precisión, la exhaustividad y la puntuación F1, junto con las puntuaciones de fuga para medir la cantidad de información sensible que no se identificó.

### Modelos a testear

De entre todos los proyectos presentados por los 18 equipos participantes en ambas categorías, seleccionamos 4 proyectos presentados por las siguientes instituciones:

1. **Bosch Center for Artificial Intelligence (Alemania)**: Centro de investigación que se especializa en el desarrollo de tecnologías de inteligencia artificial. Su enfoque abarca una variedad de aplicaciones, incluyendo el sector de la salud, donde buscan innovar en el análisis de datos y el procesamiento de lenguaje natural para mejorar la atención médica y la eficiencia de los sistemas de salud.

2. **University of Pennsylvania (Estados Unidos)**:

3. **Carlos III University of Madrid (España)**:

4. **LSI2 - Universidad Nacional de Educación a Distancia (España)**:  

### DataSet

A los datos brindados por el corpus MEDDOCAN sumaremos historias médicas desarrolladas por el generador de pacientes sintéticos open source Synthea que proporciona datos de alta calidad, realistas, pero no reales, de pacientes; en una variedad de formatos, con diferentes niveles de complejidad, cubriendo todos los aspectos de la atención médica. Los datos obtenidos no tienen coste, ni privacidad ni restricciones de seguridad.

[1301 casos clínicos oncológicos:](https://huggingface.co/datasets/PlanTL-GOB-ES/cantemist-ner) This dataset was designed for the CANcer TExt MIning Shared Task, sponsored by Plan-TL.

[3219 reportes radiológicos:](https://huggingface.co/datasets/chizhikchi/CARES) Corpus of Anonymised Radiological Evidences in Spanish

[850 reportes](https://huggingface.co/datasets/IIC/livingner1) LivingNER: Named entity recognition, normalization & classification of species, pathogens and food

[500 reportes](https://huggingface.co/datasets/bigbio/meddocan)

### Planificación
* *Semana 1*: Definir los objetivos del proyecto y el alcance del testeo de los modelos. 

* *Semana 2*: Preparar el conjunto de datos. Adaptar modelos existentes para transformar los datos de entrada y salida a los formatos requeridos por los nuevos modelos a testear.

* *Semana 3*: Configurar los modelos a testear y desarrollar el modelo de evaluación para comparar los resultados obtenidos.

* *Semana 4*: Ejecutar los modelos y realizar una evaluación preliminar de los resultados obtenidos.

* *Semana 5*: Analizar los resultados obtenidos en detalle y redactar el informe final, incluyendo conclusiones y recomendaciones para futuros trabajos.

# Preparación del conjunto de datos

- Se generaron datos sintéticos utilizando la herramienta Synthea y se transformaron a un formato compatible con el corpus de MEDDOCAN, en formato .brat.
- Disponemos actualmente del corpus de MEDDOCAN.
- Se realizó una investigación exhaustiva sobre los modelos asociados a MEDDOCAN y se identificó que algunos modelos inicialmente previstos para estudio no contaban con código accesible. Por ello, se procedió a realizar ajustes en los objetivos iniciales del proyecto.
- Actualmente, estamos considerando los siguientes modelos:
  - [lukas.lange](https://huggingface.co/llange/xlm-roberta-large-spanish-clinical)
  - [gauku](https://github.com/ionur/MEDDOCAN-Medical-Document-Anonymization)
  - [ccolon](https://github.com/ccolonruiz/MEDDOCAN)
  - [lsi2](https://github.com/alicialara/lsi2_uned_at_MEDDOCAN2019)

### Referencias
- ChengXiang Zhai & Sean Masung (2016). Text Data Management and Analysis: A Practical Introduction to Information Retrieval and Text Mining. ACM Books.
- [Automatic De-Identification of Medical Texts in Spanish: the MEDDOCAN Track, Corpus, Guidelines, Methods and Evaluation of Results.](https://ceur-ws.org/Vol-2421/MEDDOCAN_overview.pdf)
- [Synthea: Synthetic Health Data Engine.](https://github.com/synthetichealth/synthea/wiki)
- [BRAT: a Web-based Tool for NLP-Assisted Text Annotation.](https://aclanthology.org/E12-2021.pdf)
- [Guías de anotación de información de salud protegida.](https://temu.bsc.es/meddocan/wp-content/uploads/2019/02/gu%C3%ADas-de-anotaci%C3%B3n-de-informaci%C3%B3n-de-salud-protegida.pdf)
- [GDPR: General Data Protection Regulation.](https://gdpr.eu/)
- [HIPAA: Health Insurance Portability and Accountability Act.](https://www.hhs.gov/hipaa/index.html)
- [i2b2: Informatics for Integrating Biology and the Bedside.](https://www.i2b2.org/)

