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
- **Evaluar el rendimiento de modelos ya existentes**: los datos usados para el desarrollo, entrenamiento y validación de los modelos de anonimización de textos médicos en español, provienen de una fuente específica, desarrollados por un equipo en particular, con una estructura y contenido bien definidos, lo que puede afectar su rendimiento en contextos diferentes de la vida real, porque un problema común en el aprendizaje automático es que los modelos no generalizan bien a datos que no han visto antes (overfitting). Por ello, buscamos conocer su verdadera eficacia en la tarea de anonimización de textos, donde la notación de textos médicos pueden variar en su estructura y contenido, analizando la eficacia de los modelos seleccionados, con datos que no fueron utilizados en su entrenamiento y creados a partir de diferentes fuentes.
- **Identificar la dependencia de datos y optimizar estrategias de anonimización**: buscaremos determinar en qué puntos los modelos son dependientes de los datos de entrenamiento y cómo esta dependencia afecta su rendimiento. Propondremos mejoras y estrategias basadas en los hallazgos del análisis, con el fin de optimizar el rendimiento de los modelos en la anonimización de datos sensibles.
- **Fomentar la colaboración interdisciplinaria**: el trabajo será en parte guiado por las críticas de nuestros compañeros, por lo que habrá un intercambio de conocimientos dado por la colaboración entre los diferentes equipos.

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
      <tr> <th>Type</th>                             <th>Train</th><th>Dev</th><th>Test</th><th>Total</th> </tr>
    </thead>
    <tbody>
      <tr> <td>TERRITORIO</td>                       <td>1875</td> <td>987</td> <td>956</td> <td>3818</td> </tr>
      <tr> <td>FECHAS</td>                           <td>1231</td> <td>724</td> <td>611</td> <td>2566</td> </tr>
      <tr> <td>EDAD_SUJETO_ASISTENCIA</td>           <td>1035</td> <td>521</td> <td>518</td> <td>2074</td> </tr>
      <tr> <td>NOMBRE_SUJETO_ASISTENCIA</td>         <td>1009</td> <td>503</td> <td>502</td> <td>2014</td> </tr>
      <tr> <td>NOMBRE_PERSONAL_SANITARIO</td>        <td>1000</td> <td>497</td> <td>501</td> <td>1998</td> </tr>
      <tr> <td>SEXO_SUJETO_ASISTENCIA</td>           <td>925</td>  <td>455</td> <td>461</td> <td>1841</td> </tr>
      <tr> <td>CALLE</td>                            <td>862</td>  <td>434</td> <td>413</td> <td>1709</td> </tr>
      <tr> <td>PAIS</td>                             <td>713</td>  <td>347</td> <td>363</td> <td>1423</td> </tr>
      <tr> <td>ID_SUJETO_ASISTENCIA</td>             <td>567</td>  <td>292</td> <td>283</td> <td>1142</td> </tr>
      <tr> <td>CORREO_ELECTRONICO</td>               <td>469</td>  <td>241</td> <td>249</td> <td>959</td>  </tr>
      <tr> <td>ID_TITULACION_PERSONAL_SANITARIO</td> <td>471</td>  <td>226</td> <td>234</td> <td>931</td>  </tr>
      <tr> <td>ID_ASEGURAMIENTO</td>                 <td>391</td>  <td>194</td> <td>198</td> <td>783</td>  </tr>
      <tr> <td>HOSPITAL</td>                         <td>255</td>  <td>140</td> <td>130</td> <td>525</td>  </tr>
      <tr> <td>FAMILIARES_SUJETO_ASISTENCIA</td>     <td>243</td>  <td>92</td>  <td>81</td>  <td>416</td>  </tr>
      <tr> <td>INSTITUCION</td>                      <td>98</td>   <td>72</td>  <td>67</td>  <td>237</td>  </tr>
      <tr> <td>ID_CONTACTO_ASISTENCIAL</td>          <td>77</td>   <td>32</td>  <td>39</td>  <td>148</td>  </tr>
      <tr> <td>NUMERO_TELEFONO</td>                  <td>58</td>   <td>25</td>  <td>26</td>  <td>109</td>  </tr>
      <tr> <td>PROFESION</td>                        <td>24</td>   <td>4</td>   <td>9</td>   <td>37</td>   </tr>
      <tr> <td>NUMERO_FAX</td>                       <td>15</td>   <td>6</td>   <td>7</td>   <td>28</td>   </tr>
      <tr> <td>OTROS_SUJETO_ASISTENCIA</td>          <td>9</td>    <td>6</td>   <td>7</td>   <td>22</td>   </tr>
      <tr> <td>CENTRO_SALUD</td>                     <td>6</td>    <td>2</td>   <td>6</td>   <td>14</td>   </tr>
      <tr> <td>ID_EMPLEO_PERSONAL_SANITARIO</td>     <td>0</td>    <td>1</td>   <td>0</td>   <td>1</td>    </tr>
      <tr> <td>IDENTIF_VEHICULOS_NRSERIE_PLACAS</td> <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
      <tr> <td>IDENTIF_DISPOSITIVOS_NRSERIE</td>     <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
      <tr> <td>NUMERO_BENEF_PLAN_SALUD</td>          <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
      <tr> <td>URL_WEB</td>                          <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
      <tr> <td>DIREC_PROT_INTERNET</td>              <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
      <tr> <td>IDENTF_BIOMETRICOS</td>               <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
      <tr> <td>OTRO_NUMERO_IDENTIF</td>              <td>0</td>    <td>0</td>   <td>0</td>   <td>0</td>    </tr>
    </tbody>
  </table>
  <figcaption>Tabla 1: Distribución del tipo de entidad entre los conjuntos de datos.</figcaption>
</div>

Durante el proyecto se desarrollaron diversos modelos entrenados con 500 casos clínicos del corpus, seleccionados aleatoriamente, siendo los casos restantes utilizados para desarrollo y prueba. Para evaluar el rendimiento de los sistemas de desidentificación, el proyecto definió dos subtareas o categorías principales que constribuyeron al objetivo general del proyecto:

1. *NER offset and entity type classification*: Reconocimiento de entidades nombradas, posición y clasificación del tipo de entidad. Esta subtarea se centró en la identificación y clasificación precisa de la información sensible dentro de los textos médicos. Los sistemas participantes debían identificar las entidades de PHI, como nombres de pacientes, teléfonos, direcciones, etc., y clasificarlas correctamente según su tipo específico (NOMBRE_SUJETO_ASISTENCIA, ID ASEGURAMIENTO, ...).

2. *Sensitive span detection*: Detección de segmento sensible. Esta subtarea se enfocó en la detección general de tramos de texto sensibles, sin necesidad de clasificarlos por tipo. El objetivo principal era identificar y enmascarar los datos confidenciales para su divulgación segura. La subtarea 2 se dividió en dos escenarios, que reflejan diferentes niveles de granularidad en la identificación:
    * *Sensitive token detection with strict spans:*  Detección de tokens sensibles con tramos estrictos. Los sistemas debían identificar los tramos de texto sensibles con precisión, sin fusionar tramos adyacentes.
    * *Sensitive token detection with merged spans:* Detección de tokens sensibles con tramos fusionados. Se permitía la fusión de tramos de PHI conectados por caracteres no alfanuméricos, lo que refleja un escenario más práctico de desidentificación. 

Cada equipo podía enviar hasta cinco archivos de predicción (ejecuciones) en un formato predefinido (BRAT o i2b2). La evaluación se basó en la precisión, la exhaustividad y la puntuación F1, junto con las puntuaciones de fuga para medir la cantidad de información sensible que no se identificó.

### Modelos a testear

De entre todos los proyectos presentados por los 18 equipos participantes en ambas categorías, seleccionamos 4 proyectos presentados por las siguientes instituciones:

- **CLIN-X (Bosch Center for Artificial Intelligence - Spoken Language Systems, Saarland University, Alemania)**: El modelo CLIN-X es un modelo preentrenado basado en el transformador multilingüe XLM-R, entrenado en 100 idiomas. Para adaptarlo al dominio clínico español, se entrenó con un corpus de 790MB de documentos del archivo Scielo y los recursos MeSpEn, utilizando el objetivo de modelado de lenguaje enmascarado (MLM). Aunque el modelo se adaptó al español, sigue siendo multilingüe y se puede aplicar a tareas en inglés.

Notar que el comando para ejecutar train_our_model_architecture.py (python3 train_our_model_architecture.py --data_path split_files/ --train_files random_split_1.txt,random_split_2.txt,random_split_3.txt,random_split_4.txt --dev_file random_split_5.txt --model xlm-roberta-large-spanish-clinical --name model_name --storage_path models) tiene un error: hay que reemplazar --model xlm-roberta-large-spanish-clinical por --model llange/xlm-roberta-large-spanish-clinical. 
Para usar la CPU (y no GPU de NVIDIA) hay que hacer ciertas modificaciones en el archivo create_data_splits.py:
model = AutoModel.from_pretrained(args.model_path).cuda() -> model = AutoModel.from_pretrained(args.model_path)
input_ids = torch.stack(input_ids).long().cuda() -> input_ids = torch.stack(input_ids).long()
doc_vec = torch.zeros(1024).cuda() -> doc_vec = torch.zeros(1024)

- **BiLSTM-CRF (University of Pennsylvania, Estados Unidos)**: Este modelo es una red neuronal que combina Bi-LSTM para capturar el contexto y CRF para la decodificación de etiquetas. Bi-LSTM (memoria bidireccional a largo y corto plazo) es un tipo de red neuronal recurrente (RNN) que procesa datos secuenciales tanto en dirección directa como inversa. Combina la potencia de LSTM con el procesamiento bidireccional, lo que permite que el modelo capture el contexto pasado y futuro de la secuencia de entrada. Por otro lado, CRF (campo aleatorio condicional) es un modelo estocástico utilizado para modelar secuencias de etiquetas en problemas de aprendizaje supervisado.

- **NeuroNer (Carlos III University of Madrid, España)**: Este modelo se basa en la herramienta NeuroNer, que utiliza una arquitectura de aprendizaje profundo (deep learning) que combina dos capas bidireccionales de memoria a corto plazo (BiLSTM) con una capa final de campos aleatorios condicionales (CRF). Se utiliza una segunda capa BiLSTM para obtener la secuencia de probabilidades para cada token de pertenecer a una etiqueta específica. La capa CRF se utiliza para modelar las dependencias entre las etiquetas y obtener la secuencia de etiquetas más probable.
Para poder correr el modelo es necesario modificar la segunda linea de código del archivo Train.ipynb:
from keras.backend.tensorflow_backend import set_session -> from tensorflow.compat.v1.keras.backend import set_session

- **Few-shot NER (Universidad Nacional de Educación a Distancia, España)**: Este modelo preentrenado se basa en la implementación de Keras de "Few-shot Learning for Named Entity Recognition in Medical Text", un modelo LSTM y CNN bidireccional similar a Chiu y Nichols (2016) para los datos de noticias de CoNLL 2003. Es un híbrido Bi-LSTM y CNN que combina la capacidad de Bi-LSTM para capturar el contexto y la capacidad de la red neuronal convolucional (CNN) para capturar características locales. La CNN utiliza filtros convolucionales para extraer características locales de los datos de entrada, permitiendo que el modelo aprenda patrones espaciales en los datos.

### DataSet

Para llevar a cabo la comparación de los modelos seleccionados, se utilizarán los siguientes conjuntos de datos:

- **CANTEMIST NER**: Colección de 1301 informes clínicos oncológicos escritos en español, con menciones de morfología tumoral anotadas manualmente y mapeadas por expertos clínicos a una terminología controlada. Este conjunto de datos fue diseñado para la tarea compartida de minería de texto sobre cáncer, patrocinada por Plan-TL.

- **CARES**: Corpus de Anonimización de Evidencias Radiológicas en Español, que contiene 3219 informes radiológicos en español, anotados con entidades médicas y personales de manera manual y revisados por radiólogos.

- **LivingNER**: Conjunto de datos de entidades nombradas, normalización y clasificación de especies, patógenos y alimentos, que contiene 850 informes. Este conjunto de datos fue diseñado para la tarea compartida de minería de texto sobre reconocimiento de entidades nombradas, normalización y clasificación de especies, patógenos y alimentos, patrocinada por Plan-TL.

- **Synthea**: Conjunto de datos sintéticos generados por la herramienta Synthea, que simula datos de pacientes siguiendo un modelo de datos realista y estructurado.

### Preparación del conjunto de datos

### Referencias
- ChengXiang Zhai & Sean Masung (2016). Text Data Management and Analysis: A Practical Introduction to Information Retrieval and Text Mining. ACM Books.
- [Automatic De-Identification of Medical Texts in Spanish: the MEDDOCAN Track, Corpus, Guidelines, Methods and Evaluation of Results.](https://ceur-ws.org/Vol-2421/MEDDOCAN_overview.pdf)
- [Synthea: Synthetic Health Data Engine.](https://github.com/synthetichealth/synthea/wiki)
- [BRAT: a Web-based Tool for NLP-Assisted Text Annotation.](https://aclanthology.org/E12-2021.pdf)
- [Guías de anotación de información de salud protegida.](https://temu.bsc.es/meddocan/wp-content/uploads/2019/02/gu%C3%ADas-de-anotaci%C3%B3n-de-informaci%C3%B3n-de-salud-protegida.pdf)
- [GDPR: General Data Protection Regulation.](https://gdpr.eu/)
- [HIPAA: Health Insurance Portability and Accountability Act.](https://www.hhs.gov/hipaa/index.html)
- [i2b2: Informatics for Integrating Biology and the Bedside.](https://www.i2b2.org/)
- [Bhavna Saluja, Gaurav Kumar, João Sedoc, Chris Callison-Burch. Anonymization of Sensitive Information in Medical Health Records.](https://ceur-ws.org/Vol-2421/MEDDOCAN_paper_2.pdf)
- [Lukas Lange, Heike Adel, Jannik Strötgen. NLNDE: The Neither-Language-Nor-Domain-Experts’ Way of Spanish Medical Document De-Identification.](https://ceur-ws.org/Vol-2421/MEDDOCAN_paper_5.pdf)
- [Cristóbal Colón-Ruiz, Isabel Segura-Dedmar. Protected Health Information Recognition by BiLSTM-CRF.](https://ceur-ws.org/Vol-2421/MEDDOCAN_paper_6.pdf)
- [Alicia Lara-Clares, Ana Garcia-Serrano. Key Phrases Annotation in Medical Documents: MEDDOCAN 2019 Anonymization Task.](https://ceur-ws.org/Vol-2421/MEDDOCAN_paper_15.pdf)

