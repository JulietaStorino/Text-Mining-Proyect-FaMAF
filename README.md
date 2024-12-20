# Trabajo de presentación para el curso "Minería de Datos para Textos" (Text Mining)
## Impacto de la Dependencia de Datos en Modelos de Fine-Tuning para la anonimización de textos de ciencias de la vida
**Profesores**: Alemany Laura Alonso, Bordone Carranza Matias Eduardo, Teruel Milagro 

**Alumnos**: Brunello Florencia, Storino Julieta Paola, Troiano Santiago

### Índice
1. [Introducción](#introducción)
2. [Motivación](#motivación)
3. [Objetivos](#objetivos)
4. [Trabajos Previos: Campaña MEDDOCAN](#trabajos-previos-campaña-meddocan)
6. [Modelos a evaluar](#modelos-a-evaluar)
7. [Conjunto de datos](#conjunto-de-datos)
8. [Ejecución de los modelos](#ejecución-de-los-modelos)
9. [Análisis y comparación de los resultados](#análisis-y-comparación-de-los-resultados)
10. [Conclusiones](#conclusiones)
11. [Referencias](#referencias)

### Introducción

El texto (lenguaje natural) es el camino más natural para codificar el conocimiento humano. Como resultado, la mayor parte del conocimiento humano se codifica en forma de datos de texto. Por ejemplo, el conocimiento científico existe casi exclusivamente en la literatura científica, mientras que los manuales técnicos contienen explicaciones detalladas sobre cómo operar los dispositivos. El texto es, con diferencia, el tipo de información más común con el que se encuentran las personas. De hecho, la mayor parte de la información que una persona produce y consume a diario está en forma de texto. [ChengXiang Zhai & Sean Masung 2016] 

### Motivación

En el ámbito del cuidado de la salud se genera una gran cantidad de datos valiosos para la detección y caracterización de enfermedades. El acceso a estos datos se ha facilitado gracias a la digitalización del sistema de salud, permitiendo el procesamiento de grandes cantidades de datos (big data), mediante técnicas como el procesamiento del lenguaje natural (NLP) y el aprendizaje automático (ML), que permite analizar información sobre la salud de las poblaciones y descubrir factores sociales que influyen en los resultados de salud. Sin embargo, estos textos médicos están compuestos, en gran medida, por datos relativos al paciente y requieren un manejo cuidadoso para evitar la identificación de las personas, ya que pueden aludir a información que históricamente ha dado lugar a diversos tipos de discriminación.

Esto plantea retos en cuanto a la protección de la privacidad de los datos de los pacientes, en la cual resalta el delicado equilibrio entre la utilidad de los datos y la privacidad, ya que, se pueden anonimizar datos hasta un punto en que no brinden información beneficiosa, volviéndose inútil o, por lo contrario, se puede tener datos que brinden información sustancial para ser vinculados a las personas a las que refieren. Por estas razones, investigadores de todo el mundo han desarrollado técnicas y algoritmos avanzados para anonimizar datos, permitiendo su uso para los fines solicitados mientras se mantiene el anonimato del paciente. Lo que nos lleva a la necesidad de evaluar la eficacia de estos algoritmos y técnicas de anonimización de datos, con el objetivo de identificar áreas de mejora y diseñar estrategias para optimizar su rendimiento.

### Objetivos

El presente informe tiene como objetivo principal comparar la dependencia de datos presente en diferentes aproximaciones para la anonimización de textos médicos en español. Para ello, se establecen los siguientes objetivos específicos:
- **Evaluar el rendimiento de modelos ya existentes**: los datos usados para el desarrollo, entrenamiento y validación de los modelos de anonimización de textos médicos en español provienen de una fuente específica, desarrollados por un equipo en particular, con una estructura y contenido bien definidos, lo que puede afectar su rendimiento en contextos diferentes de la vida real, porque un problema común en el aprendizaje automático es que los modelos no generalizan bien a datos que no han visto antes (overfitting). Por ello, buscamos conocer su verdadera eficacia en la tarea de anonimización de textos, donde la notación de textos médicos pueden variar en su estructura y contenido, analizando la eficacia de los modelos seleccionados con datos que no fueron utilizados en su entrenamiento.
- **Comparar el rendimiento de los modelos**: buscaremos determinar en qué aspectos los modelos se destacan o presentan deficiencias al realizar la tarea de anonimización en el conjunto de datos.
- **Fomentar la colaboración interdisciplinaria**: el trabajo será en parte guiado por las críticas de nuestros compañeros, por lo que habrá un intercambio de conocimientos dado por la colaboración entre los diferentes equipos.

### Trabajos Previos: Campaña MEDDOCAN

Este trabajo se centra en la campaña "Medical Document Anonymization (MEDDOCAN)" desarrollada por el Plan de Impulso de las Tecnologías del Lenguaje del Gobierno de España. Donde se propuso la tarea de anonimización de textos médicos en español, con el objetivo de desarrollar y evaluar la eficacia de modelos de aprendizaje automático para la identificación y desidentificación de información de salud protegida (PHI) en textos médicos en español. La campaña se llevó a cabo en el marco de la conferencia internacional IberLEF 2021, que se celebró en línea en septiembre de 2021.

#### Corpus MEDDOCAN

En el marco de esta campaña, se creó un corpus sintético de 1000 registros clínicos con información de salud protegida (PHI, por sus siglas en inglés), denominado [Corpus MEDDOCAN](./SPACCC_MEDDOCAN/corpus/). Este corpus se distribuyó en texto plano con codificación UTF-8, donde cada caso clínico se almacena como un único archivo de texto y las anotaciones PHI se publican en formato BRAT e i2b2 para su visualización y evaluación.

```plaintext
Datos del paciente.
Nombre:  Pedro.
Apellidos: De Miguel Rivera.
NHC: 2569870.
Domicilio: Calle Carmen Romero, 23, 1D.
Localidad/ Provincia: Madrid.
CP: 28035.
Datos asistenciales.
Fecha de nacimiento: 10/10/1963.
País: España.
Edad: 53 años Sexo: H.
Fecha de Ingreso: 17/06/2016.
Médico: Estefanía Romero Selas  NºCol: 28 28 20943.
Informe clínico del paciente: varón de 53 años sin antecedentes de interés que ingresa procedente de urgencias con un cuadro de tromboembolismo pulmonar.
Ante la sospecha de neoplasia oculta y la presencia de hematuria no evidenciada previamente se procede a la realización de Ecografía abdominal donde se evidencia una masa renal derecha y se completa el estudio con TAC y RM..
Ambos estudios confirman la presencia de una tumoración heterogénea que infiltra los dos tercios inferiores del riñón derecho de aproximadamente 10x10 cms. con afectación del seno e hilio renal objetivándose también trombosis tumoral de la vena renal derecha y cava infrahepática. No se evidenciaban adenopatías ni metástasis.
Es intervenido quirúrgicamente realizándosele por vía anterior, una nefrectomía radical con cavotomía para la exéresis del trombo y una extensa linfadenectomía aorto-cava derecha.
El resultado anatomo-patológico fue de carcinoma de células renales claras grado 2 de Fuhrman de 9 cm. con invasión de hilio renal, grasa perinéfrica y vena renal, sin afectación metastásica de ganglios ni de los bordes de dicha grasa ni del hilio, así como uréter libres. (Estadio III, T3N0M0). El paciente fue dado de alta hospitalaria al sexto día.
A los 3 meses de la intervención el paciente refiere leve dolor e induración en el pene de reciente aparición. A la palpación se objetiva una masa indurada.
Se le realiza RM de pelvis que nos informa de la existencia de una masa que ocupa y expande el cuerpo cavernoso izquierdo, compatible con metástasis de carcinoma renal previo..
Se toma biopsia de dicha lesión, cuyo resultado nos confirma la sospecha evidenciándose en los cortes histológicos nidos aislados de células tumorales compatibles con metástasis de carcinoma de células claras.
Ante este diagnóstico, nos planteamos cuál sería la mejor actitud terapéutica para el paciente, y tuvimos en cuenta para ello, el incremento progresivo del dolor local, la edad del paciente y su buen estado general. Por ello, se optó por una penectomía total hasta confirmar intraoperatoriamente un borde quirúrgico libre de enfermedad.
Una semana después del alta ingresa en el servicio de Oncología con un cuadro de obnubilación y alteración motora y sensitiva, y presenta en el TC craneal lesiones en cerebelo y hemisferio cerebral derecho compatible con metástasis. Se realiza un TC torácico y aparecen también múltiples nódulos pulmonares y microadenopatías paratraqueales bilaterales en relación con metástasis.
El paciente fallece a los nueve meses de la primera intervención de su carcinoma renal, es decir seis meses después del diagnóstico de las metástasis en pene.
Remitido por: Dra. Estefanía Romero Selas. Email: eromeroselas@yahoo.es
```
<div align="center">
  <figcaption> Figura 1: Ejemplo de texto médico del corpus MEDDOCAN.</figcaption>
</div>

``` ann
T1	CORREO_ELECTRONICO 3042 3063	eromeroselas@yahoo.es
T2	SEXO_SUJETO_ASISTENCIA 363 368	varón
T3	EDAD_SUJETO_ASISTENCIA 372 379	53 años
T4	ID_TITULACION_PERSONAL_SANITARIO 320 331	28 28 20943
T5	NOMBRE_PERSONAL_SANITARIO 3011 3033	Estefanía Romero Selas
T6	NOMBRE_PERSONAL_SANITARIO 289 311	Estefanía Romero Selas
T7	FECHAS 269 279	17/06/2016
T8	PAIS 220 226	España
T9	SEXO_SUJETO_ASISTENCIA 248 249	H
T10	EDAD_SUJETO_ASISTENCIA 234 241	53 años
T11	FECHAS 202 212	10/10/1963
T12	TERRITORIO 153 158	28035
T13	TERRITORIO 141 147	Madrid
T14	CALLE 90 117	Calle Carmen Romero, 23, 1D
T15	ID_SUJETO_ASISTENCIA 70 77	2569870
T16	NOMBRE_SUJETO_ASISTENCIA 47 63	De Miguel Rivera
T17	NOMBRE_SUJETO_ASISTENCIA 29 34	Pedro
```
<div align="center">
  <figcaption> Figura 1: Ejemplo de notación de entidades de un texto médico del corpus MEDDOCAN en BRAT (.ann) </figcaption>
</div>

``` xml
<?xml version='1.0' encoding='UTF-8'?>
<MEDDOCAN>
  <TEXT><![CDATA[Datos del paciente.
Nombre:  Pedro.
Apellidos: De Miguel Rivera.
NHC: 2569870.
Domicilio: Calle Carmen Romero, 23, 1D.
Localidad/ Provincia: Madrid.
CP: 28035.
Datos asistenciales.
Fecha de nacimiento: 10/10/1963.
País: España.
Edad: 53 años Sexo: H.
Fecha de Ingreso: 17/06/2016.
Médico: Estefanía Romero Selas  NºCol: 28 28 20943.
Informe clínico del paciente: varón de 53 años sin antecedentes de interés que ingresa procedente de urgencias con un cuadro de tromboembolismo pulmonar.
Ante la sospecha de neoplasia oculta y la presencia de hematuria no evidenciada previamente se procede a la realización de Ecografía abdominal donde se evidencia una masa renal derecha y se completa el estudio con TAC y RM..
Ambos estudios confirman la presencia de una tumoración heterogénea que infiltra los dos tercios inferiores del riñón derecho de aproximadamente 10x10 cms. con afectación del seno e hilio renal objetivándose también trombosis tumoral de la vena renal derecha y cava infrahepática. No se evidenciaban adenopatías ni metástasis.
Es intervenido quirúrgicamente realizándosele por vía anterior, una nefrectomía radical con cavotomía para la exéresis del trombo y una extensa linfadenectomía aorto-cava derecha.
El resultado anatomo-patológico fue de carcinoma de células renales claras grado 2 de Fuhrman de 9 cm. con invasión de hilio renal, grasa perinéfrica y vena renal, sin afectación metastásica de ganglios ni de los bordes de dicha grasa ni del hilio, así como uréter libres. (Estadio III, T3N0M0). El paciente fue dado de alta hospitalaria al sexto día.
A los 3 meses de la intervención el paciente refiere leve dolor e induración en el pene de reciente aparición. A la palpación se objetiva una masa indurada.
Se le realiza RM de pelvis que nos informa de la existencia de una masa que ocupa y expande el cuerpo cavernoso izquierdo, compatible con metástasis de carcinoma renal previo..
Se toma biopsia de dicha lesión, cuyo resultado nos confirma la sospecha evidenciándose en los cortes histológicos nidos aislados de células tumorales compatibles con metástasis de carcinoma de células claras.
Ante este diagnóstico, nos planteamos cuál sería la mejor actitud terapéutica para el paciente, y tuvimos en cuenta para ello, el incremento progresivo del dolor local, la edad del paciente y su buen estado general. Por ello, se optó por una penectomía total hasta confirmar intraoperatoriamente un borde quirúrgico libre de enfermedad.
Una semana después del alta ingresa en el servicio de Oncología con un cuadro de obnubilación y alteración motora y sensitiva, y presenta en el TC craneal lesiones en cerebelo y hemisferio cerebral derecho compatible con metástasis. Se realiza un TC torácico y aparecen también múltiples nódulos pulmonares y microadenopatías paratraqueales bilaterales en relación con metástasis.
El paciente fallece a los nueve meses de la primera intervención de su carcinoma renal, es decir seis meses después del diagnóstico de las metástasis en pene.
Remitido por: Dra. Estefanía Romero Selas. Email: eromeroselas@yahoo.es
]]></TEXT>
  <TAGS>
    <NAME id="T17" start="29" end="34" text="Pedro" TYPE="NOMBRE_SUJETO_ASISTENCIA" comment=""/>
    <NAME id="T16" start="47" end="63" text="De Miguel Rivera" TYPE="NOMBRE_SUJETO_ASISTENCIA" comment=""/>
    <ID id="T15" start="70" end="77" text="2569870" TYPE="ID_SUJETO_ASISTENCIA" comment=""/>
    <LOCATION id="T14" start="90" end="117" text="Calle Carmen Romero, 23, 1D" TYPE="CALLE" comment=""/>
    <LOCATION id="T13" start="141" end="147" text="Madrid" TYPE="TERRITORIO" comment=""/>
    <LOCATION id="T12" start="153" end="158" text="28035" TYPE="TERRITORIO" comment=""/>
    <DATE id="T11" start="202" end="212" text="10/10/1963" TYPE="FECHAS" comment=""/>
    <LOCATION id="T8" start="220" end="226" text="España" TYPE="PAIS" comment=""/>
    <AGE id="T10" start="234" end="241" text="53 años" TYPE="EDAD_SUJETO_ASISTENCIA" comment=""/>
    <OTHER id="T9" start="248" end="249" text="H" TYPE="SEXO_SUJETO_ASISTENCIA" comment=""/>
    <DATE id="T7" start="269" end="279" text="17/06/2016" TYPE="FECHAS" comment=""/>
    <NAME id="T6" start="289" end="311" text="Estefanía Romero Selas" TYPE="NOMBRE_PERSONAL_SANITARIO" comment=""/>
    <ID id="T4" start="320" end="331" text="28 28 20943" TYPE="ID_TITULACION_PERSONAL_SANITARIO" comment=""/>
    <OTHER id="T2" start="363" end="368" text="varón" TYPE="SEXO_SUJETO_ASISTENCIA" comment=""/>
    <AGE id="T3" start="372" end="379" text="53 años" TYPE="EDAD_SUJETO_ASISTENCIA" comment=""/>
    <NAME id="T5" start="3011" end="3033" text="Estefanía Romero Selas" TYPE="NOMBRE_PERSONAL_SANITARIO" comment=""/>
    <CONTACT id="T1" start="3042" end="3063" text="eromeroselas@yahoo.es" TYPE="CORREO_ELECTRONICO" comment=""/>
  </TAGS>
</MEDDOCAN>
```
<div align="center">
  <figcaption> Figura 1: Ejemplo de notación de entidades de un texto médico del corpus MEDDOCAN en i2b2 (.xml) </figcaption>
</div>

#### Guía de anotación de entidades sensibles

Como parte de la campaña MEDDOCAN y para llevar a cabo la anotación manual, se construyeron las primeras [pautas públicas para PHI en español](./SPACCC_MEDDOCAN/guidelines/guías-de-anotación-de-información-de-salud-protegida.pdf), a manera de proporcionar un conjunto claro y conciso de reglas y directrices para la identificación y clasificación de información sensible dentro de los historiales clínicos en español, asegurando la coherencia y la precisión en la anotación de PHI en todo el corpus MEDDOCAN. La guía se basó en la especificaciones del Reglamento General de Protección de Datos (GDPR) de la Unión Europea para garantizar que el corpus anotado cumpliera con las normas de privacidad de datos. Se consideraron las directrices de anotación de los proyectos de desidentificación i2b2, basados en la Ley de Portabilidad y Responsabilidad del Seguro Médico (HIPAA) de los Estados Unidos. Esto permitió la adaptación de sistemas y enfoques utilizados para textos en inglés al contexto del español. Las pautas de anotación de MEDDOCAN definieron un total de 29 tipos de entidades.

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

#### Evaluación de los modelos 

Durante el proyecto se desarrollaron diversos modelos entrenados con 500 casos clínicos del corpus, seleccionados aleatoriamente, siendo los casos restantes utilizados para desarrollo y prueba. Para evaluar el rendimiento de los sistemas de desidentificación, el proyecto definió dos subtareas o categorías principales que constribuyeron al objetivo general del proyecto:

1. NER offset and entity type classification: reconocimiento de entidades nombradas, posición y clasificación del tipo de entidad. Esta subtarea se centró en la identificación y clasificación precisa de la información sensible dentro de los textos médicos. Los sistemas participantes debían identificar las entidades de PHI, como nombres de pacientes, teléfonos, direcciones, etc., y clasificarlas correctamente según su tipo específico (NOMBRE_SUJETO_ASISTENCIA, ID ASEGURAMIENTO, ...).

2. Sensitive span detection: detección de segmento sensible. Esta subtarea se enfocó en la detección general de tramos de texto sensibles, sin necesidad de clasificarlos por tipo. El objetivo principal era identificar y enmascarar los datos confidenciales para su divulgación segura. La subtarea 2 se dividió en dos escenarios, que reflejan diferentes niveles de granularidad en la identificación:
    * Sensitive token detection with strict spans:  Detección de tokens sensibles con tramos estrictos. Los sistemas debían identificar los tramos de texto sensibles con precisión, sin fusionar tramos adyacentes.
    * Sensitive token detection with merged spans: Detección de tokens sensibles con tramos fusionados. Se permitía la fusión de tramos de PHI conectados por caracteres no alfanuméricos, lo que refleja un escenario más práctico de desidentificación. 

Cada equipo podía enviar hasta cinco archivos de predicción (ejecuciones) en un formato predefinido (BRAT o i2b2). La evaluación se basó en la precisión, la exhaustividad y la puntuación F1, junto con las puntuaciones de fuga para medir la cantidad de información sensible que no se identificó.

### Modelos a evaluar

De entre todos los proyectos presentados por los 18 equipos participantes en ambas categorías, seleccionamos 3 proyectos de código abierto que presentaron un rendimiento destacado en la tarea de anonimización de textos médicos en español. Los modelos seleccionados son los siguientes:

- **CLIN-X (Bosch Center for Artificial Intelligence - Spoken Language Systems, Saarland University, Alemania)**: El modelo CLIN-X es un modelo preentrenado basado en el transformador multilingüe XLM-R, entrenado en 100 idiomas. Para adaptarlo al dominio clínico español, se entrenó con un corpus de 790MB de documentos del archivo Scielo y los recursos MeSpEn, utilizando el objetivo de modelado de lenguaje enmascarado (MLM). Aunque el modelo se adaptó al español, sigue siendo multilingüe y se puede aplicar a tareas en inglés.

- **BiLSTM-CRF (University of Pennsylvania, Estados Unidos)**: Este modelo es una red neuronal que combina Bi-LSTM para capturar el contexto y CRF para la decodificación de etiquetas. Bi-LSTM (memoria bidireccional a largo y corto plazo) es un tipo de red neuronal recurrente (RNN) que procesa datos secuenciales tanto en dirección directa como inversa. Combina la potencia de LSTM con el procesamiento bidireccional, lo que permite que el modelo capture el contexto pasado y futuro de la secuencia de entrada. Por otro lado, CRF (campo aleatorio condicional) es un modelo estocástico utilizado para modelar secuencias de etiquetas en problemas de aprendizaje supervisado.

- **NeuroNer (Carlos III University of Madrid, España)**: Este modelo se basa en la herramienta NeuroNer, que utiliza una arquitectura de aprendizaje profundo (deep learning) que combina dos capas bidireccionales de memoria a corto plazo (BiLSTM) con una capa final de campos aleatorios condicionales (CRF). Se utiliza una segunda capa BiLSTM para obtener la secuencia de probabilidades para cada token de pertenecer a una etiqueta específica. La capa CRF se utiliza para modelar las dependencias entre las etiquetas y obtener la secuencia de etiquetas más probable.
Para poder correr el modelo es necesario modificar la segunda linea de código del archivo Train.ipynb:
from keras.backend.tensorflow_backend import set_session -> from tensorflow.compat.v1.keras.backend import set_session

### Conjunto de datos

Para llevar a cabo la comparación de los modelos seleccionados se generaron nuevos [conjuntos de datos de textos](./data/) médicos en español, obtenidos a partir del módulo [Synthetic Patient Generator](./Synthetic-Patient-Generation/). Este módulo permite la generación de datos sintéticos de pacientes, incluyendo información de salud protegida (PHI), como nombres, direcciones, fechas de nacimiento, números de teléfono, entre otros. Los datos generados se almacenan en formato txt, ann y xml.

Luego, para generar el informe que contiene los datos sintéticos de los pacientes, se utilizó GPT-4o, un transformador generativo multimodal y multilingüe preentrenado, diseñado por OpenAI, que destaca por su capacidad para generar texto coherente y relevante en diferentes idiomas. De esta manera, se obtubieron informes médicos con texto sensible de manera semi-esctructurada.

Seguido de esto, se agregaron errores ortográficos y gramaticales a los informes médicos generados, como **Agusitn acheccco De la Fuente** en lugar de **Agustín Pachecco De la Fuente** y **smtp:/LSCNBarbadosYOrtega-biotechcom** en lugar de **smtp://LSCNBarbadosYOrtega-biotech.com**, para simular la variabilidad y complejidad de los datos reales. A demás se incluyeron notaciones especificadas en la guía de anotación de MEDDOCAN, como abreviaturas e iniciales de nombres , apodos, titulos nobiliarios, menciones seguidas del género, distintas notaciones de fechas, entre otros.

Por último, se modificaron las notaciones de los datos dados por el generador para que coincidieran con las modificaciones realizadas de manera manual, calculando las posiciones de las entidades con un [programa en Python](./data/brat/procesar-texto.py).

De esta manera, la cantidad de ocurrencias de cada tipo de entidad en el conjunto de datos generado es la siguiente:
* TERRITORIO: 26
* FECHAS: 16
* EDAD_SUJETO_ASISTENCIA: 9
* NOMBRE_SUJETO_ASISTENCIA: 12
* NOMBRE_PERSONAL_SANITARIO: 6
* SEXO_SUJETO_ASISTENCIA: 7
* CALLE: 8
* PAIS: 3
* ID_SUJETO_ASISTENCIA: 10
* CORREO_ELECTRONICO: 5
* ID_TITULACION_PERSONAL_SANITARIO: 4
* ID_ASEGURAMIENTO: 5
* HOSPITAL: 3
* FAMILIARES_SUJETO_ASISTENCIA: 6
* INSTITUCION: 5
* ID_CONTACTO_ASISTENCIAL: 4
* NUMERO_TELEFONO: 7
* PROFESION: 1
* NUMERO_FAX: 3
* OTROS_SUJETO_ASISTENCIA: 5
* CENTRO_SALUD: 3
* ID_EMPLEO_PERSONAL_SANITARIO: 10
* IDENTIF_VEHICULOS_NRSERIE_PLACAS: 4
* IDENTIF_DISPOSITIVOS_NRSERIE: 10
* NUMERO_BENEF_PLAN_SALUD: 0
* URL_WEB: 3
* DIREC_PROT_INTERNET: 7
* IDENTF_BIOMETRICOS: 8
* OTRO_NUMERO_IDENTIF: 5

Siendo un total de 195 entidades en el conjunto de datos generado, incluídas entidades de tipos que no habían sido utilizados en los datos del corpus original, como ID_EMPLEO_PERSONAL_SANITARIO, IDENTIF_VEHICULOS_NRSERIE_PLACAS, IDENTIF_DISPOSITIVOS_NRSERIE, URL_WEB, DIREC_PROT_INTERNET, IDENTF_BIOMETRICOS y OTRO_NUMERO_IDENTIF.
Cabe aclarar que no se agregaron entidades de tipo NUMERO_BENEF_PLAN_SALUD ya que, a pesar de que se encuentra listada como una entidad sensible en la presentación del proyecto, no se detalla en la guía de anotación utilizada.

### Ejecución de los modelos

Dado que los modelos presentados fueron desarrollador a fines del 2019 y no se han actualizado desde entonces, se presentaron problemas de compatibilidad con las versiones de las bibliotecas utilizadas y requerimientos particulares. Por lo tanto, se realizaron modificaciones en los modelos para poder ejecutarlos con éxito en el conjunto de datos generado. A continuación, se detallan los pasos necesarios para ejecutar cada modelo:

1. Descargar e inicializar los repositorios de los modelos:
``` bash
git clone --recurse-submodules "https://github.com/JulietaStorino/Text-Mining-Proyect-FaMAF.git"
cd Text-Mining-Proyect-FaMAF/SPACCC_MEDDOCAN
git clone "https://github.com/PlanTL-GOB-ES/SPACCC_MEDDOCAN.git" .
cd ..
```

#### BiLSTM-CRF:

2. Descomprimir y reestructurar los datos de entrenamiento
``` bash
cd models/BiLSTM-CRF
unzip data.zip
mv data/train .
mv data/dev .
mkdir output
mv data/test output
```
3. Crear el entorno virtual
``` bash
python3 -m venv .env
source .env/bin/activate
```
4. Instalar las dependencias
``` bash
pip install tensorflow==1.14.0 numpy==1.16.4 scipy==1.3.0 cython wheel spacy==2.3.2 nltk matplotlib gast==0.2.2 scikit-learn
```
5. Descargar pipeline en español optimizado para CPU
``` bash
python3 -m spacy download es_core_news_sm
```
6. Reemplaza todas las instancias de spacy.load('es') por spacy.load('es_core_news_sm')
``` bash
cd code
sed -i "s/spacy.load('es')/spacy.load('es_core_news_sm')/g" preprocessing.py
```
7. Descargar los recursos necesarios
``` bash
python3 ../../../requirements.py
```
8. Preprocesar los datos
``` bash
python3 preprocessing.py --dataDir ../train/gold --train
python3 preprocessing.py --dataDir ../dev/gold --dev
python3 preprocessing.py --dataDir ../test/gold --test
```
9. Descarga el embedding de palabras necesario
``` bash
cd Extension2
wget -O wiki.es.vec https://dl.fbaipublicfiles.com/fasttext/vectors-wiki/wiki.es.vec
```
10. Genera los archivos necesarios por el modelo neuronal
``` bash
python3 create_vocabs.py --trainpickle ../train_word_ner_startidx_dict.pickle --devpickle ../dev_word_ner_startidx_dict.pickle --testpickle ../test_word_ner_startidx_dict.pickle --embfile wiki.es.vec --vocabEmbFile vocab_embeddings.npz
```
11. Crear directorios para guardar el modelo y los confusion plots
``` bash
mkdir models
mkdir plots
cd Code
```
12. Correr el modelo
``` bash
python3 train.py
```
13. Evaluar el modelo con los datos del MEDDOCAN
``` bash
cd ../..
python3 evaluate.py brat ner ../output/test/gold ../output/test/system
```
14. Preprocesar los datos nuevos
``` bash
python3 preprocessing.py --dataDir ../../../data/brat/gold --test
```
15. Genera los archivos necesarios por el modelo neuronal
``` bash
cd Extension2
python3 create_vocabs.py --trainpickle ../train_word_ner_startidx_dict.pickle --devpickle ../dev_word_ner_startidx_dict.pickle --testpickle ../test_word_ner_startidx_dict.pickle --embfile wiki.es.vec --vocabEmbFile vocab_embeddings.npz
```
16. Modificar el código con la ubicación de los nuevos modelos
``` bash
cd Code
sed -i 's|\.\./\.\./\.\./output/test/system/|\.\./\.\./\.\./\.\./\.\./data/brat/system/|g' train.py
```
17. Correr el modelo
``` bash
python3 train.py
```
18. Evaluar los datos nuevos
``` bash
cd ../..
python3 evaluate.py brat ner ../../../data/brat/gold ../../../data/brat/system
```

#### CLIN-X:

2. Crear el entorno virtual
``` bash
cd models/CLIN-X
conda create -n clin-x python==3.8.5
conda activate clin-x
```
3. Instalar las dependencias
``` bash
pip3 install flair==0.8 transformers==4.6.1 torch==1.8.1 scikit-learn==0.23.1 scipy==1.6.3 numpy nltk tqdm seaborn matplotlib
```
4. Crear las carpetas *bio_files*, *split_files*, *models*, *predictions*, *models/clin_x_experiment* y *predictions/clin_x_experiment*
``` bash
mkdir bio_files
mkdir bio_files/train
mkdir bio_files/dev
mkdir bio_files/test_m
mkdir bio_files/test_o
mkdir split_files
mkdir models
mkdir predictions
mkdir models/clin_x_experiment
mkdir predictions/old_data
mkdir predictions/new_data
```
5. En caso de no tener una GPU disponible, modificar el archivo *create_data_splits.py*
``` bash
sed -i 's/.cuda()/.to("cpu")/g' create_data_splits.py
```
6. Tokenizar los archivos de entrennamiento, desarrollo y prueba (del corpus MEDDOCAN y nuevos), convirtiéndolos en archivos BIO con el formato requerido
``` bash
python tokenize_files.py --input_path ../../SPACCC_MEDDOCAN/corpus/train/brat/ --output_path bio_files/train/
python tokenize_files.py --input_path ../../SPACCC_MEDDOCAN/corpus/dev/brat/ --output_path bio_files/dev/
python tokenize_files.py --input_path ../../SPACCC_MEDDOCAN/corpus/test/brat/ --output_path bio_files/test_m/
python tokenize_files.py --input_path ../../data/brat/gold/ --output_path bio_files/test_o/
```
7. Dividir los datos tokenizados de entranamiento y desarrollo en varias partes de manera aleatoria
``` bash
python create_data_splits.py --train_files bio_files/train/ --dev_files bio_files/dev/ --method random --output_dir split_files/
```
8. Entrenar el modelo
``` bash
python train_our_model_architecture.py --data_path split_files/ --train_files random_split_1.txt,random_split_2.txt,random_split_3.txt,random_split_4.txt --dev_file random_split_5.txt --model llange/xlm-roberta-large-spanish-clinical --name clin_x_experiment --storage_path models
```
9. Realizar las predicciones
``` bash
python get_test_predictions.py --name models/clin_x_experiment/ --conll_path bio_files/test_m/ --out_path predictions/old_data/
python get_test_predictions.py --name models/clin_x_experiment/ --conll_path bio_files/test_o/ --out_path predictions/new_data/
```
10. Evaluar los datos

#### NeuroNer:

Para hacer más fácil la ejecución de este modelo, ya que originalmente se corre con una notebook de Jupyter con los parámetros ya definidos, se convirtió el archivo Train.ipynb a un archivo .py para poder ejecutarlo en la terminal, agregando a demás distinción entre los archivos de salida old_data y new_data, para poder evaluar los datos del MEDDOCAN y los datos generados respectivamente. 

2. Crear el entorno virtual
``` bash
cd models/NeuroNer/
python3 -m venv .env
source .env/bin/activate
```
3. Instalar las dependencias
``` bash
pip install cython tensorflow==1.14.0 keras==2.2.4 spacy pyneuroner
```
4. Mover el archivo NeuroNer.py a la carpeta del modelo
``` bash
mv ../../NeuroNer.py .
```
5. Crear los directorios necesarios
``` bash
mkdir output/old_data
mkdir output/new_data
```
6. Correr el modelo
``` bash
python NeuroNer.py
```
7. Evaluar los datos

### Análisis y comparación de los resultados

#### BiLSTM-CRF:
* Datos de prueba (MEDDOCAN):
``` bash
Report (SYSTEM: system):
------------------------------------------------------------
SubTrack 1 [NER]                   Measure        Micro               
------------------------------------------------------------
Total (156 docs)                   Precision      0.8918              
                                   Recall         0.8337              
                                   F1             0.8618              
------------------------------------------------------------
```
* Datos de prueba (generados):
``` bash                                                                 
Report (SYSTEM: system):
------------------------------------------------------------
SubTrack 1 [NER]                   Measure        Micro               
------------------------------------------------------------
Total (5 docs)                     Precision      0.3871              
                                   Recall         0.3077              
                                   F1             0.3429              
------------------------------------------------------------
```

#### CLIN-X:

#### NeuroNer:
  
### Comparación y análisis de los resultados

En relación con el primer modelo, BiLSTM-CRF, se observa que su desempeño en la detección de entidades y posiciones es inconsistente. Aunque logra identificar correctamente las entidades más comunes y típicas del conjunto de datos original, presenta dificultades al tratar casos más complejos. Por ejemplo, al añadir números aleatorios a las direcciones, siguiendo las guías de notación, el modelo no logra reconocerlos adecuadamente, lo que afecta su rendimiento. Asimismo, algunos números de teléfono no son detectados, lo que impide su anonimización y, por ende, podría resultar en la identificación de las personas a partir de esos datos.

Además, el modelo mostró un deterioro significativo en su desempeño al evaluar datos no vistos, con una precisión del 0.3871 y una recuperación de 0.3077. Esto sugiere que ha sufrido de sobreajuste al conjunto de entrenamiento, lo que limita su capacidad de generalización a nuevos contextos. Aunque fue capaz de identificar correctamente las entidades comunes en el conjunto original, tuvo dificultades con información menos frecuente y con errores ortográficos intencionales que se introdujeron en los datos generados. Esta falta de reconocimiento de variantes en nombres, direcciones y números de teléfono indica que el modelo necesita una mayor robustez para adaptarse a la variabilidad natural del lenguaje.

Si comparamos los cinco casos clínicos, podemos observar las siguientes diferencias entre el texto anonimizado generado por el modelo y la versión correcta, siguiendo las guías de notación:

* Fragmentación incorrecta de entidades: T8	TERRITORIO 204 212	Castilla T9	TERRITORIO 213 214	y T10	TERRITORIO 215 219	León
* Etiquetado incorrecto de datos: T9	CALLE 221 234	+34 941 83 66
* Desajuste en las posiciones de las etiquetas: T25	NOMBRE_PERSONAL_SANITARIO 3687 3713	José Ángel Guerrero Rivera en vez de T31	NOMBRE_PERSONAL_SANITARIO 1641 1667	José Ángel Guerrero Rivera
* Inclusión de información irrelevante: T20	ID_SUJETO_ASISTENCIA 1142 1215	acompañado por sus dos hijas de 20 años y 35 años debido a la reciente de
* Omisión de etiquetas: T36	IDENTF_BIOMETRICOS 3028 3039	RG545EMR302 la etiqueta para los identificadores biométricos no fue reconocida ni aplicada.
* Errores al combinar información de diferentes categorías en una única etiqueta: T6	CALLE 140 195	Paseo de Zorrilla, 1, quinto a la izquierda, Valladolid o T11	ID_ASEGURAMIENTO 232 269	Fijo 983 82 84 99, Móvil 683 78 02 65
* Inconsistencia en el etiquetado de datos: T19	ID_SUJETO_ASISTENCIA 586 640	asociada 249.36.206.164 - MAC asociada 3D9.E0F.A68.E2A
* Formato incorrecto al etiquetar datos: T32	ID_SUJETO_ASISTENCIA 3191 3221	906382828.

### Conclusiones

Si tuviéramos la oportunidad o responsabilidad de seguir con el proyecto dentro de un año, contando con un equipo de cinco empleados, podríamos optimizar diversos aspectos clave. En primer lugar, sería fundamental obtener un dataset más amplio y representativo, lo que permitiría mejorar la precisión y fiabilidad de los resultados obtenidos por los modelos. 

Asimismo, sería necesario mejorar los algoritmos utilizados, ya sea optimizando los existentes o explorando nuevas arquitecturas que se adapten mejor a la naturaleza específica de los datos. Esto incluiría mejorar algunos de los modelos existentes para garantizar que se ajusten adecuadamente a las particularidades del lenguaje y los contextos clínicos representados.

Por último, se podría seleccionar una mayor cantidad de modelos y comparar su desempeño con respecto a la anonimización de los textos médicos.
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
