import json
import os

# Ruta del archivo que contiene el archivo JSON
archivo = "CARES_DATASET.json"
salida_carpeta = "./"

# Asegurarse de que la carpeta de salida existe si no lo crea
os.makedirs(salida_carpeta, exist_ok=True)

# Función para cargar el JSON del archivo
def cargar_datos(archivo):
    with open(archivo, 'r', encoding='utf8') as file:
        return json.load(file)

# Procesar cada archivo en la carpeta
def procesar_datos():
    archivo_path = os.path.join("./", archivo)
        
    # Cargar datos del archivo
    data = cargar_datos(archivo_path)

    for entry in data['rows']:
        # Crear un archivo por cada entrada
        file_name = os.path.join(salida_carpeta, f"{entry['row_idx']}.txt")
        with open(file_name, 'w', encoding='utf8') as file:

            # Guarda la edad del paciente
            file.write(entry['row']['full_text'])

procesar_datos()
