from huggingface_hub import login
from datasets import load_dataset
from random import randint
import json
import os

def download_and_save_dataset():
    """
    Descarga el dataset de testeo del programa CARES.
    """

    # Inicia sesión en Hugging Face
    token = input("Introduce tu token de Hugging Face: ")
    login(token=token)

    # Carga el conjunto de datos
    ds = load_dataset("chizhikchi/CARES")

    # Selecciona la división del conjunto de datos 'train'
    data = ds['train']

    # Convierte los datos a una lista de diccionarios
    data_list = [dict(item) for item in data]

    # Especifica el nombre del archivo
    output_file = 'cares.json'

    # Guarda los datos en un archivo JSON
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data_list, f, ensure_ascii=False, indent=4)

def process_dataset():
    """
    Procesa el dataset de testeo del programa CARES.
    """

    # Carga el archivo JSON
    with open('cares.json', 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Procesa los datos
    for item in data:
        id = randint(1000000, 9999999)

        # Añade un ID único a cada entrada
        file_name = os.path.join("./data/", f"{id}.json")

        # Crea la carpeta si no existe
        os.makedirs("./data/", exist_ok=True)

        # Guarda los datos en un archivo JSON
        with open(file_name, 'w', encoding='utf8') as file:
            text = item['full_text']
            file.write('{ "Historia clínica": "')
            file.write(text)
            file.write('"}')

if __name__ == '__main__':
    print("Descargando el dataset CARES...")
    download_and_save_dataset()
    
    if os.path.exists('cares.json'):
        print("Descarga completada.")
        print("Procesando el dataset...")
        
        process_dataset()

        print("Procesamiento completado.")

    else:
        print("Error al descargar el dataset.")
        exit()