from huggingface_hub import login
from datasets import load_dataset as hf_load_dataset
from random import randint
import json
import os

def download_and_save_dataset():
    """
    Descarga el dataset de testeo del programa CANTEMIST.
    """

    # Inicia sesión en Hugging Face
    token = input("Introduce tu token de Hugging Face: ")
    login(token=token)

    # Cargar el conjunto de datos
    ds = hf_load_dataset("PlanTL-GOB-ES/cantemist-ner")

    # Selecciona la división del conjunto de datos que deseas guardar, por ejemplo, 'train'
    data = ds['train']

    # Convierte los datos a una lista de diccionarios
    data_list = [dict(item) for item in data]

    # Especifica el nombre del archivo
    output_file = 'cantemist.json'

    # Guarda los datos en un archivo JSON
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data_list, f, ensure_ascii=False, indent=4)

def process_dataset():
    """
    Procesa el dataset de testeo del programa CANTEMIST.
    """

    # Carga el archivo JSON
    with open('cantemist.json', 'r', encoding='utf-8') as f:
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
            tokens = item['tokens']
            file.write('{ "Historia clínica": "')
            for token in tokens:
                if token != "," and token != ".":
                    file.write(" ")
                file.write(token)
            file.write('" }')

if __name__ == '__main__':
    print("Descargando el dataset de testeo del programa CANTEMIST...")
    download_and_save_dataset()
    
    if os.path.exists('cantemist.json'):
        print("Descarga completada.\n")
        print("Procesando el dataset...")
        
        process_dataset()

    else:
        print("Error al descargar el dataset.")
        exit()