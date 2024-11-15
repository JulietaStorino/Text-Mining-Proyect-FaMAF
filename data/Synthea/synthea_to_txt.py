import json
import os
import re

def cargar_datos(archivo):
    """Carga los datos de un archivo JSON."""
    with open(archivo, 'r') as file:
        return json.load(file)

def obtener_id_paciente(data):
    """Extrae el ID del paciente de los datos."""
    for entry in data.get("entry", []):
        if entry.get("resource", {}).get("resourceType") == "Patient":
            return entry["resource"].get("id")
    return None

def obtener_lugar_nacimiento(data):
    """Extrae el lugar de nacimiento del paciente de los datos."""
    for entry in data.get("entry", []):
        if entry.get("resource", {}).get("resourceType") == "Patient":
            extensions = entry["resource"].get("extension", [])
            for ext in extensions:
                if ext.get("url") == "http://hl7.org/fhir/StructureDefinition/patient-birthPlace":
                    return ext.get("valueAddress")
    return None

def procesar_archivos(carpeta):
    for archivo_nombre in os.listdir(carpeta):
        if archivo_nombre.endswith(".txt"):
            archivo_path = os.path.join(carpeta, archivo_nombre)
            
            with open(archivo_path, 'r', encoding='utf8') as file:
                contenido = file.read()
            
            # Buscar y reemplazar números en el nombre
            contenido_modificado = re.sub(r'Nombre: (Mrs\.|Mr\.|Ms\.) (([\wáéíóúÁÉÍÓÚñÑ]+)\d* ){1,2}([\wáéíóúÁÉÍÓÚñÑ]+)\d*', r'Nombre: \1 \2 \4', contenido)
                                    
            # Guardar el contenido modificado en el mismo archivo
            with open(archivo_path, 'w', encoding='utf8') as file:
                file.write(contenido_modificado)

def obtener_datos_paciente(data):
    """Extrae los datos principales del paciente de los datos."""
    for entry in data.get("entry", []):
        if entry.get("resource", {}).get("resourceType") == "Patient":
            paciente = entry["resource"]
            # Extraer nombre completo
            nombre = f"{paciente['name'][0].get('prefix', [''])[0]} {paciente['name'][0].get('given', [''])[0]} {paciente['name'][0].get('family', '')}"
            
            # Extraer dirección
            address_data = paciente.get("address", [{}])[0]
            direccion = f"{address_data.get('line', [''])[0]}, {address_data.get('city', '')}, {address_data.get('state', '')} {address_data.get('postalCode', '')}, {address_data.get('country', '')}"
            
            # Extraer teléfono
            telefono = paciente.get("telecom", [{}])[0].get("value", "Número de teléfono no disponible")
            
            # Extraer género
            genero = paciente.get("gender", "Género no disponible")
            
            return nombre, direccion, telefono, genero
    return None, None, None, None



if (__name__ == '__main__'):
    print("Ejecutando el script synthea_to_txt.py")
    print("Este script convierte los archivos JSON generados por Synthea en archivos de texto con información de los pacientes.")
    print("Para usar este script, debes tener descargado Synthea y haber corrido el archivo 'generate_samples.sh'.")
    print("Luego, modifica la variable 'carpeta' con el path a tu carpeta de samples.")
    print("El script guardará los archivos de texto en la carpeta 'datos_synthea'.")
    print("Para más información, revisa el código del script.")
    print("¡Gracias por usar el script!")

    # Ruta de la carpeta que contiene los archivos JSON
    carpeta = "synthea/samples/" # 
    salida_carpeta = "data/"

    # Asegurarse de que la carpeta de salida existe
    os.makedirs(salida_carpeta, exist_ok=True)

    # Procesar cada archivo en la carpeta
    for archivo_nombre in os.listdir(carpeta):
        if archivo_nombre.endswith(".json"):
            archivo_path = os.path.join(carpeta, archivo_nombre)
            
            # Cargar datos del archivo
            data = cargar_datos(archivo_path)

            # Crear contenido de información del paciente
            patient_id = obtener_id_paciente(data)
            patient_info = f"ID del Paciente: {patient_id}\n" if patient_id else "No se encontró un recurso de tipo 'Patient'.\n"

            # Obtener y agregar el lugar de nacimiento
            birthplace = obtener_lugar_nacimiento(data)
            if birthplace:
                ciudad = birthplace.get("city", "Ciudad no disponible")
                estado = birthplace.get("state", "Estado no disponible")
                pais = birthplace.get("country", "País no disponible")
                patient_info += "Lugar de nacimiento:\n"
                patient_info += f"Ciudad: {ciudad}\n"
                patient_info += f"Estado: {estado}\n"
                patient_info += f"País: {pais}\n"
            else:
                patient_info += "No se encontró 'valueAddress' en el recurso 'Patient'.\n"

            # Obtener y agregar los datos del paciente
            nombre, direccion, telefono, genero = obtener_datos_paciente(data)
            patient_info += "\nDatos del Paciente:\n"
            patient_info += f"Nombre: {nombre}\n"
            patient_info += f"Dirección: {direccion}\n"
            patient_info += f"Teléfono: {telefono}\n"

            # Traducir el género
            if genero == "male":
                genero = "hombre"
            elif genero == "female":
                genero = "mujer"
            patient_info += f"Género: {genero}\n"

            # Guardar la información en un archivo de texto con el mismo nombre que el JSON
            output_file = os.path.join(salida_carpeta, f"{os.path.splitext(archivo_nombre)[0]}.txt")
            with open(output_file, "w") as file:
                file.write(patient_info)

            # Saca los numeros de los nombres
            carpeta = "./data/"
            procesar_archivos(carpeta)

    print(f"La información de los pacientes se ha guardado en la carpeta '{salida_carpeta}'.")
