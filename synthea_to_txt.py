import json
import os

# Para usar este archivo, debes tener descargado synthea y haber corrido el archivo "generate_samples.sh"
# Luego, modificar la variable "carpeta" con el path a tu carpeta de samples

# Ruta de la carpeta que contiene los archivos JSON
carpeta = "synthea/samples/" # 
salida_carpeta = "datos_synthea"

# Asegurarse de que la carpeta de salida existe
os.makedirs(salida_carpeta, exist_ok=True)

# Función para cargar el JSON del archivo
def cargar_datos(archivo):
    with open(archivo, 'r') as file:
        return json.load(file)

# Función para extraer el ID del paciente
def obtener_id_paciente(data):
    for entry in data.get("entry", []):
        if entry.get("resource", {}).get("resourceType") == "Patient":
            return entry["resource"].get("id")
    return None

# Función para obtener el lugar de nacimiento
def obtener_lugar_nacimiento(data):
    for entry in data.get("entry", []):
        if entry.get("resource", {}).get("resourceType") == "Patient":
            extensions = entry["resource"].get("extension", [])
            for ext in extensions:
                if ext.get("url") == "http://hl7.org/fhir/StructureDefinition/patient-birthPlace":
                    return ext.get("valueAddress")
    return None

# Función para obtener los datos principales del paciente
def obtener_datos_paciente(data):
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

print(f"La información de los pacientes se ha guardado en la carpeta '{salida_carpeta}'.")
