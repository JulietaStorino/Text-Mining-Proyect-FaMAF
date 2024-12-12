import re

def contar_menciones():
    # Inicializar contadores
    menciones = {
        "TERRITORIO": 0,
        "FECHAS": 0,
        "EDAD_SUJETO_ASISTENCIA": 0,
        "NOMBRE_SUJETO_ASISTENCIA": 0,
        "NOMBRE_PERSONAL_SANITARIO": 0,
        "SEXO_SUJETO_ASISTENCIA": 0,
        "CALLE": 0,
        "PAIS": 0,
        "ID_SUJETO_ASISTENCIA": 0,
        "CORREO_ELECTRONICO": 0,
        "ID_TITULACION_PERSONAL_SANITARIO": 0,
        "ID_ASEGURAMIENTO": 0,
        "HOSPITAL": 0,
        "FAMILIARES_SUJETO_ASISTENCIA": 0,
        "INSTITUCION": 0,
        "ID_CONTACTO_ASISTENCIAL": 0,
        "NUMERO_TELEFONO": 0,
        "PROFESION": 0,
        "NUMERO_FAX": 0,
        "OTROS_SUJETO_ASISTENCIA": 0,
        "CENTRO_SALUD": 0,
        "ID_EMPLEO_PERSONAL_SANITARIO": 0,
        "IDENTIF_VEHICULOS_NRSERIE_PLACAS": 0,
        "IDENTIF_DISPOSITIVOS_NRSERIE": 0,
        "NUMERO_BENEF_PLAN_SALUD": 0,
        "URL_WEB": 0,
        "DIREC_PROT_INTERNET": 0,
        "IDENTF_BIOMETRICOS": 0,
        "OTRO_NUMERO_IDENTIF": 0,
    }

    archivos ={
        "brat/067253898.ann",
        "brat/149891773.ann",
        "brat/387802160.ann",
        "brat/486642396.ann",
        "brat/572018553.ann",
    }

    for archivo in archivos:

        # Leer el archivo
        with open(archivo, 'r', encoding='utf-8') as f:
            contenido = f.read()

        # Contar menciones
        for entidad in menciones.keys():
            menciones[entidad] += len(re.findall(entidad, contenido))

    # Imprimir resultados
    for entidad, cuenta in menciones.items():
        print(f"{entidad}: {cuenta}")

    # Imprimir total
    total = sum(menciones.values())
    print(f"\n\nTotal: {total}")

if __name__ == "__main__":
    contar_menciones()
