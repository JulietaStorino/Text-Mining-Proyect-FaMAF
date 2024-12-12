def procesar_caracteres(archivo):
    """
    Imprime en pantalla una lista con las posiciones de los caracteres '×' en el archivo de texto en pares
    junto con el texto entre ellos.
    """
    
    with open(archivo, 'r') as f:
        texto = f.read()

    pos = []
    j = 0
    for i, char in enumerate(texto):
        if char == '×':
            pos.append(j)
        else:
            j += 1

    # Elimina los '×' del archivo
    texto_sin_x = texto.replace('×', '')
    with open(archivo, 'w') as f:
        f.write(texto_sin_x)

    # Imprimir posiciones de a pares
    for i in range(0, len(pos), 2):
        num = int(i / 2) + 1
        entidad = texto_sin_x[pos[i]:pos[i + 1]]
        if num < 10:
            print(f"T{num}  {pos[i]} {pos[i + 1]}\t{entidad}")
        else:
            print(f"T{num} {pos[i]} {pos[i + 1]}\t{entidad}")