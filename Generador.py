import random
import math
import sys

def generar_peptidos():
    # Genera entre 1 y 5 péptidos aleatorios de 5 caracteres únicos
    amino_acidos = 'ACDEFGHIKLMNPQRSTVWY'
    num_peptidos = random.randint(1, 5)
    peptidos = set()
    while len(peptidos) < num_peptidos:
        peptido = ''.join(random.choice(amino_acidos) for _ in range(5))
        peptidos.add(peptido)
    return list(peptidos)

def generar_celulas(n, d):
    # Genera n células con coordenadas y péptidos aleatorios
    celulas = []
    for i in range(1, n + 1):
        x = random.randint(0, d * 10)  # Coordenadas en un rango proporcional a d
        y = random.randint(0, d * 10)
        peptidos = generar_peptidos()
        celulas.append((i, x, y, peptidos))
    return celulas

def escribir_caso_prueba(n, d, celulas):
    # Escribe un caso de prueba en el formato requerido
    resultado = []
    resultado.append(f"{n} {d}")
    for (id, x, y, peptidos) in celulas:
        linea = f"{id} {x} {y} " + ' '.join(peptidos)
        resultado.append(linea)
    return '\n'.join(resultado)

def generar_entrada(num_casos, n_min, n_max, d_max, filename):
    casos = []
    for _ in range(num_casos):
        n = random.randint(n_min, n_max)
        d = random.randint(1, d_max)
        celulas = generar_celulas(n, d)
        caso = escribir_caso_prueba(n, d, celulas)
        casos.append(caso)

    with open(filename, 'w') as f:
        f.write(f"{len(casos)}\n")
        for caso in casos:
            f.write(caso + "\n")

def main():
    # Configuración inicial
    try:
        num_casos = int(input("Número de casos de prueba: "))
        n_min = int(input("Número mínimo de células: "))
        n_max = int(input("Número máximo de células: "))
        d_max = int(input("Distancia máxima para conexión: "))
        filename = input("Nombre del archivo de salida (ejemplo: casos_prueba.txt): ")
    except ValueError:
        print("Por favor, ingresa valores válidos para cada parámetro.")
        sys.exit(1)

    # Genera los escenarios de prueba
    generar_entrada(num_casos, n_min, n_max, d_max, filename)
    print(f"Escenarios de prueba generados en el archivo '{filename}'.")

if __name__ == "__main__":
    main()
