from math import sqrt
from collections import defaultdict
import math

def grid_neighbors(coords, d):
    """
    Función optimizada para encontrar vecinos dentro de una distancia `d` usando una estructura de cuadrícula.
    """
    grid = defaultdict(list)
    cell_size = d  # Tamaño de la celda en función de la distancia d

    # Asigna cada coordenada a una celda en la cuadrícula
    for i, (x, y) in enumerate(coords):
        cell_x, cell_y = int(x // cell_size), int(y // cell_size)
        grid[(cell_x, cell_y)].append(i)

    neighbors = defaultdict(list)

    # Para cada célula, buscar vecinos en celdas adyacentes
    for i, (x, y) in enumerate(coords):
        cell_x, cell_y = int(x // cell_size), int(y // cell_size)

        # Recorrer solo las celdas adyacentes y la propia celda
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                neighbor_cell = (cell_x + dx, cell_y + dy)

                # Verifica cada punto en la celda vecina
                for j in grid.get(neighbor_cell, []):
                    if i != j:
                        # Verifica la distancia real entre puntos para confirmar que está en el rango
                        dist = math.sqrt((coords[i][0] - coords[j][0]) ** 2 + (coords[i][1] - coords[j][1]) ** 2)
                        if dist <= d:
                            neighbors[i].append(j)
    
    return neighbors

def construir_grafo(celulas, d):
    grafo = defaultdict(list)
    n = len(celulas)
    
    coords = []
    for i in range(n):
        id_celula, x, y, peptidos = celulas[i]
        coords.append((x, y))  # Guardamos solo las coordenadas

    # Usamos la función grid_neighbors para encontrar los vecinos
    vecinos = grid_neighbors(coords, d)

    for i in range(n):
        id1, x1, y1, peptidos1 = celulas[i]
        for j in vecinos.get(i, []):
            id2, x2, y2, peptidos2 = celulas[j]
            if peptidos1 & peptidos2:  # Verifica si tienen péptidos comunes
                grafo[id1].append(id2)
                grafo[id2].append(id1)

    return grafo

def componentes_conexas(grafo):
    visitado = set()
    componentes = []

    def bfs(nodo):
        cola = deque([nodo])
        componente = []
        while cola:
            actual = cola.popleft()
            if actual not in visitado:
                visitado.add(actual)
                componente.append(actual)
                cola.extend(grafo[actual])
        return componente

    for nodo in grafo:
        if nodo not in visitado:
            componente = bfs(nodo)
            componentes.append(componente)
    
    return componentes

def resolver_caso(caso):
    n, d = caso["n"], caso["d"]
    celulas = []
    for entrada in caso["celulas"]:
        id_celula = entrada[0]
        x, y = entrada[1], entrada[2]
        peptidos = set(entrada[3:])
        celulas.append((id_celula, x, y, peptidos))

    grafo = construir_grafo(celulas, d)
    grupos = componentes_conexas(grafo)

    resultado = {}
    for i, grupo in enumerate(grupos, start=1):
        for celula in grupo:
            resultado[celula] = i

    return resultado

def main():
    # Casos de prueba proporcionados
    casos = [
        {
            "n": 7,
            "d": 1,
            "celulas": [
                [1, 0, 0, "AETQT", "DFTYA", "PHLYT"],
                [2, 0, 2, "DSQTS", "IYHLK", "LHGPS", "LTLLS"],
                [3, 1, 0, "AETQT", "DFTYA", "HGCYS", "LSVGG", "SRFNH"],
                [4, 1, 1, "DFTYA", "HGCYS", "IYHLK", "SRFNH"],
                [5, 1, 2, "DSQTS", "IYHLK", "LSVGG", "LTLLS", "TTVTG"],
                [6, 2, 1, "AETQT", "HGCYS", "IYHLK", "LSVGG", "LTLLS"],
                [7, 2, 2, "HGCYS", "SRFNH", "TTVTG"],
            ]
        },
        {
            "n": 7,
            "d": 2,
            "celulas": [
                [1, 0, 0, "AETQT", "DFTYA", "PHLYT"],
                [2, 0, 2, "DSQTS", "IYHLK", "LHGPS", "LTLLS"],
                [3, 1, 0, "AETQT", "DFTYA", "HGCYS", "LSVGG", "SRFNH"],
                [4, 1, 1, "DFTYA", "HGCYS", "IYHLK", "SRFNH"],
                [5, 1, 2, "DSQTS", "IYHLK", "LSVGG", "LTLLS", "TTVTG"],
                [6, 2, 1, "AETQT", "HGCYS", "IYHLK", "LSVGG", "LTLLS"],
                [7, 2, 2, "HGCYS", "SRFNH", "TTVTG"],
            ]
        },
    ]

    for i, caso in enumerate(casos, start=1):
        resultado = resolver_caso(caso)
        print(f"Caso {i}:")
        for id_celula in sorted(resultado):
            print(id_celula, resultado[id_celula])

if __name__ == "__main__":
    main()
