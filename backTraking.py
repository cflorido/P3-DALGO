from math import sqrt
from collections import defaultdict
import math

def grid_neighbors(coords, d):
    """
    Función optimizada para encontrar vecinos dentro de una distancia d usando una estructura de cuadrícula.
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

from collections import defaultdict
import math

def es_clique(grafo, clique, nodo):
    """
    Verifica si al agregar nodo al clique sigue siendo un clique válido.
    """
    for v in clique:
        if nodo not in grafo[v]:
            return False
    return True

def backtracking_cliques(grafo, cliques, nodo_actual, solucion):
    """
    Algoritmo de backtracking para encontrar la partición en cliques.
    """
    if nodo_actual == len(grafo):
        return True  # Todos los nodos han sido procesados

    nodo = list(grafo.keys())[nodo_actual]
    for i, clique in enumerate(cliques):
        if es_clique(grafo, clique, nodo):
            clique.append(nodo)
            solucion[nodo] = i + 1  # Asignar el clique al nodo
            if backtracking_cliques(grafo, cliques, nodo_actual + 1, solucion):
                return True
            clique.pop()  # Backtracking

    # Crear un nuevo clique con este nodo
    cliques.append([nodo])
    solucion[nodo] = len(cliques)
    if backtracking_cliques(grafo, cliques, nodo_actual + 1, solucion):
        return True
    cliques.pop()  # Backtracking

    return False

def resolver_caso_backtracking(caso):
    """
    Resuelve el problema usando Backtracking.
    """
    n, d = caso["n"], caso["d"]
    celulas = []
    for entrada in caso["celulas"]:
        id_celula = entrada[0]
        x, y = entrada[1], entrada[2]
        peptidos = set(entrada[3:])
        celulas.append((id_celula, x, y, peptidos))

    grafo = construir_grafo(celulas, d)
    
    # Backtracking para encontrar el mínimo número de cliques
    solucion = {}
    cliques = []
    if backtracking_cliques(grafo, cliques, 0, solucion):
        return solucion

    return None  # No se encontró solución (en teoría, no debería ocurrir)

def construir_grafo(celulas, d):
    grafo = defaultdict(list)
    n = len(celulas)
    
    # Extraemos las coordenadas de las células para pasarlas a la función grid_neighbors
    coords = [(x, y) for _, x, y, _ in celulas]

    # Usamos la función grid_neighbors para encontrar los vecinos de las células
    neighbors = grid_neighbors(coords, d)
    
    for i in range(n):
        id1, _, _, peptidos1 = celulas[i]
        for j in neighbors[i]:
            id2, _, _, peptidos2 = celulas[j]
            if peptidos1 & peptidos2:  # Si tienen péptidos en común
                grafo[id1].append(id2)
                grafo[id2].append(id1)

    return grafo

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
        resultado = resolver_caso_backtracking(caso)
        print(f"Caso {i}:")
        for id_celula in sorted(resultado):
            print(id_celula, resultado[id_celula])

if __name__ == "__main__":
    main()
