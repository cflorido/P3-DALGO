import numpy as np
import random
from collections import defaultdict
from math import sqrt
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

    return grafo

def light_backtrack(adj_mat, cliques, v=0, best=(math.inf, None)):
    n = adj_mat.shape[0]

    if v == n:
        if is_solution(cliques, adj_mat):
            num_cliques = len(set(cliques))
            if num_cliques < best[0]:
                best = (num_cliques, cliques[:])
    else:
        for i in range(1, v + 2):
            cliques[v] = i
            if is_partial_solution(cliques, adj_mat, v):
                best = light_backtrack(adj_mat, cliques, v + 1, best)

    return best

def is_partial_solution(cliques, adj_mat, v):
    for i in range(v):
        if cliques[i] == cliques[v] and adj_mat[i, v] == 0:
            return False
    return True

def is_solution(cliques, adj_mat):
    n = len(cliques)
    for i in range(n):
        for j in range(i + 1, n):
            if cliques[i] == cliques[j] and adj_mat[i, j] == 0:
                return False
    return True

def componentes_clique(grafo):
    ids = list(grafo.keys())
    n = len(ids)

    adj_mat = np.zeros((n, n), dtype=int)
    id_to_index = {id_: idx for idx, id_ in enumerate(ids)}

    for id1, vecinos in grafo.items():
        for id2 in vecinos:
            adj_mat[id_to_index[id1], id_to_index[id2]] = 1
            adj_mat[id_to_index[id2], id_to_index[id1]] = 1

    cliques = [0] * n
    _, assignment = light_backtrack(adj_mat, cliques)

    clique_assignment = {ids[i]: assignment[i] for i in range(n)}
    return clique_assignment

def resolver_caso(caso):
    n, d = caso["n"], caso["d"]
    celulas = []
    for entrada in caso["celulas"]:
        id_celula = entrada[0]
        x, y = entrada[1], entrada[2]
        peptidos = set(entrada[3:])
        celulas.append((id_celula, x, y, peptidos))

    grafo = construir_grafo(celulas, d)
    resultado = componentes_clique(grafo)

    return resultado

def main():
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
        {
            "n": 4,
            "d": 1,
            "celulas": [
                [1, 0, 0, "AETQT", "DFTYA"],
                [2, 0, 1, "AETQT", "HGCYS"],
                [3, 1, 0, "DFTYA", "IYHLK"],
                [4, 1, 1, "HGCYS", "IYHLK"],

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
