import math
from collections import defaultdict
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

import math
import random
from collections import defaultdict

def construir_grafo(celulas, d):
    """
    Construye un grafo en el que las células son nodos y hay una arista entre dos células
    si están dentro de una distancia `d` y tienen al menos un péptido en común.
    """
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

def greedy_clique_cover(grafo):
    """
    Algoritmo greedy para encontrar una cobertura mínima de cliques en el grafo.
    """
    # Inicializamos los cliques
    cliques = []
    unassigned = set(grafo.keys())  # Conjunto de células no asignadas a ningún clique

    while unassigned:
        # Elegimos un nodo arbitrario de las células no asignadas
        v = unassigned.pop()

        # Creamos un nuevo clique con el nodo v
        new_clique = {v}
        to_remove = {v}

        # Buscamos los vecinos de v que pueden ser parte del mismo clique
        for u in list(new_clique):
            for neighbor in grafo[u]:
                if neighbor not in to_remove:
                    # Si el vecino puede unirse al clique, lo agregamos
                    if all(neighbor in grafo[n] for n in new_clique):
                        new_clique.add(neighbor)
                        to_remove.add(neighbor)

        # Añadimos el clique encontrado a la lista de cliques
        cliques.append(new_clique)

        # Actualizamos el conjunto de células no asignadas
        unassigned -= to_remove

    return cliques

def resolver_caso(caso):
    """
    Resuelve el caso dado construyendo el grafo y aplicando el algoritmo greedy.
    """
    n, d = caso["n"], caso["d"]
    celulas = []
    for entrada in caso["celulas"]:
        id_celula = entrada[0]
        x, y = entrada[1], entrada[2]
        peptidos = set(entrada[3:])
        celulas.append((id_celula, x, y, peptidos))

    # Construimos el grafo
    grafo = construir_grafo(celulas, d)

    # Aplicamos el algoritmo greedy para encontrar el clique cover
    cliques = greedy_clique_cover(grafo)

    # Asignamos un número de clique a cada célula
    resultado = {}
    for i, clique in enumerate(cliques, start=1):
        for celula in clique:
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
