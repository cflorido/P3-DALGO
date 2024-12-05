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

def greedy_componentes_clique(grafo):
    """
    Asigna greedy las células a cliques minimizando el número de cliques necesarios.
    """
    cliques = {}
    
    for celula in sorted(grafo.keys()):
        # Encuentra los cliques vecinos ya asignados
        cliques_vecinos = {cliques[vecino] for vecino in grafo[celula] if vecino in cliques}
        
        # Asigna el menor número de clique disponible
        clique_asignado = 1
        while clique_asignado in cliques_vecinos:
            clique_asignado += 1
        cliques[celula] = clique_asignado
    
    return cliques

def resolver_caso(caso):
    """
    Procesa un caso específico.
    """
    n, d = caso["n"], caso["d"]
    celulas = []
    for entrada in caso["celulas"]:
        id_celula = entrada[0]
        x, y = entrada[1], entrada[2]
        peptidos = set(entrada[3:])
        celulas.append((id_celula, x, y, peptidos))

    # Construir grafo basado en las restricciones
    grafo = construir_grafo(celulas, d)

    resultado = greedy_componentes_clique(grafo)

    return resultado

def encontrar_clique(grafo, disponibles):
    """
    Encuentra un clique máximo posible usando una aproximación greedy.
    """
    clique = []
    while disponibles:

        v = disponibles.pop()
        clique.append(v)
        disponibles = {u for u in disponibles if all(u in grafo[n] for n in clique)}
    
    return clique


def minimum_clique_cover(grafo):
    """
    Encuentra una aproximación greedy para el Minimum Clique Cover.
    """
    cubiertos = set()
    nodos = set(grafo.keys())
    cliques = []

    while nodos != cubiertos:
        disponibles = nodos - cubiertos
    
        clique = encontrar_clique(grafo, disponibles.copy())
        cliques.append(clique)

      
        cubiertos.update(clique)

    return cliques


def componentes_conexas(grafo):
    """
    Encuentra las componentes conexas del grafo usando DFS iterativo.
    """
    visitados = set()
    componentes = []

    for nodo in grafo:
        if nodo not in visitados:
            componente = []
        
            pila = [nodo]
            while pila:
                actual = pila.pop()
                if actual not in visitados:
                    visitados.add(actual)
                    componente.append(actual)
                 
                    for vecino in grafo[actual]:
                        if vecino not in visitados:
                            pila.append(vecino)
            componentes.append(componente)

    return componentes



# Integrar todo en el flujo principal
def resolver_caso(n, d, celulas):
    celulas = [(int(data[0]), int(data[1]), int(data[2]), set(data[4:])) for data in celulas]
    grafo = construir_grafo(celulas, d)
    componentes = componentes_conexas(grafo)
    resultado = {}
    clique_id = 1

    for componente in componentes:
        subgrafo = {nodo: [vec for vec in grafo[nodo] if vec in componente] for nodo in componente}
        cliques = minimum_clique_cover(subgrafo)

        for clique in cliques:
            for nodo in clique:
                resultado[nodo] = clique_id
            clique_id += 1
    return resultado


def main():
    import sys
    input = sys.stdin.read
    data = input().splitlines()
    
    t = int(data[0])  # Número de casos
    index = 1
    resultados = []
    
    for _ in range(t):
        n, d = map(int, data[index].split())
        index += 1
        celulas = [data[index + i].split() for i in range(n)]
        index += n
        
        # Resolver cada caso
        resultado = resolver_caso(n, d, celulas)
        resultados.append(resultado)
    
    # Imprimir resultados
    for resultado in resultados:
        for id_celula in sorted(resultado):
            print(id_celula, resultado[id_celula])

if __name__ == "__main__":
    main()