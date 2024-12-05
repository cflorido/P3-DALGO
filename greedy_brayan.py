import sys
from collections import defaultdict
import math

def grid_neighbors(coords, d):
    grid = defaultdict(list)
    cell_size = d

    for i, (x, y) in enumerate(coords):
        cell_x, cell_y = int(x // cell_size), int(y // cell_size)
        grid[(cell_x, cell_y)].append(i)

    neighbors = defaultdict(list)

    for i, (x, y) in enumerate(coords):
        cell_x, cell_y = int(x // cell_size), int(y // cell_size)

        for dx in range(-1, 2):
            for dy in range(-1, 2):
                neighbor_cell = (cell_x + dx, cell_y + dy)
                for j in grid.get(neighbor_cell, []):
                    if i != j:
                        dist = math.sqrt((coords[i][0] - coords[j][0]) ** 2 + (coords[i][1] - coords[j][1]) ** 2)
                        if dist <= d:
                            neighbors[i].append(j)

    return neighbors

def construir_grafo(celulas, d):
    grafo = defaultdict(list)
    n = len(celulas)

    coords = [(x, y) for _, x, y, _ in celulas]
    neighbors = grid_neighbors(coords, d)

    for i in range(n):
        id1, _, _, peptidos1 = celulas[i]
        grafo[id1]=[]
        for j in neighbors[i]:
            id2, _, _, peptidos2 = celulas[j]
            if peptidos1 & peptidos2:
                grafo[id1].append(id2)

    return grafo

def find_approximate_cliques(adj_list, n):
    cliques = []
    visited = [False] * (n + 1)

    for node in range(1, n + 1):
        if not visited[node]:
            clique = set()
            clique.add(node)
            visited[node] = True

            for neighbor in adj_list[node]:
                if not visited[neighbor]:
                    can_join = True
                    for member in clique:
                        if neighbor not in adj_list[member]:
                            can_join = False
                            break
                    if can_join:
                        clique.add(neighbor)
                        visited[neighbor] = True
            cliques.append(clique)

    return cliques

def componentes_clique(grafo):
    ids = list(grafo.keys())
    adj_list = {id_: set(vecinos) for id_, vecinos in grafo.items()}
    cliques = find_approximate_cliques(adj_list, len(ids))

    clique_assignment = {}
    for i, clique in enumerate(cliques, start=1):
        for node in clique:
            clique_assignment[node] = i

    return clique_assignment

def resolver_caso(n, d, celulas):
    grafo = construir_grafo(celulas, d)
    resultado = componentes_clique(grafo)
    return resultado

def main():
    input = sys.stdin.read
    data = input().strip().split("\n")
    
    t = int(data[0])  # Número de casos de prueba
    index = 1
    results = []
    
    for _ in range(t):
        n, d = map(int, data[index].split())
        index += 1
        celulas = []
        
        for _ in range(n):
            line = data[index].split()
            index += 1
            id_celula = int(line[0])
            x, y = map(int, line[1:3])
            peptidos = set(line[3:])
            celulas.append((id_celula, x, y, peptidos))
        
        resultado = resolver_caso(n, d, celulas)
        
        for id_celula in sorted(resultado):
            results.append(f"{id_celula} {resultado[id_celula]}")
    
    sys.stdout.write("\n".join(results) + "\n")

if __name__ == "__main__":
    main()
