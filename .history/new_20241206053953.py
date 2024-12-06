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

def clique_aproximation(grafo):
    processed_nodes = set()
    groups = []
    adjacency_list = {node: set(neighbors) for node, neighbors in grafo.items()}

    # Ordenar los nodos por el grado (número de vecinos) de manera descendente
    sorted_nodes = sorted(adjacency_list, key=lambda x: len(adjacency_list[x]), reverse=True)

    def expand_clique(group):
        """ Expande el clique agregando nodos que sean adyacentes a todos los miembros del clique """
        potential_members = set(adjacency_list[group[-1]])  # Comienza con los vecinos del último miembro
        for node in group:
            potential_members &= adjacency_list[node]  # Solo vecinos comunes

        # Solo agregamos nodos que son adyacentes a todos los miembros del clique
        return group + [node for node in potential_members if node not in group]

    for current_node in sorted_nodes:
        if current_node not in processed_nodes:
            group = [current_node]
            processed_nodes.add(current_node)

            # Expandimos el clique mientras se pueda
            expanded_group = expand_clique(group)
            groups.append(expanded_group)

            # Marcamos todos los miembros del grupo como procesados
            processed_nodes.update(expanded_group)

    # Asignar un grupo único a cada nodo
    return {
        node: group_index + 1
        for group_index, group in enumerate(groups)
        for node in group
    }



def resolver_caso(n, d, celulas):
    grafo = construir_grafo(celulas, d)
    resultado = clique_aproximation(grafo)
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
