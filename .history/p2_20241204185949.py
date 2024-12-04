from collections import defaultdict, deque
import numpy as np
import time
import math
class Graph:
    def __init__(self):
        self.graph = defaultdict(list)
        self.capacity = defaultdict(int)
        self.flow_per_node = defaultdict(int)
        self.level = {}
        self.ptr = {}

    def add_edge(self, u, v, cap):
        self.graph[u].append(v)
        self.graph[v].append(u)
        self.capacity[(u, v)] += cap
        self.capacity[(v, u)] += 0

    def bfs(self, source, sink, blocked):
        self.level = {source: 0}
        queue = deque([source])

        while queue:
            u = queue.popleft()
            for v in self.graph[u]:
                if v not in self.level and self.capacity[(u, v)] > 0 and u != blocked and v != blocked:
                    self.level[v] = self.level[u] + 1
                    queue.append(v)
                    if v == sink:
                        return True
        return False

    def dfs(self, u, sink, flow, blocked, original_caps):
        if u == sink:
            return flow
        while self.ptr[u] < len(self.graph[u]):
            v = self.graph[u][self.ptr[u]]
            if self.level.get(v, -1) == self.level[u] + 1 and self.capacity[(u, v)] > 0 and u != blocked and v != blocked:
                cur_flow = min(flow, self.capacity[(u, v)])
                result = self.dfs(v, sink, cur_flow, blocked, original_caps)
                if result > 0:
                    # Restore original capacity if blocking node is removed
                    if (u, v) not in original_caps:
                        original_caps[(u, v)] = self.capacity[(u, v)]
                        original_caps[(v, u)] = self.capacity[(v, u)]

                    self.capacity[(u, v)] -= result
                    self.capacity[(v, u)] += result
                    self.flow_per_node[u] += result
                    self.flow_per_node[v] += result
                    return result
            self.ptr[u] += 1
        return 0

    def dinic(self, source, sink, blocked=None):
        max_flow = 0
        original_caps = {}

        while self.bfs(source, sink, blocked):
            self.ptr = {u: 0 for u in self.graph}
            while True:
                flow = self.dfs(source, sink, float('Inf'), blocked, original_caps)
                if flow == 0:
                    break
                max_flow += flow

        for u, v in original_caps:
            self.capacity[(u, v)] = original_caps[(u, v)]
        return max_flow

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

def calculate_cells(cases):
    results = []
    for case in cases:
        n, d, cells = case
        graph = Graph()
        source = 'source'
        sink = 'sink'

        calculators = []
        max_capacity = float("inf")
        peptides_cache = {}
        coords = []
        id_map = {}

        for i in range(n):
            id1, x1, y1, type1, *peptides1 = cells[i]
            coords.append((x1, y1))
            id_map[(x1, y1)] = (id1, type1, peptides1)

            if type1 == 1:
                graph.add_edge(source, id1, max_capacity)
            elif type1 == 2:
                calculators.append(id1)
            elif type1 == 3:
                graph.add_edge(id1, sink, max_capacity)

        neighbors = grid_neighbors(coords, d)

        for i in range(n):
            id1, x1, y1, type1, *peptides1 = cells[i]
            peptides_cache[id1] = set(peptides1)

            for j in neighbors[i]:
                id2, x2, y2, type2, *peptides2 = cells[j]

                if id2 not in peptides_cache:
                    peptides_cache[id2] = set(peptides2)

                shared_peptides = len(peptides_cache[id1] & peptides_cache[id2])
                if shared_peptides > 0:
                    if type1 == 1 and type2 == 2:
                        graph.add_edge(id1, id2, shared_peptides)
                    elif type1 == 2 and type2 == 2:
                        graph.add_edge(id1, id2, shared_peptides)
                    elif type1 == 2 and type2 == 3:
                        graph.add_edge(id1, id2, shared_peptides)

        total_flow = graph.dinic(source, sink)

        flow_values = [graph.flow_per_node[calc] for calc in calculators]
        threshold = np.percentile(flow_values, 80)
        top_calculators = [calc for calc in calculators if graph.flow_per_node[calc] >= threshold]

        max_reduction = total_flow
        blocked_cell = -1

        for calc in top_calculators:
            reduced_flow = graph.dinic(source, sink, blocked=calc)
            if reduced_flow < max_reduction:
                max_reduction = reduced_flow
                blocked_cell = calc
            elif reduced_flow == max_reduction and calc > blocked_cell:
                max_reduction = reduced_flow
                blocked_cell = calc

        results.append((blocked_cell, total_flow, max_reduction))

    return results

if __name__ == "__main__":
    cases = []
    t = int(input().strip())
    for _ in range(t):
        n, d = map(int, input().strip().split())
        cells = []
        for _ in range(n):
            data = input().strip().split()
            cell_id = int(data[0])
            x, y = int(data[1]), int(data[2])
            cell_type = int(data[3])
            peptides = data[4:]
            cells.append((cell_id, x, y, cell_type, *peptides))
        cases.append((n, d, cells))

    start_time = time.time()
    results = calculate_cells(cases)
    end_time = time.time()

    print(f"Execution time: {end_time - start_time:.6f} seconds")

    for result in results:
        print(f"{result[0]} {result[1]} {result[2]}")
