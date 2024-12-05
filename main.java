
// Framework: Java Standard Library
import java.util.*;
import java.util.stream.Collectors;

public class Main {
    public static Map<Integer, List<Integer>> gridNeighbors(List<Tuple> coordinates, double distance) {
        /*
         * Función optimizada para encontrar vecinos dentro de una distancia `d` usando
         * una estructura de cuadrícula.
         */
        Map<Tuple, List<Integer>> grid = new HashMap<>();
        double cellSize = distance; // Tamaño de la celda en función de la distancia d

        // Asigna cada coordenada a una celda en la cuadrícula
        for (int i = 0; i < coordinates.size(); i++) {
            Tuple coordinate = coordinates.get(i);
            int cellX = (int) (coordinate.x / cellSize);
            int cellY = (int) (coordinate.y / cellSize);
            Tuple cell = new Tuple(cellX, cellY);
            grid.computeIfAbsent(cell, k -> new ArrayList<>()).add(i);
        }

        Map<Integer, List<Integer>> neighbors = new HashMap<>();

        // Para cada célula, buscar vecinos en celdas adyacentes
        for (int i = 0; i < coordinates.size(); i++) {
            Tuple coordinate = coordinates.get(i);
            int cellX = (int) (coordinate.x / cellSize);
            int cellY = (int) (coordinate.y / cellSize);

            // Recorrer solo las celdas adyacentes y la propia celda
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    Tuple neighborCell = new Tuple(cellX + dx, cellY + dy);

                    // Verifica cada punto en la celda vecina
                    for (int j : grid.getOrDefault(neighborCell, new ArrayList<>())) {
                        if (i != j) {
                            // Verifica la distancia real entre puntos para confirmar que está en el rango
                            double dist = Math.sqrt(Math.pow(coordinates.get(i).x - coordinates.get(j).x, 2) +
                                    Math.pow(coordinates.get(i).y - coordinates.get(j).y, 2));
                            if (dist <= distance) {
                                neighbors.computeIfAbsent(i, k -> new ArrayList<>()).add(j);
                            }
                        }
                    }
                }
            }
        }

        return neighbors;
    }

    public static Map<Integer, List<Integer>> buildGraph(List<Tuple> cells, double distance) {
        Map<Integer, List<Integer>> graph = new HashMap<>();
        int n = cells.size();

        // Extraemos las coordenadas de las células para pasarlas a la función
        // gridNeighbors
        List<Tuple> coordinates = new ArrayList<>();
        for (Tuple cell : cells) {
            coordinates.add(new Tuple(cell.id, cell.x, cell.y));
        }

        // Usamos la función gridNeighbors para encontrar los vecinos de las células
        Map<Integer, List<Integer>> neighbors = gridNeighbors(coordinates, distance);

        for (int i = 0; i < n; i++) {
            int id1 = cells.get(i).id;
            Set<Integer> peptides1 = cells.get(i).peptides;
            for (int j : neighbors.getOrDefault(i, new ArrayList<>())) {
                int id2 = cells.get(j).id;
                Set<Integer> peptides2 = cells.get(j).peptides;
                if (!Collections.disjoint(peptides1, peptides2)) { // Si tienen péptidos en común
                    graph.computeIfAbsent(id1, k -> new ArrayList<>()).add(id2);
                }
            }
        }

        return graph;
    }

    public static Map<Integer, Integer> greedyCliqueComponents(Map<Integer, List<Integer>> graph) {
        /*
         * Asigna greedy las células a cliques minimizando el número de cliques
         * necesarios.
         */
        Map<Integer, Integer> cliques = new HashMap<>();

        for (Integer cell : graph.keySet().stream().sorted().collect(Collectors.toList())) {
            // Encuentra los cliques vecinos ya asignados
            Set<Integer> neighborCliques = new HashSet<>();
            for (Integer neighbor : graph.get(cell)) {
                if (cliques.containsKey(neighbor)) {
                    neighborCliques.add(cliques.get(neighbor));
                }
            }

            // Asigna el menor número de clique disponible
            int assignedClique = 1;
            while (neighborCliques.contains(assignedClique)) {
                assignedClique++;
            }
            cliques.put(cell, assignedClique);
        }

        return cliques;
    }

    public static Map<Integer, Integer> resolveCase(Map<String, Object> caseData) {
        /*
         * Procesa un caso específico.
         */
        int n = (int) caseData.get("n");
        double distance = (double) caseData.get("d");
        List<Tuple> cells = new ArrayList<>();
        List<List<Object>> cellEntries = (List<List<Object>>) caseData.get("cells");
        for (List<Object> entry : cellEntries) {
            int cellId = (int) entry.get(0);
            double x = (double) entry.get(1);
            double y = (double) entry.get(2);
            Set<Integer> peptides = new HashSet<>();
            for (int i = 3; i < entry.size(); i++) {
                peptides.add((int) entry.get(i));
            }
            cells.add(new Tuple(cellId, x, y, peptides));
        }

        // Construir grafo basado en las restricciones
        Map<Integer, List<Integer>> graph = buildGraph(cells, distance);
        Map<Integer, Integer> result = greedyCliqueComponents(graph);

        return result;
    }

    public static List<Integer> findClique(Map<Integer, List<Integer>> graph, Set<Integer> available) {
        /*
         * Encuentra un clique máximo posible usando una aproximación greedy.
         */
        List<Integer> clique = new ArrayList<>();
        while (!available.isEmpty()) {
            Integer vertex = available.iterator().next();
            available.remove(vertex);
            clique.add(vertex);
            available.removeIf(u -> !clique.stream().allMatch(n -> graph.get(n).contains(u)));
        }

        return clique;
    }

    public static List<List<Integer>> minimumCliqueCover(Map<Integer, List<Integer>> graph) {
        /*
         * Encuentra una aproximación greedy para el Minimum Clique Cover.
         */
        Set<Integer> covered = new HashSet<>();
        Set<Integer> nodes = new HashSet<>(graph.keySet());
        List<List<Integer>> cliques = new ArrayList<>();

        while (!nodes.equals(covered)) {
            Set<Integer> available = new HashSet<>(nodes);
            available.removeAll(covered);

            List<Integer> clique = findClique(graph, available);
            cliques.add(clique);
            covered.addAll(clique);
        }

        return cliques;
    }

    public static List<List<Integer>> connectedComponents(Map<Integer, List<Integer>> graph) {
        /*
         * Encuentra las componentes conexas del grafo usando DFS iterativo.
         */
        Set<Integer> visited = new HashSet<>();
        List<List<Integer>> components = new ArrayList<>();

        for (Integer node : graph.keySet()) {
            if (!visited.contains(node)) {
                List<Integer> component = new ArrayList<>();
                Stack<Integer> stack = new Stack<>();
                stack.push(node);
                while (!stack.isEmpty()) {
                    Integer current = stack.pop();
                    if (!visited.contains(current)) {
                        visited.add(current);
                        component.add(current);
                        for (Integer neighbor : graph.get(current)) {
                            if (!visited.contains(neighbor)) {
                                stack.push(neighbor);
                            }
                        }
                    }
                }
                components.add(component);
            }
        }

        return components;
    }

    public static Map<Integer, Integer> resolveCase(int n, double distance, List<List<Object>> cells) {
        List<Tuple> cellList = new ArrayList<>();
        for (List<Object> data : cells) {
            int cellId = (int) data.get(0);
            double x = (double) data.get(1);
            double y = (double) data.get(2);
            Set<Integer> peptides = new HashSet<>();
            for (int i = 4; i < data.size(); i++) {
                peptides.add((int) data.get(i));
            }
            cellList.add(new Tuple(cellId, x, y, peptides));
        }
        Map<Integer, List<Integer>> graph = buildGraph(cellList, distance);
        List<List<Integer>> components = connectedComponents(graph);
        Map<Integer, Integer> result = new HashMap<>();
        int cliqueId = 1;

        for (List<Integer> component : components) {
            Map<Integer, List<Integer>> subgraph = new HashMap<>();
            for (Integer node : component) {
                subgraph.put(node, new ArrayList<>());
                for (Integer neighbor : graph.get(node)) {
                    if (component.contains(neighbor)) {
                        subgraph.get(node).add(neighbor);
                    }
                }
            }
            List<List<Integer>> cliques = minimumCliqueCover(subgraph);

            for (List<Integer> clique : cliques) {
                for (Integer node : clique) {
                    result.put(node, cliqueId);
                }
                cliqueId++;
            }
        }
        return result;
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int t = Integer.parseInt(scanner.nextLine()); // Número de casos
        List<Map<Integer, Integer>> results = new ArrayList<>();

        for (int caseIndex = 0; caseIndex < t; caseIndex++) {
            String[] caseParams = scanner.nextLine().split(" ");
            int n = Integer.parseInt(caseParams[0]);
            double d = Double.parseDouble(caseParams[1]);
            List<List<Object>> cells = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                String[] cellData = scanner.nextLine().split(" ");
                List<Object> cell = new ArrayList<>();
                for (String data : cellData) {
                    cell.add(Integer.parseInt(data));
                }
                cells.add(cell);
            }

            // Resolver cada caso
            Map<Integer, Integer> result = resolveCase(n, d, cells);
            results.add(result);
        }

        // Imprimir resultados
        for (Map<Integer, Integer> result : results) {
            for (Integer cellId : new TreeSet<>(result.keySet())) {
                System.out.println(cellId + " " + result.get(cellId));
            }
        }
    }
}

class Tuple {
    int id;
    double x;
    double y;
    Set<Integer> peptides;

    Tuple(int id, double x, double y) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.peptides = new HashSet<>();
    }

    Tuple(int id, double x, double y, Set<Integer> peptides) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.peptides = peptides;
    }

    Tuple(double x, double y) {
        this.id = 0; // Asigna un id por defecto
        this.x = x;
        this.y = y;
        this.peptides = new HashSet<>();
    }
}