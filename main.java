import java.util.*;
import java.io.*;

public class Main {
    static class Cell {
        int id;
        int x;
        int y;
        Set<String> peptides;

        Cell(int id, int x, int y, Set<String> peptides) {
            this.id = id;
            this.x = x;
            this.y = y;
            this.peptides = peptides;
        }
    }

    static Map<Integer, List<Integer>> gridNeighbors(List<Cell> cells, double d) {
        Map<Pair<Integer, Integer>, List<Integer>> grid = new HashMap<>();
        double cellSize = d;

        for (int i = 0; i < cells.size(); i++) {
            Cell cell = cells.get(i);
            int cellX = (int) (cell.x / cellSize);
            int cellY = (int) (cell.y / cellSize);
            grid.computeIfAbsent(new Pair<>(cellX, cellY), k -> new ArrayList<>()).add(i);
        }

        Map<Integer, List<Integer>> neighbors = new HashMap<>();

        for (int i = 0; i < cells.size(); i++) {
            Cell cell = cells.get(i);
            int cellX = (int) (cell.x / cellSize);
            int cellY = (int) (cell.y / cellSize);

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    Pair<Integer, Integer> neighborCell = new Pair<>(cellX + dx, cellY + dy);
                    for (int j : grid.getOrDefault(neighborCell, Collections.emptyList())) {
                        if (i != j) {
                            double dist = Math.sqrt(Math.pow(cells.get(i).x - cells.get(j).x, 2)
                                    + Math.pow(cells.get(i).y - cells.get(j).y, 2));
                            if (dist <= d) {
                                neighbors.computeIfAbsent(i, k -> new ArrayList<>()).add(j);
                            }
                        }
                    }
                }
            }
        }

        return neighbors;
    }

    static Map<Integer, List<Integer>> construirGrafo(List<Cell> cells, double d) {
        Map<Integer, List<Integer>> graph = new HashMap<>();
        int n = cells.size();

        Map<Integer, List<Integer>> neighbors = gridNeighbors(cells, d);

        for (int i = 0; i < n; i++) {
            Cell cell1 = cells.get(i);
            graph.put(cell1.id, new ArrayList<>());
            for (int j : neighbors.getOrDefault(i, Collections.emptyList())) {
                Cell cell2 = cells.get(j);
                if (!Collections.disjoint(cell1.peptides, cell2.peptides)) {
                    graph.get(cell1.id).add(cell2.id);
                }
            }
        }

        return graph;
    }

    static List<Set<Integer>> findApproximateCliques(Map<Integer, Set<Integer>> adjList, int n) {
        List<Set<Integer>> cliques = new ArrayList<>();
        boolean[] visited = new boolean[n + 1];

        for (int node = 1; node <= n; node++) {
            if (!visited[node]) {
                Set<Integer> clique = new HashSet<>();
                clique.add(node);
                visited[node] = true;

                for (int neighbor : adjList.getOrDefault(node, Collections.emptySet())) {
                    if (!visited[neighbor]) {
                        boolean canJoin = true;
                        for (int member : clique) {
                            if (!adjList.getOrDefault(member, Collections.emptySet()).contains(neighbor)) {
                                canJoin = false;
                                break;
                            }
                        }
                        if (canJoin) {
                            clique.add(neighbor);
                            visited[neighbor] = true;
                        }
                    }
                }
                cliques.add(clique);
            }
        }

        return cliques;
    }

    static Map<Integer, Integer> componentesClique(Map<Integer, List<Integer>> graph) {
        List<Integer> ids = new ArrayList<>(graph.keySet());
        Map<Integer, Set<Integer>> adjList = new HashMap<>();
        for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
            adjList.put(entry.getKey(), new HashSet<>(entry.getValue()));
        }
        List<Set<Integer>> cliques = findApproximateCliques(adjList, ids.size());

        Map<Integer, Integer> cliqueAssignment = new HashMap<>();
        for (int i = 0; i < cliques.size(); i++) {
            for (int node : cliques.get(i)) {
                cliqueAssignment.put(node, i + 1);
            }
        }

        return cliqueAssignment;
    }

    static Map<Integer, Integer> resolverCaso(int n, double d, List<Cell> cells) {
        Map<Integer, List<Integer>> graph = construirGrafo(cells, d);
        return componentesClique(graph);
    }

    public static void main(String[] args) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        int t = Integer.parseInt(br.readLine().trim());
        List<String> results = new ArrayList<>();

        for (int i = 0; i < t; i++) {
            String[] line = br.readLine().trim().split(" ");
            int n = Integer.parseInt(line[0]);
            double d = Double.parseDouble(line[1]);
            List<Cell> cells = new ArrayList<>();

            for (int j = 0; j < n; j++) {
                line = br.readLine().trim().split(" ");
                int idCelula = Integer.parseInt(line[0]);
                int x = Integer.parseInt(line[1]);
                int y = Integer.parseInt(line[2]);
                Set<String> peptides = new HashSet<>(Arrays.asList(Arrays.copyOfRange(line, 3, line.length)));
                cells.add(new Cell(idCelula, x, y, peptides));
            }

            Map<Integer, Integer> resultado = resolverCaso(n, d, cells);

            for (int idCelula : new TreeSet<>(resultado.keySet())) {
                results.add(idCelula + " " + resultado.get(idCelula));
            }
        }

        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
        bw.write(String.join("\n", results));
        bw.newLine();
        bw.flush();
    }

    static class Pair<K, V> {
        K first;
        V second;

        Pair(K first, V second) {
            this.first = first;
            this.second = second;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o)
                return true;
            if (o == null || getClass() != o.getClass())
                return false;
            Pair<?, ?> pair = (Pair<?, ?>) o;
            return Objects.equals(first, pair.first) && Objects.equals(second, pair.second);
        }

        @Override
        public int hashCode() {
            return Objects.hash(first, second);
        }
    }
}
