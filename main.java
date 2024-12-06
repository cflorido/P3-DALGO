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

    static Map<Integer, List<Integer>> gridNeighbors(List<int[]> coords, double d) {
        Map<Pair<Integer, Integer>, List<Integer>> grid = new HashMap<>();
        double cellSize = d;

        for (int i = 0; i < coords.size(); i++) {
            int[] coord = coords.get(i);
            int cellX = (int) (coord[0] / cellSize);
            int cellY = (int) (coord[1] / cellSize);
            Pair<Integer, Integer> cell = new Pair<>(cellX, cellY);
            grid.computeIfAbsent(cell, k -> new ArrayList<>()).add(i);
        }

        Map<Integer, List<Integer>> neighbors = new HashMap<>();

        for (int i = 0; i < coords.size(); i++) {
            int[] coord = coords.get(i);
            int cellX = (int) (coord[0] / cellSize);
            int cellY = (int) (coord[1] / cellSize);

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    Pair<Integer, Integer> neighborCell = new Pair<>(cellX + dx, cellY + dy);
                    for (int j : grid.getOrDefault(neighborCell, Collections.emptyList())) {
                        if (i != j) {
                            double dist = Math.sqrt(Math.pow(coords.get(i)[0] - coords.get(j)[0], 2) +
                                    Math.pow(coords.get(i)[1] - coords.get(j)[1], 2));
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

        List<int[]> coords = new ArrayList<>();
        for (Cell cell : cells) {
            coords.add(new int[] { cell.x, cell.y });
        }
        Map<Integer, List<Integer>> neighbors = gridNeighbors(coords, d);

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

    static Map<Integer, Integer> cliqueApproximation(Map<Integer, List<Integer>> graph) {
        Set<Integer> processedNodes = new HashSet<>();
        List<Set<Integer>> groups = new ArrayList<>();
        Map<Integer, Set<Integer>> adjacencyList = new HashMap<>();

        for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
            adjacencyList.put(entry.getKey(), new HashSet<>(entry.getValue()));
        }

        for (int currentNode : adjacencyList.keySet()) {
            if (!processedNodes.contains(currentNode)) {
                Set<Integer> group = new HashSet<>();
                group.add(currentNode);
                processedNodes.add(currentNode);

                List<Integer> potentialMembers = new ArrayList<>(adjacencyList.get(currentNode));
                potentialMembers.sort((a, b) -> Integer.compare(
                        Sets.intersection(adjacencyList.get(b), group).size(),
                        Sets.intersection(adjacencyList.get(a), group).size()));

                for (int neighbor : potentialMembers) {
                    if (!processedNodes.contains(neighbor)) {
                        if (adjacencyList.get(neighbor).containsAll(group)) {
                            group.add(neighbor);
                            processedNodes.add(neighbor);
                        }
                    }
                }

                groups.add(group);
            }
        }

        Map<Integer, Integer> result = new HashMap<>();
        for (int i = 0; i < groups.size(); i++) {
            for (int node : groups.get(i)) {
                result.put(node, i + 1);
            }
        }

        return result;
    }

    static Map<Integer, Integer> resolverCaso(int n, double d, List<Cell> cells) {
        Map<Integer, List<Integer>> graph = construirGrafo(cells, d);
        return cliqueApproximation(graph);
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
                int id = Integer.parseInt(line[0]);
                int x = Integer.parseInt(line[1]);
                int y = Integer.parseInt(line[2]);
                Set<String> peptides = new HashSet<>(Arrays.asList(Arrays.copyOfRange(line, 3, line.length)));
                cells.add(new Cell(id, x, y, peptides));
            }

            Map<Integer, Integer> resultado = resolverCaso(n, d, cells);

            for (int cellId : new TreeSet<>(resultado.keySet())) {
                results.add(cellId + " " + resultado.get(cellId));
            }
        }

        for (String result : results) {
            System.out.println(result);
        }
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

    static class Sets {
        static <T> Set<T> intersection(Set<T> set1, Set<T> set2) {
            Set<T> result = new HashSet<>(set1);
            result.retainAll(set2);
            return result;
        }
    }
}
