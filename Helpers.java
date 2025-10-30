import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

public class Helpers {

    public static Classes.Instance readCSV(String filepath) throws IOException {
        List<Classes.Node> nodes = new ArrayList<>();
        List<String> lines = Files.readAllLines(Paths.get(filepath));

        for (int idx = 0; idx < lines.size(); idx++) {
            String[] parts = lines.get(idx).trim().split(";");
            if (parts.length != 3)
                throw new IllegalArgumentException("Invalid data format in line " + (idx + 1));

            int x = Integer.parseInt(parts[0]);
            int y = Integer.parseInt(parts[1]);
            int cost = Integer.parseInt(parts[2]);
            nodes.add(new Classes.Node(idx, x, y, cost));
        }
        return new Classes.Instance(nodes);
    }

    public static int evaluate(List<Integer> path, Classes.Instance inst, boolean formCycle) {
        if (path.size() < 2)
            return path.stream().mapToInt(i -> inst.nodes.get(i).cost).sum();

        int dist = 0;
        if (formCycle) {
            for (int i = 0; i < path.size(); i++) {
                int a = path.get(i);
                int b = path.get((i + 1) % path.size()); //includes edge between last and first node
                dist += inst.distMat[a][b];
            }
        } else {
            for (int i = 0; i < path.size() - 1; i++) {
                int a = path.get(i);
                int b = path.get(i + 1);
                dist += inst.distMat[a][b];
            }
        }

        int cost = path.stream().mapToInt(i -> inst.nodes.get(i).cost).sum();
        return dist + cost;
    }

    public static Set<Integer> getRemaining(int n, List<Integer> selected) {
        Set<Integer> remaining = IntStream.range(0, n).boxed().collect(Collectors.toSet());
        remaining.removeAll(selected);
        return remaining;
    }

    public static <T> void reverseSublist(List<T> list, int from, int to) {
        int i = from;
        int j = to;
        while (i < j) {
            Collections.swap(list, i, j);
            i++;
            j--;
        }
    }

    // ---- Initializations ----
    //Random 
    public static List<Integer> getRandomSolution(Classes.Instance inst, Random rand) {
        List<Integer> allNodes = IntStream.range(0, inst.n).boxed().collect(Collectors.toList());
        Collections.shuffle(allNodes, rand);
        int size = (inst.n + 1) /2;
        return new ArrayList<>(allNodes.subList(0, size));
    }

    // ---- Calculating deltas ----
    // Inter-route Node: switching node from the cycle with the one from outside
    public static int calculateDeltaInterNode(List<Integer> selected, int i_idx, int k_node, Classes.Instance inst) {
        int n = selected.size();

        int i_node = selected.get(i_idx);
        int prev_node = selected.get((i_idx - 1 + n) % n);
        int next_node = selected.get((i_idx + 1) % n);

        int remove_dist = inst.distMat[prev_node][i_node] + inst.distMat[i_node][next_node];
        int add_dist = inst.distMat[prev_node][k_node] + inst.distMat[k_node][next_node];

        int remove_cost = inst.nodes.get(i_node).cost;
        int add_cost = inst.nodes.get(k_node).cost;

        return (add_dist - remove_dist) + (add_cost - remove_cost);
    }

    //Intra-route Edge: switch 2 edges
    public static int calculateDeltaIntraEdge(List<Integer> selected, int i_idx, int j_idx, Classes.Instance inst) {
        int n = selected.size();
        
        if (i_idx == j_idx) {return 0;}

        // Ensure i_idx is smaller than j_idx
        if (i_idx > j_idx) {
            int temp = i_idx;
            i_idx = j_idx;
            j_idx = temp;
        }

        // Skip adjecent edges
        if (j_idx == i_idx + 1 || (i_idx == 0 && j_idx == n - 1)) {
            return 0;
        }

        int node_i = selected.get(i_idx);
        int next_i = selected.get(i_idx + 1);
        int node_j = selected.get(j_idx);
        int next_j = selected.get((j_idx + 1) % n);

        int remove_dist = inst.distMat[node_i][next_i] + inst.distMat[node_j][next_j];
        int add_dist = inst.distMat[node_i][node_j] + inst.distMat[next_i][next_j];

        return add_dist - remove_dist;
    }
}
