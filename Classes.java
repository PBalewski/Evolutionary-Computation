import java.util.*;
import java.lang.Math;

public class Classes {

    public static class Node {
        public int id;
        public int x, y, cost;

        public Node(int id, int x, int y, int cost) {
            this.id = id;
            this.x = x;
            this.y = y;
            this.cost = cost;
        }
    }

    public static class Instance {
        public List<Node> nodes;
        public int n;
        public int[][] distMat;

        public Instance(List<Node> nodes) {
            this.nodes = nodes;
            this.n = nodes.size();
            this.distMat = computeDistanceMatrix();
        }

        private int[][] computeDistanceMatrix() {
            int[][] matrix = new int[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    int dx = nodes.get(i).x - nodes.get(j).x;
                    int dy = nodes.get(i).y - nodes.get(j).y;
                    int dij = (int) Math.round(Math.sqrt(dx * dx + dy * dy));
                    matrix[i][j] = dij;
                    matrix[j][i] = dij;
                }
            }
            return matrix;
        }
    }

    public static class Pair<A, B> {
        public A first;
        public B second;

        public Pair(A first, B second) {
            this.first = first;
            this.second = second;
        }
    }

    public static class Move {
        public String type;      // "inter" or "intra_edge"
        public int A;        // inter: prev_node, intra_edge: A_node
        public int B;        // inter: next_node, intra_edge: B_node
        public int C;        // inter: selected_node, intra_edge: C_node
        public int D;        // inter: remaining_node, intra_edge: D_node
        public int delta;

        public Move(String type, int A, int B, int C, int D, int delta) {
            this.type = type;
            this.A = A;
            this.B = B;
            this.C = C;
            this.D = D;
            this.delta = delta;
        }
    }
}
