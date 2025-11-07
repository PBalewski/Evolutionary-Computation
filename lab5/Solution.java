import java.util.*;
import java.util.stream.*;

public class Solution {
    
    private static List<Classes.Move> steepestSearch_LM = new ArrayList<>();
    // ---- Initialization ----
    //Random 
    public static List<Integer> getRandomSolution(Classes.Instance inst, Random rand) {
        List<Integer> allNodes = IntStream.range(0, inst.n).boxed().collect(Collectors.toList());
        Collections.shuffle(allNodes, rand);
        int size = (inst.n + 1) /2;
        return new ArrayList<>(allNodes.subList(0, size));
    }


    // ---- Calculating deltas ----
    // Inter-route Node: switching node from the cycle with the one from outside
    private static int calculateDeltaInterNode(
                List<Integer> selected, int node_i, int node_k, Classes.Instance inst) {

        int i_idx = selected.indexOf(node_i);
        if (i_idx == -1) return Integer.MAX_VALUE; // node not in cycle

        int n = selected.size();
        int prev_node = selected.get((i_idx - 1 + n) % n);
        int next_node = selected.get((i_idx + 1) % n);

        int remove_dist = inst.distMat[prev_node][node_i] + inst.distMat[node_i][next_node];
        int add_dist = inst.distMat[prev_node][node_k] + inst.distMat[node_k][next_node];

        int remove_cost = inst.nodes.get(node_i).cost;
        int add_cost = inst.nodes.get(node_k).cost;

        return (add_dist - remove_dist) + (add_cost - remove_cost);
    }

    
    //Intra-route Edge: switch 2 edges
    private static int calculateDeltaIntraEdge(
                List<Integer> selected, int node_a, int node_c, Classes.Instance inst) {

        int i_idx = selected.indexOf(node_a);
        int j_idx = selected.indexOf(node_c);
        if (i_idx == -1 || j_idx == -1) return Integer.MAX_VALUE;

        int n = selected.size();
        // Ensure i_idx is smaller than j_idx for consistent indexing relative to the list
        if (i_idx > j_idx) {
            int tmp = i_idx; i_idx = j_idx; j_idx = tmp;
        }

        // Skip adjacent edges or same nodes
        if (i_idx == j_idx || j_idx == i_idx + 1 || (i_idx == 0 && j_idx == n - 1)) return 0;

        int prev_a = selected.get((i_idx - 1 + n) % n); // Node before A
        int next_a = selected.get((i_idx + 1) % n);    // Node after A (B)
        int next_c = selected.get((j_idx + 1) % n);    // Node after C (D)

        int remove_dist = inst.distMat[node_a][next_a] + inst.distMat[node_c][next_c]; // Edges (A,B) and (C,D)
        int add_dist = inst.distMat[node_a][node_c] + inst.distMat[next_a][next_c]; // Edges (A,C) and (B,D)

        return add_dist - remove_dist;
    }


    // ---- Local Search ----

    public static Classes.Pair<List<Integer>, Integer> localSearch(
                List<Integer> initial_solution,
                Classes.Instance inst,
                String ls_type,
                String move_type,
                Random rand)
            {
        List<Integer> current_solution = new ArrayList<>(initial_solution);
        Set<Integer> remaining = Helpers.getRemaining(inst.n, current_solution);

        boolean improved = true;
        // int iteration = 0; // Unused variable
        
        steepestSearch_LM.clear();
        while (improved) {
            // iteration++; // Unused increment
            improved = steepestSearch(current_solution, remaining, inst, move_type);
        }
        int final_cost = Helpers.evaluate(current_solution, inst, true);
        return new Classes.Pair<>(current_solution, final_cost);
    }


    // steepest local search with reuse of previous deltas (LM)
    private static boolean steepestSearch(List<Integer> selected, Set<Integer> remaining, Classes.Instance inst, String move_type) {
        List<Classes.Move> LM = steepestSearch_LM;
        int n = selected.size();

        if (LM.isEmpty()) {
            // inter-route moves
            for (int i_idx = 0; i_idx < n; i_idx++) { // FIX 1: Iterate over indices
                int i_node = selected.get(i_idx);
                for (int k_node : remaining) {
                    int delta = calculateDeltaInterNode(selected, i_node, k_node, inst);
                    if (delta < 0) {
                        int prev_node = selected.get((i_idx - 1 + n) % n); // Use i_idx
                        int next_node = selected.get((i_idx + 1) % n);     // Use i_idx
                        // A=prev, B=next, C=i_node, D=k_node
                        LM.add(new Classes.Move("inter", prev_node, next_node, i_node, k_node, delta));
                    }
                }
            }

            // intra-route moves
            for (int i_idx = 0; i_idx < n; i_idx++) {
                for (int j_idx = i_idx + 1; j_idx < n; j_idx++) {
                    int a = selected.get(i_idx); // Start of first edge (A)
                    int next_a = selected.get((i_idx + 1) % n); // End of first edge (B)
                    int c = selected.get(j_idx); // Start of second edge (C)
                    int next_c = selected.get((j_idx + 1) % n); // End of second edge (D)
                    
                    if (next_c == a) continue; // Avoid adjacent edges

                    int delta = calculateDeltaIntraEdge(selected, a, c, inst);
                    if (delta < 0) {
                        // A=A, B=nextA, C=C, D=nextC
                        LM.add(new Classes.Move("intra_edge", a, next_a, c, next_c, delta));
                    }
                }
            }
        }

        Map<Integer,Integer> succMap = Helpers.buildSuccMap(selected);
        // int[] posArr = Helpers.buildPosArray(selected, inst.n); // posArr not used here

        // ---- 1. Recheck moves in LM ----
        LM.sort(Comparator.comparingInt(m -> m.delta)); // best first
        Iterator<Classes.Move> it = LM.iterator();
        while (it.hasNext()) {
            Classes.Move m = it.next();

            // validate move
            Classes.Pair<Boolean, Boolean> val = validate(succMap, remaining, m);
            boolean applyMove = val.first;
            boolean keepMove = val.second;

            if (!keepMove) {
                it.remove();
                continue;
            }
            
            if (!applyMove) {
                continue;
            }

            applyMove(selected, remaining, m);
            it.remove();
            succMap = Helpers.buildSuccMap(selected);
            // posArr = Helpers.buildPosArray(selected, inst.n); // posArr not used here
            
            // identify changed nodes
            // A=prev, B=next, C=i_node, D=k_node (inter)
            // A=A, B=nextA, C=C, D=nextC (intra)
            Set<Integer> changed = new HashSet<>(Arrays.asList(
                        m.A, m.B, m.C, m.D));

            // add new local improving moves
            addNewLocalMovesToLM(LM, selected, remaining, inst, changed, m.type);
            steepestSearch_LM = LM;
            return true; // Found and applied the best move
        }
        return false; // No improving move found
    }


    // check if a move is still valid
    private static Classes.Pair<Boolean, Boolean> validate(Map<Integer,Integer> succMap, Set<Integer> remaining, Classes.Move m) {

        switch (m.type) {
            case "inter":
                boolean acb = Helpers.hasEdge(succMap, m.A, m.C) && Helpers.hasEdge(succMap, m.C, m.B);
                boolean bca = Helpers.hasEdge(succMap, m.B, m.C) && Helpers.hasEdge(succMap, m.C, m.A);
                boolean d = remaining.contains(m.D);
                if ((acb || bca) && d) {
                    return new Classes.Pair<>(true, true); 
                }
                return new Classes.Pair<>(false, false);

            case "intra_edge":
                boolean ab = Helpers.hasEdge(succMap, m.A, m.B);
                boolean ba = Helpers.hasEdge(succMap, m.A, m.B);
                boolean cd = Helpers.hasEdge(succMap, m.C, m.D);
                boolean dc = Helpers.hasEdge(succMap, m.C, m.D);
                if ((!ab && !ba) || (!cd && !dc)) {
                    return new Classes.Pair<>(false, false);
                }
                if ((ab && cd) || (ba && dc)) {
                    return new Classes.Pair<>(true, true);
                }
                return new Classes.Pair<>(false, true);
        }
        return new Classes.Pair<>(false, false);
    }

    // actually apply a move to the current solution
    private static void applyMove(List<Integer> selected, Set<Integer> remaining, Classes.Move m) {
        switch (m.type) {
            case "inter": {
                // Swap C (in cycle) with D (remaining)
                int idx = selected.indexOf(m.C);
                if (idx >= 0) {
                    selected.set(idx, m.D);
                    remaining.remove(m.D);
                    remaining.add(m.C);
                }
                break;
            }
            case "intra_edge": {
                // m.A=A, m.C=C. We reverse the segment from B (index of A + 1) to C (index of C).
                int i_idx = selected.indexOf(m.A); // Index of A
                int j_idx = selected.indexOf(m.C); // Index of C

                if (i_idx >= 0 && j_idx >= 0) {
                    if (i_idx > j_idx) {                         
                        int tmp = i_idx; i_idx = j_idx; j_idx = tmp;
                    }
                    
                    // Segment to reverse: [i_idx + 1, j_idx] -> [B, C]
                    int startReverseIdx = i_idx + 1; // Index of B (next_a)
                    int endReverseIdx = j_idx;       // Index of C
                    
                    // Standard reverseSublist usually takes (start inclusive, end inclusive)
                    if (startReverseIdx < endReverseIdx) {
                        Helpers.reverseSublist(selected, startReverseIdx, endReverseIdx);
                    }
                }
                break;
            }
        }
    }

   // Finalized addNewLocalMovesToLM
    private static void addNewLocalMovesToLM(
        List<Classes.Move> LM,
        List<Integer> selected,
        Set<Integer> remaining,
        Classes.Instance inst,
        Set<Integer> changedNodes,
        String move_type)
    {
        int n = selected.size();
        
        // --- 1. Inter-route Re-evaluation (No change from last proposal, as it was sound) ---
        Set<Integer> indices_to_check_inter = new HashSet<>();
        for (int node_id : changedNodes) {
            int idx = selected.indexOf(node_id);
            if (idx != -1) {
                // Check position I and its neighbors
                indices_to_check_inter.add((idx - 1 + n) % n); 
                indices_to_check_inter.add(idx);               
                indices_to_check_inter.add((idx + 1) % n);     
            }
        }

        // A. Re-evaluate neighbors of affected positions with all remaining nodes (K)
        for (int i_idx : indices_to_check_inter) {
            int node_i = selected.get(i_idx);
            for (int node_k : remaining) {
                int delta = calculateDeltaInterNode(selected, node_i, node_k, inst);
                if (delta < 0) {
                    int prev = selected.get((i_idx - 1 + n) % n);
                    int next = selected.get((i_idx + 1) % n);
                    LM.add(new Classes.Move("inter", prev, next, node_i, node_k, delta));
                }
            }
        }
        
        // B. Re-evaluate all cycle nodes (I) with the newly available node (C), if it was an "inter" move
        if (move_type.equals("inter")) {
            int node_c_old = changedNodes.stream()
                .filter(i -> remaining.contains(i))
                .findFirst().orElse(-1);

            if (node_c_old != -1) {
                for (int i_idx = 0; i_idx < n; i_idx++) {
                    int node_i = selected.get(i_idx);
                    int delta = calculateDeltaInterNode(selected, node_i, node_c_old, inst);
                    if (delta < 0) {
                        int prev = selected.get((i_idx - 1 + n) % n);
                        int next = selected.get((i_idx + 1) % n);
                        LM.add(new Classes.Move("inter", prev, next, node_i, node_c_old, delta));
                    }
                }
            }
        }


        // --- 2. Intra-route Re-evaluation (Wider Check) ---
        // Nodes A, B, C, D were involved in the previous move.
        // Check all pairs of edges (A', B') and (C', D') where A' or C' is one of the affected nodes (A, B, C, D)

        Set<Integer> nodes_to_check = new HashSet<>(changedNodes);
        
        for (int a_node : nodes_to_check) { // Node A'
            int i_idx = selected.indexOf(a_node);
            if (i_idx == -1) continue;
            
            // Iterate through all other possible start nodes C' for the second edge
            for (int j_idx = 0; j_idx < n; j_idx++) { // Index of C'
                if (i_idx == j_idx) continue;
                
                // Ensure i_idx < j_idx for consistent indexing to match delta logic
                int current_i = i_idx;
                int current_j = j_idx;
                if (current_i > current_j) {
                    int tmp = current_i; current_i = current_j; current_j = tmp;
                }
                
                // Skip adjacent edges
                if (current_j == current_i + 1 || (current_i == 0 && current_j == n - 1)) continue;

                int a = selected.get(current_i);
                int c = selected.get(current_j);

                // Recompute delta
                int delta = calculateDeltaIntraEdge(selected, a, c, inst);
                
                if (delta < 0) {
                    int next_a = selected.get((current_i + 1) % n); 
                    int next_c = selected.get((current_j + 1) % n); 
                    LM.add(new Classes.Move("intra_edge", a, next_a, c, next_c, delta));
                }
            }
        }

        // ---- Sort LM by delta (best first) ----
        LM.sort(Comparator.comparingInt(m -> m.delta));
    }

    // ---------- Main ----------
    public static void main(String[] args) {
        try {
            String instance = "TSPB";
            String filepath = "./data/" + instance + ".csv";
            Classes.Instance inst = Helpers.readCSV(filepath);
            System.out.println("Loaded instance with " + inst.n + " nodes");

            final int NUM_RUNS = inst.n;
            Random rand = new Random(42);

            Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
            Map<String, List<Long>> timings = new LinkedHashMap<>();

            String methodName = "LS_steepest_with_LM";
            results.put(methodName, new ArrayList<>());
            timings.put(methodName, new ArrayList<>());

            // --- Search loop ---
            for (int i = 0; i < NUM_RUNS; i++) {
                
                if (i % 25 == 0) {
                    System.out.println("Processing instance " + i + " / " + NUM_RUNS);
                }

                long t_start_ls = System.nanoTime();

                List<Integer> start_sol = getRandomSolution(inst, rand);

                Classes.Pair<List<Integer>, Integer> final_sol = localSearch(
                        start_sol, inst, "steepest", "edges", rand
                );

                long t_end_ls = System.nanoTime();

                results.get(methodName).add(final_sol);
                timings.get(methodName).add(t_end_ls - t_start_ls);
            }

            
            // --- Report solutions
            System.out.println("\n---------- Results (Objective Function) ----------");
            System.out.printf("%-30s | %10s | %10s | %10s%n", "Method", "Avg", "Min", "Max");
            System.out.println(new String(new char[66]).replace("\0", "-"));

            Map<String, Classes.Pair<List<Integer>, Integer>> bestSolutions = new LinkedHashMap<>();

            for (var entry : results.entrySet()) {
                String name = entry.getKey();
                List<Classes.Pair<List<Integer>, Integer>> sols = entry.getValue();

                int min = sols.stream().mapToInt(p -> p.second).min().orElse(0);
                int max = sols.stream().mapToInt(p -> p.second).max().orElse(0);
                double avg = sols.stream().mapToInt(p -> p.second).average().orElse(0.0);

                System.out.printf("%-30s | %10.0f | %10d | %10d%n", name, avg, min, max);

                bestSolutions.put(name, sols.stream().min(Comparator.comparingInt(p -> p.second)).orElse(null));
            }
            
            // --- Report running times ---
            System.out.println("\n---------- Results (Running Time [ms]) ----------");
            System.out.printf("%-30s | %10s | %10s | %10s%n", "Method", "Avg (ms)", "Min (ms)", "Max (ms)");
            System.out.println(new String(new char[66]).replace("\0", "-"));

             for (var entry : timings.entrySet()) {
                String name = entry.getKey();
                List<Long> times_ns = entry.getValue();

                // Convert nanoseconds to milliseconds
                double min = times_ns.stream().mapToLong(l -> l).min().orElse(0) / 1_000_000.0;
                double max = times_ns.stream().mapToLong(l -> l).max().orElse(0) / 1_000_000.0;
                double avg = times_ns.stream().mapToLong(l -> l).average().orElse(0.0) / 1_000_000.0;

                System.out.printf("%-30s | %10.2f | %10.2f | %10.2f%n", name, avg, min, max);
            }


            // --- Save visualizations ---
            System.out.println("\nSaving best solution plots...");
            for (var entry : bestSolutions.entrySet()) {
                if (entry.getValue() == null) continue;
                String method = entry.getKey();
                Classes.Pair<List<Integer>, Integer> best = entry.getValue();

                // Create subdirectory if it doesn't exist
                String dirPath = "./lab5/visualizations/";
                new java.io.File(dirPath).mkdirs();
                String plotPath = dirPath + "/" + instance + "_" + method;

                System.out.println("Plotting best for " + method);
                System.out.println("Best value: " + best.second);
                System.out.println("Best path: " + best.first);
                Plotter.showPlot(inst.nodes, best.first, plotPath);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
