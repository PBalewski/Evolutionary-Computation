import java.util.*;
import java.util.stream.*;

public class Solution {
    
    // ---- Initializations ----
    //Random 
    public static List<Integer> getRandomSolution(Classes.Instance inst, Random rand) {
        List<Integer> allNodes = IntStream.range(0, inst.n).boxed().collect(Collectors.toList());
        Collections.shuffle(allNodes, rand);
        int size = (inst.n + 1) /2;
        return new ArrayList<>(allNodes.subList(0, size));
    }

    //regretNearestAny(0.5)
    public static Classes.Pair<List<Integer>, Integer> regretNearestAny(Classes.Instance inst, int start, int k, double regretWeight) {

        List<Integer> selected = new ArrayList<>();
        selected.add(start);
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.remove(start);

        while (selected.size() < (inst.n + 1) / 2) {
            int bestNode = -1;
            int bestPos = -1;
            double bestScore = Double.NEGATIVE_INFINITY;

            // Iterate over each remaining node
            for (int j : remaining) {
                List<Double> deltas = new ArrayList<>();

                // For each possible insertion position in the current path
                for (int pos = 0; pos <= selected.size(); pos++) {
                    double delta;
                    if (pos == 0) {
                        // insert before the first node
                        int kNode = selected.get(0);
                        delta = inst.distMat[j][kNode] + inst.nodes.get(j).cost;
                    } else if (pos == selected.size()) {
                        // insert after the last node
                        int i = selected.get(selected.size() - 1);
                        delta = inst.distMat[i][j] + inst.nodes.get(j).cost;
                    } else {
                        // insert between i and k
                        int i = selected.get(pos - 1);
                        int kNode = selected.get(pos);
                        delta = inst.distMat[i][j] + inst.distMat[j][kNode] - inst.distMat[i][kNode] + inst.nodes.get(j).cost;
                    }
                    deltas.add(delta);
                }

                // Sort deltas ascending to compute regret
                List<Double> sorted = new ArrayList<>(deltas);
                Collections.sort(sorted);

                double bestDelta = sorted.get(0);
                double regret = 0.0;
                for (int m = 1; m < k; m++) {
                    regret += (sorted.get(m) - bestDelta);
                }

                double score = regretWeight*regret - (1-regretWeight)*bestDelta;

                if (score > bestScore) {
                    bestScore = score;
                    bestNode = j;
                    bestPos = deltas.indexOf(bestDelta);
                }
            }

            // Insert node at best position
            selected.add(bestPos, bestNode);
            remaining.remove(bestNode);
        }

        int val = Helpers.evaluate(selected, inst, true);
        return new Classes.Pair<>(selected, val);
    }


    // ---- Calculating deltas ----
    // Inter-route Node: switching node from the cycle with the one from outside
    private static int calculateDeltaInterNode(List<Integer> selected, int i_idx, int k_node, Classes.Instance inst) {
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
    
    // Intra-route Node: switching places of two nodes from the cycle
    private static int calculateDeltaIntraNode(List<Integer> selected, int i_idx, int j_idx, Classes.Instance inst) {
        int n = selected.size();
        // swapping in cycles smaller than 3 has no effect
        if (n < 3) return 0;

        // Ensure i_idx is smaller than j_idx
        if (i_idx > j_idx) {
            int temp = i_idx;
            i_idx = j_idx;
            j_idx = temp;
        }

        int i_node = selected.get(i_idx);
        int j_node = selected.get(j_idx);

        // Check if nodes are adjacent (including cycle wrap-around)
        if (j_idx == i_idx + 1 || (i_idx == 0 && j_idx == n - 1)) {
            // For adjacent nodes, we need special handling
            int before, after;
            if (i_idx == 0 && j_idx == n - 1) {
                // Special case: positions 0 and n-1 are adjacent in cycle
                // Cycle: ... [n-2] -> [n-1] -> [0] -> [1] ...
                // We need the actual neighbors, not the positions we're swapping
                int prev_j = selected.get(n - 2);
                int next_i = selected.get(1);

                before = inst.distMat[prev_j][j_node] + inst.distMat[j_node][i_node] + inst.distMat[i_node][next_i];
                after = inst.distMat[prev_j][i_node] + inst.distMat[i_node][j_node] + inst.distMat[j_node][next_i];
            } else {
                // Regular adjacent case: positions i and i+1
                // Cycle: ... [i-1] -> [i] -> [i+1] -> [i+2] ...
                int prev_i = selected.get((i_idx - 1 + n) % n);
                int next_j = selected.get((j_idx + 1) % n);
                
                before = inst.distMat[prev_i][i_node] + inst.distMat[i_node][j_node] + inst.distMat[j_node][next_j];
                after = inst.distMat[prev_i][j_node] + inst.distMat[j_node][i_node] + inst.distMat[i_node][next_j];
            }

            // For adjacent nodes we consider 3 edges that change:
            // BEFORE: prev_i -> i_node -> j_node -> next_j
            // AFTER:  prev_i -> j_node -> i_node -> next_j
            return after - before;
        } else {
            // For non-adjacent nodes we need to consider 4 edges that change
            int prev_i = selected.get((i_idx - 1 + n) % n);
            int next_i = selected.get((i_idx + 1) % n);
            int prev_j = selected.get((j_idx - 1 + n) % n);
            int next_j = selected.get((j_idx + 1) % n);

            // BEFORE: prev_i -> i_node -> next_i  and  prev_j -> j_node -> next_j
            // AFTER:  prev_i -> j_node -> next_i  and  prev_j -> i_node -> next_j
            int before = inst.distMat[prev_i][i_node] + inst.distMat[i_node][next_i]
                       + inst.distMat[prev_j][j_node] + inst.distMat[j_node][next_j];
            int after = inst.distMat[prev_i][j_node] + inst.distMat[j_node][next_i]
                      + inst.distMat[prev_j][i_node] + inst.distMat[i_node][next_j];
            return after - before;
        }
    }

    //Intra-route Edge: switch 2 edges
    private static int calculateDeltaIntraEdge(List<Integer> selected, int i_idx, int j_idx, Classes.Instance inst) {
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

        // For segment reversal [i..j], we:
        // - break edges (i-1)->i and j->(j+1)
        // - add edges (i-1)->j and i->(j+1)
        // This effectively reverses the segment [i..j]
        int prev_i = selected.get((i_idx - 1 + n) % n);
        int node_i = selected.get(i_idx);
        int node_j = selected.get(j_idx);
        int next_j = selected.get((j_idx + 1) % n);

        // Calculate total change in cycle length
        int remove_dist = inst.distMat[prev_i][node_i] + inst.distMat[node_j][next_j];
        int add_dist = inst.distMat[prev_i][node_j] + inst.distMat[node_i][next_j];

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
        int iteration = 0;

        // while(improved && iteration < max_iterations) {

        while (improved) {
            iteration++;

            if (ls_type.equals("steepest")) {
                improved = steepestSearch(current_solution, remaining, inst, move_type);
            } else {
                improved = greedySearch(current_solution, remaining, inst, move_type, rand);
            }
        }

    int final_cost = Helpers.evaluate(current_solution, inst, true);
    return new Classes.Pair<>(current_solution, final_cost);
}

    //steepest search
    private static boolean steepestSearch(List<Integer> selected, Set<Integer> remaining, Classes.Instance inst, String move_type) {
        
        int best_delta = 0;
        Runnable best_move = null;
        int n = selected.size();

        // Inter-route neighbourhood
        for (int i_idx = 0; i_idx < n; i_idx++) {
            for (int k_node: remaining) {
                int delta = calculateDeltaInterNode(selected, i_idx, k_node, inst);
                if (delta < best_delta) {
                    best_delta = delta;
                    final int final_i_idx = i_idx;
                    final int final_k_node = k_node;
                    best_move = () -> {
                        int i_node_to_remove = selected.set(final_i_idx, final_k_node);
                        remaining.remove(final_k_node);
                        remaining.add(i_node_to_remove);
                    };
                }
            }
        }

        // Intra-route neighbourhood: nodes && edges
        for (int i_idx = 0; i_idx < n; i_idx++) {
            for (int j_idx = i_idx + 1; j_idx < n; j_idx++) {
                int delta;
                Runnable current_action;

                if (move_type.equals("nodes")) {
                    delta = calculateDeltaIntraNode(selected, i_idx, j_idx, inst);
                    final int final_i_idx = i_idx;
                    final int final_j_idx = j_idx;
                    current_action = () -> {
                        Collections.swap(selected, final_i_idx, final_j_idx);
                    };
                } else {
                    delta = calculateDeltaIntraEdge(selected, i_idx, j_idx, inst);
                    final int final_i_idx = i_idx;
                    final int final_j_idx = j_idx;
                    current_action = () -> {
                        // Reversal inclusive i..j matches the delta calculation
                        Helpers.reverseSublist(selected, final_i_idx, final_j_idx);
                    };
                }

                if (delta < best_delta) {
                    best_delta = delta;
                    best_move = current_action;
                }
            }
        }

        // Perform the best move
        if (best_move != null) {
            best_move.run();
            return true;
        }
        return false;
    }

    //greedy search
    private static boolean greedySearch(List<Classes.Pair<Integer, Integer>> inter_moves, List<Integer> selected, Set<Integer> remaining, Classes.Instance inst, String move_type, Random rand) {
        int n = selected.size();

        // Create the whole neighbourhood
        List<Classes.Pair<Integer, Integer>> inter_moves = new ArrayList<>(); //(i_idx, k_node)
        for (int i_idx = 0; i_idx < n; i_idx++) {
            for (int k_node: remaining) {
                inter_moves.add(new Classes.Pair<>(i_idx, k_node));
            }
        }

        List<Classes.Pair<Integer, Integer>> intra_moves = new ArrayList<>(); // (i_idx, j_idx)
         for (int i_idx = 0; i_idx < n; i_idx++) {
            for (int j_idx = i_idx + 1; j_idx < n; j_idx++) {
                intra_moves.add(new Classes.Pair<>(i_idx, j_idx));
            }
        }

        // Shuffle lists
        Collections.shuffle(inter_moves, rand);
        Collections.shuffle(intra_moves, rand);

        // Random pick which list to look through first
        boolean checkInterFirst = rand.nextBoolean();
        if (checkInterFirst) {
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, move_type, inter_moves, true)) return true;
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, move_type, intra_moves, false)) return true;
        } else {
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, move_type, intra_moves, false)) return true;
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, move_type, inter_moves, true)) return true;
        }
        return false;
    }

    //helper for greedy
    private static boolean findAndApplyFirstImprovingMove(
            List<Integer> selected,
            Set<Integer> remaining,
            Classes.Instance inst,
            String intraMoveType,
            List<Classes.Pair<Integer, Integer>> moves,
            boolean isInterMoveList)
    {
        for (Classes.Pair<Integer, Integer> move : moves) {
            int delta;
            int i = move.first;
            int j_or_k = move.second;

            if (isInterMoveList) {
                // Inter (i = i_idx, j_or_k = k_node)
                delta = calculateDeltaInterNode(selected, i, j_or_k, inst);
                if (delta < 0) {
                    int i_node_to_remove = selected.set(i, j_or_k);
                    remaining.remove(j_or_k);
                    remaining.add(i_node_to_remove);
                    return true;
                }
            } else {
                // Intra move (i = i_idx, j_or_k = j_idx)
                if (intraMoveType.equals("nodes")) {
                    delta = calculateDeltaIntraNode(selected, i, j_or_k, inst);
                    if (delta < 0) {
                        Collections.swap(selected, i, j_or_k);
                        return true;
                    }
                } else { // "edges" (2-opt move)
                    delta = calculateDeltaIntraEdge(selected, i, j_or_k, inst);
                     if (delta < 0) {
                        // Reversal inclusive i..j matches the delta calculation
                        Helpers.reverseSublist(selected, i, j_or_k);
                        return true;
                    }
                }
            }
        }
        return false;
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

            // experiment variants
            String[] ls_types = {"steepest", "greedy"};
            String[] intra_types = {"nodes", "edges"};
            String[] start_types = {"NNregret", "random"};

            timings.put("NNregret_sol_creation", new ArrayList<>());
            timings.put("random_sol_creation", new ArrayList<>());

            // create result/timing buckets
            for (String ls_type : ls_types) {
                for (String intra_type : intra_types) {
                    for (String start_type : start_types) {
                        String methodName = String.format("%s_%s_%s", ls_type, intra_type, start_type);
                        results.put(methodName, new ArrayList<>());
                        timings.put(methodName, new ArrayList<>());
                    }
                }
            }

            // --- Search loop ---
            for (int i = 0; i < NUM_RUNS; i++) {
                
                if (i % 25 == 0) {
                    System.out.println("Processing instance " + i + " / " + NUM_RUNS);
                }

                // iterate over experiment variants
                for (String start_type : start_types) {
                    List<Integer> start_sol;

                    long t_start = System.nanoTime();
                    if (start_type.equals("NNregret")) {
                        Classes.Pair<List<Integer>, Integer> start_pair = regretNearestAny(inst, i, 2, 0.5);
                        start_sol = start_pair.first;
                        timings.get("NNregret_sol_creation").add(System.nanoTime() - t_start);
                    } else {
                        start_sol = getRandomSolution(inst, rand);
                        timings.get("random_sol_creation").add(System.nanoTime() - t_start);
                    }

                    // run all LS variants on the chosen start solution
                    for (String ls_type : ls_types) {
                        for (String intra_type : intra_types) {
                            String method_name = String.format("%s_%s_%s", ls_type, intra_type, start_type);

                            // int max_ls_iterations = 5000;

                            long t_start_ls = System.nanoTime();

                            Classes.Pair<List<Integer>, Integer> final_sol = localSearch(
                                    start_sol, inst, ls_type, intra_type, rand
                            );

                            long t_end_ls = System.nanoTime();

                            results.get(method_name).add(final_sol);
                            timings.get(method_name).add(t_end_ls - t_start_ls);
                        }
                    }
                }
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
                String dirPath = "./lab3/visualizations/";
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