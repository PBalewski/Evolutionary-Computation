import java.util.*;
import java.util.stream.*;

public class Solution {

    // Candidate list
    public static int[][] getCandidateList(Classes.Instance inst, int top) {
        int n = inst.n;
        int[][] candidate_list = new int[n][top];

        for (int i=0; i<n; i++) {
            List<Classes.Pair<Integer, Double>> node_scores = new ArrayList<>();
            for (int j=0; j<n; j++) {
                if (i == j) continue;

                double score = inst.distMat[i][j] + inst.nodes.get(j).cost;
                node_scores.add(new Classes.Pair<>(j, score));
            }

            //sort & get top
            node_scores.sort(Comparator.comparingDouble(p -> p.second));

            for(int m=0; m<top; m++) {
                candidate_list[i][m] = node_scores.get(m).first;
            }
        }
        return candidate_list;
    }

    public static Classes.Pair<List<Integer>, Integer> localSearch(
        List<Integer> initial_solution,
        Classes.Instance inst,
        int[][] candidate_list)
    {
        List<Integer> current_solution = new ArrayList<>(initial_solution);
        Set<Integer> remaining = Helpers.getRemaining(inst.n, current_solution);
        

        boolean improved = true;
        while(improved) {
            improved = steepestSearchCandidate(current_solution, remaining, inst, candidate_list);
        }

        int final_cost = Helpers.evaluate(current_solution, inst, true);
        return new Classes.Pair<>(current_solution, final_cost);
    }

    private static boolean steepestSearchCandidate(
        List<Integer> selected,
        Set<Integer> remaining,
        Classes.Instance inst,
        int[][] candidate_list)
    {
        int best_delta = 0;
        Runnable best_move = null;
        int n = selected.size();
        Set<Long> checkedMoves = new HashSet<>(); //avoid re-evaluation moves


        // Inter-route neighbourhood
        // swap(i, k) is a candidate if (p_i, k) or (k, n_i) is candidate edge
        for (int i_idx=0; i_idx<n; i_idx++) {
            int prev_node = selected.get((i_idx-1+n) %n);
            int next_node = selected.get((i_idx+1) %n);

            // Check if candidate edge: prev_node -> k_node
            for(int k_node : candidate_list[prev_node]) {
                if (remaining.contains(k_node)) {
                    // long moveKey = ((long)i_node << 32) | k_node;
                    // if (checkedMoves.add(moveKey)) { // Check if move already evaluated
                    
                    int delta = Helpers.calculateDeltaInterNode(selected, i_idx, k_node, inst);
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
            // Check if candidate edge: k_node -> next_node
            for(int k_node : candidate_list[next_node]) {
                if (remaining.contains(k_node)) {
                    // long moveKey = ((long)i_node << 32) | k_node;
                    // if (checkedMoves.add(moveKey)) { // Check if move already evaluated

                    int delta = Helpers.calculateDeltaInterNode(selected, i_idx, k_node, inst);
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

       
        }
        
        checkedMoves.clear(); // Clear for intra-route moves
        
        // Intra-route neighbourhood (edges)
        // reverse(i+1, j) is a candidate if (i,j) or (i+1, j+1) is candidate edge
        for (int i_idx=0; i_idx<n; i_idx++) {
            int i_node = selected.get(i_idx);
            int i_next = selected.get((i_idx+1) %n);

            // Check if candidate edge: i -> j
            for (int j_node : candidate_list[i_node]) {
                int j_idx = selected.indexOf(j_node);
                if (j_idx == -1) continue; //j not in the cycle
                // long moveKey = ((long)Math.min(i_idx, j_idx) << 32) | Math.max(i_idx, j_idx);
                // if (checkedMoves.add(moveKey)) {

                int delta = Helpers.calculateDeltaIntraEdge(selected, i_idx, j_idx, inst);
                if (delta < best_delta) {
                        best_delta = delta;
                        final int final_i_idx = (i_idx < j_idx) ? i_idx : j_idx;
                        final int final_j_idx = (i_idx < j_idx) ? j_idx : i_idx;
                        best_move = () -> Helpers.reverseSublist(selected, final_i_idx + 1, final_j_idx);
                }
            }

            // Check if candidate edge: i+1 -> j+1
            for (int j_next : candidate_list[i_next]) {
                int j_next_idx = selected.indexOf(j_next);
                if (j_next_idx == -1) continue;

                // Find original j_idx from j_next_idx
                int j_idx = (j_next_idx - 1 + n) % n;

                // long moveKey = ((long)Math.min(i_idx, j_idx) << 32) | Math.max(i_idx, j_idx);
                // if (checkedMoves.add(moveKey)) {

                int delta = Helpers.calculateDeltaIntraEdge(selected, i_idx, j_idx, inst);
                if (delta < best_delta) {
                        best_delta = delta;
                        final int final_i_idx = (i_idx < j_idx) ? i_idx : j_idx;
                        final int final_j_idx = (i_idx < j_idx) ? j_idx : i_idx;
                        best_move = () -> Helpers.reverseSublist(selected, final_i_idx + 1, final_j_idx);
                }
            }
        }

        if (best_move != null) {
            best_move.run();
            return true;
        }

        return false;
    }

    // ---------- Main ----------
    public static void main(String[] args) {
        try{
            String instance = "TSPA"; // Make sure to test all instances
            String filepath = "./data/" + instance + ".csv";
            Classes.Instance inst = Helpers.readCSV(filepath);
            System.out.println("Loaded instance with " + inst.n + " nodes");

            final int NUM_RUNS = inst.n;
            final int CANDIDATE_K = 10;
            Random rand = new Random(42);
            
            Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
            Map<String, List<Long>> timings = new LinkedHashMap<>();

            results.put("steepestCandidate", new ArrayList<>());
            timings.put("steepestCandidate", new ArrayList<>());
            timings.put("random_sol_creation", new ArrayList<>());
            timings.put("candidate_list_creation", new ArrayList<>());

            //candidate list
            System.out.println("Creating candidate list (k=" + CANDIDATE_K + ")");
            long t_start_candidates = System.nanoTime();
            int[][] candidateList = getCandidateList(inst, CANDIDATE_K);
            timings.get("candidate_list_creation").add(System.nanoTime() - t_start_candidates);


            // --- Main Experiment Loop ---
            for (int i = 0; i < NUM_RUNS; i++) {
                
                if ((i + 1) % 25 == 0) {
                    System.out.println("Processing run: " + (i + 1) + " / " + NUM_RUNS);
                }

                long t_start_sol = System.nanoTime();
                List<Integer> start_sol = Helpers.getRandomSolution(inst, rand);
                timings.get("random_sol_creation").add(System.nanoTime() - t_start_sol);

                // local search
                long t_start_cand = System.nanoTime();
                Classes.Pair<List<Integer>, Integer> final_sol_cand = localSearch(
                        start_sol, inst, candidateList
                );
                long t_end_cand = System.nanoTime();
                results.get("steepestCandidate").add(final_sol_cand);
                timings.get("steepestCandidate").add(t_end_cand - t_start_cand);
                
            }

            
            // --- Report solutions ---
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
            // Create subdirectory if it doesn't exist
            String dirPath = "./lab4/visualizations/";
            new java.io.File(dirPath).mkdirs();

            for (var entry : bestSolutions.entrySet()) {
                if (entry.getValue() == null) continue;
                String method = entry.getKey();
                Classes.Pair<List<Integer>, Integer> best = entry.getValue();
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


