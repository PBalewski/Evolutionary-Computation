//lab6
import java.util.*;

public class Solution {
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



    // Record to hold ILS specific results
    record ILSResult(List<Integer> solution, int cost, int lsRuns) {}

    // --- LOCAL SEARCH ---
    public static Classes.Pair<List<Integer>, Integer> localSearch(
        List<Integer> initial_solution,
        Classes.Instance inst)
    {
        List<Integer> current_solution = new ArrayList<>(initial_solution);
        Set<Integer> remaining = Helpers.getRemaining(inst.n, current_solution);
        

        boolean improved = true;
        while(improved) {
            improved = steepestSearch(current_solution, remaining, inst);
        }

        int final_cost = Helpers.evaluate(current_solution, inst, true);
        return new Classes.Pair<>(current_solution, final_cost);
    }

    //steepest search
    private static boolean steepestSearch(
        List<Integer> selected, 
        Set<Integer> remaining, 
        Classes.Instance inst) {
        
        int best_delta = 0;
        Runnable best_move = null;
        int n = selected.size();

        // 1. Inter-route: node exchange
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

        // 2. Intra-route: edges exchange
        for (int i_idx = 0; i_idx < n; i_idx++) {
            for (int j_idx = i_idx + 1; j_idx < n; j_idx++) {
                int delta;
                Runnable current_action;

                delta = calculateDeltaIntraEdge(selected, i_idx, j_idx, inst);
                final int final_i_idx = i_idx;
                final int final_j_idx = j_idx;
                current_action = () -> {
                    // Reversal inclusive i..j matches the delta calculation
                    Helpers.reverseSublist(selected, final_i_idx, final_j_idx);
                };

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

    /**
     * Segment Exchange:
     * Divide the cycle into 4 parts (A, B, C, D) wap the middle segments -> (A, C, B, D)
     * It's an intra-route move, in the worst case its a node swap
     */
    public static void perturb(
            List<Integer> solution,
            Classes.Instance inst,
            Random rand) {
        
        int n = solution.size();

        // 1. Get random 3 separating points (p1 < p2 < p3) and that each segment has length of at least 1
        int p1 = 1 + rand.nextInt(n - 3);
        int p2 = p1 + 1 + rand.nextInt(n - 2 - p1);
        int p3 = p2 + 1 + rand.nextInt(n - 1 - p2);

        // 2. Get the segments
        List<Integer> A = new ArrayList<>(solution.subList(0, p1));
        List<Integer> B = new ArrayList<>(solution.subList(p1, p2));
        List<Integer> C = new ArrayList<>(solution.subList(p2, p3));
        List<Integer> D = new ArrayList<>(solution.subList(p3, n));

        // 3. Swap B and C
        solution.clear();
        solution.addAll(A);
        solution.addAll(C);
        solution.addAll(B);
        solution.addAll(D);
    }

    // 1. Multiple Start Local Search (MSLS)
    public static Classes.Pair<List<Integer>, Integer> runMSLS(
        Classes.Instance inst,
        int iterations,
        Random rand) {
        
        Classes.Pair<List<Integer>, Integer> best_global = null;

        for(int i=0; i<iterations; i++){
            List<Integer> initial_solution = Helpers.getRandomSolution(inst, rand);
            Classes.Pair<List<Integer>, Integer> result = localSearch(initial_solution, inst);

            if( best_global == null || result.second < best_global.second) {
                best_global = result;
            }
            
        }
        return best_global;
    };

    // 2. Iterated Local Search (ILS)
    public static ILSResult runILS(
        Classes.Instance inst,
        long time_limit,
        Random rand) {

        long start_time = System.nanoTime();
        int runs = 0;
        List<Integer> initial_solution = Helpers.getRandomSolution(inst, rand);

        //initial search
        Classes.Pair<List<Integer>, Integer> best_result = localSearch(initial_solution, inst);
        List<Integer> best_solution = best_result.first;
        int best_cost = best_result.second;
        runs++;

        //main loop
        while(System.nanoTime() - start_time < time_limit){

            List<Integer> perturbed_solution = new ArrayList<>(best_solution);

            //peturb
            perturb(perturbed_solution, inst, rand);

            //local search
            Classes.Pair<List<Integer>, Integer> result = localSearch(perturbed_solution, inst);
            runs++;

            if (result.second < best_cost) {
            best_cost = result.second;
            best_solution = result.first;
            // We accepted a new best, next perturbation will start from here; otherwise peturb the best_solution again
            }
        }
        
        return new ILSResult(best_solution, best_cost, runs);
    };

    
    // ---------- Main ----------   
    public static void main(String[] args) {
        try{
            String instance = "TSPB"; // Make sure to test all instances
            String filepath = "./data/" + instance + ".csv";
            Classes.Instance inst = Helpers.readCSV(filepath);
            System.out.println("Loaded instance with " + inst.n + " nodes");

            final int EXPERIMENT_RUNS = 20;
            final int MSLS_ITERATIONS = 200;
            Random rand = new Random(42);
            
            Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
            Map<String, List<Long>> timings = new LinkedHashMap<>();
            Map<String, List<Integer>> counts = new LinkedHashMap<>();

            results.put("MSLS", new ArrayList<>());
            timings.put("MSLS", new ArrayList<>());
            
            timings.put("ILS", new ArrayList<>());
            results.put("ILS", new ArrayList<>());
            counts.put("ILS", new ArrayList<>());


            System.out.println("\n--- Starting MSLS Experiment (" + EXPERIMENT_RUNS + " runs) ---");

            // Run MSLS
            long total_msls_time = 0;
            for (int i = 0; i < EXPERIMENT_RUNS; i++) {
                long start = System.nanoTime();
                Classes.Pair<List<Integer>, Integer> result = runMSLS(inst, MSLS_ITERATIONS, rand);
                long duration = System.nanoTime() - start;
                total_msls_time += duration;

                results.get("MSLS").add(result);
                timings.get("MSLS").add(duration);

                if ((i+1) % 1 == 0) System.out.println("MSLS Run " + (i+1) + " done.");
            }
            long avg_msls_time = total_msls_time / EXPERIMENT_RUNS;
            System.out.printf("\nAverage MSLS Time: %.0f ms. Starting ILS.\n", avg_msls_time / 1e6);

            // Run ILS
            for (int i = 0; i < EXPERIMENT_RUNS; i++) {
                long start = System.nanoTime();
                ILSResult result = runILS(inst, avg_msls_time, rand);
                long duration = System.nanoTime() - start; // Should be close to avgMslsTimeNs

                results.get("ILS").add(new Classes.Pair<>(result.solution, result.cost));
                timings.get("ILS").add(duration);
                counts.get("ILS").add(result.lsRuns);

                if ((i+1) % 1 == 0) System.out.println("ILS Run " + (i+1) + " done.");
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

            // --- Report number of LS runs ---
            System.out.println("\n================ RESULTS (ILS Runs) ================");
            System.out.printf("%-30s | %10s | %10s | %10s%n", "Method", "Avg", "Min", "Max");
            System.out.println(new String(new char[66]).replace("\0", "-"));

            for (var entry : counts.entrySet()) {
                String name = entry.getKey();
                List<Integer> runs = entry.getValue();

                int min = runs.stream().mapToInt(v -> v).min().orElse(0);
                int max = runs.stream().mapToInt(v -> v).max().orElse(0);
                double avg = runs.stream().mapToInt(v -> v).average().orElse(0.0);

                System.out.printf("%-30s | %10.2f | %10d | %10d%n", name, avg, min, max);
            }

            // --- Save visualizations ---
            System.out.println("\nSaving best solution plots...");
            // Create subdirectory if it doesn't exist
            String dirPath = "./lab6/visualizations/";
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
        
    };
};