import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.io.File;

public class Solution {
    
    // --- LOCAL SEARCH ---
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
    public static Classes.Pair<List<Integer>, Integer> localSearch(
            List<Integer> initial_solution,
            Classes.Instance inst,
            Random rand)
    {
        List<Integer> current_solution = new ArrayList<>(initial_solution);
        Set<Integer> remaining = Helpers.getRemaining(inst.n, current_solution);

        boolean improved = true;
        // int iteration = 0;

        // while(improved && iteration < max_iterations) {

        while (improved) {
            // iteration++;
            improved = steepestSearch(current_solution, remaining, inst);
        }
        int final_cost = Helpers.evaluate(current_solution, inst, true);
        return new Classes.Pair<>(current_solution, final_cost);
    }
    private static boolean steepestSearch(List<Integer> selected, Set<Integer> remaining, Classes.Instance inst) {
        
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
    };

    private static boolean greedySearch( 
        List<Integer> selected, 
        Set<Integer> remaining, 
        Classes.Instance inst,
        Random rand) 
    {
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
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, inter_moves, true)) return true;
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, intra_moves, false)) return true;
        } else {
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, intra_moves, false)) return true;
            if (findAndApplyFirstImprovingMove(selected, remaining, inst, inter_moves, true)) return true;
        }
        return false;
    }
    private static boolean findAndApplyFirstImprovingMove(
            List<Integer> selected,
            Set<Integer> remaining,
            Classes.Instance inst,
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
                    // System.out.println(delta);
                    return true;
                }
            } else {
                delta = calculateDeltaIntraEdge(selected, i, j_or_k, inst);
                    if (delta < 0) {
                    // Reversal inclusive i..j matches the delta calculation
                    Helpers.reverseSublist(selected, i, j_or_k);
                    // System.out.println(delta);
                    return true;
                    }
                }
            }
        return false;
    }

    // ---- repair function based on regretNearestAny (regretWeight=0.5) greedy heuristic  ----
    public static Classes.Pair<List<Integer>, Integer> repair(Classes.Instance inst, List<Integer> selected, int k, double regretWeight) {
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.removeAll(selected);
        
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


    // --- OPERATORS ---
    private static List<Integer> operator1(
        List<Integer> p1,
        List<Integer> p2,
        int length,
        Random rand) {
        //find common edges
        Set<String> edges_p1 = Helpers.getEdgesSet(p1);
        List<List<Integer>> segments = Helpers.getCommonEdges(p2, edges_p1);

        Set<Integer> used_nodes = new HashSet<>();
        for (List<Integer> seg : segments) used_nodes.addAll(seg);

        //find common nodes
        Set<Integer> common_nodes = Helpers.getCommonNodes(p1, p2);

        //add them
        for(int node : common_nodes) {
            if (!used_nodes.contains(node)){
                List<Integer> single_node = new ArrayList<>();
                single_node.add(node);
                segments.add(single_node);
                used_nodes.add(node);
            }
        }

        //fill with nodes if needed
        List<Integer> available = IntStream.range(0, length)
                .filter(i -> !used_nodes.contains(i))
                .boxed().collect(Collectors.toList());
        Collections.shuffle(available, rand);

        int needed = length - used_nodes.size();
        for (int i = 0; i < needed && i < available.size(); i++) {
            List<Integer> singleNode = new ArrayList<>();
            singleNode.add(available.get(i));
            segments.add(singleNode);
        }

        //combine randomly
        Collections.shuffle(segments, rand);
        List<Integer> child = new ArrayList<>();

        for (List<Integer> child_seg : segments) {
            if (child_seg.size() > 1 && rand.nextBoolean()) {
                Collections.reverse(child_seg);
            }
            child.addAll(child_seg);
        }
        return child;
    };

    private static List<Integer> operator2(
            List<Integer> p1, 
            List<Integer> p2, 
            Classes.Instance inst) {
        
        boolean[] otherPresence = Helpers.getNodePresenceArray(p2, inst.n);
        List<Integer> child = new ArrayList<>(p1.size());
        for (int node : p1) {
            if (otherPresence[node]) {
                child.add(node);
            }
        }

        return repair(inst, child, 2, 0.5).first;
    }
    

    // --- HYBRID EVOLUTIONARY --

    //intialization
    private static List<Classes.Pair<List<Integer>, Integer>> getInitialPopulation(Classes.Instance inst, int pop_size, Random rand) {
        List<Classes.Pair<List<Integer>, Integer>> population = new ArrayList<>();
        while (population.size() < pop_size) {
            List<Integer> solution = Helpers.getRandomSolution(inst, rand);
            Classes.Pair<List<Integer>, Integer> result = localSearch(solution, inst, rand);
            if (!isDuplicate(population, result.second)) {
                population.add(result);
            }
        }
        return population;
    }

    public static Classes.Pair<List<Integer>, Integer> runHEA(
            Classes.Instance inst, 
            int pop_size,
            long timeLimit, 
            int operatorType, // 1 or 2
            boolean useLS,    // true/false
            Random rand) {

        long startTime = System.nanoTime();
        
        // 1. Population init
        List<Classes.Pair<List<Integer>, Integer>> population = getInitialPopulation(inst, pop_size, rand);

        // 2. Sort ascending
        population.sort(Comparator.comparingInt(p -> p.second));
        // System.out.println("Population costs" + population.stream()
        //         .map(p -> p.second.toString())
        //         .collect(Collectors.joining(", ", "[", "]")));

        // 3. Main part
        while (System.nanoTime() - startTime < timeLimit) {
            
            // A. Select (uniform)
            int idx1 = rand.nextInt(pop_size);
            int idx2 = rand.nextInt(pop_size);
            while (idx1 == idx2) idx2 = rand.nextInt(pop_size);
            
            List<Integer> p1 = population.get(idx1).first;
            List<Integer> p2 = population.get(idx2).first;

            // B. Recombination
            List<Integer> offspring;
            if (operatorType == 1) {
                offspring = operator1(p1, p2, inst.n/2, rand);
            } else {
                // random base for operator
                if (rand.nextBoolean()) offspring = operator2(p1, p2, inst);
                else offspring = operator2(p2, p1, inst);
            }

            // C. Local Search
            int offspringCost;
            offspringCost = Helpers.evaluate(offspring, inst, true);
            // System.out.println(offspringCost + ", Local search: " + useLS);
            if (useLS) {
                Classes.Pair<List<Integer>, Integer> res = localSearch(offspring, inst, rand);
                offspring = res.first;
                offspringCost = res.second;
            } else {
                offspringCost = Helpers.evaluate(offspring, inst, true);
            }
            // System.out.println(offspringCost);

            // D. Update
            Classes.Pair<List<Integer>, Integer> worst = population.get(pop_size - 1);
            
            if (offspringCost < worst.second && !isDuplicate(population, offspringCost)) {
                population.remove(pop_size - 1);
                population.add(new Classes.Pair<>(offspring, offspringCost));
                population.sort(Comparator.comparingInt(p -> p.second));
            }
        }
        return population.get(0);
    }

    // check duplicates
    private static boolean isDuplicate(List<Classes.Pair<List<Integer>, Integer>> pop, int cost) {
        for (Classes.Pair<List<Integer>, Integer> ind : pop) {
            if (ind.second == cost) return true;
        }
        return false;
    }

    public static void main(String[] args) {
        try {
            // Configuration
            String[] instances = {"TSPA", "TSPB"}; 
            Random rand = new Random(15);
            int POPULATION_SIZE = 20;
            int EXPERIMENT_RUNS = 20;
            
            // Time limit per run (based on MSLS average)
            long STOP_TIME_MS = 8903; //8903 
            long STOP_TIME_NS = STOP_TIME_MS * 1_000_000L;

            for (String instance : instances) {
                System.out.println("\n==========================================");
                System.out.println("PROCESSING INSTANCE: " + instance);
                System.out.println("==========================================");
                
                String filepath = "./data/" + instance + ".csv";
                Classes.Instance inst = Helpers.readCSV(filepath);

                // Store results and times
                Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
                Map<String, List<Long>> times = new LinkedHashMap<>();
                
                String[] methods = {"HEA_Op1_LS", "HEA_Op2_LS", "HEA_Op2_NoLS"};
                for (String m : methods) {
                    results.put(m, new ArrayList<>());
                    times.put(m, new ArrayList<>());
                }

                // --- 1. Run Experiments ---
                
                // HEA Op1 LS
                System.out.println("Running HEA (Op1 + LS)...");
                for (int i = 0; i < EXPERIMENT_RUNS; i++) {
                    long start = System.nanoTime();
                    results.get("HEA_Op1_LS").add(runHEA(inst, POPULATION_SIZE, STOP_TIME_NS, 1, true, rand));
                    long end = System.nanoTime();
                    times.get("HEA_Op1_LS").add((end - start) / 1_000_000); // Convert ns to ms
                }

                // HEA Op2 LS
                System.out.println("Running HEA (Op2 + LS)...");
                for (int i = 0; i < EXPERIMENT_RUNS; i++) {
                    long start = System.nanoTime();
                    results.get("HEA_Op2_LS").add(runHEA(inst, POPULATION_SIZE, STOP_TIME_NS, 2, true, rand));
                    long end = System.nanoTime();
                    times.get("HEA_Op2_LS").add((end - start) / 1_000_000);
                }

                // HEA Op2 NoLS
                System.out.println("Running HEA (Op2 NO LS)...");
                for (int i = 0; i < EXPERIMENT_RUNS; i++) {
                    long start = System.nanoTime();
                    results.get("HEA_Op2_NoLS").add(runHEA(inst, POPULATION_SIZE, STOP_TIME_NS, 2, false, rand));
                    long end = System.nanoTime();
                    times.get("HEA_Op2_NoLS").add((end - start) / 1_000_000);
                }

                // --- 2. Print Results Summary (Objective Function) ---
                System.out.println("\nRESULTS SUMMARY (" + instance + ") - COST");
                System.out.println("-----------------------------------------------------");
                System.out.printf("%-15s | %10s | %10s | %10s%n", "Method", "Avg", "Min", "Max");
                System.out.println("-----------------------------------------------------");

                for (String m : methods) {
                    List<Classes.Pair<List<Integer>, Integer>> all = results.get(m);
                    int min = all.stream().mapToInt(p -> p.second).min().orElse(0);
                    int max = all.stream().mapToInt(p -> p.second).max().orElse(0);
                    double avg = all.stream().mapToInt(p -> p.second).average().orElse(0.0);
                    System.out.printf("%-15s | %10.0f | %10d | %10d%n", m, avg, min, max);
                }
                System.out.println("-----------------------------------------------------");

                // --- 3. Print Time Summary (Execution Time) ---
                System.out.println("\nTIME SUMMARY (" + instance + ") - MILLISECONDS");
                System.out.println("-----------------------------------------------------");
                System.out.printf("%-15s | %10s | %10s | %10s%n", "Method", "Avg", "Min", "Max");
                System.out.println("-----------------------------------------------------");

                for (String m : methods) {
                    List<Long> tList = times.get(m);
                    long minT = tList.stream().mapToLong(l -> l).min().orElse(0);
                    long maxT = tList.stream().mapToLong(l -> l).max().orElse(0);
                    double avgT = tList.stream().mapToLong(l -> l).average().orElse(0.0);
                    System.out.printf("%-15s | %10.0f | %10d | %10d%n", m, avgT, minT, maxT);
                }
                System.out.println("-----------------------------------------------------");

                // --- 4. Print Best Paths and Save Plots ---
                System.out.println("\nBEST SOLUTIONS DETAILS (" + instance + ")");
                String plotDir = "./lab9/visualizations/";
                new File(plotDir).mkdirs();

                for (String m : methods) {
                    List<Classes.Pair<List<Integer>, Integer>> all = results.get(m);
                    Classes.Pair<List<Integer>, Integer> best = all.stream()
                        .min(Comparator.comparingInt(p -> p.second))
                        .orElseThrow();

                    System.out.println("\nMethod: " + m);
                    System.out.println("Best Cost: " + best.second);
                    System.out.println("Best Path: " + best.first);

                    // Save Plot
                    String plotPath = plotDir + instance + "_" + m + ".png";
                    Plotter plotter = new Plotter(inst.nodes, best.first);
                    Plotter.savePlotAsImage(plotter, plotPath);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
