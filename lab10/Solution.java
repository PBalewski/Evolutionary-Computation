import java.util.*;
import java.util.stream.*;


public class Solution {
    public static List<Integer> regretNearestAny(Classes.Instance inst, int start, int k, double regretWeight) {

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
        return selected;
    }



    private static List<Classes.Move> steepestSearch_LM = new ArrayList<>();

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

    public static Classes.Pair<List<Integer>, Integer> localSearch(
                List<Integer> initial_solution,
                Classes.Instance inst)
            {
        List<Integer> current_solution = new ArrayList<>(initial_solution);
        Set<Integer> remaining = Helpers.getRemaining(inst.n, current_solution);

        boolean improved = true;
        // int iteration = 0; // Unused variable
        
        steepestSearch_LM.clear();
        while (improved) {
            // iteration++; // Unused increment
            improved = steepestSearch_LM(current_solution, remaining, inst);
        }
        int final_cost = Helpers.evaluate(current_solution, inst, true);
        return new Classes.Pair<>(current_solution, final_cost);
    }

    // steepest local search with reuse of previous deltas (LM)
    private static boolean steepestSearch_LM(List<Integer> selected, Set<Integer> remaining, Classes.Instance inst) {
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


    public static List<Integer> destroy(
            List<Integer> sol,
            Classes.Instance inst,
            double removeFraction)
    {
        Random rnd = new Random();
        int n = sol.size();
        int targetRemove = (int)Math.round(n * removeFraction);

        List<Integer> remaining = new ArrayList<>(sol);
        Set<Integer> removed = new HashSet<>();

        // ----- Compute edge costs -----
        double[] edgeCost = new double[n];
        double maxCost = 0;

        for (int i = 0; i < n; i++) {
            int A = sol.get(i);
            int B = sol.get((i + 1) % n);
            double dist = inst.distMat[A][B];
            double nodesCost = inst.nodes.get(A).cost + inst.nodes.get(B).cost;
            double c = dist + nodesCost;
            edgeCost[i] = c;
            if (c > maxCost) maxCost = c;
        }

        // ----- Probabilities proportional to cost -----
        double[] prob = new double[n];
        for (int i = 0; i < n; i++) {
            // base probability 0.2 + biased term
            prob[i] = 0.2 + 0.8 * (edgeCost[i] / maxCost);
        }

        // ====== Optionally remove a whole SUBPATH (30% chance) ======
        if (rnd.nextDouble() < 0.30) {
            int start = rnd.nextInt(n);
            int min = targetRemove / 3;
            int max = targetRemove;
            int len = min + rnd.nextInt(max - min + 1);

            for (int k = 0; k < len && removed.size() < targetRemove; k++) {
                int idx = (start + k) % n;
                int node = sol.get(idx);
                if (removed.add(node)) {
                    remaining.remove(Integer.valueOf(node));
                }
            }
        }

        // ====== Probabilistic SCATTERED removal ======
        while (removed.size() < targetRemove) {
            int i = rnd.nextInt(n);   // pick a random edge index

            // biased by edge cost
            if (rnd.nextDouble() < prob[i]) {
                int node = sol.get((i + 1) % n);
                if (removed.add(node)) {
                    remaining.remove(Integer.valueOf(node));
                }
            }
        }

        return remaining;
    }


    public static List<Integer> destroyRelated(List<Integer> sol, Classes.Instance inst, double removeFraction, Random rnd) {
        int n = sol.size();
        int targetRemove = (int) Math.round(n * removeFraction);
        List<Integer> remaining = new ArrayList<>(sol);
        Set<Integer> removedSet = new HashSet<>();

        // 1. Pick a random seed node from the current solution
        int seedNode = sol.get(rnd.nextInt(n));
        removedSet.add(seedNode);
        remaining.remove(Integer.valueOf(seedNode));

        // 2. Remove nodes closest to the seed node
        List<Integer> candidates = new ArrayList<>(remaining);
        candidates.sort(Comparator.comparingDouble(node -> inst.distMat[seedNode][node]));

        for (int i = 0; i < targetRemove - 1 && i < candidates.size(); i++) {
            int toRemove = candidates.get(i);
            removedSet.add(toRemove);
            remaining.remove(Integer.valueOf(toRemove));
        }

        return remaining;
    }


    record LNS_Result(List<Integer> solution, int cost, int runs) {}

    public static LNS_Result LNS_ultimate(
            Classes.Instance inst,
            double removeFraction,
            Boolean bonus_LS_iter,
            long timeLimit,
            Random rand,
            int startNode)
    {
        long startTime = System.nanoTime();

        int runs = 0;
        List<Integer> initial_solution = regretNearestAny(inst, startNode, 2, 0.5);

        Classes.Pair<List<Integer>, Integer> xPair = localSearch(initial_solution, inst);
        List<Integer> x = xPair.first;
        int fx = xPair.second;
        steepestSearch_LM = new ArrayList<>();

        List<Integer> xBest = new ArrayList<>(x);
        int fBest = fx;

        double T = 100.0; // Initial temperature
        double coolingRate = 0.9995;
        double[] weights = {1.0, 1.0}; 
        // MAIN LNS LOOP
        while(System.nanoTime() - startTime < timeLimit){
            List<Integer> y;
            int destroyType = (rand.nextDouble() < (weights[0] / (weights[0] + weights[1]))) ? 0 : 1;
            if (destroyType == 1) {
                y = destroyRelated(x, inst, 0.3, rand); // NEW: Geographic logic
            } else {
                y = destroy(x, inst, 0.3); // Existing cost-based logic
            }


            Classes.Pair<List<Integer>, Integer> yPair = repair(inst, y, 2, 0.5);
            y = yPair.first;
            int fy = yPair.second;

            if (bonus_LS_iter) {
                yPair = localSearch(y, inst);
                y = yPair.first;
                fy = yPair.second;
                steepestSearch_LM = new ArrayList<>();
            }

            double delta = fy - fx;
            if (delta < 0 || rand.nextDouble() < Math.exp(-delta / T)) {
                x = new ArrayList<>(y);
                fx = fy;

                if (delta < 0) {
                    weights[destroyType] += 1.0;  // Small reward for improvement
                }
                
                if (fx < fBest) {
                    weights[destroyType] += 3.0; // Biger reward for new global best
                    fBest = fx;
                    xBest = new ArrayList<>(x);
                }
            }

            T *= coolingRate;
            runs++;
        }
        return new LNS_Result(xBest, fBest, runs);
    }


    // ---------- Main ----------   
    public static void main(String[] args) {
        try{
            String instance = "TSPB";
            String filepath = "./data/" + instance + ".csv";
            long avg_msls_time = 8992L * 1_000_000L; // CHANGE ACCORDING TO INSTANCE
            Classes.Instance inst = Helpers.readCSV(filepath);
            System.out.println("Loaded instance" + instance + "with " + inst.n + " nodes");

            final int EXPERIMENT_RUNS = 20;
            Random rand = new Random(42);
            double removeFraction = 0.3;
            
            Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
            Map<String, List<Long>> timings = new LinkedHashMap<>();
            Map<String, List<Integer>> counts = new LinkedHashMap<>();

            results.put("LNS_ultimate", new ArrayList<>());
            timings.put("LNS_ultimate", new ArrayList<>());
            counts.put("LNS_ultimate", new ArrayList<>());


            System.out.println("\n--- Starting LNS_ultimate Experiment (" + EXPERIMENT_RUNS + " runs) ---");
            
            // Run largeNS_LS
            for (int i = 0; i < EXPERIMENT_RUNS; i++) {
                long start = System.nanoTime();
                LNS_Result result = LNS_ultimate(inst, removeFraction, true, avg_msls_time, rand, rnd.nextInt(inst.n));
                long duration = System.nanoTime() - start;

                results.get("LNS_ultimate").add(new Classes.Pair<>(result.solution, result.cost));
                timings.get("LNS_ultimate").add(duration);
                counts.get("LNS_ultimate").add(result.runs);

                System.out.println("LNS_ultimate Run " + (i+1) + " done.");
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

            // --- Report number of main loop runs ---
            System.out.println("\n================ RESULTS (main loop runs) ================");
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
            String dirPath = "./lab10/visualizations/";
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