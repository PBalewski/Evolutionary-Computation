import java.util.*;
import java.io.File;

public class Solution {

    // --- GREEDY LOCAL SEARCH (from lab3, intra: edges) ---
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
            improved = greedySearch(current_solution, remaining, inst, rand);
        }
        int final_cost = Helpers.evaluate(current_solution, inst, true);
        return new Classes.Pair<>(current_solution, final_cost);
    }

    //greedy search
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

    //helper for greedy
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


    // --- MAIN EXPERIMENT ---

    public static void main(String[] args) {
        try {
            String[] instances = {"TSPA", "TSPB"};
            Random rand = new Random(42);
            int NUM_OPTIMA = 1000;

            for (String instanceName : instances) {
                System.out.println("\nProcessing instance: " + instanceName);
                String filepath = "./data/" + instanceName + ".csv";
                Classes.Instance inst = Helpers.readCSV(filepath);

                // 1. Generate 1000 Random Local Optima
                System.out.println("Generating " + NUM_OPTIMA + " random local optima (Greedy LS)...");
                List<List<Integer>> solutions = new ArrayList<>();
                List<Integer> costs = new ArrayList<>();
                
                for (int i = 0; i < NUM_OPTIMA; i++) {
                    if (i % 25 == 0){
                    System.out.println("Processing instance " + i + " / " + NUM_OPTIMA);
                }
                    List<Integer> start_solution = Helpers.getRandomSolution(inst, rand);
                    Classes.Pair<List<Integer>, Integer> result = localSearch(start_solution, inst, rand);
                    solutions.add(result.first);
                    costs.add(result.second);
                }

                // 2. Identify Reference Solutions
                
                // Ref A: Best of the 1000
                int best_idx = 0;
                for (int i = 1; i < costs.size(); i++) {
                    if (costs.get(i) < costs.get(best_idx)) {
                        best_idx = i;
                    }
                }
                List<Integer> best_1000 = solutions.get(best_idx);
                int best_1000_cost = costs.get(best_idx);

                // Ref B: Best Known (from ILS)
                List<Integer> best_ils;
                int best_ils_cost;

                if (instanceName.equals("TSPA")) {
                    best_ils = Arrays.asList(9, 148, 124, 94, 63, 79, 80, 176, 137, 23, 186, 89, 183, 143, 0, 117, 93, 140, 108, 18, 22, 159, 193, 41, 139, 68, 46, 115, 59, 118, 51, 151, 133, 162, 149, 131, 65, 116, 43, 42, 181, 34, 160, 48, 54, 177, 10, 190, 184, 35, 84, 4, 112, 123, 127, 70, 135, 154, 180, 53, 100, 26, 86, 75, 101, 1, 97, 152, 2, 120, 44, 25, 16, 171, 175, 113, 56, 31, 78, 145, 196, 81, 90, 165, 119, 40, 185, 179, 92, 129, 57, 55, 52, 106, 178, 49, 14, 144, 102, 62); 
                    best_ils_cost = 69176;
                } else {
                    best_ils = Arrays.asList(21, 82, 111, 144, 33, 160, 29, 0, 109, 35, 143, 106, 124, 62, 18, 55, 34, 170, 152, 183, 140, 4, 149, 28, 20, 60, 148, 47, 94, 66, 179, 22, 99, 130, 95, 185, 86, 166, 194, 176, 113, 114, 137, 127, 89, 103, 163, 187, 153, 81, 77, 141, 91, 61, 36, 177, 5, 78, 175, 142, 45, 80, 190, 136, 73, 54, 31, 193, 117, 198, 156, 1, 38, 63, 40, 107, 133, 122, 135, 131, 121, 51, 90, 147, 6, 188, 169, 132, 70, 3, 15, 145, 13, 195, 168, 139, 11, 138, 104, 8);
                    best_ils_cost = 43535;
                }

                System.out.println("Best of 1000 Cost: " + best_1000_cost);
                System.out.println("Best Known (ILS) Cost: " + best_ils_cost);


                // 3. Calculate Correlations and Prepare Plots (12 Charts Logic)
                // Measures: [Nodes, Edges] x Methods: [Avg, BestOf1000, BestKnown] x Instance [TSPA, TSPB]

                String plotDir = "./lab8/visualizations/";
                new File(plotDir).mkdirs();

                // --- DATA COLLECTION ---
                
                // Structure: Data for [Method][Measure] -> List of {Value, Similarity}
                // Methods: 0=Avg, 1=BestOf1000, 2=BestKnown
                // Measures: 0=Nodes, 1=Edges
                List<double[]>[][] chartData = new ArrayList[3][2];
                for(int i=0; i<3; i++) for(int j=0; j<2; j++) chartData[i][j] = new ArrayList<>();

                System.out.println("Calculating similarities...");
                
                for (int i = 0; i < NUM_OPTIMA; i++) {
                    List<Integer> current = solutions.get(i);
                    double objective = costs.get(i); 

                    // --- Measure 1: Common NODES ---
                    
                    // 1.1 Similarity to best_1000
                    // Exclude the best solution itself to avoid outlier (100% similarity)
                    if (i != best_idx) { 
                        double nodes_sim_1000 = Helpers.similarityNodes(current, best_1000);
                        chartData[1][0].add(new double[]{objective, nodes_sim_1000});
                    }
                    
                    // 1.2 Similarity to best_ils
                    double nodes_sim_ils = Helpers.similarityNodes(current, best_ils);
                    chartData[2][0].add(new double[]{objective, nodes_sim_ils});

                    // 1.3 Average Similarity (Nodes)
                    double nodes_sim_sum = 0;
                    for (int j = 0; j < NUM_OPTIMA; j++) {
                        if (i == j) continue;
                        nodes_sim_sum += Helpers.similarityNodes(current, solutions.get(j));
                    }
                    double nodes_sim_avg = nodes_sim_sum / (NUM_OPTIMA - 1);
                    chartData[0][0].add(new double[]{objective, nodes_sim_avg});


                    // --- Measure 2: Common EDGES ---

                    // 2.1 Similarity to best_1000
                    if (i != best_idx) {
                        double edges_sim_1000 = Helpers.similarityEdges(current, best_1000);
                        chartData[1][1].add(new double[]{objective, edges_sim_1000});
                    }

                    // 2.2 Similarity to BestKnown
                    double edges_sim_ils = Helpers.similarityEdges(current, best_ils);
                    chartData[2][1].add(new double[]{objective, edges_sim_ils});

                    // 2.3 Average Similarity (Edges)
                    double edges_sim_sum = 0;
                    for (int j = 0; j < NUM_OPTIMA; j++) {
                        if (i == j) continue;
                        edges_sim_sum += Helpers.similarityEdges(current, solutions.get(j));
                    }
                    double edges_sim_avg = edges_sim_sum / (NUM_OPTIMA - 1);
                    chartData[0][1].add(new double[]{objective, edges_sim_avg});
                }
                
                // --- GLOBAL SCALE FOR PLOTS ---
                double minCost = costs.stream().min(Integer::compare).get();
                double maxCost = costs.stream().max(Integer::compare).get();


                // --- PLOTTING ---
                String[] methods = {"Average", "Best_1000", "Best_ILS"};
                String[] measures = {"Nodes", "Edges"};

                for (int m = 0; m < 3; m++) {     // Method
                    for (int t = 0; t < 2; t++) { // Measure
                        String title = String.format("%s - %s Sim (%s)", instanceName, measures[t], methods[m]);
                        String filename = String.format("%s%s_%s_%s.png", plotDir, instanceName, measures[t], methods[m]);
                        
                        Plotter.saveScatterPlot(
                            chartData[m][t],
                            title,
                            "Objective Function (Cost)",
                            "Similarity (" + measures[t] + ")",
                            filename,
                            minCost, maxCost, // Fixed X-axis scale
                            0.2, 1              // Fixed Y-axis scale (normalized)
                        );
                    }
                }
            }
            System.out.println("\nAll charts generated successfully.");

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}