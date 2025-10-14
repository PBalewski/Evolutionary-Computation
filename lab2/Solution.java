import java.util.*;
import java.util.stream.*;

public class Solution {

    public static Classes.Pair<List<Integer>, Integer> regretNearestAny(Classes.Instance inst, int start, int k) {

        List<Integer> selected = new ArrayList<>();
        selected.add(start);
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.remove(start);

        while (selected.size() < (inst.n + 1) / 2) {
            int bestNode = -1;
            int bestPos = -1;
            double bestRegret = Double.NEGATIVE_INFINITY;

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

                if (regret > bestRegret) {
                    bestRegret = regret;
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

    public static Classes.Pair<List<Integer>, Integer> regretGreedyCycle(Classes.Instance inst, int start, int k) {
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.remove(start);

        // ---- initialize with best pair ----
        int bestSecond = -1, bestVal = Integer.MAX_VALUE;
        for (int j : remaining) {
            int val = inst.distMat[start][j] + inst.nodes.get(j).cost;
            if (val < bestVal) {
                bestVal = val;
                bestSecond = j;
            }
        }

        List<Integer> selected = new ArrayList<>(Arrays.asList(start, bestSecond));
        remaining.remove(bestSecond);

        while (selected.size() < (inst.n + 1) / 2) {
            int bestNode = -1, bestPos = -1;
            double bestRegret = Double.NEGATIVE_INFINITY;

            // evaluate regret for each remaining node
            for (int j : remaining) {
                List<Double> deltas = new ArrayList<>();

                // compute insertion delta for all edges (cycle!)
                for (int pos = 0; pos < selected.size(); pos++) {
                    int i = selected.get(pos);
                    int kNode = selected.get((pos + 1) % selected.size()); // pairs last with first node
                    double delta = inst.distMat[i][j] + inst.distMat[j][kNode] - inst.distMat[i][kNode] + inst.nodes.get(j).cost;
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

                if (regret > bestRegret) {
                    bestRegret = regret;
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


    // ---------- Main ----------

    public static void main(String[] args) {
        try {
            String instance = "TSPA";
            String filepath = ".\\data\\" + instance + ".csv";
            Classes.Instance inst = Helpers.readCSV(filepath);
            System.out.println("Loaded instance with " + inst.n + " nodes");

            Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
            results.put("regretNearestAny", new ArrayList<>());
            results.put("regretGreedyCycle", new ArrayList<>());

            for (int start = 0; start < inst.n; start++) {
                if (start % 20 == 0) System.out.println("Processing start node: " + start);

                results.get("regretNearestAny").add(regretNearestAny(inst, start, 2));
                results.get("regretGreedyCycle").add(regretGreedyCycle(inst, start, 2));
            }


            Map<String, Classes.Pair<List<Integer>, Integer>> bestSolutions = new LinkedHashMap<>();

            for (var entry : results.entrySet()) {
                String name = entry.getKey();
                List<Classes.Pair<List<Integer>, Integer>> sols = entry.getValue();

                int min = sols.stream().mapToInt(p -> p.second).min().orElse(0);
                int max = sols.stream().mapToInt(p -> p.second).max().orElse(0);
                double avg = sols.stream().mapToInt(p -> p.second).average().orElse(0.0);

                System.out.printf("%-15s -> min: %d  max: %d  avg: %.2f%n", name, min, max, avg);

                Classes.Pair<List<Integer>, Integer> best = sols.stream()
                        .min(Comparator.comparingInt(p -> p.second))
                        .orElse(null);

                bestSolutions.put(name, best);
            }

            for (var entry : bestSolutions.entrySet()) {
                String method = entry.getKey();
                Classes.Pair<List<Integer>, Integer> best = entry.getValue();

                System.out.println("\nMethod: " + method);
                System.out.println("Best value: " + best.second);
                System.out.println("Best path: " + best.first);

                Plotter.showPlot(inst.nodes, best.first, "lab2/visualizations/" + instance + "_" + method);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

