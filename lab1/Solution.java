import java.util.*;
import java.util.stream.*;

public class Solution {

    public static Classes.Pair<List<Integer>, Integer> randomSolution(Classes.Instance inst) {
        int count = (inst.n + 1) / 2;
        List<Integer> nodes = IntStream.range(0, inst.n).boxed().collect(Collectors.toList());
        Collections.shuffle(nodes, new Random());
        List<Integer> selected = nodes.subList(0, count);
        int val = Helpers.evaluate(selected, inst, true);
        return new Classes.Pair<>(new ArrayList<>(selected), val);
    }

    public static Classes.Pair<List<Integer>, Integer> nearestEnd(Classes.Instance inst, int start) {
        List<Integer> selected = new ArrayList<>();
        selected.add(start);
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.remove(start);

        while (selected.size() < (inst.n + 1) / 2) {
            int last = selected.get(selected.size() - 1);
            int bestNode = -1, bestVal = Integer.MAX_VALUE;

            for (int j : remaining) {
                int val = inst.distMat[last][j] + inst.nodes.get(j).cost;
                if (val < bestVal) {
                    bestVal = val;
                    bestNode = j;
                }
            }
            selected.add(bestNode);
            remaining.remove(bestNode);
        }
        int val = Helpers.evaluate(selected, inst, true);
        return new Classes.Pair<>(selected, val);
    }

    public static Classes.Pair<List<Integer>, Integer> nearestAny(Classes.Instance inst, int start) {
        List<Integer> selected = new ArrayList<>();
        selected.add(start);
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.remove(start);

        while (selected.size() < (inst.n + 1) / 2) {
            int bestNode = -1, bestPos = -1, bestVal = Integer.MAX_VALUE;

            for (int j : remaining) {
                for (int pos = 0; pos <= selected.size(); pos++) {
                    List<Integer> trial = new ArrayList<>(selected);
                    trial.add(pos, j);
                    int val = Helpers.evaluate(trial, inst, false);
                    if (val < bestVal) {
                        bestVal = val;
                        bestNode = j;
                        bestPos = pos;
                    }
                }
            }
            selected.add(bestPos, bestNode);
            remaining.remove(bestNode);
        }
        int val = Helpers.evaluate(selected, inst, true);
        return new Classes.Pair<>(selected, val);
    }

    public static Classes.Pair<List<Integer>, Integer> greedyCycle(Classes.Instance inst, int start) {
        Set<Integer> remaining = IntStream.range(0, inst.n).boxed().collect(Collectors.toSet());
        remaining.remove(start);

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
            int bestNode = -1, bestPos = -1, bestEval = Integer.MAX_VALUE;

            for (int j : remaining) {
                for (int pos = 0; pos < selected.size(); pos++) {
                    List<Integer> trial = new ArrayList<>(selected);
                    trial.add(pos + 1, j);
                    int val = Helpers.evaluate(trial, inst, true);
                    if (val < bestEval) {
                        bestEval = val;
                        bestNode = j;
                        bestPos = pos + 1;
                    }
                }
            }
            selected.add(bestPos, bestNode);
            remaining.remove(bestNode);
        }
        int val = Helpers.evaluate(selected, inst, true);
        return new Classes.Pair<>(selected, val);
    }

    // ---------- Main ----------

    public static void main(String[] args) {
        try {
            String filepath = "TSPA.csv";
            Classes.Instance inst = Helpers.readCSV(filepath);
            System.out.println("Loaded instance with " + inst.n + " nodes");

            Map<String, List<Classes.Pair<List<Integer>, Integer>>> results = new LinkedHashMap<>();
            results.put("randomSolution", new ArrayList<>());
            results.put("nearestEnd", new ArrayList<>());
            results.put("nearestAny", new ArrayList<>());
            results.put("greedyCycle", new ArrayList<>());

            Random rnd = new Random(42);
            int runs = 200;

            for (int i = 0; i < runs; i++) {
                results.get("randomSolution").add(randomSolution(inst));
            }

            for (int start = 0; start < 5; start++) {
                if (start % 20 == 0) System.out.println("Processing start node: " + start);

                results.get("nearestEnd").add(nearestEnd(inst, start));
                results.get("nearestAny").add(nearestAny(inst, start));
                results.get("greedyCycle").add(greedyCycle(inst, start));
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

                Plotter.showPlot(inst.nodes, best.first, method + " Best Solution");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
