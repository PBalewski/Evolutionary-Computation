import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

public class Helpers {

    public static Classes.Instance readCSV(String filepath) throws IOException {
        List<Classes.Node> nodes = new ArrayList<>();
        List<String> lines = Files.readAllLines(Paths.get(filepath));

        for (int idx = 0; idx < lines.size(); idx++) {
            String[] parts = lines.get(idx).trim().split(";");
            if (parts.length != 3)
                throw new IllegalArgumentException("Invalid data format in line " + (idx + 1));

            int x = Integer.parseInt(parts[0]);
            int y = Integer.parseInt(parts[1]);
            int cost = Integer.parseInt(parts[2]);
            nodes.add(new Classes.Node(idx, x, y, cost));
        }
        return new Classes.Instance(nodes);
    }

    public static int evaluate(List<Integer> path, Classes.Instance inst, boolean formCycle) {
        if (path.size() < 2)
            return path.stream().mapToInt(i -> inst.nodes.get(i).cost).sum();

        int dist = 0;
        if (formCycle) {
            for (int i = 0; i < path.size(); i++) {
                int a = path.get(i);
                int b = path.get((i + 1) % path.size()); //includes edge between last and first node
                dist += inst.distMat[a][b];
            }
        } else {
            for (int i = 0; i < path.size() - 1; i++) {
                int a = path.get(i);
                int b = path.get(i + 1);
                dist += inst.distMat[a][b];
            }
        }

        int cost = path.stream().mapToInt(i -> inst.nodes.get(i).cost).sum();
        return dist + cost;
    }
}
