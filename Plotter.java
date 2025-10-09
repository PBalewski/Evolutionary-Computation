import javax.swing.*;
import java.awt.*;
import java.util.List;

public class Plotter extends JPanel {

    private final List<Classes.Node> nodes;
    private final List<Integer> path;

    public Plotter(List<Classes.Node> nodes, List<Integer> path) {
        this.nodes = nodes;
        this.path = path;
        setPreferredSize(new Dimension(800, 800));
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        Graphics2D g2 = (Graphics2D) g;
        g2.setStroke(new BasicStroke(2));

        // Find max/min coordinates to scale to panel
        int margin = 200;
        int minX = nodes.stream().mapToInt(n -> n.x).min().orElse(0) - margin;
        int maxX = nodes.stream().mapToInt(n -> n.x).max().orElse(1) + margin;
        int minY = nodes.stream().mapToInt(n -> n.y).min().orElse(0) - margin;
        int maxY = nodes.stream().mapToInt(n -> n.y).max().orElse(1) + margin;

        double scaleX = getWidth() / (double) (maxX - minX + 1);
        double scaleY = getHeight() / (double) (maxY - minY + 1);

        // Draw edges (cycle)
        g2.setColor(Color.BLACK);
        for (int i = 0; i < path.size(); i++) {
            Classes.Node a = nodes.get(path.get(i));
            Classes.Node b = nodes.get(path.get((i + 1) % path.size()));
            int x1 = (int) ((a.x - minX) * scaleX);
            int y1 = (int) ((a.y - minY) * scaleY);
            int x2 = (int) ((b.x - minX) * scaleX);
            int y2 = (int) ((b.y - minY) * scaleY);
            g2.drawLine(x1, y1, x2, y2);
        }

        // Find min/max cost to map to color or size
        int minCost = nodes.stream().mapToInt(n -> n.cost).min().orElse(0);
        int maxCost = nodes.stream().mapToInt(n -> n.cost).max().orElse(1);

        // Draw nodes
        for (int idx = 0; idx < nodes.size(); idx++) {
            Classes.Node n = nodes.get(idx);

            // Map cost to color (blue = low, red = high)
            float ratio = (n.cost - minCost) / (float) (maxCost - minCost);
            Color color = new Color(ratio, 0f, 1f - ratio); // RGB gradient
            g2.setColor(color);

            int x = (int) ((n.x - minX) * scaleX);
            int y = (int) ((n.y - minY) * scaleY);

            int radius = 10 + (int) (10 * ratio); // larger radius for higher cost
            g2.fillOval(x - radius / 2, y - radius / 2, radius, radius);
        }
    }

    // ---------- Utility to display ----------
    public static void showPlot(List<Classes.Node> nodes, List<Integer> path, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.add(new Plotter(nodes, path));
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}
