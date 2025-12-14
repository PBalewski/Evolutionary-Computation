import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;
import javax.imageio.ImageIO;
import javax.swing.*;

public class Plotter extends JPanel {
    private final List<Classes.Node> nodes;
    private final List<Integer> path;

    public Plotter(List<Classes.Node> nodes, List<Integer> path) {
        this.nodes = nodes;
        this.path = path;
        setPreferredSize(new Dimension(800, 800));
    }

    // ====================================
    // ===       TSP DRAWING            ===
    // ====================================

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setStroke(new BasicStroke(2));

        if (nodes == null || nodes.isEmpty()) return;

        // --- 1. Calculate Scaling (Old Style with Margin 200) ---
        int margin = 200;
        int minX = nodes.stream().mapToInt(n -> n.x).min().orElse(0) - margin;
        int maxX = nodes.stream().mapToInt(n -> n.x).max().orElse(1) + margin;
        int minY = nodes.stream().mapToInt(n -> n.y).min().orElse(0) - margin;
        int maxY = nodes.stream().mapToInt(n -> n.y).max().orElse(1) + margin;

        double scaleX = getWidth() / (double) (maxX - minX + 1);
        double scaleY = getHeight() / (double) (maxY - minY + 1);

        // --- 2. Draw Edges (Black) ---
        g2.setColor(Color.BLACK);
        if (path != null && path.size() > 1) {
            for (int i = 0; i < path.size(); i++) {
                // Handle cycle connection safely
                if (i == path.size() - 1 && path.size() < 2) continue;
                
                Classes.Node a = nodes.get(path.get(i));
                Classes.Node b = nodes.get(path.get((i + 1) % path.size()));
                
                int x1 = (int) ((a.x - minX) * scaleX);
                int y1 = (int) ((a.y - minY) * scaleY);
                int x2 = (int) ((b.x - minX) * scaleX);
                int y2 = (int) ((b.y - minY) * scaleY);
                g2.drawLine(x1, y1, x2, y2);
            }
        }

        // --- 3. Draw Nodes (Gradient Blue->Red, Variable Size) ---
        // Find min/max cost to map to color/size
        int minCost = nodes.stream().mapToInt(n -> n.cost).min().orElse(0);
        int maxCost = nodes.stream().mapToInt(n -> n.cost).max().orElse(1);

        for (Classes.Node n : nodes) {
            // Map cost to color (blue = low, red = high)
            float ratio = (float) (n.cost - minCost) / (maxCost - minCost);
            // Clamp ratio to 0-1 just in case
            ratio = Math.max(0, Math.min(1, ratio));
            
            Color color = new Color(ratio, 0f, 1f - ratio); // RGB gradient
            g2.setColor(color);

            int x = (int) ((n.x - minX) * scaleX);
            int y = (int) ((n.y - minY) * scaleY);

            int radius = 10 + (int) (10 * ratio); // larger radius for higher cost
            g2.fillOval(x - radius / 2, y - radius / 2, radius, radius);
        }
    }

    // =========================================================
    // ===           UTILITY METHODS (Display/Save)          ===
    // =========================================================

    public static void showPlot(List<Classes.Node> nodes, List<Integer> path, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        Plotter panel = new Plotter(nodes, path);
        frame.add(panel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
        
        // Auto-save when shown (optional, consistent with old code)
        savePlotAsImage(panel, title + ".png"); 
    }

    public static void savePlotAsImage(Plotter plotter, String filepath) {
        int width = 800;
        int height = 800;
        
        // Match Old Style: ARGB for transparency support
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2 = image.createGraphics();
        
        // Ensure the plotter has a size for layout calculation
        plotter.setSize(width, height);
        
        // Optional: Fill white background if you don't want transparency
        // g2.setColor(Color.WHITE);
        // g2.fillRect(0, 0, width, height);
        
        plotter.paintComponent(g2);
        g2.dispose();

        try {
            File file = new File(filepath);
            file.getParentFile().mkdirs(); 
            ImageIO.write(image, "png", file);
            System.out.println("Saved TSP plot: " + filepath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static final int SCATTER_WIDTH = 800;
    private static final int SCATTER_HEIGHT = 600;
    private static final int SCATTER_PADDING = 60;
    private static final int LABEL_PADDING = 30;
    private static final int POINT_SIZE = 6;

    /**
     * Scatter Plot and save the image
     */
    public static void saveScatterPlot(List<double[]> data, String title, String xAxis, String yAxis, String filepath) {
        double minX = Double.MAX_VALUE, maxX = Double.MIN_VALUE;
        double minY = Double.MAX_VALUE, maxY = Double.MIN_VALUE;

        if (data.isEmpty()) {
            minX = 0; maxX = 10; minY = 0; maxY = 10;
        } else {
            for (double[] point : data) {
                if (point[0] < minX) minX = point[0];
                if (point[0] > maxX) maxX = point[0];
                if (point[1] < minY) minY = point[1];
                if (point[1] > maxY) maxY = point[1];
            }
        }
        saveScatterPlot(data, title, xAxis, yAxis, filepath, minX, maxX, minY, maxY);
    }

    public static void saveScatterPlot(List<double[]> data, String title, String xAxis, String yAxis, String filepath, double minX, double maxX, double minY, double maxY) {
        BufferedImage image = new BufferedImage(SCATTER_WIDTH, SCATTER_HEIGHT, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = image.createGraphics();

        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // White background for scatter plots
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, SCATTER_WIDTH, SCATTER_HEIGHT);

        double xScale = (double) (SCATTER_WIDTH - 2 * SCATTER_PADDING - LABEL_PADDING) / (maxX - minX);
        double yScale = (double) (SCATTER_HEIGHT - 2 * SCATTER_PADDING - LABEL_PADDING) / (maxY - minY);

        // Axes
        g2.setColor(Color.BLACK);
        g2.drawLine(SCATTER_PADDING + LABEL_PADDING, SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING, SCATTER_PADDING + LABEL_PADDING, SCATTER_PADDING);
        g2.drawLine(SCATTER_PADDING + LABEL_PADDING, SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING, SCATTER_WIDTH - SCATTER_PADDING, SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING);

        // Ticks and Labels
        g2.setFont(new Font("Arial", Font.PLAIN, 12));
        FontMetrics metrics = g2.getFontMetrics();
        int numTicks = 5;

        // X Axis
        for (int i = 0; i <= numTicks; i++) {
            double value = minX + (maxX - minX) * i / numTicks;
            int xPos = (int) ((value - minX) * xScale + SCATTER_PADDING + LABEL_PADDING);
            int yPos = SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING;
            g2.drawLine(xPos, yPos, xPos, yPos + 5);
            String label = String.format("%.1f", value);
            int labelWidth = metrics.stringWidth(label);
            g2.drawString(label, xPos - labelWidth / 2, yPos + 20);
        }

        // Y Axis
        for (int i = 0; i <= numTicks; i++) {
            double value = minY + (maxY - minY) * i / numTicks;
            int xPos = SCATTER_PADDING + LABEL_PADDING;
            int yPos = (int) ((maxY - value) * yScale + SCATTER_PADDING);
            g2.drawLine(xPos, yPos, xPos - 5, yPos);
            String label = String.format("%.2f", value);
            int labelWidth = metrics.stringWidth(label);
            g2.drawString(label, xPos - labelWidth - 10, yPos + (metrics.getAscent() / 2));
        }

        // Titles
        g2.setFont(new Font("Arial", Font.PLAIN, 14));
        g2.drawString(xAxis, SCATTER_WIDTH / 2, SCATTER_HEIGHT - SCATTER_PADDING / 2);
        
        g2.rotate(-Math.PI / 2);
        g2.drawString(yAxis, -SCATTER_HEIGHT / 2, SCATTER_PADDING - 10);
        g2.rotate(Math.PI / 2);

        g2.setFont(new Font("Arial", Font.BOLD, 18));
        g2.drawString(title, SCATTER_PADDING, SCATTER_PADDING / 2 + 5);

        // Correlation
        double correlation = calculateCorrelation(data);
        g2.setFont(new Font("Arial", Font.PLAIN, 14));
        String corrString = String.format("Correlation: %.4f", correlation);
        int corrStringWidth = g2.getFontMetrics().stringWidth(corrString);
        g2.drawString(corrString, SCATTER_WIDTH - SCATTER_PADDING - corrStringWidth, SCATTER_PADDING);
        
        // Points
        g2.setColor(new Color(44, 102, 230, 150));
        for (double[] point : data) {
            double xVal = point[0];
            double yVal = point[1];
            int x = (int) ((xVal - minX) * xScale + SCATTER_PADDING + LABEL_PADDING);
            int y = (int) ((maxY - yVal) * yScale + SCATTER_PADDING);
            Shape circle = new Ellipse2D.Double(x - POINT_SIZE / 2.0, y - POINT_SIZE / 2.0, POINT_SIZE, POINT_SIZE);
            g2.fill(circle);
        }

        g2.dispose();

        try {
            File file = new File(filepath);
            file.getParentFile().mkdirs();
            ImageIO.write(image, "png", file);
            System.out.println("Saved Scatter Plot: " + filepath);
        } catch (IOException e) {
            System.err.println("Error saving scatter plot: " + e.getMessage());
        }
    }

    private static double calculateCorrelation(List<double[]> data) {
        if (data.isEmpty()) return 0.0;
        double sumX = 0.0, sumY = 0.0, sumXY = 0.0;
        double sumX2 = 0.0, sumY2 = 0.0;
        int n = data.size();

        for (double[] point : data) {
            sumX += point[0];
            sumY += point[1];
            sumXY += point[0] * point[1];
            sumX2 += point[0] * point[0];
            sumY2 += point[1] * point[1];
        }

        double numerator = n * sumXY - sumX * sumY;
        double denominator = Math.sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));

        return (denominator == 0) ? 0 : numerator / denominator;
    }
}