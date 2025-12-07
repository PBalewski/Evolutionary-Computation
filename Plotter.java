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

    // Stałe dla rysowania tras TSP
    private static final int PADDING = 20;
    private static final int NODE_SIZE = 8;

    public Plotter(List<Classes.Node> nodes, List<Integer> path) {
        this.nodes = nodes;
        this.path = path;
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        if (nodes == null || nodes.isEmpty()) return;

        // Skalowanie dla trasy TSP
        int minX = Integer.MAX_VALUE, maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE, maxY = Integer.MIN_VALUE;

        for (Classes.Node node : nodes) {
            if (node.x < minX) minX = node.x;
            if (node.x > maxX) maxX = node.x;
            if (node.y < minY) minY = node.y;
            if (node.y > maxY) maxY = node.y;
        }

        double scaleX = (double) (getWidth() - 2 * PADDING) / (maxX - minX);
        double scaleY = (double) (getHeight() - 2 * PADDING) / (maxY - minY);

        // Rysowanie krawędzi
        g2.setColor(Color.LIGHT_GRAY);
        if (path != null && path.size() > 1) {
            for (int i = 0; i < path.size() - 1; i++) {
                Classes.Node n1 = nodes.get(path.get(i));
                Classes.Node n2 = nodes.get(path.get(i + 1));
                drawEdge(g2, n1, n2, minX, minY, scaleX, scaleY);
            }
            // Zamknięcie cyklu
            Classes.Node last = nodes.get(path.get(path.size() - 1));
            Classes.Node first = nodes.get(path.get(0));
            drawEdge(g2, last, first, minX, minY, scaleX, scaleY);
        }

        // Rysowanie węzłów
        for (int i = 0; i < nodes.size(); i++) {
            Classes.Node node = nodes.get(i);
            int x = PADDING + (int) ((node.x - minX) * scaleX);
            int y = PADDING + (int) ((node.y - minY) * scaleY);

            if (path != null && path.contains(i)) {
                g2.setColor(Color.RED);
                g2.fillOval(x - NODE_SIZE / 2, y - NODE_SIZE / 2, NODE_SIZE, NODE_SIZE);
            } else {
                g2.setColor(Color.GRAY);
                g2.fillOval(x - NODE_SIZE / 4, y - NODE_SIZE / 4, NODE_SIZE / 2, NODE_SIZE / 2);
            }
        }
    }

    private void drawEdge(Graphics2D g2, Classes.Node n1, Classes.Node n2, int minX, int minY, double scaleX, double scaleY) {
        int x1 = PADDING + (int) ((n1.x - minX) * scaleX);
        int y1 = PADDING + (int) ((n1.y - minY) * scaleY);
        int x2 = PADDING + (int) ((n2.x - minX) * scaleX);
        int y2 = PADDING + (int) ((n2.y - minY) * scaleY);
        g2.drawLine(x1, y1, x2, y2);
    }

    // --- STARE METODY (DLA TSP) ---

    public static void showPlot(List<Classes.Node> nodes, List<Integer> path, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(800, 800);
        frame.add(new Plotter(nodes, path));
        frame.setVisible(true);
    }

    public static void savePlotAsImage(Plotter plotter, String filepath) {
        BufferedImage image = new BufferedImage(800, 800, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = image.createGraphics();
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, 800, 800);
        plotter.setSize(800, 800);
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

    // =========================================================
    // === NOWE METODY (DLA ZADANIA GLOBAL CONVEXITY / SCATTER PLOT) ===
    // =========================================================

    private static final int SCATTER_WIDTH = 800;
    private static final int SCATTER_HEIGHT = 600;
    private static final int SCATTER_PADDING = 60;
    private static final int LABEL_PADDING = 30;
    private static final int POINT_SIZE = 6;

    /**
     * Rysuje wykres punktowy (Scatter Plot) i zapisuje do pliku.
     * Używane do wizualizacji korelacji Fitness-Distance.
     * * @param data      Lista tablic double[], gdzie [0]=Wartość funkcji celu, [1]=Podobieństwo
     * @param title     Tytuł wykresu
     * @param xAxis     Opis osi X
     * @param yAxis     Opis osi Y
     * @param filepath  Ścieżka do pliku wyjściowego
     */
    public static void saveScatterPlot(List<double[]> data, String title, String xAxis, String yAxis, String filepath) {
        // Znajdź min/max dla skali
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

        // Antyaliasing
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // Białe tło
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, SCATTER_WIDTH, SCATTER_HEIGHT);

        double xScale = (double) (SCATTER_WIDTH - 2 * SCATTER_PADDING - LABEL_PADDING) / (maxX - minX);
        double yScale = (double) (SCATTER_HEIGHT - 2 * SCATTER_PADDING - LABEL_PADDING) / (maxY - minY);

        // Rysowanie osi
        g2.setColor(Color.BLACK);
        // Oś Y
        g2.drawLine(SCATTER_PADDING + LABEL_PADDING, SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING, SCATTER_PADDING + LABEL_PADDING, SCATTER_PADDING);
        // Oś X
        g2.drawLine(SCATTER_PADDING + LABEL_PADDING, SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING, SCATTER_WIDTH - SCATTER_PADDING, SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING);

        // Rysowanie skali i etykiet na osiach
        g2.setFont(new Font("Arial", Font.PLAIN, 12));
        FontMetrics metrics = g2.getFontMetrics();
        int numTicks = 5;

        // Skala osi X
        for (int i = 0; i <= numTicks; i++) {
            double value = minX + (maxX - minX) * i / numTicks;
            int xPos = (int) ((value - minX) * xScale + SCATTER_PADDING + LABEL_PADDING);
            int yPos = SCATTER_HEIGHT - SCATTER_PADDING - LABEL_PADDING;
            g2.drawLine(xPos, yPos, xPos, yPos + 5); // Tick mark
            String label = String.format("%.1f", value);
            int labelWidth = metrics.stringWidth(label);
            g2.drawString(label, xPos - labelWidth / 2, yPos + 20);
        }

        // Skala osi Y
        for (int i = 0; i <= numTicks; i++) {
            double value = minY + (maxY - minY) * i / numTicks;
            int xPos = SCATTER_PADDING + LABEL_PADDING;
            int yPos = (int) ((maxY - value) * yScale + SCATTER_PADDING);
            g2.drawLine(xPos, yPos, xPos - 5, yPos); // Tick mark
            String label = String.format("%.2f", value);
            int labelWidth = metrics.stringWidth(label);
            g2.drawString(label, xPos - labelWidth - 10, yPos + (metrics.getAscent() / 2));
        }


        // Opisy osi
        g2.setFont(new Font("Arial", Font.PLAIN, 14));
        g2.drawString(xAxis, SCATTER_WIDTH / 2, SCATTER_HEIGHT - SCATTER_PADDING / 2);
        
        // Obrócony opis osi Y
        g2.rotate(-Math.PI / 2);
        g2.drawString(yAxis, -SCATTER_HEIGHT / 2, SCATTER_PADDING - 10);
        g2.rotate(Math.PI / 2);

        // Tytuł
        g2.setFont(new Font("Arial", Font.BOLD, 18));
        g2.drawString(title, SCATTER_PADDING, SCATTER_PADDING / 2 + 5);

        // Korelacja
        double correlation = calculateCorrelation(data);
        g2.setFont(new Font("Arial", Font.PLAIN, 14));
        String corrString = String.format("Correlation: %.4f", correlation);
        int corrStringWidth = g2.getFontMetrics().stringWidth(corrString);
        g2.drawString(corrString, SCATTER_WIDTH - SCATTER_PADDING - corrStringWidth, SCATTER_PADDING);
        
        // Rysowanie punktów
        g2.setColor(new Color(44, 102, 230, 150)); // Niebieski, półprzezroczysty
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

    /**
     * Pomocnicza metoda do obliczania korelacji Pearsona.
     */
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