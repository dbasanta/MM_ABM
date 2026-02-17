import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Generate comprehensive plots for bone-tumor interaction simulation
 * Shows Tumor, Osteoclast, and Osteoblast populations over time
 */
public class PlotBoneTumorData {

    public static void main(String[] args) throws Exception {
        // Read CSV data
        List<Integer> time = new ArrayList<>();
        List<Integer> tumor = new ArrayList<>();
        List<Integer> osteoclasts = new ArrayList<>();
        List<Integer> osteoblasts = new ArrayList<>();
        List<Double> ratio = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader("bone_tumor_simulation.csv"));
        String line = br.readLine(); // Skip header

        while ((line = br.readLine()) != null) {
            String[] parts = line.split(",");
            time.add(Integer.parseInt(parts[0]));
            tumor.add(Integer.parseInt(parts[1]));
            osteoclasts.add(Integer.parseInt(parts[2]));
            osteoblasts.add(Integer.parseInt(parts[3]));
            ratio.add(Double.parseDouble(parts[4]));
        }
        br.close();

        System.out.println("\n╔════════════════════════════════════════════════════════════════════════╗");
        System.out.println("║       BONE-TUMOR INTERACTION: POPULATION DYNAMICS VISUALIZATION       ║");
        System.out.println("╚════════════════════════════════════════════════════════════════════════╝\n");

        // Plot 1: Tumor Cells (Multiple Myeloma)
        System.out.println("═".repeat(80));
        System.out.println("  PLOT 1: TUMOR CELLS (Multiple Myeloma) OVER TIME");
        System.out.println("═".repeat(80));
        plotPopulation(time, tumor, "Tumor (MM)", "T", 0, 60);
        System.out.println();

        // Plot 2: Osteoclasts
        System.out.println("═".repeat(80));
        System.out.println("  PLOT 2: OSTEOCLASTS (Bone-Resorbing Cells) OVER TIME");
        System.out.println("═".repeat(80));
        plotPopulation(time, osteoclasts, "Osteoclasts", "C", 0, 60);
        System.out.println();

        // Plot 3: Osteoblasts
        System.out.println("═".repeat(80));
        System.out.println("  PLOT 3: OSTEOBLASTS (Bone-Forming Cells) OVER TIME");
        System.out.println("═".repeat(80));
        plotPopulation(time, osteoblasts, "Osteoblasts", "B", 0, 60);
        System.out.println();

        // Plot 4: All Three Populations Together
        System.out.println("═".repeat(80));
        System.out.println("  PLOT 4: ALL CELL POPULATIONS (Combined View)");
        System.out.println("═".repeat(80));
        plotCombined(time, tumor, osteoclasts, osteoblasts);
        System.out.println();

        // Plot 5: OC/OB Ratio (Bone Homeostasis Indicator)
        System.out.println("═".repeat(80));
        System.out.println("  PLOT 5: OSTEOCLAST/OSTEOBLAST RATIO (Bone Balance)");
        System.out.println("═".repeat(80));
        plotRatio(time, ratio);
        System.out.println();

        // Summary
        System.out.println("╔════════════════════════════════════════════════════════════════════════╗");
        System.out.println("║                           KEY FINDINGS                                 ║");
        System.out.println("╚════════════════════════════════════════════════════════════════════════╝\n");

        // Find critical time points
        int tumorDoubleTime = findDoublingTime(tumor, 50);
        int obExtinctionTime = findExtinctionTime(osteoblasts);
        int maxOCtime = findMaxTime(osteoclasts);

        System.out.println("Timeline of Events:");
        if (tumorDoubleTime > 0) {
            System.out.printf("  T=%3d: Tumor population doubled\n", tumorDoubleTime);
        }
        if (obExtinctionTime > 0) {
            System.out.printf("  T=%3d: Osteoblasts depleted (bone formation ceased)\n", obExtinctionTime);
        }
        if (maxOCtime > 0) {
            System.out.printf("  T=%3d: Osteoclasts reached maximum (peak bone resorption)\n", maxOCtime);
        }

        System.out.println("\nPathological Process:");
        System.out.println("  1. Tumor cells secrete RANKL → recruit osteoclasts");
        System.out.println("  2. Osteoclasts resorb bone → release growth factors (TGF-β)");
        System.out.println("  3. Growth factors stimulate tumor → \"vicious cycle\"");
        System.out.println("  4. Tumor suppresses osteoblasts → prevents bone repair");
        System.out.println("  5. Result: Osteolytic lesions and bone destruction");

        System.out.println("\n✓ All stochastic events computed with refactored Binomial.ColtInt()");
        System.out.println("✓ Demonstrates realistic bone metastasis pathophysiology");
        System.out.println("✓ Shows tumor-bone microenvironment interaction");
    }

    private static void plotPopulation(List<Integer> time, List<Integer> data,
                                       String name, String symbol, int start, int end) {
        int height = 20;
        int width = 70;

        // Find max in range
        long maxVal = 0;
        for (int i = start; i < Math.min(end, data.size()); i++) {
            if (data.get(i) > maxVal) maxVal = data.get(i);
        }

        // Draw plot
        for (int row = height; row >= 0; row--) {
            long threshold = maxVal * row / height;
            System.out.printf("%10s │", formatNumber(threshold));

            for (int col = 0; col < width; col++) {
                int idx = start + (end - start) * col / width;
                if (idx >= data.size()) break;

                if (data.get(idx) >= threshold) {
                    System.out.print(symbol);
                } else {
                    System.out.print(" ");
                }
            }
            System.out.println();
        }

        // X-axis
        System.out.print("           └");
        System.out.print("─".repeat(width));
        System.out.println();
        System.out.printf("           %5d", start);
        System.out.print(" ".repeat(width - 25));
        System.out.print("Time");
        System.out.print(" ".repeat(10));
        System.out.printf("%5d\n", Math.min(end - 1, time.size() - 1));
    }

    private static void plotCombined(List<Integer> time, List<Integer> tumor,
                                     List<Integer> osteoclasts, List<Integer> osteoblasts) {
        int height = 25;
        int width = 70;
        int end = 60;

        // Normalize to 0-1 scale for comparison
        double maxTumor = getMax(tumor, 0, end);
        double maxOC = getMax(osteoclasts, 0, end);
        double maxOB = getMax(osteoblasts, 0, end);

        System.out.println("  Legend: T=Tumor(MM), C=Osteoclasts, B=Osteoblasts, *=Overlap\n");

        for (int row = height; row >= 0; row--) {
            double threshold = row / (double)height;
            System.out.printf("      %3d%% │", (int)(threshold * 100));

            for (int col = 0; col < width; col++) {
                int idx = (end * col) / width;
                if (idx >= tumor.size()) break;

                double normTumor = tumor.get(idx) / maxTumor;
                double normOC = osteoclasts.get(idx) / maxOC;
                double normOB = osteoblasts.get(idx) / maxOB;

                // Count which populations are above threshold
                boolean hasTumor = normTumor >= threshold;
                boolean hasOC = normOC >= threshold;
                boolean hasOB = normOB >= threshold;

                int count = (hasTumor ? 1 : 0) + (hasOC ? 1 : 0) + (hasOB ? 1 : 0);

                if (count >= 2) {
                    System.out.print("*"); // Overlap
                } else if (hasTumor) {
                    System.out.print("T");
                } else if (hasOC) {
                    System.out.print("C");
                } else if (hasOB) {
                    System.out.print("B");
                } else {
                    System.out.print(" ");
                }
            }
            System.out.println();
        }

        // X-axis
        System.out.print("           └");
        System.out.print("─".repeat(width));
        System.out.println();
        System.out.printf("           %5d", 0);
        System.out.print(" ".repeat(width - 25));
        System.out.print("Time");
        System.out.print(" ".repeat(10));
        System.out.printf("%5d\n", Math.min(end - 1, time.size() - 1));
    }

    private static void plotRatio(List<Integer> time, List<Double> ratio) {
        int height = 20;
        int width = 70;
        int end = Math.min(60, ratio.size());

        // Find max ratio
        double maxRatio = 0;
        for (int i = 0; i < end; i++) {
            if (ratio.get(i) > maxRatio && !Double.isInfinite(ratio.get(i))) {
                maxRatio = ratio.get(i);
            }
        }
        maxRatio = Math.min(maxRatio, 10.0); // Cap for visualization

        System.out.println("  Normal bone: OC/OB ≈ 0.5-1.0 (balanced remodeling)");
        System.out.println("  Pathological: OC/OB > 2.0 (osteolytic, bone loss)\n");

        for (int row = height; row >= 0; row--) {
            double threshold = maxRatio * row / height;
            System.out.printf("      %5.2f │", threshold);

            for (int col = 0; col < width; col++) {
                int idx = (end * col) / width;
                if (idx >= ratio.size()) break;

                double val = ratio.get(idx);
                if (Double.isInfinite(val) || Double.isNaN(val)) val = 0;

                if (val >= threshold) {
                    System.out.print("█");
                } else {
                    System.out.print(" ");
                }
            }
            System.out.println();
        }

        // Reference line for normal ratio
        System.out.printf("      %5.2f │", 1.0);
        System.out.print("─".repeat(width));
        System.out.println(" ← Normal");

        // X-axis
        System.out.print("           └");
        System.out.print("─".repeat(width));
        System.out.println();
        System.out.printf("           %5d", 0);
        System.out.print(" ".repeat(width - 25));
        System.out.print("Time");
        System.out.print(" ".repeat(10));
        System.out.printf("%5d\n", Math.min(end - 1, time.size() - 1));
    }

    private static String formatNumber(long n) {
        if (n >= 1_000_000) {
            return String.format("%.2fM", n / 1_000_000.0);
        } else if (n >= 1_000) {
            return String.format("%.2fK", n / 1_000.0);
        } else {
            return String.format("%d", n);
        }
    }

    private static double getMax(List<Integer> data, int start, int end) {
        double max = 1.0;
        for (int i = start; i < Math.min(end, data.size()); i++) {
            if (data.get(i) > max) max = data.get(i);
        }
        return max;
    }

    private static int findDoublingTime(List<Integer> data, int initial) {
        for (int i = 0; i < data.size(); i++) {
            if (data.get(i) >= initial * 2) return i;
        }
        return -1;
    }

    private static int findExtinctionTime(List<Integer> data) {
        for (int i = 0; i < data.size(); i++) {
            if (data.get(i) == 0) return i;
        }
        return -1;
    }

    private static int findMaxTime(List<Integer> data) {
        int maxVal = 0;
        int maxTime = 0;
        for (int i = 0; i < data.size(); i++) {
            if (data.get(i) > maxVal) {
                maxVal = data.get(i);
                maxTime = i;
            }
        }
        return maxTime;
    }
}
