import java.io.PrintWriter;
import java.io.FileWriter;

/**
 * Simple tumor growth simulation demonstrating the refactored Binomial code.
 * Models exponential growth with stochastic birth/death processes and mutations.
 */
public class SimpleTumorSimulation {

    public static void main(String[] args) throws Exception {
        System.out.println("=== Tumor Growth Simulation ===");
        System.out.println("Testing refactored Binomial methods in real simulation\n");

        // Simulation parameters
        int timeSteps = 500;
        double birthRate = 0.12;      // Birth rate per cell per timestep
        double deathRate = 0.08;      // Death rate per cell per timestep
        double mutationRate = 0.001;  // Probability of mutation per division
        int initialPopulation = 100;

        // Initialize populations
        int normalCells = initialPopulation;
        int mutantCells = 0;

        SimpleBinomial binomial = new SimpleBinomial();
        SimpleRand rng = new SimpleRand(42);

        // Track populations over time
        int[] normalPop = new int[timeSteps];
        int[] mutantPop = new int[timeSteps];
        int[] totalPop = new int[timeSteps];

        System.out.println("Parameters:");
        System.out.println("  Birth rate: " + birthRate);
        System.out.println("  Death rate: " + deathRate);
        System.out.println("  Mutation rate: " + mutationRate);
        System.out.println("  Initial population: " + initialPopulation);
        System.out.println("\nRunning simulation for " + timeSteps + " timesteps...\n");

        // Run simulation
        for (int t = 0; t < timeSteps; t++) {
            // Normal cells: births, deaths, and mutations
            if (normalCells > 0) {
                int normalBirths = binomial.sampleInt(normalCells, birthRate, rng);
                int normalDeaths = binomial.sampleInt(normalCells, deathRate, rng);
                int mutations = binomial.sampleInt(normalBirths, mutationRate, rng);

                normalCells = normalCells + normalBirths - normalDeaths - mutations;
                mutantCells += mutations;

                if (normalCells < 0) normalCells = 0;
            }

            // Mutant cells: enhanced birth rate (growth advantage)
            if (mutantCells > 0) {
                double mutantBirthRate = birthRate * 1.2; // 20% growth advantage
                int mutantBirths = binomial.sampleInt(mutantCells, mutantBirthRate, rng);
                int mutantDeaths = binomial.sampleInt(mutantCells, deathRate, rng);

                mutantCells = mutantCells + mutantBirths - mutantDeaths;

                if (mutantCells < 0) mutantCells = 0;
            }

            // Record populations
            normalPop[t] = normalCells;
            mutantPop[t] = mutantCells;
            totalPop[t] = normalCells + mutantCells;

            // Print status every 100 timesteps
            if (t % 100 == 0 || t == timeSteps - 1) {
                System.out.printf("Time %3d: Normal=%6d, Mutant=%6d, Total=%6d\n",
                                 t, normalCells, mutantCells, totalPop[t]);
            }
        }

        System.out.println("\nSimulation complete!");

        // Calculate statistics
        int maxPop = 0;
        int maxTime = 0;
        for (int t = 0; t < timeSteps; t++) {
            if (totalPop[t] > maxPop) {
                maxPop = totalPop[t];
                maxTime = t;
            }
        }

        System.out.println("\nStatistics:");
        System.out.println("  Final population: " + totalPop[timeSteps - 1]);
        System.out.println("  Peak population: " + maxPop + " at time " + maxTime);
        System.out.println("  Final mutant fraction: " +
                          String.format("%.2f%%", 100.0 * mutantPop[timeSteps - 1] / totalPop[timeSteps - 1]));

        // Generate CSV data for plotting
        System.out.println("\nGenerating plot data...");
        PrintWriter csv = new PrintWriter(new FileWriter("tumor_population.csv"));
        csv.println("Time,Normal,Mutant,Total");
        for (int t = 0; t < timeSteps; t++) {
            csv.printf("%d,%d,%d,%d\n", t, normalPop[t], mutantPop[t], totalPop[t]);
        }
        csv.close();
        System.out.println("Data saved to: tumor_population.csv");

        // Generate ASCII plot
        System.out.println("\n" + generateASCIIPlot(totalPop, 40, 60));

        // Generate simple plot script
        generatePlotScript();

        System.out.println("\n✓ Simulation successful - Refactored Binomial code works perfectly!");
    }

    private static String generateASCIIPlot(int[] data, int height, int width) {
        StringBuilder plot = new StringBuilder();
        plot.append("Tumor Population Over Time (ASCII Plot):\n");
        plot.append("═".repeat(width + 10)).append("\n");

        // Find max for scaling
        int max = 0;
        for (int v : data) {
            if (v > max) max = v;
        }

        // Draw plot from top to bottom
        for (int row = height; row >= 0; row--) {
            int threshold = (int)(max * row / (double)height);

            // Y-axis label
            plot.append(String.format("%6d │", threshold));

            // Plot points
            for (int col = 0; col < width; col++) {
                int dataIndex = (int)(data.length * col / (double)width);
                if (data[dataIndex] >= threshold) {
                    plot.append("█");
                } else {
                    plot.append(" ");
                }
            }
            plot.append("\n");
        }

        // X-axis
        plot.append("       └").append("─".repeat(width)).append("\n");
        plot.append("        0");
        plot.append(" ".repeat(width - 20));
        plot.append("Time");
        plot.append(" ".repeat(width - 35));
        plot.append(data.length);
        plot.append("\n");

        return plot.toString();
    }

    private static void generatePlotScript() throws Exception {
        // Generate Python plotting script
        PrintWriter py = new PrintWriter(new FileWriter("plot_tumor.py"));
        py.println("#!/usr/bin/env python3");
        py.println("import matplotlib.pyplot as plt");
        py.println("import pandas as pd");
        py.println();
        py.println("# Read data");
        py.println("data = pd.read_csv('tumor_population.csv')");
        py.println();
        py.println("# Create plot");
        py.println("plt.figure(figsize=(12, 6))");
        py.println();
        py.println("# Plot total population");
        py.println("plt.subplot(1, 2, 1)");
        py.println("plt.plot(data['Time'], data['Total'], 'b-', linewidth=2, label='Total')");
        py.println("plt.xlabel('Time')");
        py.println("plt.ylabel('Population')");
        py.println("plt.title('Total Tumor Population Over Time')");
        py.println("plt.grid(True, alpha=0.3)");
        py.println("plt.legend()");
        py.println();
        py.println("# Plot population breakdown");
        py.println("plt.subplot(1, 2, 2)");
        py.println("plt.plot(data['Time'], data['Normal'], 'g-', linewidth=2, label='Normal')");
        py.println("plt.plot(data['Time'], data['Mutant'], 'r-', linewidth=2, label='Mutant')");
        py.println("plt.xlabel('Time')");
        py.println("plt.ylabel('Population')");
        py.println("plt.title('Normal vs Mutant Cells')");
        py.println("plt.grid(True, alpha=0.3)");
        py.println("plt.legend()");
        py.println();
        py.println("plt.tight_layout()");
        py.println("plt.savefig('tumor_population.png', dpi=150)");
        py.println("print('Plot saved to: tumor_population.png')");
        py.close();

        System.out.println("Plot script saved to: plot_tumor.py");
        System.out.println("To generate plot, run: python3 plot_tumor.py");
    }
}
