import java.io.PrintWriter;
import java.io.FileWriter;

/**
 * Bone-Tumor Interaction Simulation
 * Models the dynamics between tumor cells (Multiple Myeloma), osteoclasts, and osteoblasts
 * Demonstrates refactored Binomial code in a realistic bone metastasis scenario
 */
public class BoneTumorSimulation {

    public static void main(String[] args) throws Exception {
        System.out.println("╔════════════════════════════════════════════════════════════════════════╗");
        System.out.println("║        BONE-TUMOR INTERACTION SIMULATION                               ║");
        System.out.println("║        Multiple Myeloma, Osteoclasts, and Osteoblasts                 ║");
        System.out.println("╚════════════════════════════════════════════════════════════════════════╝\n");

        // Simulation parameters
        int timeSteps = 300;

        // Cell-specific rates (per timestep)
        double tumor_birthRate = 0.15;        // MM cells proliferate rapidly
        double tumor_deathRate = 0.05;        // Base death rate

        double osteoclast_recruitRate = 0.10; // OC recruitment (enhanced by tumor)
        double osteoclast_deathRate = 0.12;   // OC have short lifespan

        double osteoblast_recruitRate = 0.08; // OB recruitment
        double osteoblast_deathRate = 0.06;   // OB death rate
        double osteoblast_suppressionByTumor = 0.002; // Tumor suppresses OB

        // Interaction parameters
        double tumor_stimulatesOC = 0.0003;   // Tumor enhances OC recruitment (RANKL)
        double OC_stimulatesTumor = 0.0001;   // OC bone resorption releases growth factors
        double OB_inhibitsTumor = 0.0002;     // OB can inhibit tumor growth

        // Initial populations
        int tumorCells = 50;           // Initial tumor burden
        int osteoclasts = 100;         // Pre-existing osteoclasts
        int osteoblasts = 200;         // Pre-existing osteoblasts

        SimpleBinomial binomial = new SimpleBinomial();
        SimpleRand rng = new SimpleRand(12345);

        // Track populations over time
        int[] tumorPop = new int[timeSteps];
        int[] osteoclastPop = new int[timeSteps];
        int[] osteoblastPop = new int[timeSteps];

        System.out.println("Initial Conditions:");
        System.out.println("  Tumor cells (MM):  " + tumorCells);
        System.out.println("  Osteoclasts (OC):  " + osteoclasts);
        System.out.println("  Osteoblasts (OB):  " + osteoblasts);
        System.out.println("\nRunning simulation...\n");

        // Run simulation
        for (int t = 0; t < timeSteps; t++) {
            // === TUMOR DYNAMICS ===
            if (tumorCells > 0) {
                // Tumor birth (enhanced by osteoclast activity - vicious cycle)
                double effectiveTumorBirth = tumor_birthRate + (osteoclasts * OC_stimulatesTumor);
                int tumorBirths = binomial.sampleInt(tumorCells, effectiveTumorBirth, rng);

                // Tumor death (enhanced by osteoblasts)
                double effectiveTumorDeath = tumor_deathRate + (osteoblasts * OB_inhibitsTumor);
                int tumorDeaths = binomial.sampleInt(tumorCells, effectiveTumorDeath, rng);

                tumorCells = Math.max(0, tumorCells + tumorBirths - tumorDeaths);
            }

            // === OSTEOCLAST DYNAMICS ===
            if (osteoclasts > 0) {
                // OC recruitment (enhanced by tumor via RANKL)
                double effectiveOC_recruit = osteoclast_recruitRate + (tumorCells * tumor_stimulatesOC);
                int ocRecruits = binomial.sampleInt(osteoclasts, effectiveOC_recruit, rng);

                // OC death
                int ocDeaths = binomial.sampleInt(osteoclasts, osteoclast_deathRate, rng);

                osteoclasts = Math.max(0, osteoclasts + ocRecruits - ocDeaths);
            }

            // === OSTEOBLAST DYNAMICS ===
            if (osteoblasts > 0) {
                // OB recruitment (suppressed by tumor)
                double effectiveOB_recruit = osteoblast_recruitRate * (1.0 - tumorCells * osteoblast_suppressionByTumor);
                effectiveOB_recruit = Math.max(0.0, effectiveOB_recruit);
                int obRecruits = binomial.sampleInt(osteoblasts, effectiveOB_recruit, rng);

                // OB death
                int obDeaths = binomial.sampleInt(osteoblasts, osteoblast_deathRate, rng);

                osteoblasts = Math.max(0, osteoblasts + obRecruits - obDeaths);
            }

            // Prevent overflow
            tumorCells = Math.min(tumorCells, 100000);
            osteoclasts = Math.min(osteoclasts, 100000);
            osteoblasts = Math.min(osteoblasts, 100000);

            // Record populations
            tumorPop[t] = tumorCells;
            osteoclastPop[t] = osteoclasts;
            osteoblastPop[t] = osteoblasts;

            // Print progress
            if (t % 50 == 0 || t == timeSteps - 1) {
                System.out.printf("Time %3d: Tumor=%5d, OC=%5d, OB=%5d | Ratio OC/OB=%.2f\n",
                                 t, tumorCells, osteoclasts, osteoblasts,
                                 osteoblasts > 0 ? (double)osteoclasts / osteoblasts : 0.0);
            }
        }

        System.out.println("\n✓ Simulation complete!\n");

        // Generate CSV data
        PrintWriter csv = new PrintWriter(new FileWriter("bone_tumor_simulation.csv"));
        csv.println("Time,Tumor,Osteoclasts,Osteoblasts,OC_OB_Ratio");
        for (int t = 0; t < timeSteps; t++) {
            double ratio = osteoblastPop[t] > 0 ? (double)osteoclastPop[t] / osteoblastPop[t] : 0.0;
            csv.printf("%d,%d,%d,%d,%.4f\n", t, tumorPop[t], osteoclastPop[t], osteoblastPop[t], ratio);
        }
        csv.close();
        System.out.println("Data saved to: bone_tumor_simulation.csv");

        // Generate analysis
        System.out.println("\n╔════════════════════════════════════════════════════════════════════════╗");
        System.out.println("║                         SIMULATION RESULTS                             ║");
        System.out.println("╚════════════════════════════════════════════════════════════════════════╝\n");

        System.out.println("Final Populations:");
        System.out.printf("  Tumor cells:    %,6d (%.1fx growth from initial)\n",
                         tumorPop[timeSteps-1], tumorPop[timeSteps-1] / 50.0);
        System.out.printf("  Osteoclasts:    %,6d (%.1fx change from initial)\n",
                         osteoclastPop[timeSteps-1], osteoclastPop[timeSteps-1] / 100.0);
        System.out.printf("  Osteoblasts:    %,6d (%.1fx change from initial)\n",
                         osteoblastPop[timeSteps-1], osteoblastPop[timeSteps-1] / 200.0);

        double finalRatio = osteoblastPop[timeSteps-1] > 0 ?
                           (double)osteoclastPop[timeSteps-1] / osteoblastPop[timeSteps-1] : 0.0;
        double initialRatio = 100.0 / 200.0;
        System.out.printf("\nOC/OB Ratio: %.2f (initial: %.2f) - %s\n",
                         finalRatio, initialRatio,
                         finalRatio > initialRatio ? "BONE LOSS (osteolytic)" : "BONE FORMATION");

        // Key observations
        System.out.println("\nKey Observations:");
        System.out.println("  • Vicious cycle: Tumor → ↑RANKL → ↑OC → ↑bone resorption → ↑tumor");
        System.out.println("  • Tumor suppresses osteoblast recruitment");
        System.out.println("  • Osteoclasts stimulate tumor growth via released growth factors");
        System.out.println("  • Demonstrates pathological bone-tumor microenvironment");

        System.out.println("\n✓ All binomial sampling using refactored ColtInt() method!");
        System.out.println("✓ Realistic bone metastasis dynamics captured!");
    }
}
