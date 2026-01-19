import HAL.Rand;
import HAL.Tools.Internal.Binomial;

/**
 * Standalone test to verify refactored Binomial methods produce correct results.
 */
public class BinomialRefactorTest {

    public static void main(String[] args) {
        System.out.println("Testing refactored Binomial class...");

        boolean allPassed = true;

        // Test 1: Small n, various p values
        allPassed &= testBinomial(10, 0.5, 100000, "Small n, p=0.5");
        allPassed &= testBinomial(10, 0.1, 100000, "Small n, p=0.1");
        allPassed &= testBinomial(10, 0.9, 100000, "Small n, p=0.9");

        // Test 2: Medium n, various p values
        allPassed &= testBinomial(100, 0.5, 100000, "Medium n, p=0.5");
        allPassed &= testBinomial(100, 0.01, 100000, "Medium n, p=0.01");
        allPassed &= testBinomial(100, 0.99, 100000, "Medium n, p=0.99");

        // Test 3: Large n, various p values
        allPassed &= testBinomial(10000, 0.5, 100000, "Large n, p=0.5");
        allPassed &= testBinomial(10000, 0.001, 100000, "Large n, p=0.001");
        allPassed &= testBinomial(10000, 0.999, 100000, "Large n, p=0.999");

        // Test 4: Very large n (long version)
        allPassed &= testBinomialLong(1000000L, 0.5, 100000, "Very large n, p=0.5");
        allPassed &= testBinomialLong(1000000L, 0.0001, 100000, "Very large n, p=0.0001");

        // Test 5: Edge cases
        allPassed &= testBinomial(5, 0.0, 10000, "Edge case: p=0");
        allPassed &= testBinomial(5, 1.0, 10000, "Edge case: p=1");

        System.out.println("\n" + (allPassed ? "ALL TESTS PASSED!" : "SOME TESTS FAILED!"));
        System.exit(allPassed ? 0 : 1);
    }

    private static boolean testBinomial(int n, double p, int trials, String testName) {
        Rand rng = new Rand(12345); // Fixed seed for reproducibility

        // Collect samples
        long sum = 0;
        long sumSquares = 0;
        for (int i = 0; i < trials; i++) {
            int sample = rng.Binomial(n, p);
            sum += sample;
            sumSquares += (long)sample * sample;
        }

        // Calculate statistics
        double mean = (double)sum / trials;
        double variance = ((double)sumSquares / trials) - (mean * mean);
        double stdDev = Math.sqrt(variance);

        // Expected values
        double expectedMean = n * p;
        double expectedStdDev = Math.sqrt(n * p * (1 - p));

        // Tolerance (within 3 standard errors)
        double meanTolerance = expectedStdDev / Math.sqrt(trials) * 3;
        double stdDevTolerance = expectedStdDev * 0.05; // 5% tolerance for stddev

        boolean meanOk = Math.abs(mean - expectedMean) <= meanTolerance;
        boolean stdDevOk = Math.abs(stdDev - expectedStdDev) <= stdDevTolerance;

        System.out.printf("%-30s: ", testName);
        System.out.printf("Mean %.4f (expected %.4f, diff %.4f) %s | ",
            mean, expectedMean, Math.abs(mean - expectedMean), meanOk ? "✓" : "✗");
        System.out.printf("StdDev %.4f (expected %.4f, diff %.4f) %s\n",
            stdDev, expectedStdDev, Math.abs(stdDev - expectedStdDev), stdDevOk ? "✓" : "✗");

        return meanOk && stdDevOk;
    }

    private static boolean testBinomialLong(long n, double p, int trials, String testName) {
        Rand rng = new Rand(12345); // Fixed seed for reproducibility

        // Collect samples
        double sum = 0;
        double sumSquares = 0;
        for (int i = 0; i < trials; i++) {
            long sample = rng.Binomial(n, p);
            sum += sample;
            sumSquares += (double)sample * sample;
        }

        // Calculate statistics
        double mean = sum / trials;
        double variance = (sumSquares / trials) - (mean * mean);
        double stdDev = Math.sqrt(variance);

        // Expected values
        double expectedMean = n * p;
        double expectedStdDev = Math.sqrt(n * p * (1 - p));

        // Tolerance (within 3 standard errors)
        double meanTolerance = expectedStdDev / Math.sqrt(trials) * 3;
        double stdDevTolerance = expectedStdDev * 0.05; // 5% tolerance for stddev

        boolean meanOk = Math.abs(mean - expectedMean) <= meanTolerance;
        boolean stdDevOk = Math.abs(stdDev - expectedStdDev) <= stdDevTolerance;

        System.out.printf("%-30s: ", testName);
        System.out.printf("Mean %.4f (expected %.4f, diff %.4f) %s | ",
            mean, expectedMean, Math.abs(mean - expectedMean), meanOk ? "✓" : "✗");
        System.out.printf("StdDev %.4f (expected %.4f, diff %.4f) %s\n",
            stdDev, expectedStdDev, Math.abs(stdDev - expectedStdDev), stdDevOk ? "✓" : "✗");

        return meanOk && stdDevOk;
    }
}
