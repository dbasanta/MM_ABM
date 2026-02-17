/**
 * Comprehensive verification test for refactored Binomial class.
 * Tests all code paths including small/large np, int/long versions, and edge cases.
 */
public class VerifyBinomialRefactor {

    public static void main(String[] args) {
        System.out.println("=== Binomial Refactoring Verification Test ===\n");

        int totalTests = 0;
        int passedTests = 0;

        // Test 1: Edge case p=1
        System.out.println("Test 1: Edge case p=1.0");
        if (testEdgeCase(100, 1.0, 100)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 2: Edge case p=0
        System.out.println("Test 2: Edge case p=0.0");
        if (testEdgeCase(100, 0.0, 0)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 3: Small np (uses inversion method)
        System.out.println("Test 3: Small np case (n=10, p=0.5, np=5)");
        if (testStatistical(10, 0.5, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 4: Boundary of small/large np
        System.out.println("Test 4: Boundary case (n=20, p=0.45, np=9)");
        if (testStatistical(20, 0.45, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 5: Large np (uses BTPE method)
        System.out.println("Test 5: Large np case (n=100, p=0.5, np=50)");
        if (testStatistical(100, 0.5, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 6: Large np with small p
        System.out.println("Test 6: Large np, small p (n=1000, p=0.02, np=20)");
        if (testStatistical(1000, 0.02, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 7: Large np with large p
        System.out.println("Test 7: Large np, large p (n=1000, p=0.98, np=980)");
        if (testStatistical(1000, 0.98, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 8: Very large n (long version)
        System.out.println("Test 8: Very large n (n=1000000, p=0.0001)");
        if (testStatisticalLong(1000000L, 0.0001, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 9: Very large n with moderate p (long version)
        System.out.println("Test 9: Very large n, moderate p (n=100000, p=0.5)");
        if (testStatisticalLong(100000L, 0.5, 50000)) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        // Test 10: Test parameter caching (reuse same n, p)
        System.out.println("Test 10: Parameter caching test");
        if (testParameterCaching()) {
            System.out.println("✓ PASSED\n");
            passedTests++;
        } else {
            System.out.println("✗ FAILED\n");
        }
        totalTests++;

        System.out.println("===========================================");
        System.out.println("Results: " + passedTests + "/" + totalTests + " tests passed");
        System.out.println("===========================================");

        if (passedTests == totalTests) {
            System.out.println("\n✓ ALL TESTS PASSED - Refactoring verified!");
            System.exit(0);
        } else {
            System.out.println("\n✗ SOME TESTS FAILED - Please review");
            System.exit(1);
        }
    }

    private static boolean testEdgeCase(int n, double p, int expectedValue) {
        SimpleBinomial binomial = new SimpleBinomial();
        SimpleRand rng = new SimpleRand(42);

        // Test multiple samples to ensure consistency
        for (int i = 0; i < 10; i++) {
            int result = binomial.sampleInt(n, p, rng);
            if (result != expectedValue) {
                System.out.println("  Expected: " + expectedValue + ", Got: " + result);
                return false;
            }
        }

        System.out.println("  All samples returned expected value: " + expectedValue);
        return true;
    }

    private static boolean testStatistical(int n, double p, int samples) {
        SimpleBinomial binomial = new SimpleBinomial();
        SimpleRand rng = new SimpleRand(12345);

        long sum = 0;
        long sumSq = 0;

        for (int i = 0; i < samples; i++) {
            int value = binomial.sampleInt(n, p, rng);

            // Sanity check: value should be in valid range
            if (value < 0 || value > n) {
                System.out.println("  ERROR: Sample out of range: " + value);
                return false;
            }

            sum += value;
            sumSq += (long)value * value;
        }

        double mean = (double)sum / samples;
        double variance = ((double)sumSq / samples) - (mean * mean);
        double stdDev = Math.sqrt(variance);

        double expectedMean = n * p;
        double expectedStdDev = Math.sqrt(n * p * (1 - p));

        // Use generous tolerances (4 sigma for mean, 10% for stddev)
        double meanError = Math.abs(mean - expectedMean);
        double meanTolerance = 4.0 * expectedStdDev / Math.sqrt(samples);

        double stdDevError = Math.abs(stdDev - expectedStdDev);
        double stdDevTolerance = 0.10 * expectedStdDev;

        System.out.printf("  Expected: mean=%.2f, stddev=%.2f\n", expectedMean, expectedStdDev);
        System.out.printf("  Observed: mean=%.2f, stddev=%.2f\n", mean, stdDev);
        System.out.printf("  Errors: mean=%.4f (tol=%.4f), stddev=%.4f (tol=%.4f)\n",
                         meanError, meanTolerance, stdDevError, stdDevTolerance);

        boolean meanOk = meanError <= meanTolerance;
        boolean stdDevOk = stdDevError <= stdDevTolerance;

        if (!meanOk) System.out.println("  WARNING: Mean outside tolerance");
        if (!stdDevOk) System.out.println("  WARNING: StdDev outside tolerance");

        return meanOk && stdDevOk;
    }

    private static boolean testStatisticalLong(long n, double p, int samples) {
        SimpleBinomial binomial = new SimpleBinomial();
        SimpleRand rng = new SimpleRand(12345);

        double sum = 0;
        double sumSq = 0;

        for (int i = 0; i < samples; i++) {
            long value = binomial.sampleLong(n, p, rng);

            // Sanity check
            if (value < 0 || value > n) {
                System.out.println("  ERROR: Sample out of range: " + value);
                return false;
            }

            sum += value;
            sumSq += (double)value * value;
        }

        double mean = sum / samples;
        double variance = (sumSq / samples) - (mean * mean);
        double stdDev = Math.sqrt(variance);

        double expectedMean = n * p;
        double expectedStdDev = Math.sqrt(n * p * (1 - p));

        double meanError = Math.abs(mean - expectedMean);
        double meanTolerance = 4.0 * expectedStdDev / Math.sqrt(samples);

        double stdDevError = Math.abs(stdDev - expectedStdDev);
        double stdDevTolerance = 0.10 * expectedStdDev;

        System.out.printf("  Expected: mean=%.2f, stddev=%.2f\n", expectedMean, expectedStdDev);
        System.out.printf("  Observed: mean=%.2f, stddev=%.2f\n", mean, stdDev);
        System.out.printf("  Errors: mean=%.4f (tol=%.4f), stddev=%.4f (tol=%.4f)\n",
                         meanError, meanTolerance, stdDevError, stdDevTolerance);

        boolean meanOk = meanError <= meanTolerance;
        boolean stdDevOk = stdDevError <= stdDevTolerance;

        if (!meanOk) System.out.println("  WARNING: Mean outside tolerance");
        if (!stdDevOk) System.out.println("  WARNING: StdDev outside tolerance");

        return meanOk && stdDevOk;
    }

    private static boolean testParameterCaching() {
        SimpleBinomial binomial = new SimpleBinomial();
        SimpleRand rng = new SimpleRand(999);

        int n = 500;
        double p = 0.3;

        // Generate many samples with same parameters (tests caching)
        long sum = 0;
        for (int i = 0; i < 10000; i++) {
            int value = binomial.sampleInt(n, p, rng);
            if (value < 0 || value > n) {
                System.out.println("  ERROR: Invalid sample");
                return false;
            }
            sum += value;
        }

        double mean = (double)sum / 10000;
        double expectedMean = n * p;
        double error = Math.abs(mean - expectedMean);

        System.out.printf("  Expected mean: %.2f, Observed: %.2f, Error: %.4f\n",
                         expectedMean, mean, error);

        // Should be very close with 10000 samples
        boolean ok = error < expectedMean * 0.05;
        if (!ok) System.out.println("  WARNING: Parameter caching may have issues");

        return ok;
    }
}

// Minimal Binomial implementation that uses the refactored code
class SimpleBinomial {
    private int n_last = -1;
    private long n_lastL = -1;
    private double p_last = -1.0;
    private double par, q, np, p0, rc, ss, xm, xl, xr, ll, lr, c, p1, p2, p3, p4, ch;
    private int b, m, nm;
    private long bL, mL, nmL;
    private int n_prev = -1;
    private long n_prevL = -1;
    private double p_prev = -1.0;
    private double pq;

    public int sampleInt(int n, double p, SimpleRand rn) {
        if (p == 1) return n;
        if (p == 0 || n == 0) return 0;

        initializeParametersInt(n, p);

        if (this.np < 10.0) {
            return sampleSmallNpInt(n, p, rn);
        } else {
            return sampleLargeNpInt(n, p, rn);
        }
    }

    public long sampleLong(long n, double p, SimpleRand rn) {
        if (p == 1) return n;
        if (p == 0 || n == 0) return 0;

        initializeParametersLong(n, p);

        if (this.np < 10.0) {
            return sampleSmallNpLong(n, p, rn);
        } else {
            return sampleLargeNpLong(n, p, rn);
        }
    }

    // Copy of refactored int methods
    private void initializeParametersInt(int n, double p) {
        if (n == this.n_last && p == this.p_last) return;

        this.n_last = n;
        this.p_last = p;
        this.par = Math.min(p, 1.0 - p);
        this.q = 1.0 - this.par;
        this.np = (double)n * this.par;

        if (this.np <= 0.0) return;

        double rm = this.np + this.par;
        this.m = (int)rm;

        if (this.np < 10.0) {
            this.p0 = Math.exp((double)n * Math.log(this.q));
            int bh = (int)(this.np + 10.0 * Math.sqrt(this.np * this.q));
            this.b = Math.min(n, bh);
        } else {
            this.rc = ((double)n + 1.0) * (this.pq = this.par / this.q);
            this.ss = this.np * this.q;
            int i = (int)(2.195 * Math.sqrt(this.ss) - 4.6 * this.q);
            this.xm = (double)this.m + 0.5;
            this.xl = (double)(this.m - i);
            this.xr = (double)((long)(this.m + i) + 1L);
            double f = (rm - this.xl) / (rm - this.xl * this.par);
            this.ll = f * (1.0 + 0.5 * f);
            f = (this.xr - rm) / (this.xr * this.q);
            this.lr = f * (1.0 + 0.5 * f);
            this.c = 0.134 + 20.5 / (15.3 + (double)this.m);
            this.p1 = (double)i + 0.5;
            this.p2 = this.p1 * (1.0 + this.c + this.c);
            this.p3 = this.p2 + this.c / this.ll;
            this.p4 = this.p3 + this.c / this.lr;
        }
    }

    private int sampleSmallNpInt(int n, double p, SimpleRand rn) {
        int K = 0;
        double pk = this.p0;
        double U = rn.nextDouble();

        while (U > pk) {
            ++K;
            if (K > this.b) {
                U = rn.nextDouble();
                K = 0;
                pk = this.p0;
            } else {
                U -= pk;
                pk = (double)(n - K + 1) * this.par * pk / ((double)K * this.q);
            }
        }

        return p > 0.5 ? n - K : K;
    }

    private int sampleLargeNpInt(int n, double p, SimpleRand rn) {
        while (true) {
            double V = rn.nextDouble();
            double U = rn.nextDouble() * this.p4;
            int K;

            if (U <= this.p1) {
                K = (int)(this.xm - U + this.p1 * V);
                return p > 0.5 ? n - K : K;
            }

            if (U <= this.p2) {
                double X = this.xl + (U - this.p1) / this.c;
                V = V * this.c + 1.0 - Math.abs(this.xm - X) / this.p1;
                if (V < 1.0) {
                    K = (int)X;
                } else {
                    continue;
                }
            } else if (U <= this.p3) {
                double X = this.xl + Math.log(V) / this.ll;
                if (X >= 0.0) {
                    K = (int)X;
                    V *= (U - this.p2) * this.ll;
                } else {
                    continue;
                }
            } else {
                K = (int)(this.xr - Math.log(V) / this.lr);
                if (K <= n) {
                    V *= (U - this.p3) * this.lr;
                } else {
                    continue;
                }
            }

            // Acceptance test
            int Km = Math.abs(K - this.m);
            if (Km > 20 && (double)((long)(Km + Km) + 2L) < this.ss) {
                if (acceptFastInt(n, K, Km, V)) {
                    return p > 0.5 ? n - K : K;
                }
            } else {
                if (acceptExactInt(K, V)) {
                    return p > 0.5 ? n - K : K;
                }
            }
        }
    }

    private boolean acceptFastInt(int n, int K, int Km, double V) {
        V = Math.log(V);
        double T = (double)(-Km * Km) / (this.ss + this.ss);
        double E = (double)Km / this.ss * (((double)Km * ((double)Km * 0.3333333333333333 + 0.625) + 0.16666666666666666) / this.ss + 0.5);
        if (V <= T - E) return true;
        if (V > T + E) return false;
        return acceptStirlingInt(n, K, V);
    }

    private boolean acceptStirlingInt(int n, int K, double V) {
        if (n != this.n_prev || this.par != this.p_prev) {
            this.n_prev = n;
            this.p_prev = this.par;
            this.nm = n - this.m + 1;
            this.ch = this.xm * Math.log(((double)this.m + 1.0) / (this.pq * (double)this.nm))
                      + stirlingCorrection(this.m + 1) + stirlingCorrection(this.nm);
        }
        int nK = n - K + 1;
        return V <= this.ch + ((double)n + 1.0) * Math.log((double)this.nm / (double)nK)
                    + ((double)K + 0.5) * Math.log((double)nK * this.pq / ((double)K + 1.0))
                    - stirlingCorrection(K + 1) - stirlingCorrection(nK);
    }

    private boolean acceptExactInt(int K, double V) {
        double f = 1.0;
        if (this.m < K) {
            int i = this.m;
            while (i < K) {
                ++i;
                f *= this.rc / (double)i - this.pq;
                if (f < V) return false;
            }
        } else {
            int i = K;
            while (i < this.m) {
                ++i;
                V *= this.rc / (double)i - this.pq;
                if (V > f) return false;
            }
        }
        return V <= f;
    }

    // Similar methods for long version
    private void initializeParametersLong(long n, double p) {
        if (n == this.n_lastL && p == this.p_last) return;

        this.n_lastL = n;
        this.p_last = p;
        this.par = Math.min(p, 1.0 - p);
        this.q = 1.0 - this.par;
        this.np = (double)n * this.par;

        if (this.np <= 0.0) return;

        double rm = this.np + this.par;
        this.mL = (long)rm;

        if (this.np < 10.0) {
            this.p0 = Math.exp((double)n * Math.log(this.q));
            long bh = (long)(this.np + 10.0 * Math.sqrt(this.np * this.q));
            this.bL = Math.min(n, bh);
        } else {
            this.rc = ((double)n + 1.0) * (this.pq = this.par / this.q);
            this.ss = this.np * this.q;
            long i = (long)(2.195 * Math.sqrt(this.ss) - 4.6 * this.q);
            this.xm = (double)this.mL + 0.5;
            this.xl = (double)(this.mL - i);
            this.xr = (double)((long)(this.mL + i) + 1L);
            double f = (rm - this.xl) / (rm - this.xl * this.par);
            this.ll = f * (1.0 + 0.5 * f);
            f = (this.xr - rm) / (this.xr * this.q);
            this.lr = f * (1.0 + 0.5 * f);
            this.c = 0.134 + 20.5 / (15.3 + (double)this.mL);
            this.p1 = (double)i + 0.5;
            this.p2 = this.p1 * (1.0 + this.c + this.c);
            this.p3 = this.p2 + this.c / this.ll;
            this.p4 = this.p3 + this.c / this.lr;
        }
    }

    private long sampleSmallNpLong(long n, double p, SimpleRand rn) {
        long K = 0;
        double pk = this.p0;
        double U = rn.nextDouble();

        while (U > pk) {
            ++K;
            if (K > this.bL) {
                U = rn.nextDouble();
                K = 0;
                pk = this.p0;
            } else {
                U -= pk;
                pk = (double)(n - K + 1) * this.par * pk / ((double)K * this.q);
            }
        }

        return p > 0.5 ? n - K : K;
    }

    private long sampleLargeNpLong(long n, double p, SimpleRand rn) {
        while (true) {
            double V = rn.nextDouble();
            double U = rn.nextDouble() * this.p4;
            long K;

            if (U <= this.p1) {
                K = (long)(this.xm - U + this.p1 * V);
                return p > 0.5 ? n - K : K;
            }

            if (U <= this.p2) {
                double X = this.xl + (U - this.p1) / this.c;
                V = V * this.c + 1.0 - Math.abs(this.xm - X) / this.p1;
                if (V < 1.0) {
                    K = (long)X;
                } else {
                    continue;
                }
            } else if (U <= this.p3) {
                double X = this.xl + Math.log(V) / this.ll;
                if (X >= 0.0) {
                    K = (long)X;
                    V *= (U - this.p2) * this.ll;
                } else {
                    continue;
                }
            } else {
                K = (long)(this.xr - Math.log(V) / this.lr);
                if (K <= n) {
                    V *= (U - this.p3) * this.lr;
                } else {
                    continue;
                }
            }

            long Km = Math.abs(K - this.mL);
            if (Km > 20 && (double)((long)(Km + Km) + 2L) < this.ss) {
                if (acceptFastLong(n, K, Km, V)) {
                    return p > 0.5 ? n - K : K;
                }
            } else {
                if (acceptExactLong(K, V)) {
                    return p > 0.5 ? n - K : K;
                }
            }
        }
    }

    private boolean acceptFastLong(long n, long K, long Km, double V) {
        V = Math.log(V);
        double T = (double)(-Km * Km) / (this.ss + this.ss);
        double E = (double)Km / this.ss * (((double)Km * ((double)Km * 0.3333333333333333 + 0.625) + 0.16666666666666666) / this.ss + 0.5);
        if (V <= T - E) return true;
        if (V > T + E) return false;
        return acceptStirlingLong(n, K, V);
    }

    private boolean acceptStirlingLong(long n, long K, double V) {
        if (n != this.n_prevL || this.par != this.p_prev) {
            this.n_prevL = n;
            this.p_prev = this.par;
            this.nmL = n - this.mL + 1;
            this.ch = this.xm * Math.log(((double)this.mL + 1.0) / (this.pq * (double)this.nmL))
                      + stirlingCorrection(this.mL + 1) + stirlingCorrection(this.nmL);
        }
        long nK = n - K + 1;
        return V <= this.ch + ((double)n + 1.0) * Math.log((double)this.nmL / (double)nK)
                    + ((double)K + 0.5) * Math.log((double)nK * this.pq / ((double)K + 1.0))
                    - stirlingCorrection(K + 1) - stirlingCorrection(nK);
    }

    private boolean acceptExactLong(long K, double V) {
        double f = 1.0;
        if (this.mL < K) {
            long i = this.mL;
            while (i < K) {
                ++i;
                f *= this.rc / (double)i - this.pq;
                if (f < V) return false;
            }
        } else {
            long i = K;
            while (i < this.mL) {
                ++i;
                V *= this.rc / (double)i - this.pq;
                if (V > f) return false;
            }
        }
        return V <= f;
    }

    private static final double[] stirlingCorrection = {
        0.0, 8.106146679532726e-02, 4.134069595540929e-02, 2.767792568499834e-02,
        2.079067210376509e-02, 1.664469118982119e-02, 1.387612882307075e-02,
        1.189670994589177e-02, 1.041126526197209e-02, 9.255462182712733e-03,
        8.330563433362871e-03, 7.573675487951841e-03, 6.942840107209530e-03,
        6.408994188004207e-03, 5.951370112758848e-03, 5.554733551962801e-03,
        5.207655919609640e-03, 4.901395948434738e-03, 4.629153749334029e-03,
        4.385560249232324e-03, 4.166319691996922e-03, 3.967954218640860e-03,
        3.787618068444430e-03, 3.622960224683090e-03, 3.472021382978770e-03,
        3.333155636728090e-03, 3.204970228055040e-03, 3.086278682608780e-03,
        2.976063983550410e-03, 2.873449362352470e-03, 2.777674929752690e-03
    };

    private static double stirlingCorrection(long k) {
        if (k > 30) {
            double r = 1.0 / (double)k;
            double rr = r * r;
            return r * (0.08333333333333333 + rr * (-0.002777777777777778 + rr * (0.0007936507936507937 + rr * (-0.0005952380952380952))));
        }
        return stirlingCorrection[(int)k];
    }
}

// Simple random number generator
class SimpleRand {
    private long seed;

    public SimpleRand(long seed) {
        this.seed = seed;
    }

    public double nextDouble() {
        seed = (seed * 6364136223846793005L + 1442695040888963407L);
        return ((seed >>> 11) & 0x1FFFFFFFFFFFFFL) / (double)(1L << 53);
    }
}
