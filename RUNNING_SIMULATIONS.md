# Running Simulations Locally

Complete guide to running bone-tumor simulations on your local machine.

---

## 🎯 Quick Start (30 seconds)

```bash
cd /home/user/MM_ABM

# Run the bone-tumor simulation
javac BoneTumorSimulation.java && java BoneTumorSimulation

# View the plots
javac PlotBoneTumorData.java && java PlotBoneTumorData
```

Done! You'll see ASCII plots and a CSV file with results.

---

## 📋 Prerequisites

### Required
- **Java JDK** 11 or higher (you have: Java 21 ✓)
- **Terminal/Command line** access

### Optional (for better plots)
- **Python 3** with matplotlib (for graphical plots)
- **gnuplot** (alternative plotting tool)

---

## 🧪 Available Simulations

### 1️⃣ **BoneTumorSimulation** (Recommended)

Models tumor-bone microenvironment with three cell types.

**What it simulates:**
- Tumor cells (Multiple Myeloma)
- Osteoclasts (bone-resorbing)
- Osteoblasts (bone-forming)
- "Vicious cycle" pathology

**How to run:**
```bash
# Basic run (300 timesteps, default parameters)
java BoneTumorSimulation

# View plots
java PlotBoneTumorData
```

**Customize parameters:**
Edit `BoneTumorSimulation.java` lines 13-24:
```java
int timeSteps = 300;              // Duration
double tumor_birthRate = 0.15;    // Tumor proliferation
double tumor_deathRate = 0.05;    // Tumor death
double osteoclast_recruitRate = 0.10;  // OC recruitment
double osteoblast_recruitRate = 0.08;  // OB recruitment
// ... etc
```

**Output files:**
- `bone_tumor_simulation.csv` - Raw data
- Console: ASCII plots of all populations

**Runtime:** <1 second

---

### 2️⃣ **SimpleTumorSimulation**

Basic tumor growth with normal/mutant cell competition.

**What it simulates:**
- Normal tumor cells
- Mutant tumor cells (with growth advantage)
- Stochastic birth/death/mutation events

**How to run:**
```bash
java SimpleTumorSimulation
```

**Customize parameters:**
Edit lines 15-19:
```java
int timeSteps = 500;
double birthRate = 0.12;
double deathRate = 0.08;
double mutationRate = 0.001;
int initialPopulation = 100;
```

**Output files:**
- `tumor_population.csv`
- `plot_tumor.py` (Python plotting script)

**Runtime:** <1 second

---

### 3️⃣ **Verification Tests**

Test that the refactored Binomial code works correctly.

**How to run:**
```bash
# Comprehensive unit tests (10 tests)
java VerifyBinomialRefactor

# Expected output: "ALL TESTS PASSED!"
```

**What it tests:**
- Edge cases (p=0, p=1)
- Small np (inversion method)
- Large np (BTPE method)
- Both int and long versions
- Statistical accuracy

**Runtime:** ~10 seconds

---

## 🎨 Generating Better Plots

### Option A: Python with matplotlib

**If you have Python + matplotlib:**
```bash
# Run simulation first
java BoneTumorSimulation

# Generate plot with provided script
python3 plot_tumor.py
# Creates: tumor_population.png
```

**Install matplotlib if needed:**
```bash
pip3 install matplotlib pandas
```

### Option B: Custom Python Script

Create `plot_bone_tumor.py`:
```python
#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd

# Read data
data = pd.read_csv('bone_tumor_simulation.csv')

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: All populations
axes[0,0].plot(data['Time'], data['Tumor'], 'r-', label='Tumor', linewidth=2)
axes[0,0].plot(data['Time'], data['Osteoclasts'], 'orange', label='Osteoclasts', linewidth=2)
axes[0,0].plot(data['Time'], data['Osteoblasts'], 'b-', label='Osteoblasts', linewidth=2)
axes[0,0].set_xlabel('Time')
axes[0,0].set_ylabel('Population')
axes[0,0].set_title('All Cell Populations')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

# Plot 2: Tumor only
axes[0,1].plot(data['Time'], data['Tumor'], 'r-', linewidth=2)
axes[0,1].set_xlabel('Time')
axes[0,1].set_ylabel('Tumor Cells')
axes[0,1].set_title('Tumor Growth')
axes[0,1].grid(True, alpha=0.3)

# Plot 3: OC vs OB
axes[1,0].plot(data['Time'], data['Osteoclasts'], 'orange', label='Osteoclasts', linewidth=2)
axes[1,0].plot(data['Time'], data['Osteoblasts'], 'b-', label='Osteoblasts', linewidth=2)
axes[1,0].set_xlabel('Time')
axes[1,0].set_ylabel('Population')
axes[1,0].set_title('Bone Cell Balance')
axes[1,0].legend()
axes[1,0].grid(True, alpha=0.3)

# Plot 4: OC/OB Ratio
axes[1,1].plot(data['Time'], data['OC_OB_Ratio'], 'purple', linewidth=2)
axes[1,1].axhline(y=1.0, color='g', linestyle='--', label='Normal')
axes[1,1].axhline(y=2.0, color='r', linestyle='--', label='Pathological')
axes[1,1].set_xlabel('Time')
axes[1,1].set_ylabel('OC/OB Ratio')
axes[1,1].set_title('Bone Homeostasis')
axes[1,1].legend()
axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('bone_tumor_analysis.png', dpi=300)
print('✓ Plot saved: bone_tumor_analysis.png')
```

Run it:
```bash
python3 plot_bone_tumor.py
```

---

## 🔧 Customizing Simulations

### Change Simulation Duration

```java
// Short test (quick preview)
int timeSteps = 50;

// Default (good balance)
int timeSteps = 300;

// Long run (see full dynamics)
int timeSteps = 1000;
```

### Adjust Cell Parameters

**Make tumor more aggressive:**
```java
double tumor_birthRate = 0.20;     // Increased from 0.15
double tumor_deathRate = 0.03;     // Decreased from 0.05
```

**Increase osteoclast activity:**
```java
double osteoclast_recruitRate = 0.15;  // Increased from 0.10
```

**Simulate treatment effect:**
```java
// Add this in the simulation loop after time = 100
if (t > 100) {
    // Anti-RANKL therapy (reduce OC recruitment)
    effectiveOC_recruit *= 0.5;  // 50% reduction
}
```

### Change Initial Conditions

```java
// Start with more tumor burden
int tumorCells = 200;  // Instead of 50

// Healthier bone initially
int osteoblasts = 500;  // Instead of 200
```

### Add Random Variation

```java
// Use different random seed for different results
SimpleRand rng = new SimpleRand(System.currentTimeMillis());
// Or use a specific seed for reproducibility
SimpleRand rng = new SimpleRand(42);
```

---

## 📊 Understanding the Output

### CSV File Format

```csv
Time,Tumor,Osteoclasts,Osteoblasts,OC_OB_Ratio
0,50,100,200,0.5000
1,52,100,192,0.5208
2,54,102,188,0.5426
...
```

**Columns:**
- `Time`: Timestep number
- `Tumor`: Tumor cell count
- `Osteoclasts`: Osteoclast count
- `Osteoblasts`: Osteoblast count
- `OC_OB_Ratio`: Ratio (indicates bone balance)

### Import into Excel/Sheets

1. Open Excel/Google Sheets
2. File → Import → CSV
3. Select `bone_tumor_simulation.csv`
4. Create charts from the data

### Analyze with R

```r
# Read data
data <- read.csv("bone_tumor_simulation.csv")

# Basic plots
plot(data$Time, data$Tumor, type='l', col='red', 
     main='Tumor Growth', xlab='Time', ylab='Cells')

# Statistical summary
summary(data)
```

---

## 🔬 Advanced: Running Parameter Sweeps

Create `parameter_sweep.sh`:

```bash
#!/bin/bash
# Test multiple tumor birth rates

for rate in 0.10 0.12 0.15 0.18 0.20
do
    echo "Testing birth rate: $rate"
    
    # Modify the Java file (you'd do this programmatically)
    # Or create multiple versions with different parameters
    
    # Run simulation
    java BoneTumorSimulation > "output_rate_${rate}.txt"
    
    # Move CSV to unique name
    mv bone_tumor_simulation.csv "bone_tumor_rate_${rate}.csv"
done

echo "All parameter sweep runs complete!"
```

---

## 🐛 Troubleshooting

### "Command not found: javac"

**Solution:** Install Java JDK
```bash
# Ubuntu/Debian
sudo apt-get install openjdk-21-jdk

# macOS
brew install openjdk@21

# Check installation
javac -version
```

### "Exception in thread ... NoClassDefFoundError"

**Solution:** Make sure you're in the correct directory
```bash
pwd  # Should show .../MM_ABM
ls   # Should show BoneTumorSimulation.java
```

### Compilation errors

**Solution:** Recompile dependencies
```bash
# Compile in correct order
javac SimpleBinomial.java SimpleRand.java
javac BoneTumorSimulation.java
```

### "Cannot find symbol" errors

**Solution:** Check that all required files exist
```bash
ls -la *.java | grep -E "(Simple|Bone|Verify)"
```

Should show:
- SimpleBinomial.java (in VerifyBinomialRefactor.java)
- SimpleRand.java (in VerifyBinomialRefactor.java)
- BoneTumorSimulation.java
- PlotBoneTumorData.java

---

## 💾 Saving Your Results

### Create an experiment log

```bash
# Create experiment directory
mkdir -p experiments/exp_001

# Run simulation
java BoneTumorSimulation | tee experiments/exp_001/output.log

# Save data
cp bone_tumor_simulation.csv experiments/exp_001/
cp BoneTumorSimulation.java experiments/exp_001/parameters.java

# Add notes
echo "Experiment 1: Baseline tumor-bone dynamics" > experiments/exp_001/README.txt
echo "Date: $(date)" >> experiments/exp_001/README.txt
echo "Parameters: Default" >> experiments/exp_001/README.txt
```

---

## 🚀 Next Steps

**Beginner:**
1. Run `BoneTumorSimulation` with default parameters
2. View the ASCII plots
3. Open CSV in Excel to see the data

**Intermediate:**
4. Modify one parameter (e.g., `tumor_birthRate`)
5. Run again and compare results
6. Create Python plots for better visualization

**Advanced:**
7. Run parameter sweeps
8. Implement treatment effects
9. Add new cell types or interactions
10. Compare to published bone metastasis data

---

## 📚 Documentation

- **Code documentation:** See comments in Java files
- **Algorithm details:** REFACTORING_SUMMARY.md
- **Verification results:** VERIFICATION_COMPLETE.md
- **Bone-tumor biology:** BONE_TUMOR_RESULTS.md

---

## ✅ Quick Checklist

Before running:
- [ ] In `/home/user/MM_ABM` directory?
- [ ] Java installed? (`javac -version`)
- [ ] Files compiled? (`ls *.class`)

After running:
- [ ] CSV file created?
- [ ] Plots displayed correctly?
- [ ] Results make biological sense?

---

**Happy simulating! 🎉**

For questions or issues, see the documentation files or check the code comments.
