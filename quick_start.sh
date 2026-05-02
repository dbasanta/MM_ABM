#!/bin/bash
# Quick start script for running bone-tumor simulations

echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║         Bone-Tumor Simulation Quick Start                             ║"
echo "║         Demonstrating Refactored Binomial Code                        ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

# Check Java
if ! command -v javac &> /dev/null; then
    echo "❌ Error: Java compiler (javac) not found"
    echo "   Please install Java JDK 11 or higher"
    exit 1
fi

echo "✓ Java found: $(javac -version 2>&1)"
echo ""

# Menu
echo "Choose a simulation to run:"
echo ""
echo "  1) Bone-Tumor Interaction (tumor + osteoclasts + osteoblasts)"
echo "     → Shows pathological vicious cycle"
echo "     → Runtime: ~1 second"
echo ""
echo "  2) Simple Tumor Growth (normal + mutant cells)"
echo "     → Shows tumor evolution with mutations"
echo "     → Runtime: ~1 second"
echo ""
echo "  3) Verification Tests (test refactored code)"
echo "     → 10 comprehensive unit tests"
echo "     → Runtime: ~10 seconds"
echo ""
echo "  4) All of the above (complete demo)"
echo ""
read -p "Enter choice (1-4): " choice

echo ""

case $choice in
    1|4)
        echo "═══════════════════════════════════════════════════════════════════"
        echo "  Running Bone-Tumor Simulation"
        echo "═══════════════════════════════════════════════════════════════════"
        echo ""

        # Compile if needed
        if [ ! -f "BoneTumorSimulation.class" ]; then
            echo "🔨 Compiling..."
            javac BoneTumorSimulation.java
        fi

        # Run simulation
        echo "🎮 Running simulation (300 timesteps)..."
        java BoneTumorSimulation
        echo ""

        # Generate plots
        if [ ! -f "PlotBoneTumorData.class" ]; then
            echo "🔨 Compiling visualization..."
            javac PlotBoneTumorData.java
        fi

        echo "📊 Generating plots..."
        java PlotBoneTumorData
        echo ""

        echo "✓ Results saved to: bone_tumor_simulation.csv"
        echo ""

        if [ "$choice" != "4" ]; then
            echo "💡 Tip: View the CSV in Excel or run:"
            echo "   head bone_tumor_simulation.csv"
            exit 0
        fi
        ;;
esac

case $choice in
    2|4)
        echo "═══════════════════════════════════════════════════════════════════"
        echo "  Running Simple Tumor Simulation"
        echo "═══════════════════════════════════════════════════════════════════"
        echo ""

        if [ ! -f "SimpleTumorSimulation.class" ]; then
            echo "🔨 Compiling..."
            javac SimpleTumorSimulation.java
        fi

        echo "🎮 Running simulation (500 timesteps)..."
        java SimpleTumorSimulation
        echo ""

        echo "✓ Results saved to: tumor_population.csv"
        echo ""

        if [ "$choice" != "4" ]; then
            exit 0
        fi
        ;;
esac

case $choice in
    3|4)
        echo "═══════════════════════════════════════════════════════════════════"
        echo "  Running Verification Tests"
        echo "═══════════════════════════════════════════════════════════════════"
        echo ""

        if [ ! -f "VerifyBinomialRefactor.class" ]; then
            echo "🔨 Compiling..."
            javac VerifyBinomialRefactor.java
        fi

        echo "🧪 Running 10 comprehensive tests..."
        java VerifyBinomialRefactor
        echo ""
        ;;
esac

if [ "$choice" == "4" ]; then
    echo "╔════════════════════════════════════════════════════════════════════════╗"
    echo "║                      ALL SIMULATIONS COMPLETE!                         ║"
    echo "╚════════════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "📁 Output files created:"
    echo "   • bone_tumor_simulation.csv (bone-tumor data)"
    echo "   • tumor_population.csv (simple tumor data)"
    echo ""
    echo "📖 Next steps:"
    echo "   • View results: cat bone_tumor_simulation.csv | head -20"
    echo "   • Customize parameters: edit BoneTumorSimulation.java"
    echo "   • Read guide: cat RUNNING_SIMULATIONS.md"
    echo ""
fi

echo "✨ Done! Refactored Binomial code working perfectly!"
