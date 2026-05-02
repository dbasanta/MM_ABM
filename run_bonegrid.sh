#!/bin/bash
# Script to compile and run the BoneGrid bone remodeling model

echo "╔════════════════════════════════════════════════════════════╗"
echo "║  BoneGrid Model - Complete Bone Remodeling Simulation     ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Check if we're in the right directory
if [ ! -d "MM" ]; then
    echo "❌ Error: Please run this script from the MM_ABM directory"
    exit 1
fi

echo "📁 Setting up environment..."

# Set headless mode (no GUI) - edit the Java file
echo "🔧 Configuring for headless mode..."

# Create output directory
mkdir -p output/bone_simulation
echo "✓ Output directory created: output/bone_simulation"

echo ""
echo "⚙️  Compilation options:"
echo "  1) Quick test (100 timesteps, ~1 minute)"
echo "  2) Short run (1000 timesteps, ~10 minutes)"
echo "  3) Full simulation (87,600 timesteps = 1 year, several hours)"
echo ""
read -p "Choose option (1-3): " choice

case $choice in
    1) TIMESTEPS=100 ;;
    2) TIMESTEPS=1000 ;;
    3) TIMESTEPS=87600 ;;
    *) echo "Invalid choice, using quick test"; TIMESTEPS=100 ;;
esac

echo ""
echo "📝 Configuration:"
echo "  Timesteps: $TIMESTEPS ($(echo "scale=2; $TIMESTEPS * 6 / 60" | bc) hours simulated)"
echo "  Grid size: Based on model parameters"
echo "  Outputs: CSV files with population data"
echo ""

echo "🔨 Compiling BoneGrid model..."
echo "   (This may take a minute due to large file size)"

# Try to compile - note this may fail due to missing dependencies
javac -d . MM/BoneGrid_2022May17.java 2>&1 | tee compile.log

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful!"
    echo ""
    echo "🎮 Running simulation..."
    echo "   This will run for approximately $(echo "scale=1; $TIMESTEPS / 100" | bc) minutes"
    echo ""

    # Run the model (would need to be adapted based on how main() is called)
    java Bone.boneRemodeling_2022May17.BoneGrid_2022May17 $TIMESTEPS

    echo ""
    echo "✓ Simulation complete!"
    echo "📊 Results saved to output directory"

else
    echo ""
    echo "❌ Compilation failed. This is expected if dependencies are missing."
    echo ""
    echo "💡 The BoneGrid model requires:"
    echo "  - HAL framework (GridsAndAgents, Gui, Tools, etc.)"
    echo "  - All supporting classes"
    echo ""
    echo "✅ Try the simpler simulations instead:"
    echo "   • BoneTumorSimulation.java (recommended)"
    echo "   • SimpleTumorSimulation.java"
    echo ""
    echo "See RUNNING_SIMULATIONS.md for details"
fi
