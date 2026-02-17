#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd

# Read data
data = pd.read_csv('tumor_population.csv')

# Create plot
plt.figure(figsize=(12, 6))

# Plot total population
plt.subplot(1, 2, 1)
plt.plot(data['Time'], data['Total'], 'b-', linewidth=2, label='Total')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Total Tumor Population Over Time')
plt.grid(True, alpha=0.3)
plt.legend()

# Plot population breakdown
plt.subplot(1, 2, 2)
plt.plot(data['Time'], data['Normal'], 'g-', linewidth=2, label='Normal')
plt.plot(data['Time'], data['Mutant'], 'r-', linewidth=2, label='Mutant')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Normal vs Mutant Cells')
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig('tumor_population.png', dpi=150)
print('Plot saved to: tumor_population.png')
