import re
import matplotlib.pyplot as plt
import numpy as np

# Path to the data file
data_file = "./data.txt"

# Read the data file
with open(data_file, 'r') as f:
    data = f.read()

# Extract size and frequency using regex
pattern = r"Size (\d+): (\d+)"
matches = re.findall(pattern, data)

# Convert matches to lists of integers
sizes = [int(match[0]) for match in matches]
frequencies = [int(match[1]) for match in matches]

# Create the histogram with improved spacing
plt.figure(figsize=(16, 8))
bars = plt.bar(sizes, frequencies, width=0.7, edgecolor='black', alpha=0.85)
plt.xlabel('Clique Size', fontsize=12, labelpad=10)
plt.ylabel('Frequency (Number of Cliques)', fontsize=12, labelpad=10)
plt.title('Distribution of Maximal Cliques by Size', fontsize=16, pad=20)

# Improve x-axis spacing
plt.xticks(np.arange(min(sizes), max(sizes)+1, 2), rotation=45, ha='right')
plt.margins(x=0.01)

# Add grid for better readability
plt.grid(axis='y', linestyle='--', alpha=0.4)

# Add data labels only for significant values
threshold = max(frequencies) * 0.2
for bar, size, freq in zip(bars, sizes, frequencies):
    if freq > threshold:  # Only label significant bars to avoid clutter
        plt.text(bar.get_x() + bar.get_width()/2, freq + max(frequencies)*0.015, 
                str(freq), ha='center', va='bottom', fontsize=9)

# Save the figure with better spacing
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)  # Add space at the bottom for rotated labels
plt.savefig('./graph.png', dpi=300, bbox_inches='tight')
print("Improved histogram saved as docs/graphs/graph.png")