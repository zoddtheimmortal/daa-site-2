import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter

# Define the algorithm names and datasets
algorithms = ["Chiba", "Tomita", "Bron"]
datasets = ["wiki", "email", "skitter"]

# Execution times in milliseconds for each algorithm on each dataset
execution_times = {
    "Chiba": {
        "wiki": 419446,
        "email": 558357,
        "skitter": 1.5368e10  # Very large number
    },
    "Tomita": {
        "wiki": 1215,
        "email": 1727,
        "skitter": 1820871
    },
    "Bron": {
        "wiki": 7546,
        "email": 6208,
        "skitter": 3810505
    }
}

# Graph sizes for each dataset (number of nodes)
dataset_sizes = {
    "wiki": 7115,
    "email": 36692,
    "skitter": 1696415
}

# Create a figure for bar chart comparison
fig, ax = plt.subplots(figsize=(14, 8))

# Set up the bar positions
bar_width = 0.25
r1 = np.arange(len(datasets))
r2 = [x + bar_width for x in r1]
r3 = [x + bar_width for x in r2]

# Create the grouped bar chart
for i, algorithm in enumerate(algorithms):
    data = [execution_times[algorithm][dataset] for dataset in datasets]
    positions = [r1, r2, r3][i]
    ax.bar(positions, data, width=bar_width, label=algorithm)

# Add labels and title
ax.set_xlabel('Dataset', fontsize=14)
ax.set_ylabel('Execution Time (ms)', fontsize=14)
ax.set_title('Algorithm Performance on Different Datasets', fontsize=16)
ax.set_xticks([r + bar_width for r in range(len(datasets))])
ax.set_xticklabels([f"{d.capitalize()}\n({dataset_sizes[d]} nodes)" for d in datasets])

# Use logarithmic scale for better visualization
ax.set_yscale('log')

# Add data labels
for i, algorithm in enumerate(algorithms):
    positions = [r1, r2, r3][i]
    data = [execution_times[algorithm][dataset] for dataset in datasets]
    for x, y in zip(positions, data):
        # Format large numbers appropriately
        if y > 1000000:
            label = f"{y/1000000:.1f}M"
        elif y > 1000:
            label = f"{y/1000:.1f}K"
        else:
            label = f"{y:.0f}"
            
        ax.text(x, y*1.1, label, ha='center', va='bottom', fontsize=9)

# Add legend
ax.legend(fontsize=12)

# Add grid for better readability
ax.grid(True, axis='y', linestyle='--', alpha=0.7)

# Adjust layout
plt.tight_layout()

# Create a second visualization: horizontal bar chart for better comparison
fig2, ax2 = plt.subplots(figsize=(14, 8))

# Prepare data for horizontal bars
y_pos = np.arange(len(algorithms))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

# Setup width for dataset groups
dataset_spacing = 0.3
total_width = 0.8

# For each dataset, create a group of horizontal bars
for i, dataset in enumerate(datasets):
    data = [execution_times[algo][dataset] for algo in algorithms]
    positions = y_pos - total_width/2 + (i+0.5)*total_width/len(datasets)
    
    ax2.barh(positions, data, height=total_width/len(datasets), 
            label=f"{dataset.capitalize()} ({dataset_sizes[dataset]} nodes)",
            color=colors[i])

# Add labels and formatting
ax2.set_yticks(y_pos)
ax2.set_yticklabels(algorithms, fontsize=12)
ax2.set_xlabel('Execution Time (ms)', fontsize=14)
ax2.set_title('Algorithm Performance Comparison', fontsize=16)

# Use logarithmic scale for better visualization
ax2.set_xscale('log')

# Add grid for better readability
ax2.grid(True, axis='x', linestyle='--', alpha=0.7)

# Add legend
ax2.legend(fontsize=10, loc='upper right')

# Tight layout
plt.tight_layout()

# Create a third visualization: line graph showing algorithm performance vs graph size
fig3, ax3 = plt.subplots(figsize=(14, 8))

# Sort datasets by size for proper line graph
sorted_datasets = sorted(datasets, key=lambda d: dataset_sizes[d])
sorted_sizes = [dataset_sizes[d] for d in sorted_datasets]

# Line styles and markers for different algorithms
line_styles = ['-o', '--s', '-.^']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

# Plot line for each algorithm
for i, algorithm in enumerate(algorithms):
    # Get execution times in the same order as sorted datasets
    times = [execution_times[algorithm][d] for d in sorted_datasets]
    
    # Plot the line
    ax3.plot(sorted_sizes, times, line_styles[i], linewidth=2, 
             label=algorithm, color=colors[i], markersize=8)
    
    # Add data point labels
    for x, y in zip(sorted_sizes, times):
        if y > 1000000:
            label = f"{y/1000000:.1f}M"
        elif y > 1000:
            label = f"{y/1000:.1f}K"
        else:
            label = f"{y:.0f}"
            
        ax3.annotate(label, (x, y), textcoords="offset points", 
                    xytext=(0,10), ha='center', fontsize=9)

# Set labels and title
ax3.set_xlabel('Graph Size (nodes)', fontsize=14)
ax3.set_ylabel('Execution Time (ms)', fontsize=14)
ax3.set_title('Algorithm Performance vs Graph Size', fontsize=16)

# Use logarithmic scales for both axes
ax3.set_xscale('log')
ax3.set_yscale('log')

# Add grid for better readability
ax3.grid(True, linestyle='--', alpha=0.7)

# Add legend
ax3.legend(fontsize=12, loc='upper left')

# Format axis tick labels
ax3.xaxis.set_major_formatter(ScalarFormatter())
ax3.yaxis.set_major_formatter(ScalarFormatter())

# Add dataset names at each data point on x-axis
for i, (size, dataset) in enumerate(zip(sorted_sizes, sorted_datasets)):
    ax3.annotate(dataset.capitalize(), (size, ax3.get_ylim()[0]), 
                textcoords="offset points", xytext=(0,-20), 
                ha='center', fontsize=10, rotation=0)

# Adjust layout
plt.tight_layout()

# Save the figures
fig.savefig('algorithm_performance_bar.png', dpi=300, bbox_inches='tight')
fig2.savefig('algorithm_performance_horizontal.png', dpi=300, bbox_inches='tight')
fig3.savefig('algorithm_performance_line.png', dpi=300, bbox_inches='tight')

# Show the figures
plt.show()

def save_data_to_csv():
    """Save the execution time data to a CSV file for reference"""
    data = []
    
    for algo in algorithms:
        for dataset in datasets:
            data.append({
                'Algorithm': algo,
                'Dataset': dataset,
                'Graph Size (nodes)': dataset_sizes[dataset],
                'Execution Time (ms)': execution_times[algo][dataset]
            })
    
    df = pd.DataFrame(data)
    df.to_csv('algorithm_performance_data.csv', index=False)
    print("Data saved to algorithm_performance_data.csv")

# Uncomment to save data to CSV
# save_data_to_csv()