import matplotlib.pyplot as plt
import numpy as np

# --- Simulation Data ---
# Parameters: simtime = 1000, nx = 2000, hybrid mpi + omp (4 ranks/node, 8 threads/rank)

# Data structure: 
# [Node Configuration Name, T_init, T_compute, T_communicate]
data = [
    ['8 Nodes (32 Ranks)', 0.0462541431, 916.703613, 60.3360596],
    ['16 Nodes (64 Ranks)', 0.0747045502, 504.719849, 72.6496201],
    ['32 Nodes (128 Ranks)', 0.0443070084, 804.914673, 572.242981]
]

# Separate the data into components
node_labels = [item[0] for item in data]
T_init = np.array([item[1] for item in data])
T_compute = np.array([item[2] for item in data])
T_communicate = np.array([item[3] for item in data])

# Calculate the total time for labeling
T_total = T_init + T_compute + T_communicate

# --- Plotting Setup ---
fig, ax = plt.subplots(figsize=(8, 6))
bar_width = 0.5  # Width of the stacked bars
x = np.arange(len(node_labels)) # the label locations

# 1. Plot T_init (The base of the stack)
ax.bar(x, T_init, bar_width, label='T_init (Initialization)', color='#4c72b0')

# 2. Plot T_compute (Stacked on top of T_init)
# The bottom of this bar is T_init
rects_compute = ax.bar(x, T_compute, bar_width, bottom=T_init, label='T_compute (Main Calculation)', color='#dd8452')

# 3. Plot T_communicate (Stacked on top of T_init + T_compute)
# The bottom of this bar is T_init + T_compute
bottom_comm = T_init + T_compute
rects_comm = ax.bar(x, T_communicate, bar_width, bottom=bottom_comm, label='T_communicate (MPI/OMP Overhead)', color='#55a868')


# --- Customizing and Labeling ---

# Add labels for the total time above the stacks
for i in range(len(x)):
    total_height = T_total[i]
    ax.text(
        x[i], 
        total_height, 
        f'{total_height:.1f} s', 
        ha='center', 
        va='bottom', 
        fontsize=10, 
        fontweight='bold'
    )
    
# Add component labels (especially for T_compute/T_communicate)
def label_component(rects, base_offset, precision=1):
    for i, rect in enumerate(rects):
        height = rect.get_height()
        if height > 50: # Only label components that are large enough to be seen clearly
            mid_point = rect.get_y() + rect.get_height() / 2
            ax.text(
                rect.get_x() + rect.get_width() / 2, 
                mid_point, 
                f'{height:.{precision}f}', 
                ha='center', 
                va='center', 
                color='white',
                fontsize=8
            )
        elif height > 0:
            # Optionally label small components (like T_init) outside the bar
            pass

# Label the T_compute components
label_component(rects_compute, T_init)
# Label the T_communicate components
label_component(rects_comm, T_init + T_compute)


# Set titles and labels
ax.set_ylabel('Time (Seconds)', fontsize=12)
ax.set_title(
    'Total Simulation Time Breakdown by Node Count\n(simtime=1000, nx=2000, Hybrid MPI+OpenMP)',
    fontsize=14,
    fontweight='bold'
)
ax.set_xticks(x)
ax.set_xticklabels(node_labels, rotation=15, ha='right')
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Adjust layout to make room for the rotated x-axis labels and legend
fig.tight_layout() 

# Display the plot
plt.show()

print("Stacked plot generation complete. The total time for each configuration is labeled above the bar.")
