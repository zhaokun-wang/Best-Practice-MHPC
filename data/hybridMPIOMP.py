import matplotlib.pyplot as plt
import numpy as np

# --- Simulation Data and Recalculation ---
# Parameters: simtime = 1000, nx = 2000, hybrid mpi + omp (4 ranks/node, 8 threads/rank)

# Data structure: 
# [Node Configuration Name, T_init_raw, T_compute_raw, T_communicate]
raw_data = [
    # New 4 Node Data
    ['4 Nodes (16 Ranks)', 0.0649059862, 1893.17725, 102.034714],
    # Existing Data
    ['8 Nodes (32 Ranks)', 0.0462541431, 916.703613, 60.3360596],
    ['16 Nodes (64 Ranks)', 0.0747045502, 504.719849, 72.6496201],
    ['32 Nodes (128 Ranks)', 0.0443070084, 804.914673, 572.242981]
]

# Separate and process the data
node_labels = [item[0] for item in raw_data]

# 1. T_init (Bottom layer)
T_init = np.array([item[1] for item in raw_data])
# 3. T_communicate (Top layer)
T_communicate = np.array([item[3] for item in raw_data])
# 2. T_Net_Compute (Middle layer = Raw Compute - Communication)
T_compute_raw = np.array([item[2] for item in raw_data])
T_Net_Compute = T_compute_raw - T_communicate

# Calculate the total time (should equal T_init + T_Net_Compute + T_communicate)
T_total = T_init + T_Net_Compute + T_communicate

# --- Plotting Setup ---
fig, ax = plt.subplots(figsize=(10, 6))
bar_width = 0.6  # Width of the stacked bars
x = np.arange(len(node_labels)) # the label locations

# 1. Plot T_init (The base of the stack)
ax.bar(x, T_init, bar_width, label='T_init (Initialization)', color='#4c72b0')

# 2. Plot T_Net_Compute (Stacked on top of T_init)
bottom_net_compute = T_init
rects_net_compute = ax.bar(
    x, T_Net_Compute, bar_width, 
    bottom=bottom_net_compute, 
    label='T_Net_Compute (Raw Compute - Communication)', 
    color='#dd8452'
)

# 3. Plot T_communicate (Stacked on top of T_init + T_Net_Compute)
# The bottom of this bar is T_init + T_Net_Compute.
bottom_comm = T_init + T_Net_Compute
rects_comm = ax.bar(
    x, T_communicate, bar_width, 
    bottom=bottom_comm, 
    label='T_communicate (MPI/OMP Overhead)', 
    color='#55a868'
)


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
    
# Label components that are large enough to be clearly visible
def label_component(rects, data_array, color='white', precision=1):
    for i, rect in enumerate(rects):
        height = rect.get_height()
        # Only label components that are large enough (e.g., > 10% of max height)
        if data_array[i] > 100: 
            mid_point = rect.get_y() + rect.get_height() / 2
            ax.text(
                rect.get_x() + rect.get_width() / 2, 
                mid_point, 
                f'{data_array[i]:.{precision}f}', 
                ha='center', 
                va='center', 
                color=color,
                fontsize=8
            )

# Label the T_Net_Compute components
label_component(rects_net_compute, T_Net_Compute, color='black', precision=0)
# Label the T_communicate components
label_component(rects_comm, T_communicate, color='white', precision=0)


# Set titles and labels
ax.set_ylabel('Time (Seconds)', fontsize=12)
ax.set_title(
    'Total Simulation Time Breakdown by Node Count (T_Net_Compute)\n(simtime=1000, nx=2000, Hybrid MPI+OpenMP)',
    fontsize=14,
    fontweight='bold'
)
ax.set_xticks(x)
ax.set_xticklabels(node_labels, rotation=15, ha='right')
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), frameon=True)
ax.grid(axis='y', linestyle='--', alpha=0.7)
ax.set_ylim(0, max(T_total) * 1.1) # Set y-limit to give space for total time label

# Adjust layout to make room for the rotated x-axis labels and legend
fig.tight_layout() 

# Display the plot
plt.show()

print("Updated stacked plot generation complete. The total time and the three components (Init, Net Compute, Communication) are shown for 4, 8, 16, and 32 nodes.")
