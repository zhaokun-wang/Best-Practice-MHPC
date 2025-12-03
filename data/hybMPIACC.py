import matplotlib.pyplot as plt
import numpy as np

# --- Simulation Data and Recalculation (OpenACC-MPI) ---
# Configuration: 4 GPUs / 4 Ranks per node, 1 GPU per task
# T_Net_Compute = T_compute_raw - T_communicate

# Data structure: 
# [Node Configuration Name, T_init_raw, T_compute_raw, T_communicate]
raw_data = [
    # 4 Node Data
    ['4 Nodes (16 GPUs)', 0.054568682, 381.5820, 48.62727],
    # 8 Node Data (T_output component has been removed)
    ['8 Nodes (32 GPUs)', 0.055355541, 374.7395, 56.22041],
    # 16 Node Data (Revised)
    ['16 Nodes (64 GPUs)', 0.088977821, 322.2292, 81.76420],
    # 32 Node Data
    ['32 Nodes (128 GPUs)', 0.083296478, 304.7112, 94.49528]
]

# Separate and process the data
node_labels = [item[0] for item in raw_data]

# 1. T_init (Bottom layer)
T_init = np.array([item[1] for item in raw_data])
# 3. T_communicate (Top layer of the stack)
T_communicate = np.array([item[3] for item in raw_data])

# Calculate T_Net_Compute (Middle layer = Raw Compute - Communication)
T_compute_raw = np.array([item[2] for item in raw_data])
T_Net_Compute = T_compute_raw - T_communicate

# Calculate the total time 
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
    label='T_Net_Compute (GPU Kernel Execution)', 
    color='#dd8452'
)

# 3. Plot T_communicate (Stacked on top of T_init + T_Net_Compute)
bottom_comm = T_init + T_Net_Compute
rects_comm = ax.bar(
    x, T_communicate, bar_width, 
    bottom=bottom_comm, 
    label='T_communicate (MPI/Data Transfer Overhead)', 
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
        # Only label components that are large enough (e.g., > 100 seconds)
        if data_array[i] > 50: 
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
    'OpenACC-MPI (GPU-Accelerated) Time Breakdown\n(4 GPUs/Node, 1 GPU/Task, simtime=1000, nx=2000)',
    fontsize=14,
    fontweight='bold'
)
ax.set_xticks(x)
ax.set_xticklabels(node_labels, rotation=15, ha='right')
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), frameon=True)
ax.grid(axis='y', linestyle='--', alpha=0.7)
ax.set_ylim(0, max(T_total) * 1.1) 

# Adjust layout to make room for the rotated x-axis labels and legend
fig.tight_layout() 

# Display the plot
plt.show()

print("Revised OpenACC-MPI stacked plot generation complete. The 16-node data has been updated, and the T_output component has been removed.")
