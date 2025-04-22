#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Configuration
###############################################################################
# NOTE: These values are specific to the AMD Ryzen 7 5700U CPU.
# You may need to adjust them for your specific CPU architecture.
# ===========================================================================
# ===========================================================================
# Peak performance in GFLOP/s
peak_gflops = 419114.90 / 1000.0  # e.g. 419.1149 GFLOP/s

# Memory/cache bandwidths in GB/s
bandwidths = [
    # OLD STUFF
    # ("L1", 949004.17 / 1000.0, "red"),      # ~949.0 GB/s
    # ("L2", 534038.39 / 1000.0, "green"),    # ~534.0 GB/s
    # ("L3", 304131.64 / 1000.0, "blue"),     # ~304.1 GB/s
    # ("RAM",35591.18  / 1000.0, "magenta"),  # ~35.6 GB/s

    # ("L1_Multi",   1501047.02 / 1000.0, "orange"),
    # ("L2_Multi",   908374.20 / 1000.0, "purple"),
    # ("L3_Multi",   784738.97 / 1000.0, "brown"),
    # ("RAM_Multi",  34631.31 / 1000.0, "cyan"),

    ("L1", 250063.13 / 1000.0, "red"),
    ("L2", 114290.69 / 1000.0, "green"),
    ("L3", 50153.26 / 1000.0, "blue"),
    ("RAM",22817.37  / 1000.0, "magenta"),
]

######### PUT YOUR APPLICATION POINT HERE #########
# Application point
# ====== EXAMPLE ======
# op_intensity = 0.1602              # FLOPs/Byte
# app_perf     = 15105.8580 / 1000.0 # => ~15.106 GFLOP/s
# app_label    = "epsmim1"
# ====================
app_points = [
    ("epsmim1", 0.1585, 13898.3465 / 1000.0),
    ("epsmim2", 0.1588 , 13739.3066 / 1000.0)
]
######### END OF APPLICATION POINT #########

# Plot range
x_min = 1e-3
x_max = 1e3

###############################################################################
# Generate log-spaced x-values
###############################################################################
x_vals = np.logspace(np.log10(x_min), np.log10(x_max), num=500)

###############################################################################
# Plot setup
###############################################################################
plt.figure(figsize=(8,6))

# Plot each bandwidth line from left to right
for name, bw, color in bandwidths:
    # y = min(peak, bw*x)
    y_vals = np.minimum(peak_gflops, bw*x_vals)
    plt.plot(x_vals, y_vals, color=color, linewidth=2,
             label=f"{name} ({bw:.1f} GB/s)")

# Horizontal line for peak
plt.axhline(y=peak_gflops, color="black", linestyle="--",
            label=f"Peak FLOPs ({peak_gflops:.1f} GFLOP/s)")
# single core performance
single_core_gflops = 67994.91 / 1000.0
plt.axhline(y=single_core_gflops, color="orange", linestyle="--",
            label=f"Single-Core Peak ({single_core_gflops:.1f} GFLOP/s)")

# Define a list of colors for each application point
point_colors = ['red', 'blue', 'green', 'magenta', 'cyan']

# Mark the application points with different colors
for (label, intensity, perf), color in zip(app_points, point_colors):
    plt.plot(intensity, perf, "o", color=color, markersize=8)
    plt.text(intensity*1.2, perf,
             f"{label} ({intensity:.3f}, {perf:.2f})",
             color=color)

###############################################################################
# Axis settings
###############################################################################
plt.xscale("log")
plt.yscale("log")
plt.xlim(x_min, x_max)

# Set a reasonable y-limit so you can see everything
y_min = 1e-1
# compute max over all application points
max_app_perf = max(perf for _, _, perf in app_points)
y_max = max(peak_gflops*2, max_app_perf*2)
plt.ylim(y_min, y_max)

plt.xlabel("Operational Intensity [FLOP/Byte]", fontsize=14)
plt.ylabel("Performance [GFLOP/s]", fontsize=14)
plt.title("Roofline Model for AMD Ryzen 7 5700U", fontsize=16)

plt.grid(True, which="both", linestyle=":")
plt.legend(loc="best")

plt.tight_layout()
plt.savefig("roofline.png", dpi=300)
plt.show()
