import matplotlib.pyplot as plt
import numpy as np # Import numpy for tick manipulation
from matplotlib.ticker import MaxNLocator # Import MaxNLocator

# Read data from file
file_path = 'cmake-build-release/performance_real.dat'  # <-- Change this to your actual file path

x = []
y = []

with open(file_path, 'r') as file:
    for line in file:
        # Skip empty lines
        if line.strip() == '':
            continue
        parts = line.split()
        if len(parts) >= 2:
            x_value = int(parts[0])
            y_value = int(parts[1])
            x.append(x_value)
            y.append(y_value)

# Create the plot using Axes object for more control
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y, marker='+', linestyle='-', color='r')

# Add labels and title
ax.set_xlabel('Временных слоев за проход')
ax.set_ylabel('Время исполнения')
ax.set_title('Время исполнения по отношению к количеству слоев за раз')
ax.grid(True)

# --- Increase the number of ticks on the x-axis ---
# You can adjust the 'nbins' parameter to control the approximate number of ticks
ax.xaxis.set_major_locator(MaxNLocator(nbins=25)) # Request approximately 15 ticks

# Optional: Rotate tick labels if they overlap
plt.xticks(rotation=0, ha='right') # Rotate labels for better readability if needed

# Adjust layout to prevent labels from being cut off
plt.tight_layout()

# Show the plot
plt.show()
