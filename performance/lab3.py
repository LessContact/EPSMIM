import matplotlib.pyplot as plt

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

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(x, y, marker='+', linestyle='-', color='r')

# Add labels and title
plt.xlabel('X Values')
plt.ylabel('Y Values')
plt.title('Plot from File Data')
plt.grid(True)

# Show the plot
plt.show()
