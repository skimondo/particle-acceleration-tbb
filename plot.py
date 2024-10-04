import matplotlib.pyplot as plt

# Read data from the local file
filename = 'bench-random.dat'

num_threads = []
parallel_times = []
accelerations = []

with open(filename, 'r') as file:
    lines = file.readlines()

# Parsing the data
for line in lines[1:]:  # Skip the header line
    if line.strip() and not line.startswith("#"):  # Check for non-empty line and ignore comments
        parts = line.split()
        num_threads.append(int(parts[0]))
        parallel_times.append(int(parts[1]))
        accelerations.append(float(parts[2]))

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Execution time plot
axs[0].plot(num_threads, parallel_times, 'bo-', label='Parallel Time', markersize=8)
axs[0].set_title('Parallel Execution Time vs Number of CPUs')
axs[0].set_xlabel('Number of cores (ncpu)')
axs[0].set_ylabel('Execution Time (nanoseconds)')
axs[0].grid()
axs[0].legend()
axs[0].set_xticks(range(1, 21))  # Set x-axis ticks from 1 to 20

# Acceleration plot
serial_time = parallel_times[0]  # Using the first time as the serial time
accelerated_values = [serial_time / t for t in parallel_times]  # Calculate acceleration
axs[1].plot(num_threads, accelerated_values, 'bo-', label='Parallel Acceleration', markersize=8)
axs[1].set_title('Parallel Acceleration vs Number of CPUs')
axs[1].set_xlabel('Number of cores (ncpu)')
axs[1].set_ylabel('Acceleration')
axs[1].grid()
axs[1].legend()
axs[1].set_xticks(range(1, 21))  # Set x-axis ticks from 1 to 20

plt.tight_layout()

# Save the plot to a file
plt.savefig('execution_time_and_acceleration_comparison.png')
plt.show()

