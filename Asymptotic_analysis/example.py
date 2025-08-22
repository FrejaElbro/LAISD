from background import *

# Import ISD algorithms
from BJMM2 import ISD as BJMM2_ISD
from MitM_LA import ISD as MitM_LA_ISD
from MitM_LA_memC import ISD as MitM_LA_memC_ISD
from Dumer import ISD as Dumer_ISD

# Example 1: Single algorithm optimization
print("=== Single Algorithm Optimization ===")
q = 9
k = 0.5 
w = Hqi(1-k, q)

result, x = Dumer_ISD(k, w, q)
print(f"Runtime complexity: {result.fun}")
print(f"Optimal parameters: {x}")

# Get memory complexity
from Dumer import memory as Dumer_memory
print(f"Memory complexity: {Dumer_memory(x)}")

print("\n=== Using Precomputed Results ===")
# Load and use precomputed results
from rate_sweep import load_results_for_q, plot_complexity_vs_rate

q = 9
results = load_results_for_q(q)
k_values = results['k_values']
optimal_BJMM2_times = results['BJMM2']['times']
optimal_BJMM2_parameters = results['BJMM2']['x_values']

print(f"Available k values: {len(k_values)} points from {min(k_values)} to {max(k_values)}")
print(f"BJMM2 time for k=0.5: {optimal_BJMM2_times[49]}")  # k=0.5 is at index 49
print(f"Optimal BJMM2 parameters for k=0.5: {optimal_BJMM2_parameters[49]}")

# Uncomment to plot results
# plot_complexity_vs_rate(q)

print("\n=== Custom Optimization Sweep (commented out as it edits the pickles, and every available pickle is automatically included in the graphs above) ===")
# Uncomment to run custom optimization sweep
# from rate_sweep import run_optimization_sweep
# k_values_custom = [0.45, 0.46, 0.47, 0.48, 0.49, 0.50]
# run_optimization_sweep(51, k_values=k_values_custom) # Run the sweep for q = 51. Note that any existing results for q = 51 will be overwritten

print("\n=== Paper Figures (commented out) ===")
# Uncomment to recreate paper figures
# from rate_sweep import plot_worst_case_differences, plot_worst_case_memory
# plot_worst_case_differences(q_min=0, q_max=100)  # Use smaller q_max for testing
# plot_worst_case_memory(q_min=0, q_max=100) 