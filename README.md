# LAISD

Accompanying code for the article **"Decoupling Support Enumeration and Value Discovery in Non-Binary ISD"** by Freja Elbro and Paolo Santini.

## Overview

This repository contains implementations of Information Set Decoding (ISD) algorithms for analyzing runtime complexity in the finite and assymptotic regime. The code includes:

- **LA-ISD** (Linear Algebra Information Set Decoding)
- **MitM-LA** (Meet-in-the-Middle Linear Algebra) 
- Comparison algorithms from the literature
- Tools for parameter optimization and complexity analysis


## Prerequisites

- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Asymptotic analysis
In the paper, variables like `k`, `w`, `ell` represent values that scale with `n`. In the asymptotic implementation, we optimize over the ratios `(variable/n)` but maintain the original variable names for simplicity.

### Quick Start

1. **Run the example script:**
   ```python
   python example.py
   ```

2. **Or use the interactive examples below:**


### Single Algorithm Optimization

```python
from BJMM2 import ISD as BJMM2_ISD
from MitM_LA import ISD as MitM_LA_ISD
from MitM_LA_memC import ISD as MitM_LA_memC_ISD
from Dumer import ISD as Dumer_ISD
from background import Hqi # Inverse of q-ary entropy function

# Define parameters
q = 9          # Field size
k = 0.5        # Code rate
w = Hqi(1-k,q) # Weight parameter

# Run algorithm and get results
result, x = Dumer_ISD(k, w, q)
print(f"Runtime complexity: {result.fun}")
print(f"Optimal parameters: {x}")

# Get memory complexity
from Dumer import memory as Dumer_memory
print(f"Memory complexity: {Dumer_memory(x)}")
```

### Using Precomputed Results

The repository includes precomputed optimization sweeps for various field sizes:

```python
from rate_sweep import load_results_for_q, plot_complexity_vs_rate

# Load precomputed results for q=9
q = 9
results = load_results_for_q(q)

# Access results
k_values = results['k_values']                    # List of code rates
optimal_times = results['BJMM2']['times']         # Optimal runtimes
optimal_params = results['BJMM2']['x_values']     # Optimal parameters

# Available algorithms: 'BJMM2', 'MitM_LA', 'MitM_LA_memC', 'Dumer'

# Plot complexity vs rate
plot_complexity_vs_rate(q)
```

### Available Precomputed Data

Optimization sweeps are available for:
- **Field size:** q ∈ [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 31, 37, 43, 49, 59, 71, 81, 97, 149, 199, 251, 307, 359, 409, 461, 563, 673, 773, 877, 977]
- **Code rates:** k = 0.01, 0.02, ..., 0.99
- **Weight:** Given by Gilbert-Varshamov bound

### Running Custom Optimization Sweeps

```python
from rate_sweep import run_optimization_sweep

# Define custom parameters
k_values = [0.45, 0.46, 0.47, 0.48, 0.49, 0.50]
q = 51

# Run optimization (this may take a while, and be aware that it overwrites any other data saved for the selected q-value)
run_optimization_sweep(q, k_values=k_values)
```

### Recreating Paper Figures

```python
from rate_sweep import plot_worst_case_differences, plot_worst_case_memory

# Plot complexity differences
plot_worst_case_differences(q_min=0, q_max=1000)

# Plot memory usage
plot_worst_case_memory(q_min=0, q_max=1000)
```

**On the use of AI:** The code in rate_sweep.py was recreated with help from GitHub Copilot. This was done to create a clearer, more understandable version of the logic, as our original human-written implementation contained numerous optimizations and special cases that made it difficult to follow. The generated code was manually reviewed for correctness, and its results were verified to be identical to those from the original implementation.






