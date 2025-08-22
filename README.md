# LAISD

Accompanying code for the article **"Decoupling Support Enumeration and Value Discovery in Non-Binary ISD"** by Freja Elbro and Paolo Santini.

## Overview

This repository contains implementations of Information Set Decoding (ISD) algorithms for analyzing both finite regime and asymptotic runtime complexity. The code includes:

- **LA-ISD** (Linear Algebra ISD) implementations
- **MitM-LA** (Meet-in-the-Middle Linear Algebra) algorithms  
- Comparison algorithms from the literature
- Tools for parameter optimization and complexity analysis

## Repository Structure

```
LAISD/
├── Asymptotic_analysis/    # Asymptotic complexity analysis
│   ├── BJMM2.py           # BJMM2 algorithm implementation
│   ├── Dumer.py           # Dumer algorithm implementation  
│   ├── MitM_LA.py         # MitM-LA algorithm
│   ├── MitM_LA_memC.py    # MitM-LA with memory constraints
│   ├── example.py         # Ready-to-run examples
│   ├── background.py      # Utility functions
│   ├── rate_sweep.py      # Rate sweep analysis and plotting
│   └── pickles/           # Precomputed optimization results
├── Finite_regime/         # Finite regime analysis
└── README.md
```

## Prerequisites

- Python 3.x
- NumPy
- SciPy
- Matplotlib


## Asymptotic Analysis

### Quick Start

1. **Navigate to the asymptotic analysis directory:**
   ```bash
   cd Asymptotic_analysis/
   ```

2. **Run the example script:**
   ```python
   python example.py
   ```

3. **Or use the interactive examples below:**

### Usage Examples

#### Single Algorithm Optimization

```python
from BJMM2 import ISD as BJMM2_ISD
from MitM_LA import ISD as MitM_LA_ISD
from MitM_LA_memC import ISD as MitM_LA_memC_ISD
from Dumer import ISD as Dumer_ISD
from background import Hqi

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

#### Using Precomputed Results

The repository includes precomputed optimization sweeps for various field sizes:

```python
from rate_sweep_background import load_results_for_q, plot_complexity_vs_rate

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

#### Available Precomputed Data

Optimization sweeps are available for:
- **Field size:** q ∈ [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 31, 37, 43, 49, 59, 71, 81, 97, 149, 199, 251, 307, 359, 409, 461, 563, 673, 773, 877, 977]
- **Code rates:** k = 0.01, 0.02, ..., 0.99
- **Weight:** Given by Gilbert-Varshamov bound

#### Running Custom Optimization Sweeps

```python
from rate_sweep_background import run_optimization_sweep

# Define custom parameters
k_values = [0.45, 0.46, 0.47, 0.48, 0.49, 0.50]
q = 51

# Run optimization (this may take a while, and be aware that it overwrites any other data saved for the selected q-value)
run_optimization_sweep(q, k_values=k_values)
```

#### Recreating Paper Figures

```python
from rate_sweep_background import plot_worst_case_differences, plot_worst_case_memory

# Plot complexity differences
plot_worst_case_differences(q_min=0, q_max=1000)

# Plot memory usage
plot_worst_case_memory(q_min=0, q_max=1000)
```

**Note:** In the paper, variables like `k`, `w`, `ell` represent values that scale with `n`. In the asymptotic implementation, we optimize over the ratios `(variable/n)` but maintain the original variable names for simplicity.





## Finite Regime Analysis

*[Section to be added - guide for finite regime analysis]*

## Reference

### Algorithms Implemented

1. **BJMM2**: Becker-Joux-May-Meurer algorithm with two levels
2. **Dumer**: Dumer's algorithm
3. **MitM_LA**: Meet-in-the-Middle Linear Algebra ISD
4. **MitM_LA_memC**: MitM-LA with possibility to constrain memory

### Variable Conventions



## License

See [LICENSE](LICENSE) file for details.

