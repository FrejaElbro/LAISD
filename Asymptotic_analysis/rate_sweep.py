import pickle
import numpy as np
from background import *
import collections
import matplotlib.pyplot as plt
import os
import glob

# Define the namedtuple structures for each algorithm
BJMM2_vars = collections.namedtuple('BJMM2_vars', 'p d ell i')
Dumer_vars = collections.namedtuple('Dumer_vars', 'wp ell')
MitM_LA_vars = collections.namedtuple('MitM_LA_vars', 'wp d ell i')
MitM_LA_memC_vars = collections.namedtuple('MitM_LA_memC_vars', 'wr wc d ell i kr')

def namedtuple_to_dict(obj, algorithm):
    """Convert namedtuple to dictionary for pickling"""
    if obj is None:
        return None
    return obj._asdict()

def dict_to_namedtuple(obj_dict, algorithm):
    """Convert dictionary back to namedtuple after unpickling"""
    if obj_dict is None:
        return None
    
    if algorithm == 'BJMM2':
        return BJMM2_vars(**obj_dict)
    elif algorithm == 'Dumer':
        return Dumer_vars(**obj_dict)
    elif algorithm == 'MitM_LA':
        return MitM_LA_vars(**obj_dict)
    elif algorithm == 'MitM_LA_memC':
        return MitM_LA_memC_vars(**obj_dict)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")

def run_optimization_sweep(q, k_values=[i/100.0 for i in range(1, 100)] ):
    """Run ISD optimization for multiple algorithms across k values and q values"""
    
    # Import the ISD functions from different algorithms
    from BJMM2 import ISD as BJMM2_ISD
    from MitM_LA import ISD as MitM_LA_ISD
    from MitM_LA_memC import ISD as MitM_LA_memC_ISD
    from Dumer import ISD as Dumer_ISD, L0 as DL0  # Need DL0 for MitM_LA_memC
    
    
    print(f"Running optimization for q = {q} and {len(k_values)} k values")
    
    
    # Initialize storage 
    q_results = {
        'BJMM2': {'times': [], 'x_values': []},
        'MitM_LA': {'times': [], 'x_values': []},
        'MitM_LA_memC': {'times': [], 'x_values': []},
        'Dumer': {'times': [], 'x_values': []},
        'k_values': k_values,
        'q': q
    }
    
    for k in k_values:
        w = Hqi(1-k, q)  # Calculate w using the Hqi function
        
        print(f"  k = {k:.2f}")
            
        # Dumer
        try:
            print(f"    Dumer ", end=' ')
            result, x = Dumer_ISD(k, w, q)
            q_results['Dumer']['times'].append(result.fun)
            q_results['Dumer']['x_values'].append(namedtuple_to_dict(x, 'Dumer'))
            dumer_x = x  # Store for Paolo1mem use
            print("✓")
        except Exception as e:
            print(f"✗ Error: {e}")
            dumer_x = None
            q_results['Dumer']['times'].append(None)
            q_results['Dumer']['x_values'].append(None)
        
        # BJMM2
        try:
            print(f"    BJMM2 ", end=' ')
            result, x = BJMM2_ISD(k, w, q)
            q_results['BJMM2']['times'].append(result.fun)
            q_results['BJMM2']['x_values'].append(namedtuple_to_dict(x, 'BJMM2'))
            print("✓")
        except Exception as e:
            print(f"✗ Error: {e}")
            q_results['BJMM2']['times'].append(None)
            q_results['BJMM2']['x_values'].append(None)
        
        # MitM_LA
        try:
            print(f"    MitM_LA ", end=' ')
            result, x = MitM_LA_ISD(k, w, q)
            q_results['MitM_LA']['times'].append(result.fun)
            q_results['MitM_LA']['x_values'].append(namedtuple_to_dict(x, 'MitM_LA'))
            print("✓")
        except Exception as e:
            print(f"✗ Error: {e}")
            q_results['MitM_LA']['times'].append(None)
            q_results['MitM_LA']['x_values'].append(None)
        
        # MitM_LA_memC (requires Dumer L0)
        try:
            print(f"    MitM_LA_memC ", end=' ')
            if dumer_x is not None:
                result, x = MitM_LA_memC_ISD(k, w, q, memory_constraint = DL0(dumer_x))
                q_results['MitM_LA_memC']['times'].append(result.fun)
                q_results['MitM_LA_memC']['x_values'].append(namedtuple_to_dict(x, 'MitM_LA_memC'))
                print("✓")
            else:
                # If Dumer failed, we can't run MitM_LA_memC
                q_results['MitM_LA_memC']['times'].append(None)
                q_results['MitM_LA_memC']['x_values'].append(None)
                print("✗ Skipped (Dumer failed)")
        except Exception as e:
            print(f"✗ Error: {e}")
            q_results['MitM_LA_memC']['times'].append(None)
            q_results['MitM_LA_memC']['x_values'].append(None)
    
    # Save intermediate results after each q value
    with open(f'pickles/optimization_sweep_results_q{q}.pkl', 'wb') as f:
        pickle.dump(q_results, f)
    print(f"  Saved results for q = {q} to pickles/optimization_sweep_results_q{q}.pkl")
    
    
    return q_results

def load_results_for_q(q):
    """Load the optimization sweep results for a specific q value"""
    try:
        with open(f'pickles/optimization_sweep_results_q{q}.pkl', 'rb') as f:
            results = pickle.load(f)
        return results
    except FileNotFoundError:
        print(f"Results file for q={q} not found. Run the optimization first.")
        return None



def find_optimal_k_for_q(q, algorithm='BJMM2'):
    """Find the optimal k value for a given q and algorithm"""
    results = load_results_for_q(q)
    if results is None:
        return None
    
    times = results[algorithm]['times']
    k_values = results['k_values']
    
    # Find the maximum time (best result) and corresponding k
    max_time = float('-inf')
    optimal_k = None
    optimal_x = None
    
    for i, (time_val, k_val) in enumerate(zip(times, k_values)):
        if time_val is not None and time_val > max_time:
            max_time = time_val
            optimal_k = k_val
            # Convert dictionary back to namedtuple
            x_dict = results[algorithm]['x_values'][i]
            optimal_x = dict_to_namedtuple(x_dict, algorithm)
    
    return {
        'optimal_k': optimal_k,
        'max_time': max_time,
        'optimal_x': optimal_x,
        'q': q,
        'algorithm': algorithm
    }

def compare_algorithms_for_q(q):
    """Compare all algorithms for a given q value"""
    results = load_results_for_q(q)
    if results is None:
        return None
    
    comparison = {}
    for algorithm in ['Dumer', 'BJMM2', 'MitM_LA', 'MitM_LA_memC']:
        optimal = find_optimal_k_for_q(q, algorithm)
        if optimal:
            comparison[algorithm] = optimal
    
    return comparison

def plot_complexity_vs_rate(q):
    """
    Plot runtime and memory complexity for all algorithms vs code rate
    
    Args:
        q: Field size
    """
    
    # Load optimization results for this specific q
    results = load_results_for_q(q)
    if results is None:
        print(f"No optimization results found for q={q}. Run optimization first.")
        return
    
    # Get k values (rates) from results
    k_values = results['k_values']
    
    # Initialize storage for complexities
    algorithms = ['Prange', 'Dumer', 'MitM_LA_memC', 'MitM_LA', 'BJMM2']
    
    time_complexities = {alg: [] for alg in algorithms}
    memory_complexities = {alg: [] for alg in algorithms}
    
    print(f"Computing complexities for q = {q} with {len(k_values)} rate values...")
    
    for i, rate in enumerate(k_values):
        w = Hqi(1-rate, q)  # Convert rate to weight
        
        for alg in algorithms:
            if alg == 'Prange':
                # Prange complexity
                time_comp = binomH(1, w) - binomH(1-rate, w)
                memory_comp = 0  # Prange uses minimal memory
            else:
                # Get the result from optimization data
                time_val = results[alg]['times'][i]
                x_dict = results[alg]['x_values'][i]
                
                if time_val is not None and x_dict is not None:
                    # Convert dictionary back to namedtuple
                    x = dict_to_namedtuple(x_dict, alg)
                    time_comp = time_val
                    
                    # Get memory using the namedtuple - need to set globals first
                    if alg == 'BJMM2':
                        import BJMM2
                        BJMM2.k = rate
                        BJMM2.w = w
                        BJMM2.q = q
                        memory_comp = BJMM2.memory(x)
                    elif alg == 'Dumer':
                        import Dumer
                        Dumer.k = rate
                        Dumer.w = w
                        Dumer.q = q
                        memory_comp = Dumer.memory(x)
                    elif alg == 'MitM_LA':
                        import MitM_LA
                        MitM_LA.k = rate
                        MitM_LA.w = w
                        MitM_LA.q = q
                        memory_comp = MitM_LA.memory(x)
                    elif alg == 'MitM_LA_memC':
                        import MitM_LA_memC
                        MitM_LA_memC.k = rate
                        MitM_LA_memC.w = w
                        MitM_LA_memC.q = q
                        memory_comp = MitM_LA_memC.memory(x)
                else:
                    time_comp = float('inf')
                    memory_comp = float('inf')
            
            time_complexities[alg].append(time_comp)
            memory_complexities[alg].append(memory_comp)
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot time complexity
    ax1.set_title(f'Time Complexity vs Code Rate (q = {q})')
    ax1.set_xlabel('Code Rate')
    ax1.set_ylabel('Time Complexity (log₂)')
    ax1.grid(True, alpha=0.3)
    
    for i, alg in enumerate(algorithms):
        # Filter out infinite values for plotting
        valid_times = [(r, t) for r, t in zip(k_values, time_complexities[alg]) if t != float('inf')]
        if valid_times:
            linestyle = '--' if alg == 'MitM_LA_memC' else '-'
            valid_rates, valid_time_vals = zip(*valid_times)
            ax1.plot(valid_rates, valid_time_vals, label=alg, linestyle=linestyle)

    ax1.legend()
    ax1.set_ylim(bottom=0)
    
    # Plot memory complexity
    ax2.set_title(f'Memory Complexity vs Code Rate (q = {q})')
    ax2.set_xlabel('Code Rate')
    ax2.set_ylabel('Memory Complexity (log₂)')
    ax2.grid(True, alpha=0.3)
    
    for i, alg in enumerate(algorithms):
        # Filter out infinite values for plotting
        valid_memory = [(r, m) for r, m in zip(k_values, memory_complexities[alg]) if m != float('inf')]
        if valid_memory:
            linestyle = '--' if alg == 'MitM_LA_memC' else '-'
            valid_rates, valid_memory_vals = zip(*valid_memory)
            ax2.plot(valid_rates, valid_memory_vals, label=alg, linestyle=linestyle)

    ax2.legend()
    ax2.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.show()
    



def plot_worst_case_memory(q_min = 0, q_max = 1000):
    """
    Load all stored pickles and plot the worst-case memory complexity 
    for each algorithm as a function of q.
    """
    
    # Find all available pickle files
    pickle_pattern = 'pickles/optimization_sweep_results_q*.pkl'
    pickle_files = glob.glob(pickle_pattern)
    
    if not pickle_files:
        print("No pickle files found. Run optimization first.")
        return
    
    # Extract q values from filenames and sort them
    q_values = []
    for file in pickle_files:
        # Extract q value from filename like 'optimization_sweep_results_q3.pkl'
        filename = os.path.basename(file)
        q_str = filename.replace('optimization_sweep_results_q', '').replace('.pkl', '')
        try:
            q = int(q_str)
            q_values.append(q)
        except ValueError:
            print(f"Warning: Could not parse q value from {filename}")
    
    q_values.sort()
    q_values = [q for q in q_values if q_min <= q <= q_max]
    
    if not q_values:
        print(f"No q values found in range [{q_min}, {q_max}]")
        return
    
    print(f"Found results for q values in range [{q_min}, {q_max}]: {q_values}")
    
    # Initialize storage for worst-case memory complexities
    algorithms = ['Prange', 'Dumer', 'BJMM2', 'MitM_LA', 'MitM_LA_memC']
    colors = ['purple', 'blue', 'red', 'green', 'orange']
    
    # Store worst-case memory for each algorithm
    worst_case_memory = {alg: [] for alg in algorithms}
    
    #print("Computing worst-case memory complexities...")
    
    for q in q_values:
        #print(f"Processing q = {q}")
        
        # Load results for this q
        results = load_results_for_q(q)
        if results is None:
            print(f"  Warning: Could not load results for q = {q}")
            # Add None values to maintain alignment
            for alg in algorithms:
                worst_case_memory[alg].append(None)
            continue
        
        k_values = results['k_values']
        
        # Compute worst-case memory for each algorithm
        for alg in algorithms:
            if alg == 'Prange':
                # Prange uses minimal memory (essentially 0)
                worst_case_memory[alg].append(0)
            else:
                memory_complexities = []
                
                for i, rate in enumerate(k_values):
                    w = Hqi(1-rate, q)
                    x_dict = results[alg]['x_values'][i]
                    
                    if x_dict is not None:
                        # Convert dictionary back to namedtuple
                        x = dict_to_namedtuple(x_dict, alg)
                        
                        # Get memory using the namedtuple - need to set globals first
                        try:
                            if alg == 'BJMM2':
                                import BJMM2
                                BJMM2.k = rate
                                BJMM2.w = w
                                BJMM2.q = q
                                memory_comp = BJMM2.memory(x)
                            elif alg == 'Dumer':
                                import Dumer
                                Dumer.k = rate
                                Dumer.w = w
                                Dumer.q = q
                                memory_comp = Dumer.memory(x)
                            elif alg == 'MitM_LA':
                                import MitM_LA
                                MitM_LA.k = rate
                                MitM_LA.w = w
                                MitM_LA.q = q
                                memory_comp = MitM_LA.memory(x)
                            elif alg == 'MitM_LA_memC':
                                import MitM_LA_memC
                                MitM_LA_memC.k = rate
                                MitM_LA_memC.w = w
                                MitM_LA_memC.q = q
                                memory_comp = MitM_LA_memC.memory(x)
                            
                            memory_complexities.append(memory_comp)
                        except Exception as e:
                            # Skip this data point if memory calculation fails
                            pass
                
                # Find maximum memory complexity
                if memory_complexities:
                    alg_worst_memory = max(memory_complexities)
                    worst_case_memory[alg].append(alg_worst_memory)
                else:
                    worst_case_memory[alg].append(None)
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    for i, alg in enumerate(algorithms):
        # Filter out None values for plotting
        valid_data = [(q, mem) for q, mem in zip(q_values, worst_case_memory[alg]) if mem is not None]
        
        if valid_data:
            valid_q_vals, valid_memory = zip(*valid_data)
            linestyle = '--' if alg == 'MitM_LA_memC' else '-'
            plt.plot(valid_q_vals, valid_memory, 'o-', color=colors[i], label=alg, 
                    linestyle=linestyle, markersize=6, linewidth=2)
    
    plt.xlabel('Field Size (q)', fontsize=12)
    plt.ylabel('Worst-case Memory Complexity (log₂)', fontsize=12)
    plt.title('Worst-case Memory Complexity vs Field Size', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Set x-axis to show integer q values
    plt.xticks(q_values)
    plt.ylim(bottom=0)
    
    plt.tight_layout()
    plt.show()
    
    
    return {
        'q_values': q_values,
        'worst_case_memory': worst_case_memory,
        'algorithms': algorithms
    }

def plot_worst_case_differences(q_min = 0, q_max = 1000):
    """
    Load all stored pickles and plot the difference between worst-case time complexity 
    for Prange's algorithm and each of the algorithms as a function of q.
    Difference = Algorithm_worst - Prange_worst
    """
    
    # Find all available pickle files
    pickle_pattern = 'pickles/optimization_sweep_results_q*.pkl'
    pickle_files = glob.glob(pickle_pattern)
    
    if not pickle_files:
        print("No pickle files found. Run optimization first.")
        return
    
    # Extract q values from filenames and sort them
    q_values = []
    for file in pickle_files:
        # Extract q value from filename like 'optimization_sweep_results_q3.pkl'
        filename = os.path.basename(file)
        q_str = filename.replace('optimization_sweep_results_q', '').replace('.pkl', '')
        try:
            q = int(q_str)
            q_values.append(q)
        except ValueError:
            print(f"Warning: Could not parse q value from {filename}")
    
    q_values.sort()
    q_values = [q for q in q_values if q_min <= q <= q_max]
    print(f"Found results for q values: {q_values}")
    
    # Initialize storage for worst-case complexities
    algorithms = ['Dumer', 'BJMM2', 'MitM_LA', 'MitM_LA_memC']
    colors = ['blue', 'red', 'green', 'orange']
    
    # Store worst-case differences for each algorithm
    worst_case_differences = {alg: [] for alg in algorithms}
    
    #print("Computing worst-case complexities...")
    
    for q in q_values:
        #print(f"Processing q = {q}")
        
        # Load results for this q
        results = load_results_for_q(q)
        if results is None:
            print(f"  Warning: Could not load results for q = {q}")
            # Add None values to maintain alignment
            for alg in algorithms:
                worst_case_differences[alg].append(None)
            continue
        
        k_values = results['k_values']
        
        # Compute Prange worst-case over all k values
        prange_complexities = []
        for rate in k_values:
            w = Hqi(1-rate, q)
            prange_time = binomH(1, w) - binomH(1-rate, w)
            prange_complexities.append(prange_time)
        
        prange_worst = max(prange_complexities)
        #print(f"  Prange worst-case: {prange_worst:.3f}")
        
        # Compute worst-case for each algorithm
        for alg in algorithms:
            times = results[alg]['times']
            
            # Filter out None values and find maximum
            valid_times = [t for t in times if t is not None]
            
            if valid_times:
                alg_worst = max(valid_times)
                difference = alg_worst - prange_worst
                worst_case_differences[alg].append(difference)
                #print(f"  {alg} worst-case: {alg_worst:.3f}, difference: {difference:.3f}")
            else:
                worst_case_differences[alg].append(None)
                #print(f"  {alg}: No valid results")
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    for i, alg in enumerate(algorithms):
        # Filter out None values for plotting
        valid_data = [(q, diff) for q, diff in zip(q_values, worst_case_differences[alg]) if diff is not None]
        
        if valid_data:
            valid_q_vals, valid_diffs = zip(*valid_data)
            linestyle = '--' if alg == 'MitM_LA_memC' else '-'
            plt.plot(valid_q_vals, valid_diffs, 'o-', color=colors[i], label=alg, 
                    linestyle=linestyle, markersize=6, linewidth=2)
    
    plt.xlabel('Field Size (q)', fontsize=12)
    plt.ylabel('Worst-case Difference\n(Algorithm - Prange)', fontsize=12)
    plt.title('Worst-case Time Complexity Differences vs Prange Algorithm', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Add horizontal line at y=0 for reference
    plt.axhline(y=0, color='black', linestyle=':', alpha=0.5, linewidth=1)
    
    # Set x-axis to show integer q values
    plt.xticks(q_values)
    
    plt.tight_layout()
    plt.show()
    
    
    return {
        'q_values': q_values,
        'differences': worst_case_differences,
        'algorithms': algorithms
    }