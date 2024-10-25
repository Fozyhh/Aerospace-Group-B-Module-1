import re
import pandas as pd
import matplotlib.pyplot as plt

# Read the cachegrind.out file
with open('../resources/memprof_reports/cachegrind.out', 'r') as file:
    content = file.read()

# Fictionary to hold function cache data
functions = {}

lines = content.strip().split('\n')

current_file = None
current_function = None

for line in lines:
    # Skip header lines and empty lines
    if line.startswith(('version:', 'creator:', 'pid:', 'cmd:', 'part:', 'desc:', 'positions:', 'events:', 'summary:', 'totals:')) or not line.strip():
        continue

    # Check for file and function definitions
    if line.startswith('fl='):
        current_file = line[3:]
        continue
    if line.startswith('fn='):
        current_function = line[3:]
        continue
        
    # Process metrics line
    if current_file and current_function and re.match(r'^\d', line):
        parts = line.strip().split()
        if len(parts) >= 10:  # Line number + 9 metrics
            metrics = list(map(int, parts[1:]))
            function_key = f"{current_file}:{current_function}"
            if function_key not in functions:
                functions[function_key] = {
                    'Ir': metrics[0],    # Instructions
                    'I1mr': metrics[1],  # L1 instruction cache misses
                    'ILmr': metrics[2],  # LL instruction cache misses
                    'Dr': metrics[3],    # Data reads
                    'D1mr': metrics[4],  # L1 data cache read misses
                    'DLmr': metrics[5],  # LL data cache read misses
                    'Dw': metrics[6],    # Data writes
                    'D1mw': metrics[7],  # L1 data cache write misses
                    'DLmw': metrics[8],  # LL data cache write misses
                }
            else:
                # Add metrics if same function appears multiple times
                for i, key in enumerate(['Ir', 'I1mr', 'ILmr', 'Dr', 'D1mr', 'DLmr', 'Dw', 'D1mw', 'DLmw']):
                    functions[function_key][key] += metrics[i]

# Convert to DataFrame
cache_df = pd.DataFrame.from_dict(functions, orient='index')

print("\nCache Usage DataFrame:")
print(cache_df)

# Save to CSV
# cache_df.to_csv('cachegrind_analysis.csv')

# Calculate total cache misses (sum of all miss events)
cache_df['total_cache_misses'] = cache_df[['I1mr', 'ILmr', 'D1mr', 'DLmr', 'D1mw', 'DLmw']].sum(axis=1)

# Sort by total_cache_misses descending
cache_sorted = cache_df.sort_values(by='total_cache_misses', ascending=False)

# Plot all functions with cache misses
plt.figure(figsize=(12,6))
plt.barh(cache_sorted.index, cache_sorted['total_cache_misses'], color='coral')
plt.xlabel('Total Cache Misses')
plt.title('Functions Ranked by Cache Misses')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('../resources/memprof_reports/cachegrind_cache_misses.png')
plt.show()
