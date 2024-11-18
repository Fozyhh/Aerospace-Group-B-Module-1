import re
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Available metrics
METRICS = {
    'IR': ('Ir', 'Instruction Reads'),
    'L1IMISS': ('I1mr', 'L1 Instruction Cache Misses'),
    'LLIMISS': ('ILmr', 'LL Instruction Cache Misses'),
    'DR': ('Dr', 'Data Reads'),
    'L1DMISS': ('D1mr', 'L1 Data Cache Misses'),
    'LLDMISS': ('DLmr', 'LL Data Cache Misses'),
    'DW': ('Dw', 'Data Writes'),
    'L1WMISS': ('D1mw', 'L1 Data Write Misses'),
    'LLWMISS': ('DLmw', 'LL Data Write Misses'),
    'BC': ('Bc', 'Conditional Branches'),
    'BCM': ('Bcm', 'Conditional Branches Mispredicted'),
    'BI': ('Bi', 'Indirect Branches'),
    'BIM': ('Bim', 'Indirect Branches Mispredicted')
}

def print_available_metrics():
    print("Available metrics:")
    for key, (_, description) in METRICS.items():
        print(f"  {key}: {description}")

if len(sys.argv) != 2 or sys.argv[1] not in METRICS:
    print("Usage: python3 mem_profile.py METRIC")
    print_available_metrics()
    sys.exit(1)

# Select metric
selected_metric_key = sys.argv[1]
selected_metric, metric_description = METRICS[selected_metric_key]

# Read file
with open('../resources/memprof_reports/cachegrind.out', 'r') as file:
    content = file.read()

# Cache configuration from desc lines
cache_config = {}
for line in content.split('\n'):
    if line.startswith('desc:'):
        match = re.search(r'(\w+)\s+cache:\s+(\d+)\s+B', line)
        if match:
            cache_type, size = match.groups()
            cache_config[cache_type] = int(size)

data = []

# Profile data
current_file = ""
current_function = ""

for line in content.split('\n'):
    if line.startswith('fl='):
        current_file = line[3:]
    elif line.startswith('fn='):
        current_function = line[3:]
    elif re.match(r'^\d+', line):
        values = line.split()
        if len(values) >= 13: 
            data.append({
                'file': current_file,
                'function': current_function,
                'line': int(values[0]),
                'Ir': int(values[1]),    # Instruction reads
                'I1mr': int(values[2]),  # L1 instruction cache misses
                'ILmr': int(values[3]),  # LL instruction cache misses
                'Dr': int(values[4]),    # Data reads
                'D1mr': int(values[5]),  # L1 data cache misses
                'DLmr': int(values[6]),  # LL data cache misses
                'Dw': int(values[7]),    # Data writes
                'D1mw': int(values[8]),  # L1 data write misses
                'DLmw': int(values[9]),  # LL data write misses
                'Bc': int(values[10]),   # Conditional branches
                'Bcm': int(values[11]),  # Conditional branches mispredicted
                'Bi': int(values[12]),   # Indirect branches
                'Bim': int(values[13])   # Indirect branches mispredicted
            })

df = pd.DataFrame(data)

# Total metrics per function
function_metrics = df.groupby('function').sum()

top_functions = function_metrics.sort_values(selected_metric, ascending=False).head(10)

# Create visualization
plt.figure(figsize=(12, 12))
plt.bar(range(len(top_functions)), top_functions[selected_metric])
plt.xticks(range(len(top_functions)), top_functions.index, rotation=45, ha='right')
plt.title(f'Top 10 Functions by {metric_description}')
plt.xlabel('Function Name')
plt.ylabel(metric_description)
plt.tight_layout()

plt.savefig(f'../resources/memprof_reports/cachegrind_{selected_metric_key.lower()}.png')

print(f"\nTop 5 Functions by {metric_description}:")
print(top_functions[selected_metric].head().to_string())
