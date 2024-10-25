import re
import pandas as pd

# Input from gprof_analysis.txt to be parsed into DataFrame
gprof_report = '../resources/cpuprof_reports/gprof_analysis.txt'
with open(gprof_report, 'r') as f:
    lines = f.readlines() 

# index where the profile starts
header_idx = next(i for i, line in enumerate(lines) if line.strip().startswith('%'))
start_idx = header_idx + 1

# Extract relevant lines
data_lines = [line for line in lines[start_idx:] if line.strip()]

# Define regex pattern
# Regex def:
# \s* - 0 or more whitespace characters
# (\d+\.\d+) - 1 or more digits followed by a dot followed by 1 or more digits
# \s+ - 1 or more whitespace characters
# (\d+) - 1 or more digits
# (\S+) - 1 or more non-whitespace characters
# \s* - 0 or more whitespace characters
pattern = re.compile(r'\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\S+)')

# Parse lines
data = []
for line in data_lines:
    match = pattern.match(line)
    if match:
        data.append(match.groups())

# Create DataFrame
df = pd.DataFrame(data, columns=[
    '%_time', 'cumulative_seconds', 'self_seconds',
    'calls', 'self_Ts_per_call', 'total_Ts_per_call', 'function_name'
])

# Convert numeric columns
numeric_cols = ['%_time', 'cumulative_seconds', 'self_seconds', 'calls', 'self_Ts_per_call', 'total_Ts_per_call']
df[numeric_cols] = df[numeric_cols].astype(float)

print(df)
