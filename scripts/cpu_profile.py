import re
import pandas as pd

# Input from gprof_analysis.txt to be parsed into DataFrame
gprof_report = '../resources/cpuprof_reports/gprof_analysis.txt'

# Split content into data_lines
lines = open(gprof_report).readlines()

# Find the index where the flat profile starts
start_idx = lines.index('   %   cumulative   self              self     total           ') + 1

# Extract relevant lines
data_lines = lines[start_idx:]

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
