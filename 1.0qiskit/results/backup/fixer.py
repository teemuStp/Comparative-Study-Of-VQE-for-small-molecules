import pandas as pd
import ast

# Read the CSV file
filename = 'default_run.csv'
file_path = '/home/steemu/Documents/GitHub/Gradu/1.0qiskit/results/'+filename
df = pd.read_csv(file_path)

# Function to convert string representation of arrays to list of floats
def convert_to_floats(array_str):
    # Use ast.literal_eval to safely evaluate the string as a Python expression
    array_list = ast.literal_eval(array_str)
    # Extract the float values from the arrays
    float_list = [float(arr[0]) for arr in array_list]
    return float_list

# Apply the conversion function to the 'vqe_energies' column
df['vqe_energies'] = df['vqe_energies'].apply(convert_to_floats)

# Write the modified DataFrame back to the CSV file
df.to_csv(file_path, index=False)