################### Imports ###################


# General imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

################### Functions ###################


def mapping_renamer(map_name):

    """Renames the mapping to a shorter form
    
    Input: map_name (str) - the name of the mapping
    
    Output: new_name (str) - the new name of the mapping"""

    new_name=''

    if map_name == 'parity':
        new_name='PR'
    elif map_name == 'jordan_wigner':
        new_name='JW'
    elif map_name == 'bravyi_kitaev':
        new_name='BK'
    elif map_name == 'neven':
        new_name='OTT'
    else:
        new_name='Unknown'
    return new_name


def reformat_dataframe(data):

    molecules = data['molecule']
    num_qubits = data['num_qubits']
    z2_values = data['z2Symmetries']
    mapping = data['mapping']

    for mol in molecules:
        for num_q in num_qubits:
            for z2 in z2_values:
                for map in mapping:
                    # Create a name
                    abr_map = mapping_renamer(map)
                    name = mol +'-'+abr_map+'-Z2_'+str(z2)+'-'+str(num_q)+'q'
                    print(name)
    return ''


def JW_scaling(num_qubits):
    
        """Calculates the scaling of the Jordan-Wigner mapping
        
        Input: num_qubits (int) - the number of qubits
        
        Output: scaling (int) - the scaling of the Jordan-Wigner mapping"""

        return np.array(num_qubits)


def BK_scaling(num_qubits):
    
        """Calculates the scaling of the Bravyi-Kitaev mapping
        
        Input: num_qubits (int) - the number of qubits
        
        Output: scaling (int) - the scaling of the Bravyi-Kitaev mapping"""

        num_qubits = np.array(num_qubits)

        return np.log2(num_qubits)


def PR_scaling(num_qubits):
        
    """Calculates the scaling of the Parity mapping
            
            Input: num_qubits (int) - the number of qubits
            
            Output: scaling (int) - the scaling of the Parity mapping"""
    
    num_qubits = np.array(num_qubits)
    
    return num_qubits


def OTT_scaling(num_qubits):
            
    """Calculates the scaling of the Parity mapping
                
                Input: num_qubits (int) - the number of qubits
                
                Output: scaling (int) - the scaling of the OTT mapping"""
        
    num_qubits = np.array(num_qubits)
        
    return math.log((2*num_qubits),3)


# [ansatz_circuit', 'hamiltonian','avg_pauli_weight','avg_hardware_pauli_weight','num_pauli_strings','num_qubits','vqe_energies','iterations','exact_energies','exact_solution','parameters','error'


################### Main ###################

# plot the data
if __name__ == '__main__':


    # read in the data
    filename = 'vqe_results.csv'
    data = pd.read_csv('../results/'+filename)
    

    # reformat the data
    dropping = ['ansatz', 'ansatz_circuit',
       'hamiltonian','vqe_energies', 'iterations',
       'exact_energies', 'exact_solution', 'parameters', 'error']
    data = data.drop(columns=dropping)

    # Numbers of qubits 
    num_qubits = data['num_qubits']
    pauli_weights = data['avg_pauli_weight']
    hardware_pauli_weights = data['avg_hardware_pauli_weight']

    print(data[(data['molecule']=='H2O')])




