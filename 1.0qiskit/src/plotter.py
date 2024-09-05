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

    mappings = ['parity', 'jordan_wigner', 'bravyi_kitaev', 'neven']

    # Numbers of qubits 
    for map in mappings:
        print('\n'+map+'\n')

        # Filter the data
        try:
            map_data = data[data['mapping']==map]
            num_qubits = map_data['num_qubits']
            avg_pauli_weight = map_data['avg_pauli_weight']
            avg_hardware_pauli_weight = map_data['avg_hardware_pauli_weight']
            num_pauli_strings = map_data['num_pauli_strings']

            # Calculate the scaling
            plot_range = np.arange(2,(max(num_qubits)+4))

            if map == 'parity':
                scaling = PR_scaling(plot_range)
                name = 'Parity'
            elif map == 'jordan_wigner':
                scaling = JW_scaling(plot_range)
                name = 'Jordan-Wigner'
            elif map == 'bravyi_kitaev':
                scaling = BK_scaling(plot_range)
                name = 'Bravyi-Kitaev'
            elif map == 'neven':
                scaling = OTT_scaling(plot_range)
                name = 'Ternary Tree'


            

            # Create scatter plot, since avg_pauli_weight can have multiple values for the same x values
            plt.plot(num_qubits, avg_pauli_weight, 'o', label='Average Pauli Weight')

            #plt.plot(num_qubits, num_pauli_strings, 'o', label='Number of Pauli Strings')
            plt.plot(plot_range, scaling, label=map+' scaling')
            plt.legend()
            plt.xlabel('Number of qubits')
            plt.ylabel('Value')
            plt.title(name+' scaling')
            plt.grid('both',linestyle='--')
            plt.savefig('../results/Pauli_weight_'+map+'.pdf', format='pdf', dpi=1000)
            plt.show()

            # Plot the hardware weight and pauli weight 
            plt.plot(num_qubits, avg_hardware_pauli_weight, 'o', label='Rz based hadrware gates')
            plt.plot(num_qubits, avg_pauli_weight, 'o', label='Ideal qubit gates')
            plt.plot(plot_range, scaling, label='Theoretical scaling')
            plt.legend()
            plt.xlabel('Number of qubits')
            plt.ylabel('Number of qubit operations')
            plt.title(map+'scaling')
            plt.grid('both',linestyle='--')
            plt.savefig('../results/Hardware_Pauli_weight_'+map+'.pdf', format='pdf', dpi=1000)
            plt.show()


            
        except:
            # Jump to next mapping
            continue



