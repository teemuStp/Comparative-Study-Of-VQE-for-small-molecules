################### Imports ###################


# General imports
import matplotlib.pyplot as plt
import pandas as pd


def mapping_renamer(map_name):
    new_name=''

    if map_name == 'parity':
        new_name='PR'
    elif map_name == 'jordan_wigner':
        new_name='JW'
    elif map_name == 'bravyi_kitaev':
        new_name='BK'
    elif map_name == 'neven':
        new_name='OTT'
    return new_name



def reformat_dataframe(data):

    molecules = data['molecule']
    num_qubits = data['num_qubits']
    z2_values = data['z2Symmetries']
    mapping = data['mapping']
    ansatz = data['ansatz']

    for mol in molecules:
        for num_q in num_qubits:
            for z2 in z2_values:
                for map in mapping:
                    for ans in ansatz:
                        # Create a name
                        abr_map = mapping_renamer(map)
                        name = mol +'-'+abr_map+'-'+ans+'-'+str(z2)+'-'+str(num_q)+'q'


# [ansatz_circuit', 'hamiltonian','avg_pauli_weight','avg_hardware_pauli_weight','num_pauli_strings','num_qubits','vqe_energies','iterations','exact_energies','exact_solution','parameters','error'



# plot the data
if __name__ == '__main__':

    # read in the data
    filename = 'vqe_results.csv'
    data = pd.read_csv('../results/'+filename)
    print(data)

    properties = data.keys()
    print(properties)
    molecules = data['molecule']
    num_qubits = data['num_qubits']


