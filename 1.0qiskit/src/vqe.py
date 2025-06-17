###################### Imports ############################


# ignore warnings
import warnings
warnings.filterwarnings('ignore')


# Command line argument imports
import sys


# General imports
import numpy as np
import pandas as pd
import csv
import time
from pathlib import Path


# SciPy minimizer routine
from scipy.optimize import minimize
import spsa


# imort NUmpy algorithm to solve minumun energy exactly
from qiskit_algorithms import NumPyEigensolver


# runtime imports 1.0 qiskit properties
#from qiskit_ibm_runtime import QiskitRuntimeService, Session
#from qiskit_ibm_runtime import EstimatorV2 as Estimator


# Some qiskit imports
from qiskit.primitives import Estimator


# qiskit_nature imports
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo

# Import the self-built functions
import functions

#############    Define some constant methods   ############


# Error message for the command line arguments
arg_error_message = 'Please provide the correct number of arguments \n run_type=all/default/custom filename VQE=Y/N mol1,mol2,mol3 map1,map2,map3 ansatz1,ansatz2,ansatz3 mes1,mes2,mes3 z2_symmetry=True,False'


path_to_this_file = str(Path.cwd())


# All possible implementation methods
all_mappers = ['parity','bravyi_kitaev', 'jordan_wigner','neven']
all_ansatzes = ['EfficientSU2','kUpCCGSD']
all_Z2Symmetries_list = [True,False]
all_measurement_schemes = ['pauli_scheme','QWC']
all_anzats_reps = [1,2,3]
# Molecule list is below the chemistry_molecues dictionary


# Default options
default_measurement_scheme = ['pauli_scheme']
default_ansatz = ['EfficientSU2']
default_mapping = ['jordan_wigner']
default_z2 = [True]
default_ansatz_reps = [2]


ansatz_reps = [1,2,3]


optimizer = 'cobyla'
max_iter = 200
shots = 1024

#############################################################

# Create the chemistry molecules remove this later
# Molecules included in the study
# H2, LiH, H2O, NH3, He, C2H4, Cr2, O3

chemistry_molecules = {
    "H2": MoleculeInfo(
        symbols=["H","H"],
        coords=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.75)],
        charge=0,
        multiplicity=1,
    ),
    "LiH": MoleculeInfo(
        symbols=["Li","H"],
        coords=[(0.0, 0.0, 0.0),(0.0, 0.0, 1.6)],
        charge=0,
        multiplicity=1,
    ),
    "H2O": MoleculeInfo(
        symbols = ["O", "H", "H"],
        #coords=[(0.0,0.0,0.0),(0.757, 0.586, 0.0),(-0.757, 0.586, 0.0)],
        coords=[(0.0,0.0,0.116321),(0.0,0.751155,-0.465285),(0.0,-0.751155,-0.465285)],
        charge=0,
        multiplicity=1,
    ),
    "NH3": MoleculeInfo(
        symbols=["N","H","H","H"],
        coords=[
        (0.888, 1.538, 1.069),     # Nitrogen
        (0.0, 0.0, 0.0),  # First Hydrogen
        (0.0, 3.076, 0.0),# Second Hydrogen
        (2.664, 1.538, 0.0)# Third Hydrogen
        ],
        charge=0,
        multiplicity=1,
    ),
    "He": MoleculeInfo(
        symbols=["He"],
        coords=[(0.0,0.0,0.0)],
        charge=0,
        multiplicity=1,
    ),
    'C2H4': MoleculeInfo(
        symbols=["C","C","H","H","H","H"],
        coords=[
        (0.0, 0.0, 0.0),    # First Carbon
        (1.3304, 0.0, 0.0),   # Second Carbon
        (-0.56380327217722, 0.922210, 0.0), # First Hydrogen
        (-0.56380327217722, -0.922210, 0.0),# Second Hydrogen
        (1.89420327217722, 0.922210, 0.0),  # Third Hydrogen
        (1.89420327217722, -0.922210, 0.0)  # Fourth Hydrogen
        ],
        charge=0,
        multiplicity=1,
    ),
    #'Cr2': MoleculeInfo(
    #    symbols=["Cr","Cr"],
    #    coords=[(0.0,0.0,0.0),(1.0,0.0,0.0)],
    #    charge=0,
    #    multiplicity=1,
    #),
    #'O3': MoleculeInfo(
    #    symbols=["O","O","O"],
    ##    coords=[(-1.082,-0.665,0.0),(0.0,0.0,0.0),(1.082,-0.665,0.0)],
     #   charge=1,
    #    multiplicity=1,
    #),
    'O2': MoleculeInfo(
        symbols=["O","O"],
        coords=[(0.0,0.0,0.0),(0.0,0.0,1.2)],
        charge=0,
        multiplicity=1,
    )
}

all_molecule_names = list(chemistry_molecules.keys())

################### Define the functions ####################



############################################################


# Main function
if __name__=='__main__':


    # read arguments from the commandline when executed  run_type=all/default_mapping/default_ansatz/default_Z2Symmetries/default_measurement_scheme molecules=all/H2...
    try:
        arguments  = sys.argv[1:]
    except:
        print(arg_error_message)
        sys.exit()
    

    # parse the command line arguments
    filename,VQE,molecules, mappers, ansatzes, measurement_schemes, Z2Symmetries_list,ansatz_reps = command_line_parser(arguments)


    # Create a pandas DataFrame to store the Hamiltonians
    data = pd.DataFrame(columns=['molecule',         'z2Symmetries', 'mapping', 'ansatz',
                                'vqe_time',           'hamiltonian',
                                'avg_pauli_weight',  'num_pauli_strings',
                                'num_qubits',        'vqe_energies',
                                'iterations',        'parameters',
                                'error',             'exact_energies',
                                'exact_solution',    'avg_hardware_pauli_weight',
                                'max_pauli_weight',  'max_hrdwr_pauli_weight',
                                'num_parameters',    'gates',
                                'depth',             'ansatz_reps',
                                'classical_time',
                                'accuracies_shots',
                                'vqe_acc'])


    # flush the vqe_results.csv file adn write the header
    with open(path_to_this_file+'/../results/'+filename+'.csv','w') as file:
        writer = csv.writer(file)
        writer.writerow(data.columns)
        file.close()

    print('')
    # Iterate over all different parameters and store the results in the DataFrame named data
    for molecule in molecules:
        #for reps in ansatz_reps:
        for reps in [2]:
            for z2sym in Z2Symmetries_list:
                for map in mappers:
                    for ansatz in ansatzes:
                        for measurement in measurement_schemes:
                            print("Preparing: molecule:", molecule, "| z2Symmetries:", z2sym, "| mapping:", map, "| ansatz:", ansatz, "| measurement:", measurement)


                            # create the problem hamiltoniansource ~/demo/VQE/bin/activate

                            hamiltonian,num_spatial_orbitals,num_particles = prepare_hamiltonian(molecule_name=molecule, z2symmetry_reduction=z2sym, mapping=map)

                            # Sve Hamitlonian to a file
                            if map != 'neven':
                                _ = save_hamiltonian(hamiltonian,path_to_this_file+'/../hamiltonians/'+molecule+'-'+map+'-'+str(z2sym)+'.txt')

                            # retrieve and calculate some useful information
                            num_qubits = hamiltonian.num_qubits
                            num_pauli_strings = len(hamiltonian.paulis)


                            # Pauli weight information
                            pauli_weights = [pauli_weight(pauli) for pauli in hamiltonian.paulis]
                            max_pauli_weight = np.max([pauli_weight(pauli) for pauli in hamiltonian.paulis])
                            avg_pauli_weight = np.mean([pauli_weight(pauli) for pauli in hamiltonian.paulis])
                            hardware_pauli_weights = [hardware_pauli_weight(pauli) for pauli in hamiltonian.paulis]
                            max_hardware_pauli_weight = np.max([hardware_pauli_weight(pauli) for pauli in hamiltonian.paulis])
                            avg_hardware_pauli_weight = np.mean([hardware_pauli_weight(pauli) for pauli in hamiltonian.paulis])
            
                            # solve exact solution for the mapping
                            start  = time.time()
                            exact_solution = np.real(NumPyEigensolver(k=1).compute_eigenvalues(hamiltonian).eigenvalues[0])
                            print(exact_solution)
                            class_time_cost = time.time() - start


                            # spaghetti code for now
                            run_vqe = (VQE=='Y')


                            print('Preparation done!')

                            if False:
                                #if molecule=='NH3' or molecule=='C2H4':
                                #    shot = 10_000#,2000]#,5000,10_000,20_000,30_000,40_000,60_000,80_000,100_000]
                                #    vars = variance(hamiltonian,shot,ansatz_circuit)
                                #    shots = [1000,2000,5000,10_000,20_000,40_000,60_000,80_000,100_000,150_000,200_000]
                                #    accuracies = [approx_accracy(hamiltonian,shot,vars) for shot in shots]
                                #else:
                                shots = [100]#,2000,5000,10_000,20_000,30_000,40_000,60_000,80_000,100_000]
                                accuracies = [calc_accuracy(hamiltonian,shot,ansatz_circuit,True) for shot in shots]
                                accuracies_shots = (accuracies,shots)
                            else: 
                                accuracies_shots = (None,None)
            
                            # perform the vqe calculation for the given molecule if the system is small (<10 qubits)
                            if(run_vqe):    


                                print('Running molecule:',molecule,'with',num_qubits,'qubits')


                                #perform the VQE calculation
                                if (measurement == 'pauli_scheme' or measurement == 'QWC'):
                                    if measurement == 'pauli_scheme':
                                        grouping = False
                                    else:
                                        grouping = True
                                    start = time.time()
                                    vqe_results = vqe(num_qubits,ansatz_circuit,hamiltonian,grouping,shots=shots)
                                    vqe_time_cost = time.time() - start
                        
                                vqe_energies = vqe_results['cost_history']
                                iterations = vqe_results['iters']
                                parameters = vqe_results['vectors']
                                vqe_acc = vqe_results['vqe_acc']


                                # calculate the exact energies for the given parameters
                                estimator = Estimator()
                                exact_energies = []
                                for i in range(len(parameters)):
                                    exact_energies.append(estimator.run(circuits=ansatz_circuit,observables=hamiltonian,parameter_values=parameters[i]).result().values[0])



                                # calculate the error                                
                                error = [float(abs(exact_energies[i]-vqe_energies[i])) for i in range(len(vqe_energies))]
                                                            # calculate the exact solution
        


                            # if the system is too large, we do not perform the VQE calculation and return empty values
                            else:
                                vqe_energies = []
                                iterations = 0
                                parameters = []
                                exact_energies = []
                                error = []
                                vqe_time_cost = 0.0
                                class_time_cost = 0.0
                                vqe_acc  = []

            
                            # store the results in the DataFrame

                            new_row = { 
                                        'molecule': molecule,
                                        'z2Symmetries': str(z2sym),
                                        'mapping': map,
                                        'ansatz': ansatz,
                                        'vqe_time': vqe_time_cost,
                                        'hamiltonian': hamiltonian,
                                        'avg_pauli_weight': avg_pauli_weight,
                                        'num_pauli_strings': num_pauli_strings,
                                        'num_qubits': num_qubits,
                                        'vqe_energies': vqe_energies,
                                        'iterations': iterations,
                                        'parameters': parameters,
                                        'error': error,
                                        'exact_energies': exact_energies,
                                        'exact_solution': exact_solution,
                                        'avg_hardware_pauli_weight': avg_hardware_pauli_weight,
                                        'max_pauli_weight': max_pauli_weight,
                                        'max_hrdwr_pauli_weight': max_hardware_pauli_weight,
                                        'num_parameters': num_params,
                                        'gates': cx_gates,  # assuming cx_gates corresponds to 'gates' in your original DataFrame
                                        'depth': depth,
                                        'ansatz_reps': reps,
                                        'classical_time': class_time_cost,
                                        'accuracies_shots': accuracies_shots,
                                        'vqe_acc': vqe_acc}


                            # Write to csv
                            with open(path_to_this_file+'/../results/'+filename+'.csv','a') as file:
                                writer = csv.writer(file)
                                writer.writerow(new_row.values())

                            print('One system done and save. Moving to next system...\n')
        print('Moving to the next molecule...\n')

    print('Results saved successfully')