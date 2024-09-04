###################### Imports ############################


# ignore warnings
import warnings
warnings.filterwarnings('ignore')


# Command line argument imports
import sys


# General imports
import numpy as np
import pandas as pd


# Pre-defined ansatz circuit and operator class for Hamiltonian
from qiskit.circuit.library import EfficientSU2
from qiskit.quantum_info import Z2Symmetries


# SciPy minimizer routine
from scipy.optimize import minimize


# imort NUmpy algorithm to solve minumun energy exactly
from qiskit_algorithms import NumPyEigensolver


# Plotting functions
import matplotlib.pyplot as plt


# runtime imports 1.0 qiskit properties
#from qiskit_ibm_runtime import QiskitRuntimeService, Session
#from qiskit_ibm_runtime import EstimatorV2 as Estimator


# Some qiskit imports
from qiskit.primitives import Estimator,BackendEstimator
from qiskit.providers.fake_provider import GenericBackendV2 # qiskit 1.0 property
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager # qiskit 1.0 property


# qiskit_nature imports
from qiskit_nature.second_q.mappers import ParityMapper, JordanWignerMapper, BravyiKitaevMapper,QubitMapper,TaperedQubitMapper # QubitMapper is the new QubitConverter/FermionicTransformation
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo
from qiskit_nature.second_q.transformers import FreezeCoreTransformer
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD

# Import the ternary tree mapper from self made neven_mapping.py
# from neven_mapping import NevenMapper 


#############    Define some constant methods   ############


# All possible implementation methods
all_mappers = ['parity','bravyi_kitaev', 'jordan_wigner']
all_ansatzes = ['EfficientSU2','UCCSD']
all_Z2Symmetries_list = [True,False]
all_measurement_schemes = ['pauli_scheme','QWC']
# Molecule list is below the chemistry_molecues dictionary

# Default options
default_measurement_scheme = ['pauli_scheme']
default_ansatz = ['EfficientSU2']
default_mapping = ['jordan_wigner']
default_z2 = [True]


#############################################################


# Create the chemistry molecules remove this later
nh3 = 0.5
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
        coords=[(0.0,0.0,0.0),(0.757, 0.586, 0.0),(-0.757, 0.586, 0.0)],
        charge=0,
        multiplicity=1,
    ),
    "NH3": MoleculeInfo(
        symbols=["N","H","H","H"],
        coords=[
        (0.0, 0.0, 0.1),     # Nitrogen
        (0.94, 0.0, -0.27),  # First Hydrogen
        (-0.47, 0.81, -0.27),# Second Hydrogen
        (-0.47, -0.81, -0.27)# Third Hydrogen
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
        (1.34, 0.0, 0.0),   # Second Carbon
        (-0.67, 0.94, 0.0), # First Hydrogen
        (-0.67, -0.94, 0.0),# Second Hydrogen
        (2.01, 0.94, 0.0),  # Third Hydrogen
        (2.01, -0.94, 0.0)  # Fourth Hydrogen
        ],
        charge=0,
        multiplicity=1,
    ),
    'Cr2': MoleculeInfo(
        symbols=["Cr","Cr"],
        coords=[(0.0,0.0,0.0),(1.0,0.0,0.0)],
        charge=0,
        multiplicity=1,
    ),
    'O3': MoleculeInfo(
        symbols=["O","O","O"],
        coords=[(0.0,0.0,0.0),(1.0,0.0,0.0),(-1.0,0.0,0.0)],
        charge=0,
        multiplicity=1,
    ),
    'O2': MoleculeInfo(
        symbols=["O","O"],
        coords=[(0.0,0.0,0.0),(1.0,0.0,0.0)],
        charge=0,
        multiplicity=1,
    )
}

all_molecule_names = list(chemistry_molecules.keys())

################### Define the functions ####################

def prepare_hamiltonian(
    molecule_name, z2symmetry_reduction, mapping,freeze_core=True, basis='sto-3g'
):  
    """Creates the Hamiltoanin for a molecule with given mapping 
    
    Args: molecule(MoleculeInfo object): the molecule to create the Hamiltonian for
          z2symmetry_reduction(bool): whether to reduce the Hamiltonian using Z2 symmetries
          mapping(str): the mapping to use
          freeze_core(bool):  whether to freeze the core orbitals
          basis(str):  the basis set to use
    Returns: qubitOp(SparsePauliOp): the Hamiltonian as a (paulis=[],coeff=[])"""

    # Execute the PySCF driver
    molecule = chemistry_molecules[molecule_name]
    driver = PySCFDriver.from_molecule(molecule=molecule, basis=basis)

    if mapping == "parity":
        #qubit_mapping = FermionicQubitMappingType.PARITY
        qubit_mapping = ParityMapper()
    elif mapping == "bravyi_kitaev":
        #qubit_mapping = FermionicQubitMappingType.BRAVYI_KITAEV
        qubit_mapping = BravyiKitaevMapper()
    #elif mapping == "neven":
    #    qubit_mapping = FermionicQubitMappingType.NEVEN
    elif mapping == "jordan_wigner":
        #qubit_mapping = FermionicQubitMappingType.JORDAN_WIGNER
        qubit_mapping = JordanWignerMapper()
    else:
        raise ValueError("Wrong mapping")
    total_hamiltonian = driver.run()
    #total_hamiltonian = problem.hamiltonian.second_q_op()
    
    #total_number_of_orbitals = 6

    # Apply the freeze core transformation
    transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    #transformer.prepare_active_space(molecule,total_number_of_orbitals)
    reduced_hamiltonian = transformer.transform(total_hamiltonian)

    # Convert the Hamiltonian to second quantized form
    hamiltonian = reduced_hamiltonian.hamiltonian.second_q_op()
    
    # Apply the given fermion to Pauli encoding
    qubitOp = QubitMapper.map(self=qubit_mapping,
        second_q_ops=hamiltonian
        )
    
    # Apply Z2 symmetries (except for hydrogen)
    if molecule_name != 'H2':
        if z2symmetry_reduction:
            z2symmetries = Z2Symmetries.find_z2_symmetries(qubitOp)
            qubitOp = z2symmetries.taper(qubitOp)
    
    if molecule_name == 'H2':
        return qubitOp
    else:
        return qubitOp[0]

def pauli_weight(pauli_string):
    """Caclulates the Pauli weight of a given Pauli string (NUmber of Pauyli operators in a string). This method calculates the theoretical value if circut can execute
        X, Y and Z in on operation.
        
    Args:
        pauli_string (str): Pauli string
    Returns:
            Pauli Weight (int)"""
    weight = 0
    for oper in str(pauli_string):
        if oper == 'Z' or oper == 'X' or oper == 'Y':
            weight += 1
    return weight

def hardware_pauli_weight(pauli_string):
    """This calculates the weight of a Pauli string that can be executed on the RZ based Quantum circuit. Z and I do not add any Pauli weight,
        , but X adds two and Y adds 4.

        Args: pauli_string (str): Pauli string
        Returns: Pauli Weight (int)"""
    weight = 0
    for oper in str(pauli_string):
        # To measure X, we need to apply Hadamard gate, and in actaul circuit H = Rz(theta1)*Rz(theta2)
        if oper == 'X':
            weight += 2
        # To measure Y we need to apply Sdg gate, and in actual circuit SH = 
        elif oper == 'Y':
            weight += 4
    return weight

def build_ansatz(ansatz_name,mapper,num_qubits,num_particles,reps):
    """Create hartreeFock anstaz with selected parametrization

    Args: ansatz_name (str): name of the ansatz
          mapper (str): name of the mapper
          num_qubits (int): number of qubits
          num_particles (int): number of particles
          reps (int): number of repetitions
          
    Returns: ansatz (QuantumCircuit): ansatz circuit"""
    
    if mapper == 'parity':
        qubit_mapping = ParityMapper()
    elif mapper == 'bravyi_kitaev': 
        qubit_mapping = BravyiKitaevMapper()
    elif mapper == 'jordan_wigner':
        qubit_mapping = JordanWignerMapper()
    # elif mapper == 'neven':
    #     qubit_mapping = FermionicQubitMappingType.NEVEN
    else:
        raise ValueError('Invalid mapper name')

    # Spin is zero for our molecules, so there is equal number of alpha and beta particles
    alpha = int(num_particles/2)
    beta =  int(num_particles/2)

    # create HartreeFock state
    hf_state = HartreeFock(num_spatial_orbitals=int(num_qubits/2), num_particles=(alpha,beta), qubit_mapper=qubit_mapping)

    if ansatz_name=='UCCSD':
        ansatz = UCCSD(num_spatial_orbitals=int(num_qubits/2), num_particles=(alpha,beta), qubit_mapper=qubit_mapping)
    elif ansatz_name=='EfficientSU2':
        ansatz = EfficientSU2(num_qubits=int(num_qubits), entanglement='linear', reps=reps)
    else: 
        raise ValueError('Invalid ansatz name')

    return hf_state.compose(ansatz)

def cost_func(params, ansatz, hamiltonian, estimator,cost_history_dict):
    """Return estimate of energy from estimator

    Parameters:
        params (ndarray): Array of ansatz parameters
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
        estimator (EstimatorV2): Estimator primitive instance
        cost_history_dict: Dictionary for storing intermediate results

    Returns:
        float: Energy estimate
    """
    pub = (ansatz, [hamiltonian], [params])
    result = estimator.run(ansatz,hamiltonian,params).result()
    energy = result.values

    cost_history_dict["iters"] += 1
    cost_history_dict["prev_vector"] = params
    cost_history_dict["cost_history"].append(energy)
    cost_history_dict["vectors"].append(params)
    #print(f"Iters. done: {cost_history_dict['iters']} [Current cost: {energy}]")

    return energy
    #optimize the ansatz for backend

def vqe(num_qubits,ansatz,hamiltonian):
    """Run the VQE optimization algorithm. This function uses the COBYLA optimizer from SciPy. 

    Args: num_qubits (int): number of qubits
          ansatz (QuantumCircuit): ansatz circuit
          hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
    Returns: cost_history_dict (dict): dictionary containing the cost and parameter history of the optimization"""
    
    # Get the backend
    backend = GenericBackendV2(num_qubits)
    
    # Get the number of parameters in the ansatz
    num_params = ansatz.num_parameters
    
    # Optimize the ansatz for the backend
    target = backend.target
    pm = generate_preset_pass_manager(target=target, optimization_level=3)
    ansatz_isa = pm.run(ansatz)
    hamiltonian_isa = hamiltonian.apply_layout(layout=ansatz_isa.layout)

    cost_history_dict = {
        "prev_vector": None,
        "iters": 0,
        "cost_history": [],
        "vectors": [],
    }
    #x0 = 2 * np.pi * np.random.random(num_params)
    x0 = np.zeros(num_params)

    estimator = BackendEstimator(backend)
    estimator.options.default_shots = 10_000
    # set the maximum number of iterations
    options = {'maxiter': 400}

    res = minimize(
        cost_func,
        x0,
        args=(ansatz_isa, hamiltonian_isa, estimator,cost_history_dict),
        method="cobyla",
        options=options,
    )
    return cost_history_dict

def command_line_parser(arguments):


    """Parse the command line arguments and return the correct arguments for the VQE calculation
    Args: arguments (list): list of command line arguments without the first element
    Returns: mappings, ansatzes, molecules, measurement_schemes, z2_symmetry (list): list of the arguments for the VQE calculation"""

    try:
        arguments[0]    
    except:
        print("Please provide the correct number of arguments \n")
        print('run_type=all/default/custom VQE=Y/N mol1,mol2,mol3 map1,map2,map3 ansatz1,ansatz2,ansatz3 mes1,mes2,mes3 z2_symmetry=True,False')
        sys.exit()

        

    if arguments[0] == "all" or arguments[0] == "default":


        # Run all the possible combinations
        if arguments[0] == "all":
            print('Running all the molecules, mappings, ansatz, measurement schemes and z2_symmetry')
            print("!!! This will take a long time !!!")
            mappings = all_mappers
            ansatzes = all_ansatzes
            molecules = all_molecule_names
            measurement_schemes = all_measurement_schemes
            z2_symmetry = all_Z2Symmetries_list
            VQE = 'Y'
    

        # Run the default methods
        else:
            print('Running the comparison of methods respect to default methods')
            mappings = default_mapping
            ansatzes = default_ansatz
            molecules = all_molecule_names
            measurement_schemes = default_measurement_scheme
            z2_symmetry = default_z2
            VQE = 'Y'


    # create a check thst the arguments are correct
    elif len(arguments) != 7: 
        print("Please provide the correct number of arguments")
        print('run_type=all/default/custom mol1,mol2,mol3 map1,map2,map3 ansatz1,ansatz2,ansatz3 mes1,mes2,mes3 z2_symmetry=True,False')
        sys.exit()
    else: 
        print('Running a custom set of arguments')
        VQE = arguments[1]
        molecules = arguments[2].split(',')
        mappings = arguments[3].split(',')
        ansatzes = arguments[4].split(',')    
        measurement_schemes = arguments[5].split(',')
        # map the z2_symmetry to a boolean
        z2_s = arguments[6].split(',') 
        z2_symmetry=[]
        for z in z2_s:
            if z == 'True':
                z = True
            elif z == 'False':
                z = False
            z2_symmetry.append(z)

        # print all the methods
        print('VQE:',VQE)
        print('Molecules:',molecules)
        print('Mappings:',mappings)
        print('Ansatzes:',ansatzes)
        print('Measurement schemes:',measurement_schemes)
        print('Z2 Symmetry:',z2_symmetry)


    
        # checke that the arguments can be foun from list that conatin all different methods
        if not all([m in all_mappers for m in mappings]):
            print('The mapping method is not correct')
            sys.exit()
        if not all([a in all_ansatzes for a in ansatzes]):
            print('The ansatz method is not correct')
            sys.exit()
        if not all([m in all_molecule_names for m in molecules]):
            print('The molecule name is not correct')
            sys.exit()
        if not all([m in all_measurement_schemes for m in measurement_schemes]):
            print('The measurement scheme list is not correct')
            sys.exit()
        if not all([m in all_Z2Symmetries_list for m in z2_symmetry]):
            print('The z2_symmetry arg is not correct')
            sys.exit()

        print('All arguments correct!')
    return VQE,molecules, mappings, ansatzes, measurement_schemes, z2_symmetry
  

############################################################


# Main function
if __name__=='__main__':

    # read arguments from the commandline when executed  run_type=all/default_mapping/default_ansatz/default_Z2Symmetries/default_measurement_scheme molecules=all/H2...
    try:
        arguments  = sys.argv[1:]
    except:
        print('Please provide the correct number of arguments')
        print('run_type=all/default/custom VQE=Y/N mol1,mol2,mol3 map1,map2,map3 ansatz1,ansatz2,ansatz3 mes1,mes2,mes3 z2_symmetry=True,False')
        sys.exit()
    
    # parse the command line arguments
    VQE,molecules, mappers, ansatzes, measurement_schemes, Z2Symmetries_list = command_line_parser(arguments)
    

    # Create a pandas DataFrame to store the Hamiltonians
    data = pd.DataFrame(columns=['molecule','z2Symmetries', 'mapping', 'ansatz','ansatz_circuit', 'hamiltonian','avg_pauli_weight','avg_hardware_pauli_weight','num_pauli_strings','num_qubits','vqe_energies','iterations','exact_energies','exact_solution','parameters','error'])


    # Iterate over all different parameters and store the results in the DataFrame named data

    for map in mappers:
        for z2sym in Z2Symmetries_list:
            for molecule in molecules:
                for ansatz in ansatzes:
                    print("Running VQE for molecule:", molecule, "| z2Symmetries:", z2sym, "| mapping:", map, "| ansatz:", ansatz)


                    # create the problem hamiltonian
                    hamiltonian = prepare_hamiltonian(molecule_name=molecule, z2symmetry_reduction=z2sym, freeze_core=True, mapping=map)
            

                    # retrieve and calculate some useful information
                    num_qubits = hamiltonian.num_qubits
                    num_pauli_strings = len(hamiltonian.paulis)
                    avg_pauli_weight = np.mean([pauli_weight(pauli) for pauli in hamiltonian.paulis])
                    avg_hardware_pauli_weight = np.mean([hardware_pauli_weight(pauli) for pauli in hamiltonian.paulis])
            

                    # For laptops, simulating more than 10 qubits becomes too cucumbersome so we create a limit
                    run_vqe = (VQE=='Y')

            
                    # perform the vqe calculation for the given molecule if the system is small (<10 qubits)
                    if(run_vqe):    


                        # create the ansatz circuit
                        ansatz_circuit = build_ansatz(ansatz_name=ansatz, mapper=map, num_qubits=num_qubits, num_particles=2, reps=2)


                        # perform the VQE calculation
                        vqe_results = vqe(num_qubits,ansatz_circuit,hamiltonian)
                        vqe_energies = vqe_results['cost_history']
                        iterations = vqe_results['iters']
                        parameters = vqe_results['vectors']


                        # calculate the exact energies for the given parameters
                        estimator = Estimator()
                        exact_energies = []
                        for i in range(len(parameters)):
                            exact_energies.append(estimator.run(circuits=ansatz_circuit,observables=hamiltonian,parameter_values=parameters[i]).result().values[0])


                        # calculate the exact solution
                        exact_solution = np.real(NumPyEigensolver(k=1).compute_eigenvalues(hamiltonian).eigenvalues[0])
                        

                        # calculate the error
                        error = [float(abs(exact_energies[i]-vqe_energies[i])) for i in range(len(vqe_energies))]


                    # if the system is too large, we do not perform the VQE calculation and return empty values
                    else:
                        vqe_energies = []
                        iterations = 0
                        parameters = []
                        exact_energies = []
                        exact_solution = 0.0
                        error = []
                        ansatz_circuit = None

            
                    # store the results in the DataFrame
                    new_row = {'molecule':molecule, 'z2Symmetries': str(z2sym),'mapping':map, 'ansatz':ansatz,'ansatz_circuit':ansatz_circuit,
                                'hamiltonian':hamiltonian, 'avg_pauli_weight':avg_pauli_weight, 'num_pauli_strings':num_pauli_strings,
                                'num_qubits':num_qubits, 'vqe_energies':vqe_energies,'iterations':iterations,'parameters':parameters,
                                'error':error,'exact_energies':exact_energies,'exact_solution':exact_solution,'avg_hardware_pauli_weight':avg_hardware_pauli_weight}
                    data.loc[len(data)] = new_row


    # save the results in a csv file
    print('All VQE calculations done. Saving the results to a file vqe_results.csv')
    data.to_csv('../results/vqe_results.csv',index=False)
    print('Results saved successfully')
