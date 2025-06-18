

# ignore warnings
import warnings
warnings.filterwarnings('ignore')


# Command line argument imports
import sys


# General imports
import numpy as np
from collections import OrderedDict


# Pre-defined ansatz circuit and operator class for Hamiltonian
from qiskit.circuit.library import EfficientSU2
from qiskit.quantum_info import Z2Symmetries,SparsePauliOp
from qiskit.converters import circuit_to_dag, dag_to_circuit


# SciPy minimizer routine
from scipy.optimize import minimize


# runtime imports 1.0 qiskit properties
#from qiskit_ibm_runtime import QiskitRuntimeService, Session
#from qiskit_ibm_runtime import EstimatorV2 as Estimator


# Some qiskit imports
from qiskit.primitives import BackendEstimator
from qiskit.providers.fake_provider import GenericBackendV2 # qiskit 1.0 property
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager # qiskit 1.0 property
from qiskit.circuit import QuantumCircuit


# qiskit_nature imports
from qiskit_nature.second_q.mappers import ParityMapper, JordanWignerMapper, BravyiKitaevMapper,QubitMapper,OptimalGeneralTernaryTreeMapper # QubitMapper is the new QubitConverter/FermionicTransformation
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.transformers import FreezeCoreTransformer
from qiskit_nature.second_q.circuit.library import HartreeFock, PUCCSD

from constants import chemistry_molecules


##### Functions for Qubit wise commuting Pauli strings ######


def group_generator(pauli_list: SparsePauliOp):
    """Returns the generator of a Pauli list that contains QUBIT WISE COMMUTING pauli strings"""
    pauli_ops = pauli_list.paulis
    max_length = len(pauli_ops[0])
    generator_terms = ['I']*max_length
    
    Pauli_ops_by_qubits = []
    
    for qubit in range(max_length):
        Pauli_ops_by_qubits = [str(pauli)[qubit] for pauli in pauli_ops]
        
        for op in Pauli_ops_by_qubits:
            if op != 'I':
                generator_terms[qubit] = op
                # jump to next iteration of qubit
                break
                

    # Construct the generator
    generator = ''
    for op in generator_terms:
        generator += op
 
    
    return SparsePauliOp(data=[generator],coeffs=[1.0])

def pauli_exp_val(pauli: str, counts: dict):
    # indices to check for the Pauli operators
    indices = [i for i in range(len(pauli)) if pauli[i] != 'I']
    # Initialize the expectation value
    exp = 0
    shots = sum(counts.values())
    for binary in counts.keys():
        sign = 1
        for i in indices:
            if binary[i] == '1':
                sign *= -1
        exp += sign * counts[binary]/shots
    return exp, 1-exp**2   

def eval_commuting_group(group: SparsePauliOp, counts:dict):
    """Evaluates the expectation value of a commuting group
        Input:  group: SparsePauliOp                
                counts: dict
        Returns: expectation: float
                 variances: list"""
    
    # Get the Pauli operators
    paulis = group.paulis
    coeffs = group.coeffs
    
    # Initialize the expectation value
    expectation = 0
    variances = []
    # Iterate over the Pauli operators
    for i in range(len(paulis)):
        # Get the Pauli operator and its coefficient
        pauli = paulis[i]
        coeff = coeffs[i]
        # Compute the expectation value
        pauli_exp, pauli_var = pauli_exp_val(pauli, counts)
        expectation += coeff * pauli_exp
        variances.append(pauli_var)
    # Return the expectation value
    return expectation,variances
    
def apply_observable(circuit: QuantumCircuit,observable: str):
    """Applies the given observable to the circuit
    
    Args: circuit: (QuantumCircuit) the circuit to apply the observable to
          observable: (SparsePauliOp) the observable to apply
    Returns: circuit: (QuantumCircuit) the circuit with the observable applied"""
    
    # Get the list of Pauli operators
    paulis = str(observable)
    qubit = 0
    # Apply each Pauli operator to the circuit
    for pauli in paulis:
    
        if pauli == 'X':
            circuit.h(qubit)
        elif pauli == 'Y':
            circuit.sdg(qubit)
            circuit.h(qubit)
        elif pauli == 'Z' or pauli == 'I':
            pass
        else:
            raise ValueError("Wrong operator")
        qubit +=1
    return circuit 

def remove_idle_qwires(circ):


    """Remove the idle qubits from the circuit
        This functions is copied from stackexhange page: https://quantumcomputing.stackexchange.com/questions/25672/remove-inactive-qubits-from-qiskit-circuit
        
        Input: circ: QuantumCircuit
        Returns: QuantumCircuit without idle qubits"""
    

    dag = circuit_to_dag(circ)

    idle_wires = list(dag.idle_wires())
    for w in idle_wires:
        dag._remove_idle_wire(w)
        dag.qubits.remove(w)

    dag.qregs = OrderedDict()

    return dag_to_circuit(dag)


#############################################################


def retrieve_neven_mapper(filename):
    """Retrieve the Neven mapper from a file in ../hamiltoninas/name.txt

    input: filename (str) - the name of the file containing the Neven mapping
    
    Return: SparsePauliOp - the Neven (Ternary tree) mapping"""

    with open(filename,'r') as file:
        lines = file.readlines()

   
    paulis = []
    coeffs = []
   
   
    for pauli in lines:
        # The hamiltonian is formatted as 
        # -float * String   i.e. -0.5 * ZZIXY
        coeffs.append(float(pauli.split('*')[0].strip(' ')))
        paulis.append(pauli.split('*')[1].strip())


    hamiltonian = SparsePauliOp(data=paulis,coeffs=coeffs)

    return hamiltonian

def prepare_hamiltonian(
    molecule_name, z2symmetry_reduction, mapping, basis='sto-3g'
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


    # Run the driver
    total_hamiltonian = driver.run()


    # Get the information about the problem
    total__num_particles = total_hamiltonian.num_particles # Returns (alpha,beta)
    total_num_spatial_orbitals = total_hamiltonian.num_spatial_orbitals
    total_active_orbitals = total_hamiltonian.orbital_energies


    # Apply the freeze core transformation
    transformer = FreezeCoreTransformer(freeze_core=True)
    #transformer.prepare_active_space(molecule,total_number_of_orbitals)
    reduced_hamiltonian = transformer.transform(total_hamiltonian)
    

    # Convert the Hamiltonian to second quantized form
    hamiltonian = reduced_hamiltonian.hamiltonian.second_q_op()


    num_spatial_orbitals = reduced_hamiltonian.num_spatial_orbitals
    num_particles = reduced_hamiltonian.num_particles


    if mapping == "parity":
        qubit_mapping = ParityMapper(num_particles=num_particles)
    elif mapping == "bravyi_kitaev":
        qubit_mapping = BravyiKitaevMapper()
    elif mapping == "neven":
    #    filename = molecule_name +'-neven-'+str(z2symmetry_reduction)+ '.txt'
    #    qubitOp = retrieve_neven_mapper(path_to_this_file+'/../hamiltonians/'+filename)
    #   if z2symmetry_reduction and molecule_name != 'H2':
    #        reduction = 2
    #    else:
    #        reduction = 0
    #    return qubitOp,num_spatial_orbitals-reduction,num_particles
        qubit_mapping = OptimalGeneralTernaryTreeMapper()
    
    elif mapping == "jordan_wigner":
        qubit_mapping = JordanWignerMapper()
    else:
        raise ValueError("Wrong mapping")


    # Apply the given fermion to Pauli encoding
    qubitOp = QubitMapper.map(self=qubit_mapping,
        second_q_ops=hamiltonian
        )
    
    reduction = 0
    # Apply Z2 symmetries (except for hydrogen)
    if molecule_name != 'H2':
        if z2symmetry_reduction:
            z2symmetries = Z2Symmetries.find_z2_symmetries(qubitOp)
            qubitOp = z2symmetries.taper(qubitOp)
            reduction = 2


    # Z2 symmetries need some specific handling
    if molecule_name == 'H2':
        result = qubitOp,num_spatial_orbitals-reduction,num_particles
    elif z2symmetry_reduction:
        result = qubitOp[0],num_spatial_orbitals-reduction,num_particles
    else:
        result = qubitOp,num_spatial_orbitals-reduction,num_particles
    return result

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

def build_ansatz(ansatz_name,mapper,num_qubits: int,num_particles: tuple[int,int],num_spatial_orbitals: int,reps: int):
    """Create hartreeFock anstaz with selected parametrization. Neven mapper is handled as a special case.

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
    elif mapper == 'neven':
        print('I hope you are not using kUPCCGSD with Neven mapper')
        qubit_mapping = JordanWignerMapper()
    else:
        raise ValueError('Invalid mapper name')

    # create HartreeFock state
    hf_state = HartreeFock(num_spatial_orbitals=num_spatial_orbitals, num_particles=num_particles, qubit_mapper=qubit_mapping)

    if ansatz_name=='kUpCCGSD':
        ansatz = PUCCSD(num_spatial_orbitals=num_spatial_orbitals, num_particles=num_particles, qubit_mapper=qubit_mapping,reps=reps)
    elif ansatz_name=='EfficientSU2':
        ansatz = EfficientSU2(num_qubits=num_qubits, entanglement='linear', reps=reps)
    else: 
        print('You should not be able to get this message')
        raise ValueError('Invalid ansatz name')

    #return hf_state.compose(ansatz)
    return ansatz


    #optimize the ansatz for backend

def vqe(num_qubits,ansatz,hamiltonian,grouping,shots=1_000):
    
    
    """Run the VQE optimization algorithm. This function uses the COBYLA optimizer from SciPy. 

    Args: num_qubits (int): number of qubits
          ansatz (QuantumCircuit): ansatz circuit
          hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
    Returns: cost_history_dict (dict): dictionary containing the cost and parameter history of the optimization"""
    

    # Get the backend
    backend = GenericBackendV2(ansatz.num_qubits)
    backend.set_options(shots=shots)
    # https://docs.quantum.ibm.com/api/qiskit/qiskit.providers.fake_provider.GenericBackendV2

    # Get the number of parameters in the ansatz
    num_params = ansatz.num_parameters
    
    # Optimize the ansatz for the backend
    target = backend.target
    pm = generate_preset_pass_manager(target=target, optimization_level=3)
    
    
    # Remove idle wires from the ansatz
    ansatz_no_idles = remove_idle_qwires(ansatz)


    # Add here the evaluation of the number of shots


    # Optimize the circuit
    ansatz_isa = pm.run(ansatz_no_idles)
    hamiltonian_isa = hamiltonian.apply_layout(layout=ansatz_isa.layout)


    cost_history_dict = {
        "prev_vector": None,
        "iters": 0,
        "cost_history": [],
        "vectors": [],
        "vqe_acc":[]
    }


    #x0 = 2 * np.pi * np.random.random(num_params)
    x0 = np.zeros(num_params)


    estimator = BackendEstimator(backend,abelian_grouping=grouping)
    #estimator.options.default_shots = shots
    # set the maximum number of iterations
    options = {'maxiter': max_iter}

    def cost_func(params):
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
        result = estimator.run(ansatz,hamiltonian,params,shots=shots).result()
        energy = result.values[0]
        vars = result.metadata[0]['variance']
        acc = accuracy(weights=hamiltonian.coeffs,variances=vars,shots=shots)


        cost_history_dict["iters"] += 1
        cost_history_dict["prev_vector"] = params
        cost_history_dict["cost_history"].append(energy)
        cost_history_dict["vectors"].append(params.tolist())
        cost_history_dict["vqe_acc"].append(acc)
        #print(f"Iters. done: {cost_history_dict['iters']} [Current cost: {energy}]")

        return energy


    # Perform the minimaxation routine
    # COBYLA
    if optimizer == 'cobyla':
        res = minimize(
            cost_func,
            x0,
            method='cobyla',
            options=options,
        )

    # SPSA
    else:
        res = spsa.minimize(
            cost_func,
            x0,
            iterations=max_iter,
            )
    print(cost_history_dict['iters'])

    return cost_history_dict

def command_line_parser(arguments):


    """Parse the command line arguments and return the correct arguments for the VQE calculation
    Args: arguments (list): list of command line arguments without the first element
    Returns: mappings, ansatzes, molecules, measurement_schemes, z2_symmetry (list): list of the arguments for the VQE calculation"""

    try:
        arguments[0]    
    except:
        print(arg_error_message)
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
            filename = 'all_methods'
            ansatz_reps = all_anzats_reps
    

        # Run the default methods
        else:
            print('Running the comparison of methods respect to default methods')
            mappings = default_mapping
            ansatzes = default_ansatz
            molecules = all_molecule_names
            measurement_schemes = default_measurement_scheme
            z2_symmetry = default_z2
            VQE = 'Y'
            filename = 'default_methods'
            ansatz_reps = default_ansatz_reps


    # create a check thst the arguments are correct
    elif len(arguments) != 8: 
        print(arg_error_message)
        sys.exit()
    else: 
        print('Running a custom set of arguments')
        filename = arguments[1]
        VQE = arguments[2]
        molecules = arguments[3].split(',')
        mappings = arguments[4].split(',')
        ansatzes = arguments[5].split(',')    
        measurement_schemes = arguments[6].split(',')
        # map the z2_symmetry to a boolean
        z2_s = arguments[7].split(',') 
        z2_symmetry=[]
        ansatz_reps = [1,2,3]
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
        print('Results will be saved in:',filename+'.csv')


    
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


    return filename,VQE,molecules, mappings, ansatzes, measurement_schemes, z2_symmetry,ansatz_reps
  
def remove_idle_qwires(circ):


    """Remove the idle qubits from the circuit
        This functions is copied from stackexhange page: https://quantumcomputing.stackexchange.com/questions/25672/remove-inactive-qubits-from-qiskit-circuit
        
        Input: circ: QuantumCircuit
        Returns: QuantumCircuit without idle qubits"""
    

    dag = circuit_to_dag(circ)

    idle_wires = list(dag.idle_wires())
    for w in idle_wires:
        dag._remove_idle_wire(w)
        dag.qubits.remove(w)

    dag.qregs = OrderedDict()

    return dag_to_circuit(dag)

def required_shots(weight,variation,accuracy=0.001):
    weight = np.array(weight)
    variation = np.array(variation)
    return (np.abs(weight)**2 *variation) / (accuracy**2) 

# Functions to calculate the accuracy of the VQE calculation
def accuracy(weights,variances,shots):
    return np.sqrt(sum(np.abs(weights)**2*variances/shots))

def pauli_splitter(pauli_op):
    """Split a Pauli operator into its individual terms.
    Args:
        pauli_op (PauliOp): The Pauli operator to split.
    Returns:
        list: The list of individual Pauli terms.
    """
    paulis = pauli_op.paulis
    coeffs = pauli_op.coeffs
    terms = []
    for pauli,coeff in zip(paulis,coeffs):
        terms.append(SparsePauliOp(data=[pauli],coeffs=[coeff]))
    return terms

def variance(hamiltonian,shots,ansatz):
    vars = []
    splitted = pauli_splitter(hamiltonian)
    num_qubits = len(hamiltonian.paulis[0])
    backend_estimator = BackendEstimator(backend=GenericBackendV2(num_qubits),abelian_grouping=False)
    #backend_estimator = Estimator()
    parameters = [0.001] * ansatz.num_parameters # initial parameters, whit very close to zero angles
    for obs in splitted:
        job = backend_estimator.run(circuits=ansatz,observables=obs,parameter_values=parameters,shots=shots)
        vars.append(job.result().metadata[0]['variance'])

    return vars

def QWC_variance(hamiltonian,shots,ansatz):
    vars = []
    num_qubits = len(hamiltonian.paulis[0])
    backend_estimator = BackendEstimator(backend=GenericBackendV2(num_qubits),abelian_grouping=True)
    #backend_estimator = Estimator()
    parameters = [0.001] * ansatz.num_parameters # initial parameters, whit very close to zero angles
    job = backend_estimator.run(circuits=ansatz,observables=hamiltonian,parameter_values=parameters,shots=shots)
    vars.append(job.result().metadata[0]['variance'])

    return vars


def calc_accuracy(hamiltonian,shots,ansatz,grouping):
    if  grouping ==True:
        vars = QWC_variance(hamiltonian,shots,ansatz)
    else:
        vars = variance(hamiltonian,shots,ansatz)
    return accuracy(weights=hamiltonian.coeffs,variances=vars,shots=shots)

def approx_accracy(hamiltonian,shots,vars):
    return accuracy(weights=hamiltonian.coeffs,variances=vars,shots=shots)

def save_hamiltonian(hamiltonian,filename):
    with open(filename,'w') as file:
        for i in range(len(hamiltonian.paulis)):
            file.write(str(hamiltonian.coeffs[i]) + ' * ' + str(hamiltonian.paulis[i]) + '\n')
        file.close()
    return ''

def apply_observable(circuit,observable):
    """Applies the given observable to the circuit
    
    Args: circuit: (QuantumCircuit) the circuit to apply the observable to
          observable: (SparsePauliOp) the observable to apply
    Returns: circuit: (QuantumCircuit) the circuit with the observable applied"""
    
    # Get the list of Pauli operators
    paulis = str(observable)
    qubit = 0
    # Apply each Pauli operator to the circuit
    for pauli in paulis:
    
        if pauli == 'X':
            circuit.h(qubit)
        elif pauli == 'Y':
            circuit.sdg(qubit)
            circuit.h(qubit)
        elif pauli == 'Z' or pauli == 'I':
            pass
        else:
            raise ValueError("Wrong operator")
        qubit +=1
    return circuit