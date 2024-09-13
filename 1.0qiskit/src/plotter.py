################### Imports ###################


# Ignore warnings
import warnings
warnings.filterwarnings("ignore")


# General imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math


# Import curve fitting from scipy
from scipy.optimize import curve_fit


# Sympy imports for printing
from sympy import symbols, log
from sympy.printing.latex import latex


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

def BK_scaling(num_qubits,a=0,b=1,c=1):
    
        """Calculates the scaling of the Bravyi-Kitaev mapping
            O(log2(N)) = a +b*log2(c*num_qubits)
        
        Input: num_qubits (int) - the number of qubits
        
        Output: scaling (int) - the scaling of the Bravyi-Kitaev mapping"""

        num_qubits = np.array(num_qubits)

        return [float(a +b*np.log2(c*n)) for n in num_qubits]

def PR_scaling(num_qubits):
        
    """Calculates the scaling of the Parity mapping
            
            Input: num_qubits (int) - the number of qubits
            
            Output: scaling (int) - the scaling of the Parity mapping"""
    
    num_qubits = np.array(num_qubits)
    
    return num_qubits

def OTT_scaling(num_qubits,a=0,b=1,c=1):
            
    """Calculates the scaling of the Parity mapping
                
                Input: num_qubits (int) - the number of qubits
                
                Output: scaling (int) - the scaling of the OTT mapping"""
        
    num_qubits = np.array(num_qubits)
    args = c*2*num_qubits
        
    return [a+b*np.log(arg)/np.log(3) for arg in args]


################### Main ###################

# plot the data
if __name__ == '__main__':


    figsize = (15,10)

    # read in the data
    filename = 'vqe_results.csv'
    data = pd.read_csv('../results/'+filename)

    
    # Define the methods
    mappings = ['parity', 'jordan_wigner', 'bravyi_kitaev', 'neven']
    
    pauli_plots = input("Plot pauli results?: (y/n)")
    
    if(pauli_plots == 'y'):


        # reformat the data
        dropping = ['ansatz', 'ansatz_circuit',
            'hamiltonian','vqe_energies', 'iterations',
            'exact_energies', 'exact_solution', 'parameters', 'error']
        data = data.drop(columns=dropping)


        # Create a four figure sublot for the results
        fig, axs = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Pauli weights for different mappings')
        axs = axs.ravel()


        format = 'png'


        # Plot the the Pauli properties
        for ax,map in zip(axs.flat,mappings):
            print('\n'+map+'\n')


            # Filter the data
            if True:
                map_data = data[data['mapping']==map]
                num_qubits = map_data['num_qubits']
                avg_pauli_weight = map_data['avg_pauli_weight']
                avg_hardware_pauli_weight = map_data['avg_hardware_pauli_weight']
                num_pauli_strings = map_data['num_pauli_strings']
                max_pauli_weights = map_data['max_pauli_weight']
                max_hrdwwr_pauli_weights = map_data['max_hrdwr_pauli_weight']


                # Create the range for the scaling
                plot_range = np.arange(2,(max(num_qubits)+4))

                N = symbols('N')

                if map == 'parity':
                    scaling = PR_scaling(plot_range)
                    name = 'Parity'
                    scaling_name='O(N)'
                elif map == 'jordan_wigner':
                    scaling = JW_scaling(plot_range)
                    name = 'Jordan-Wigner'
                    scaling_name='O(N)'
                elif map == 'bravyi_kitaev':
                    scaling = BK_scaling(plot_range)
                    name = 'Bravyi-Kitaev'
                    scaling_name='log2(N)'

                    init_guess = [1,1,1]
                    params, covariance = curve_fit(BK_scaling,num_qubits,max_pauli_weights, p0=init_guess)
                    a,b,c = params
                    fitted_scaling = BK_scaling(plot_range,a,b,c)
                    ax.plot(plot_range, fitted_scaling, label=f"{a:.2f}+{b:.2f}log2({c:.2f}N)")


                elif map == 'neven':
                    scaling = OTT_scaling(plot_range)
                    name = 'Ternary Tree'
                    scaling_name='O(log₃(2N))'
                    init_guess = [1,1,1]
                    params, covariance = curve_fit(OTT_scaling,num_qubits,max_pauli_weights, p0=init_guess)
                    a,b,c = params
                    fitted_scaling = OTT_scaling(plot_range,a,b,c)
                    ax.plot(plot_range, fitted_scaling, label=f"{a:.2f}+{b:.2f}log3({2*c:.2f}N)")

               

                    


                # Create scatter plot, since avg_pauli_weight can have multiple values for the same x values
                ax.plot(num_qubits, avg_pauli_weight, 'o', label='Average Pauli Weight')
                ax.plot(num_qubits, max_pauli_weights, 'ko', label='Max Pauli Weights')
                ax.plot(plot_range, scaling, label=scaling_name)
                ax.set_xlabel('Number of qubits')
                ax.set_ylabel('Weight')
                ax.legend(loc='upper left')
       

                ax.set_ylim(0, 25)
                
       
                # Add the title and grid
                ax.set_title(name)
                ax.grid(True, which='both', linestyle='--', color='grey')
                


                # Plot the hardware weight and pauli weight 
                #plt.plot(num_qubits, avg_hardware_pauli_weight, 'o', label='Rz based hadrware gates')
                #plt.plot(num_qubits, avg_pauli_weight, 'o', label='Ideal qubit gates')
                #plt.plot(plot_range, scaling, label='Theoretical scaling')
                #plt.legend()
                #plt.xlabel('Number of qubits')
                #plt.ylabel('Number of qubit operations')
                #plt.title('Number of gates to measure Pauli string'+name)
                #plt.grid('both',linestyle='--')
                #plt.savefig('../results/Hardware_Pauli_weight_'+map+'.'+format, format=format, dpi=1000)
                #plt.show()
            
            #try:
            #    print('moi')

            #except:
                # Jump to next mapping
            #    print('No data or error occured for '+map+' mapping')
            #    continue
        plt.tight_layout()
        plt.savefig('../results/Pauli_weights.'+format, format=format, bbox_inches = 'tight',dpi=1000)
        plt.show()


    # Plot the VQE results
    
    # read in the data (again)
    data = pd.read_csv('../results/'+filename)
    
    
    # reformat the data
    dropping = ['avg_pauli_weight','avg_hardware_pauli_weight'
                 ,'parameters','ansatz_circuit','hamiltonian'
                 ,'num_pauli_strings'  
    ]
    
     
    data = data.drop(columns=dropping)
    plot_vqe = input('Plot vqe results?: (y/n)')
    
    
    if(plot_vqe  == 'y'):
        
        print(data.head())
        
        for map in mappings:
            
         # Plot the the Pauli properties
            print('\n'+map+'\n')

            # Filter the data
            try:
                map_data = data[data['mapping']==map]
                num_qubits = map_data['num_qubits']
                xxx = map_data['avg_pauli_weight']
                avg_hardware_pauli_weight = map_data['avg_hardware_pauli_weight']
                num_pauli_strings = map_data['num_pauli_strings']
                max_pauli_weights = map_data['max_pauli_weight']
                max_hrdwwr_pauli_weights = map_data['max_hrdwr_pauli_weight']

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


                format = 'png'

                # Create scatter plot, since avg_pauli_weight can have multiple values for the same x values
                plt.plot(num_qubits, avg_pauli_weight, 'o', label='Average Pauli Weight')
                plt.plot(num_qubits, max_pauli_weights, 'ko', label='Max Pauli Weights')
                plt.plot(plot_range, scaling, label='Theoretical scaling')
                plt.legend()
                plt.xlabel('Number of qubits')
                plt.ylabel('Value')
                plt.title('Pauli weights for '+name+' mapping')
                plt.grid('both',linestyle='--')
                plt.savefig('../results/Pauli_weight_'+map+'.'+format, format=format, dpi=1000)
                plt.show()


                # Plot the hardware weight and pauli weight 
                plt.plot(num_qubits, avg_hardware_pauli_weight, 'o', label='Rz based hadrware gates')
                plt.plot(num_qubits, avg_pauli_weight, 'o', label='Ideal qubit gates')
                plt.plot(plot_range, scaling, label='Theoretical scaling')
                plt.legend()
                plt.xlabel('Number of qubits')
                plt.ylabel('Number of qubit operations')
                plt.title('Number of gates to measure Pauli string'+name)
                plt.grid('both',linestyle='--')
                plt.savefig('../results/Hardware_Pauli_weight_'+map+'.'+format, format=format, dpi=1000)
                plt.show()


            
            except:
                # Jump to next mapping
                continue

if False:
    

    # Generate data for the plot
    x = np.linspace(0, 2 * np.pi, 100)
    y = np.sin(x)


    # Random scatter data
    scatter_x = np.random.uniform(0, 2 * np.pi, 30)
    scatter_y = np.sin(scatter_x) + np.random.normal(0, 0.1, 30)


    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    names = ['Parity', 'Jordan-Wigner', 'Bravyi-Kitaev', 'Neven']

    # Iterate over all the axes and customize each plot
    for ax,name in zip(axs.flat,names):
        # Plot the continuous curve (sin(x))
        ax.plot(x, y, label='Sin(x)')
    
        # Plot scattered data points
        ax.scatter(scatter_x, scatter_y, color='red', label='Data points')
    
        # Set title and grid
        ax.set_title(name)
        ax.grid(True, which='both', linestyle='--', color='grey')


    # Adjust layout
    plt.tight_layout()


    # Show the plot
    plt.show()
