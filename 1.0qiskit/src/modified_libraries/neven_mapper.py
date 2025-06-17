# This script needs to be run in "tomo" environment

from qiskit.chemistry.transformations import (
    FermionicTransformation,
    FermionicTransformationType,
    FermionicQubitMappingType,
)
from qiskit.chemistry import FermionicOperator


from qiskit.quantum_info import Pauli
import numpy as np


################
qubit_mapping = FermionicQubitMappingType.NEVEN


fermionic_transformation = FermionicTransformation(
        transformation=FermionicTransformationType.FULL,
        qubit_mapping=qubit_mapping,
        two_qubit_reduction=False
)

#operator = FermionicOperator("1.0 [0^2^]")
#qubitOp, _ = fermionic_transformation.transform(aux_operators=operator)

#print(qubitOp)



# Read in the bitstring from bitstring.txt
#with open("bitstring.txt", "r") as f:


#    bitstr = f.read()

bitstr = "1.0 * ( +_0 +_2 )"

#state = FermionicOperator(bitstr)
#state.mapping('neven')

#qubitOp, _ = qubit_mapping.transform(bitstr)


#print(qubitOp)

def neven_mode(n):

    def nodeindex(p,l):
        """Returns index of a node of depth l on the path p on the ternary tree of height h.
            Ternary tree is just a string of natural numbers starting from 0 ordered in a tree where
            each node has three children. For image and formula see arXiv:1910.10746, Eq. (3).


            Args:
                p (list): list of strings of integers 0,1 or 2 of length h
                l (int): depth of tree to calculate the index on, l <= h

            Returns:
                int: index of a node corresponding to a qubit number
            """
        h = len(p)
        assert l <= h, "should be l <= h, where h = len(p)"
        for idx in p:
            assert (idx == '0') or (idx == '1') or (idx == '2'), "indices should be 0,1 or 2"
           
        prefactor = (3**l - 1)/2 
        sumfactors = [3**(l-1-j)*int(p[j]) for j in range(l)]
        return int(prefactor + sum(sumfactors))

    def hfromnum_qubits(num_qubits):
            """Calculates the needed height of the tree from number of qubits.

            Args:
                num_qubits (int): number of qubits
            
            Returns:
                height (int): required height of the ternary tree.
            """
            height = np.ceil(np.log(2*num_qubits+1)/np.log(3)) # base 3 logarithm

            return int(height) # integererize the output

    def xyzstring(h):
        """Generate a list of repeating 'X','Y','Z' pattern to fill the ternary tree.

            Args:
                h (int): Height of the ternary tree

            Returns:
                list: List of strings ['null', 'X','Y','Z','X','Y,'Z','X',...]
            """
        num_idxs = int((3**h - 1) // 2) # number of indices (qubits) to add
        num_idxs_triplets = int(num_idxs / 3) # number of triplets 
        output = ['null'] # add an index for the 0th qubit
        for _ in range(num_idxs_triplets):
            output += ['X','Y','Z']

        return output 


    def paulipaths_full(h):
        """Generate all Pauli paths from a tree of height h.
            
            
            Args:
                h (int): height of the ternary tree

            Returns:
                list: List of Pauli strings
            """


        xyzs = xyzstring(h+1) # generate the xyz string for tree height + 1
        num_qubits = int((3**(h) - 1) // 2) # get number of qubits from the height
        # generate all paths by looping over the number of paths in ternary base 
        paths = []
        for i in range(3**(h)):
            paths += [np.base_repr(i,base=3).rjust(h,'0')]
            
        # generate the Pauli strings from paths by substituting I with appropriate Pauli gate
        paulistrings = []
        for path in paths:
            pstring = ['I']*num_qubits # initialize a string with I's
            # for each depth, get the index at which the IIII.. path should be substituted and 
            # idx2 at which the substitution Pauli is located
            for depth in range(h):
                idx = nodeindex(path,depth)
                idx2 = nodeindex(path,depth+1)
                pstring[idx] = xyzs[idx2]
            # add the resulting string to the list, converted from string to list
            paulistrings += ["".join(pstring)]

        return paulistrings

    def paulipaths(num_qubits):
        h = hfromnum_qubits(num_qubits) # get tree max height
         
        # generate two full trees, the larger one should 
        # accommodate all qubits
        pphm1 = paulipaths_full(h-1)
        pph = paulipaths_full(h)
           
        # get number of qubits in the smaller tree 
        # and the number of extra qubits
        n_qubits_hm1 = len(pphm1[0])
        n_extraqubits = num_qubits - n_qubits_hm1


        paulistrings = []
            # loop over each gate in the larger tree for extra qubits 
            # to add the extra paths, also truncate the added paths 
            # up to the real number of qubits,
            # then add the path, converted to a Pauli object
        for i in range(n_extraqubits * 3):
            path = pph[i][:len(pphm1[0])+n_extraqubits]
            # print(path)
            paulistrings += [Pauli.from_label(path)]
                
            # paulistrings += path
            # loop over each gate in the smaller tree, skipping over 
            # the first gates of extra qubits, also add extra 'I' gates,
            # then add the path, converted to a Pauli object
        for i in range(n_extraqubits, len(pphm1)):
            path = ''.join([pphm1[i]] + ['I'] * n_extraqubits)
            paulistrings += [Pauli.from_label(path)]
            # paulistrings += path
        
        return paulistrings

    # reshape the pauli strings from a list into a paired list, skipping the very last one
    allpaths = paulipaths(num_qubits=n)
    paths = [[allpaths[2*k],allpaths[2*k + 1]] for k in range(int((len(allpaths)-1)//2))] 

    return paths 

result = neven_mode(4)

for res in result:
    for pauli in res:
        print(pauli.to_label())

