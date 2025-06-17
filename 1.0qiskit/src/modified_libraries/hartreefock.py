# This code is part of a Qiskit project.
#
# (C) Copyright IBM 2018, 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Hartree-Fock initial state."""

from __future__ import annotations
import numpy as np
from qiskit import QuantumRegister
from qiskit.circuit.library import BlueprintCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.utils.validation import validate_min

from qiskit_nature.second_q.operators import FermionicOp

import subprocess


def hartree_fock_bitstring_mapped(
    num_spatial_orbitals: int,
    num_particles: tuple[int, int],
) -> list[bool]:
    # pylint: disable=unused-argument
    """Compute the bitstring representing the mapped Hartree-Fock state for the specified system.

    Args:
        num_spatial_orbitals: The number of spatial orbitals, has a min. value of 1.
        num_particles: The number of particles as a tuple (alpha, beta) containing the number of
            alpha- and  beta-spin electrons, respectively.
        qubit_mapper: A QubitMapper.

    Returns:
        The bitstring representing the mapped state of the Hartree-Fock state as array of bools.
    """
    # get the bitstring encoding the Hartree Fock state
    bitstr = hartree_fock_bitstring(num_spatial_orbitals, num_particles)

    # encode the bitstring as a `FermionicOp`
    bitstr_op = FermionicOp(
        {" ".join(f"+_{idx}" for idx, bit in enumerate(bitstr) if bit): 1.0},
        num_spin_orbitals=2 * num_spatial_orbitals,
    )

    # Write the bitstring to a file
    with open("bitstring.txt", "w") as f:
        f.write(str(bitstr_op))
     
    # map the `FermionicOp` to a qubit operator
    qubit_op: SparsePauliOp


    # We need to use different conda enviroment to perform the mapping
    script_name = "neven_mapper.py"
    env_name = 'tomo'

    command = f"conda run -n {env_name} python {script_name} {bitstr_op}"

    output = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    print(output.stdout)



    #qubit_op = qubit_mapper.map(bitstr_op)

    # We check the mapped operator `x` part of the paulis because we want to have particles
    # i.e. True, where the initial state introduced a creation (`+`) operator.
    bits = []
    #for bit in qubit_op.paulis.x[0]:
    #    bits.append(bit)

    return bits


def hartree_fock_bitstring(num_spatial_orbitals: int, num_particles: tuple[int, int]) -> list[bool]:
    """Compute the bitstring representing the Hartree-Fock state for the specified system.

    Args:
        num_spatial_orbitals: The number of spatial orbitals, has a min. value of 1.
        num_particles: The number of particles as a tuple storing the number of alpha- and beta-spin
                       electrons in the first and second number, respectively.

    Returns:
        The bitstring representing the state of the Hartree-Fock state as array of bools.

    Raises:
        ValueError: If the total number of particles is larger than the number of orbitals.
    """
    # validate the input
    validate_min("num_spatial_orbitals", num_spatial_orbitals, 1)
    num_alpha, num_beta = num_particles

    if any(n > num_spatial_orbitals for n in num_particles):
        raise ValueError("# of particles must be less than or equal to # of orbitals.")

    half_orbitals = num_spatial_orbitals
    bitstr = np.zeros(2 * num_spatial_orbitals, bool)
    bitstr[:num_alpha] = True
    bitstr[half_orbitals : (half_orbitals + num_beta)] = True

    return bitstr.tolist()


hartree_fock_bitstring_mapped(2, (1, 1))