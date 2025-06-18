from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo


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
        (0.0, 0.0, 0.0),           # First Hydrogen
        (0.0, 3.076, 0.0),         # Second Hydrogen
        (2.664, 1.538, 0.0)        # Third Hydrogen
        ],
        charge=0,
        multiplicity=1,
    ),
    #"He": MoleculeInfo(
    #    symbols=["He"],
    #    coords=[(0.0,0.0,0.0)],
    #    charge=0,
    #    multiplicity=1,
    #),
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
    #    coords=[(-1.082,-0.665,0.0),(0.0,0.0,0.0),(1.082,-0.665,0.0)],
    #    charge=1,
    #    multiplicity=1,
    #),
    'O2': MoleculeInfo(
        symbols=["O","O"],
        coords=[(0.0,0.0,0.0),(0.0,0.0,1.2)],
        charge=0,
        multiplicity=1,
    )
}