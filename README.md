## Data Structures

### elements_data Structure

The `elements_data` object is a dictionary containing information about each element in the calculation:

```python
elements_data = {
    'element_symbol': {
        'blocks': ['S', 'P', 'D', ...],  # List of angular momentum blocks
        'S': {  # Structure for S orbital
            'num_prim': <number_of_primitives>,
            'num_contr': <number_of_contractions>,
            'exp': [<exponents>],  # List of exponents
            'contractions': [<contraction_coefficients>],  # List of contraction coefficients
            'idx': [<indices>]  # List of indices
        },
        'P': {
            # Similar structure to 'S'
        },
        # Other angular momentum blocks (D, F, etc.)
    },
    # Other elements follow the same structure
}
```

#### Explanation of elements_data fields:
- `element_symbol`: The chemical symbol of the element (e.g., 'H', 'C', 'O')
- `blocks`: List of angular momentum blocks present for this element
- For each block ('S', 'P', 'D', etc.):
  - `num_prim`: Total number of primitive functions
  - `num_contr`: Number of contracted functions
  - `exp`: List of exponents for primitive functions
  - `contractions`: List of contraction coefficients
  - `idx`: List of indices mapping primitives to contractions

### data Structure

The `data` object is a dictionary containing general information about the calculation:

```python
data = {
    'basis': <basis_set_name>,
    'basis_path': <path_to_basis_file>,
    'gammcor_path': <path_to_gammcor>,
    'interface': <interface_type>,
    'units': <unit_type>,
    'geometry': [
        {
            'atom': <atom_symbol>,
            'atomic_number': <atomic_number>,
            'coordinates': [<x>, <y>, <z>]
        },
        # Additional atoms...
    ],
    'charge': <system_charge>,
    'atom_list': [<list_of_atom_symbols>]
}
```

#### Explanation of data fields:
- `basis`: Name of the basis set
- `basis_path`: Path to the basis set
- `gammcor_path`: Path to the GAMMCOR executable
- `interface`: Type of interface used (e.g., 'DALTON', 'MOLPRO')
- `units`: Units used for coordinates (e.g., any string starting from a or A for Angstroms, any string starting from b or B for Bohr)
- `geometry`: List of dictionaries, each containing information about an atom in the system
- `charge`: Total charge of the system
- `atom_list`: List of all atom symbols in the order they appear in the geometry

#### Detailed structure of 'geometry':
The `geometry` field is a list of dictionaries, where each dictionary represents an atom in the molecular structure. Here's a more detailed look at its structure:

```python
'geometry': [
    {
        'atom': 'C',            # Atom symbol (string)
        'atomic_number': 6,     # Atomic number (integer)
        'coordinates': [0.0, 0.0, 1.89]  # [x, y, z] coordinates (list of floats)
    },
    {
        'atom': 'H',
        'atomic_number': 1,
        'coordinates': [0.0, 1.63, 2.49]
    },
    # ... more atoms ...
]
```

Each dictionary in the `geometry` list contains:
- `atom`: A string representing the chemical symbol of the atom (e.g., 'C' for carbon, 'H' for hydrogen).
- `atomic_number`: An integer representing the atomic number of the element.
- `coordinates`: A list of three float values representing the x, y, and z coordinates of the atom in the specified units (see the `units` field in the main `data` structure).

The order of atoms in this list is significant and corresponds to the order in which they were specified in the input file. This order is also reflected in the `atom_list` field of the main `data` structure.

