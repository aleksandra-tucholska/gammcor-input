# Quantum Chemistry Input Parser

This README describes the input format and key data structures used in the quantum chemistry calculations input parser.

## Table of Contents
1. [Input File Format](#input-file-format)
2. [Key Data Structures](#key-data-structures)
   - [molecule_data](#molecule_data)
   - [basis_set_data](#basis_set_data)
   - [project_paths](#project_paths)
3. [Main Functionality](#main-functionality)

## Input File Format

The input file should contain the following sections:

### General Settings
```
basis <basis_set_name>
charge <charge>
units <unit_type>
interface <interface_type>
my_template_path <path_to_custom_template>  # Optional
```

### Geometry
```
xyz
<atom> <x> <y> <z>
...
end
```

### Interface-Specific Block (e.g., DALTON)
```
dalton
nsym <number_of_symmetries>
symmetry <symmetry>
multip <multiplicity>
inactive <inactive_orbitals>
electrons <number_of_electrons>
cas <complete_active_space>
state <state>
end
```

### Interface Executable Path
```
dalton_path <path_to_dalton_executable>
molpro_path <path_to_molpro_executable>
orca_path <path_to_orca_executable>
```

## Key Data Structures

### molecule_data

`molecule_data` is a dictionary containing information about the molecular structure and calculation parameters:

```python
molecule_data = {
    'basis': <basis_set_name>,
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

#### Detailed description of molecule_data fields:
- `basis`: String, name of the basis set used for the calculation.
- `interface`: String, type of interface used (e.g., 'DALTON', 'MOLPRO', 'ORCA').
- `units`: String, units used for coordinates (e.g., 'ANGS' for Angstroms, 'BOHR' for Bohr).
- `geometry`: List of dictionaries, each containing information about an atom in the system:
  - `atom`: String, chemical symbol of the atom (e.g., 'H', 'C', 'O').
  - `atomic_number`: Integer, atomic number of the element.
  - `coordinates`: List of three floats, representing the x, y, and z coordinates of the atom.
- `charge`: Integer or float, total charge of the system.
- `atom_list`: List of strings, chemical symbols of all atoms in the order they appear in the geometry.

### basis_set_data

`basis_set_data` is a dictionary containing detailed information about the basis set for each element:

```python
basis_set_data = {
    'element_symbol': {
        'angmoms': ['S', 'P', 'D', ...],
        'S': {
            'num_prim': <number_of_primitives>,
            'num_contr': <number_of_contractions>,
            'exp': [<exponents>],
            'contractions': [<contraction_coefficients>],
            'idx': [<indices>]
        },
        'P': {
            # Similar structure to 'S'
        },
        # Other angular momentum blocks...
    },
    # Other elements...
}
```

#### Detailed description of basis_set_data fields:
- `element_symbol`: String, chemical symbol of the element (e.g., 'H', 'C', 'O').
- `angmoms`: List of strings, angular momentum blocks present for this element (e.g., ['S', 'P', 'D']).
- For each angular momentum block ('S', 'P', 'D', etc.):
  - `num_prim`: Integer, total number of primitive functions.
  - `num_contr`: Integer, number of contracted functions.
  - `exp`: List of floats, exponents for primitive functions.
  - `contractions`: List of lists of floats, contraction coefficients for each contracted function.
  - `idx`: List of lists of integers, indices mapping primitives to contractions for each contracted function.

### project_paths

`project_paths` is a dictionary containing various file and directory paths relevant to the project:

```python
project_paths = {
    'basis_path': <path_to_basis_file>,
    'interface': <interface_type>,
    'gammcor_path': <path_to_gammcor>,
    'dalton_path': <path_to_dalton_executable>,
    'molpro_path': <path_to_molpro_executable>,
    'orca_path': <path_to_orca_executable>,
    'input_file': <path_to_input_file>,
    'my_template_path': <path_to_custom_template>,
    'input_name': <input_file_name_without_extension>,
    'input_dir': <input_file_directory>,
    'output_dir': <output_directory>,
    'dalton_dir': <dalton_files_directory>,
    '.dal': <path_to_dalton_input_file>,
    '.mol': <path_to_dalton_molecule_file>
}
```

#### Detailed description of project_paths fields:
- `basis_path`: String, path to the basis set definition file.
- `interface`: String, type of interface used (e.g., 'DALTON', 'MOLPRO', 'ORCA').
- `gammcor_path`: String, path to the GAMMCOR installation.
- `dalton_path`: String, path to the DALTON executable.
- `molpro_path`: String, path to the MOLPRO executable.
- `orca_path`: String, path to the ORCA executable.
- `input_file`: String, path to the main input file.
- `my_template_path`: String or False, path to a custom template file (if provided).
- `input_name`: String, name of the input file without extension.
- `input_dir`: String, directory containing the input file.
- `output_dir`: String, directory where output files will be saved.
- `dalton_dir`: String, directory for DALTON-specific files.
- `.dal`: String, path to the generated DALTON input file.
- `.mol`: String, path to the generated DALTON molecule file.

## Main Functionality

1. Parse input file and initialize data structures.
2. Create output directory structure.
3. Parse basis set data.
4. Generate interface-specific input files (e.g., DALTON).
5. Run interface-specific calculations.
6. Run GAMMCOR calculations.

The script supports multiple interfaces (DALTON, MOLPRO, ORCA), but currently, only DALTON is fully implemented.

This README provides an overview of the input file structure and the key data structures used in the quantum chemistry calculations input parser. 
