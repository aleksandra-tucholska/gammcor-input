#!/usr/bin/env python3
import sys
import os
import re
from atomic_number import elements_dict
from atomic_number import elements_short
from atomic_number import angular_momentum_dict
import shutil
from copy import deepcopy
from collections import defaultdict
import subprocess

molpro_path = "/molpro"
orca_path = "/orca"

debug = False

    
def error_msg(message):
    width = 80
    border_top = "╔" + "═" * (width - 2) + "╗"
    border_bottom = "╚" + "═" * (width - 2) + "╝"
    border_empty = "║" + " " * (width - 2) + "║"

    print("\033[91m")  # Rozpocznij czerwony tekst
    print(border_top)
    print("║" + "ERROR".center(width - 2) + "║")
    print(border_empty)

    words = message.split()
    line = ""
    for word in words:
        if len(line) + len(word) + 1 > width - 4:
            print("║ " + line.ljust(width - 4) + " ║")
            line = word
        else:
            line += (" " + word if line else word)
    if line:
        print("║ " + line.ljust(width - 4) + " ║")

    print(border_empty)
    print(border_bottom)
    print("\033[0m")  # Resetuj kolor tekstu

def inf_msg(message, color=None):
    width = 80
    border_top = "╔" + "═" * (width - 2) + "╗"
    border_bottom = "╚" + "═" * (width - 2) + "╝"
    border_empty = "║" + " " * (width - 2) + "║"

    color_codes = {
        'black': '\033[30m',
        'red': '\033[91m',
        'green': '\033[92m',
        'yellow': '\033[93m',
        'blue': '\033[94m',
        'magenta': '\033[95m',
        'cyan': '\033[96m',
        'white': '\033[97m'
    }

    # Jeśli kolor nie jest podany lub jest nieprawidłowy, użyj zielonego
    color_code = color_codes.get(color.lower() if color else None, '\033[92m')
    
    print(color_code)  # Rozpocznij kolorowy tekst
    print(border_top)
    print("║" + "INFORMATION".center(width - 2) + "║")
    print(border_empty)

    words = message.split()
    line = ""
    for word in words:
        if len(line) + len(word) + 1 > width - 4:
            print("║ " + line.ljust(width - 4) + " ║")
            line = word
        else:
            line += (" " + word if line else word)
    if line:
        print("║ " + line.ljust(width - 4) + " ║")

    print(border_empty)
    print(border_bottom)
    print("\033[0m")
    
# def inf_msg(message):
#     width = 80
#     border_top = "╔" + "═" * (width - 2) + "╗"
#     border_bottom = "╚" + "═" * (width - 2) + "╝"
#     border_empty = "║" + " " * (width - 2) + "║"

#     print("\033[92m")  # Rozpocznij zielony tekst
#     print(border_top)
#     print("║" + "INFORMATION".center(width - 2) + "║")
#     print(border_empty)

#     words = message.split()
#     line = ""
#     for word in words:
#         if len(line) + len(word) + 1 > width - 4:
#             print("║ " + line.ljust(width - 4) + " ║")
#             line = word
#         else:
#             line += (" " + word if line else word)
#     if line:
#         print("║ " + line.ljust(width - 4) + " ║")

#     print(border_empty)
#     print(border_bottom)
#     print("\033[0m") 

def initialize_project_structure(project_paths):

    # sciezka gdzie jest input glowny
    # sciezka do outputu
    # sciezka do outputu molpro, dalton, orki
    # sciezka do inputu gammcora

    input_name, _ = os.path.splitext(project_paths['input_file'])
    input_dir = os.getcwd()

    project_paths['input_name'] = input_name
    project_paths['input_dir'] = input_dir
    
    #------------CREATE INPUT DIRECTORY-----------------------
    if debug:
        create_output_dir_debug(project_paths)
    else:
        create_output_dir(project_paths)


    if project_paths['interface']  == 'DALTON':
        create_dalton_dir(project_paths)
    elif project_paths['interface'] == 'MOLPRO':
        create_molpro_dir(project_paths)
    elif project_paths['interface'] == 'ORCA':
        create_orca_dir(project_paths)

    

def parse_basis_content(molecule_data, project_paths):

    file_path = project_paths['basis_path']
    atom_list = molecule_data['atom_list']

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    data_section_started = False
    element_name = None
    basis_set_data = {}  
    current_block = None

    n_of_elems = 0

    for i, line in enumerate(lines):
        if "$DATA" in line:
            data_section_started = True
            continue

        if not data_section_started:
            continue

        if re.match(r'^[A-Z]+', line.strip()) and not re.search(r'\d', line):
            element_name = line.strip()
            element_name = elements_short[element_name]
            if element_name in atom_list:
                basis_set_data[element_name] = {'limits':{}, 'angmoms':[]}  
                basis_set_data[element_name]['limits']['startline'] = i+1
                n_of_elems += 1
                
                endline = None


                for j in range(i + 1, len(lines)):
                    
                    if re.match(r'^[A-Z]+', lines[j].strip()) and not re.search(r'\d', lines[j]):
                        endline = j - 1
                        break
                    elif "$END" in lines[j]:
                        endline = j - 1
                        break
                
                basis_set_data[element_name]['limits']['endline'] = endline


        
    keys = list(basis_set_data.keys())

    for l in range(len(keys)):
        elem = keys[l]

        st = basis_set_data[elem]['limits']['startline']
        nd = basis_set_data[elem]['limits']['endline']


        selected_lines = lines[st:nd]
        for i, line in enumerate(selected_lines):


            if line.startswith(('S', 'P', 'D', 'F','G','H','I','J','K','L','M','N')):

                coefficients_list = []
                angular_momentum = line[0]
                parts = line.split()
            
                if angular_momentum not in basis_set_data[elem]:
                    basis_set_data[elem]['angmoms'].append(angular_momentum)
                    basis_set_data[elem][angular_momentum] = {
                        'num_prim': int(parts[1]),     # Największa liczba w bloku
                        'num_contr': 0,           # Liczba podbloków
                        'exponents': [],            # Lista eksponentów 
                        'coefficients': []          # Lista współczynników 
                    }
                else:
                    if int(parts[1]) > basis_set_data[elem][angular_momentum]['num_prim']:
                        basis_set_data[elem][angular_momentum]['num_prim'] =  int(parts[1])

                basis_set_data[elem][angular_momentum]['num_contr'] += 1  #


                coefficients = []
                exponents = []
                k = 0
                for next_line in selected_lines[i+1:]:
                    if next_line.startswith(('S', 'P', 'D', 'F','G','H')):
                        basis_set_data[elem][angular_momentum]['exponents'].append(exponents)
                        basis_set_data[elem][angular_momentum]['coefficients'].append(coefficients)
                        break
                    elif i+st+k+2 == nd:
                        parsed_values = re.findall(r'[\d\.\-Ee\+]+', next_line)
                        if len(parsed_values) >= 3:
                            exponents.append(parsed_values[1])
                            coefficients.append(parsed_values[2])
                            k+=1
                        basis_set_data[elem][angular_momentum]['exponents'].append(exponents)
                        basis_set_data[elem][angular_momentum]['coefficients'].append(coefficients)
                        break
                    else:
                        parsed_values = re.findall(r'[\d\.\-Ee\+]+', next_line)
                        if len(parsed_values) >= 3:
                            exponents.append(parsed_values[1])
                            coefficients.append(parsed_values[2])
                            k+=1
                            


    for elem in basis_set_data:
        del basis_set_data[elem]['limits']
        for angular_momentum in basis_set_data[elem]['angmoms']:
            sorted_data= sort_and_index(basis_set_data[elem][angular_momentum]['exponents'], basis_set_data[elem][angular_momentum]['coefficients'])
            basis_set_data[elem][angular_momentum]['exp'] = sorted_data['exp']
            basis_set_data[elem][angular_momentum]['num_prim'] = len(basis_set_data[elem][angular_momentum]['exp'])
            basis_set_data[elem][angular_momentum]['contractions'] = sorted_data['contractions']
            basis_set_data[elem][angular_momentum]['idx'] = sorted_data['idx']
            del basis_set_data[elem][angular_momentum]['coefficients']
            del basis_set_data[elem][angular_momentum]['exponents']
            block = basis_set_data[elem][angular_momentum]
    
    return basis_set_data



def sort_and_index(exp, coef):

    exp_float = [[float(num) for num in sublist] for sublist in exp]
    number_of_contractions = len(coef)
    combined_list = [[val, i+1, j+1] for i, sublist in enumerate(exp_float) for j, val in enumerate(sublist)]
    sorted_list = sorted(combined_list, key=lambda x: float(x[0]), reverse=True)
    grouped_dict = defaultdict(list)

    for value, list_idx, sublist_idx in sorted_list:
        grouped_dict[float(value)].append([list_idx, sublist_idx])

    grouped_sorted_list = [(value, indices) for value, indices in grouped_dict.items()]

    grouped_sorted_list.sort(reverse=True, key=lambda x: x[0])


    this_dict = {}
    this_dict['exp'] = []
    this_dict['contractions'] = coef
    this_dict['idx'] = [[] for _ in range(number_of_contractions)]

    for k, elem in enumerate(grouped_sorted_list):        
        this_dict['exp'].append(f"{elem[0]:.6E}")
        for  l in range(0, number_of_contractions):
            for minilist in elem[1]:
                if (minilist[0]-1) == l:
                    this_dict['idx'][l].append(k)


    return this_dict
    

def format_dalton(molecule_data, basis_set_data):
    # Formats the parsed BASIS GAMESS US data into the Dalton format

    dalton_formatted = ""
    geometry = molecule_data['geometry']
    atomdict = {}
    
    for i, elem in enumerate(geometry):
        if elem['atom'] not in atomdict.keys():
            atomdict[elem['atom']] = [i]
        else:
            atomdict[elem['atom']].append(i)

    # sort geometry for dalton input
    new_geom = []
    for elem in atomdict:
        for i in range(0, len(atomdict[elem])):
            new_geom.append(geometry[atomdict[elem][i]])

    geometry = new_geom

    # copy new geometry for gammcor input for consistency
    molecule_data['geometry'] = deepcopy(geometry)
    basis_for_this_atom = []

    for atom in atomdict:
        instances = len(atomdict[atom])

    k = 0
    for dct in geometry:
        elem = dct['atom']
        instances = len(atomdict[elem])

        bl = len(basis_set_data[elem]['angmoms'])
        block_string = f"{bl}"
        for kk in range(0, bl):
            block_string += " 1"
        
        this_atom = ""

        if k==0:
            this_atom += f"Charge={dct['atomic_number']:.1f} Atoms={instances} Blocks={block_string} \n"
        
        dct['coordinates'] = [f"{coord:.6f}" for coord in dct['coordinates']]

        spaces_needed = 9 - len(dct['coordinates'][0].split('.')[0])

        mini = f"{elem}" + ' '*spaces_needed
        for coord in dct['coordinates']:
            mini += coord +' ' * 5
        this_atom += mini+'\n'
        dalton_formatted+= this_atom
        basis_for_this_atom = ""

        k+=1
        if  k == instances: 
            k = 0
            for angular_momentum in basis_set_data[elem]['angmoms']:
                mx = basis_set_data[elem][angular_momentum]['num_prim']

                bl = basis_set_data[elem][angular_momentum]['num_contr']
                dalton_formatted += f"H  {mx}    {bl}\n"
                for ln in range(0, mx):

                    exp = basis_set_data[elem][angular_momentum]['exp'][ln]
                    first_offset = " " *5
                    current_line = ""
                    line = ""
                    for x in range(0, bl):
                        coef = basis_set_data[elem][angular_momentum]['contractions'][x]
                        coef_idx = basis_set_data[elem][angular_momentum]['idx'][x]


                        if ln in coef_idx:

                            index = coef_idx.index(ln)
                            part = f"{float(coef[index]):>20.6E}"
                        else:
                            part = f"{float(0.000000):>20.6E}"

                        if len(current_line + part) > 80:
                            line = deepcopy(current_line) + "\n"
                            offset = " "* 20 + first_offset
                            current_line = offset + part
                        else:
                            current_line += part

                    line += current_line
                    line += "\n"
                    line = f"{first_offset}{float(exp):>15.6E}" + line
                    dalton_formatted += line
    

    return dalton_formatted



def angs_to_bohr(r):

    return r* 1.8897259886

def standarize_basis(basis):
    standarized = basis.lower().replace('_', '-').replace(' ', '-')
    return f"{standarized}"

def standarize_interface(interface):

    standarized = None
    if re.match(r'(?i)dal.*', interface):
        standarized = 'DALTON'        
    elif re.match(r'(?i)mol.*', interface):
        standarized = 'MOLPRO'
    elif re.match(r'(?i)orc.*', interface):
        standarized = 'ORCA'
    else:
        error_msg(f'NO SUCH INTERFACE IMPLEMENTED  {interface}')
        sys.exit(1)
        
    return standarized

def standarize_units(units):

    standarized = None
    if re.match(r'(?i)a.*', units):
        standarized = 'ANGS'
    elif re.match(r'(?i)b.*', units):
        standarized = 'BOHR'
    else:
        inf_msg(f'UNITS NOT GIVEN, DEFAULT UNITS USED [BOHR] {units}. Continue.')
        
    return standarized

def read_geometry(geometry, units, interface):

    znucl = 0
    geom_lines = geometry.split('\n')

    geometry_data = []
    coords_set = set()
    atom_count = 0

    for line in geom_lines:
 
        line_content = line.split()
        
        if len(line_content) > 4:
            error_msg(f'INCORRECT GEOMETRY in LINE:{line}. \nExiting')
            sys.exit(1)

        atom_symbol = line_content[0].upper()
        atomic_number = elements_dict[atom_symbol]
        znucl += atomic_number
        atom_count += 1

        if len(line_content) > 1:
            coord_content = " ".join(line_content[1:]).replace(",", " ")
#            coords = [f"{float(coord):.13f}" for coord in coord_content.split()]
            coords = [round(float(coord), 13) for coord in coord_content.split()]
        else:
            if atom_count == 1:
#                coords = [f"{0.0:.13f}", f"{0.0:.13f}", f"{0.0:.13f}"]
                coords = [0.0000, 0.0000, 0.0000]
            else:
                error_msg(f'ERROR: Atom {atom_symbol} missing coordinates. \nExiting')
                sys.exit(1)

        coord_tuple = tuple(coords)
        if coord_tuple in coords_set:
            error_msg(f'ERROR: Atom {atom_symbol} have duplicate coordinates {coords}. \nExiting')
            sys.exit(1)
        else:
            coords_set.add(coord_tuple)


        geometry_data.append({'atom': atom_symbol,
            'atomic_number': atomic_number, 
            'coordinates': coords})
    
    atom_list = []
    for dct in geometry_data:
        atom_list.append(dct['atom'])


    return geometry_data, atom_list

def parse_dalton(project_paths):

    input_file = project_paths['input_file']
    
    with open(input_file, 'r') as f:
        content = f.read()

    dalton_match = re.search(r'^\s*dalton\s*\n(.*?)\nend', content, re.DOTALL | re.IGNORECASE | re.MULTILINE)
    if not dalton_match:
        error_msg(f'INTERFACE DALTON BUT NO DALTON BLOCK SPECIFIED, EXITING\n')
        sys.exit(1)

    dalton_block = dalton_match.group(1).strip()

    dalton_data = {}

    nsym_match = re.search(r'nsym\s+(\d+)', dalton_block, re.IGNORECASE)
    symmetry_match = re.search(r'symmetry\s+(\d+)', dalton_block, re.IGNORECASE)
    multip_match = re.search(r'mul.*\s+(\d+)', dalton_block, re.IGNORECASE)
    inactive_match = re.search(r'inactive\s+(.+)', dalton_block, re.IGNORECASE)
    electrons_match = re.search(r'el.*\s+(\d+)', dalton_block, re.IGNORECASE)
    cas_match = re.search(r'cas\s+(.+)', dalton_block, re.IGNORECASE)
    state_match = re.search(r'state.*\s+(\d+)', dalton_block, re.IGNORECASE)

    if nsym_match:
        nsym = int(nsym_match.group(1))
    else:
        nsym = 1  

    if symmetry_match:
        symmetry = int(symmetry_match.group(1))
        if symmetry > nsym:
            error_msg(f' Dalton block: symmetry value exceeds number of symmetries, exiting\n')
            sys.exit(1)
    else:
        symmetry = 1


    
    if multip_match:
        mult = int(multip_match.group(1))
    else:
        inf_msg('MULTIPLICITY NOT SPECIFIED, ASSUMING SINGLET, MULT = 1\n')
        mult = 1

    if state_match:
        state = int(state_match.group(1))
    else:
        state = 1  


    if inactive_match:
        inactive_list = inactive_match.group(1).split()

        if len(inactive_list) != nsym:
            error_msg(f'Daltion_block: Expected {nsym} values for inactive, got {len(inactive_list)}.\nExiting')
            sys.exit(1)
        else:
            inactive = inactive_match.group(1).strip()            
            
    if cas_match:
        cas_list = cas_match.group(1).split()

        if len(cas_list) != nsym:
            error_msg(f'Dalton Block: Expected {nsym} values for cas, got {len(cas_list)}.\nExiting')
            sys.exit(1)
        else:
            cas = cas_match.group(1).strip()            

    if electrons_match:
        electrons = int(electrons_match.group(1))
    else:
        error_msg(f'Dalton Block: ELECRTONS NOT SPECIFIED, EXITING\n')
        sys.exit(1)

    dalton_data = {
        'nsym': nsym,
        'symmetry': symmetry,
        'mult': mult,
        'inactive': inactive,
        'electrons': electrons,
        'cas': cas,
        'state': state
    }

    return dalton_data

def parse_input_file(input_file):

    
    with open(input_file, 'r') as f:
        full_content = f.read()

    content = '\n'.join(line for line in full_content.split('\n') if not line.strip().startswith('!'))


    #----------BASIS--------------------------------------------------------------------#
    basis_match = re.search(r'^\s*basis\s*(\S+)', content, re.IGNORECASE | re.MULTILINE)
    if not(basis_match):
        error_msg(f'Input: NO BASIS SET SPECIFIED, EXITING \n')
        sys.exit(1)
    else:
        basis = basis_match.group(1)
        basis = standarize_basis(basis)

    #----------CHARGE--------------------------------------------------------------------#
    charge_match = re.search(r'charge\s*(\S+)', content, re.IGNORECASE)
    if not(charge_match):
        # ASSUMING CHARGE = 0
        charge = 0
    else:
        charge = charge_match.group(1)

    #----------UNITS--------------------------------------------------------------------#    
    units_match = re.search(r'units\s*(\S+)', content, re.IGNORECASE)
    if not(units_match):
        inf_msg('NO UNITS SPECIFIED, USING DEFAULT [BOHR]\n')
        units = 'BOHR'
    else:
        units = units_match.group(1)
        units = standarize_units(units)

    #----------INTERFACE----------------------------------------------------------------#    
    interface_match = re.search(r'interface\s*(\S+)', content, re.IGNORECASE)
    if not(interface_match):
        error_msg(f'NO INTERFACE SPECIFIED, EXITING\n')
        sys.exit(1)
    else:
        interface = interface_match.group(1)
        interface = standarize_interface(interface)

    #----------GAMMCOR-PATH-------------------------------------------------------------#
    gammcor_path = os.path.dirname(os.path.abspath(__file__))
    print(gammcor_path)

    
    # gammcor_match = re.search(r'gammcor_path\s+(\S+)', content, re.IGNORECASE)

    # if not(gammcor_match):
    #     error_msg(f'NO GAMMCOR PATH SPECIFIED, EXITING \n')
    #     sys.exit(1)
    # else:
    #     gammcor_path = gammcor_match.group(1)


    #----------MY_TEMPLATE-PATH-------------------------------------------------------------#                                                                                                                          
    my_template_match = re.search(r'my_template_path\s+(\S+)', content, re.IGNORECASE)

    if not(my_template_match):
        my_template_path = False
    else:
        my_template_path = my_template_match.group(1)


    #----------GEOMETRY---------------------------------------------------------------#    
    geometry_match =  re.search(r'xyz\n(.*?)\nend', content, re.DOTALL | re.IGNORECASE)
    if not(geometry_match):
        error_msg(f'NO GEOMETRY SPECIFIED, EXITING\n')
        sys.exit(1)
    else:
        geometry = geometry_match.group(1).strip()    
        geometry_data, atom_list = read_geometry(geometry, units, interface)

    #------------------------------------------------------------------------------#

    molecule_data = {
        'basis': basis,
        'interface': interface,
        'units': units,
        'geometry': geometry_data,
        'charge': charge,
        'atom_list': atom_list
    }


    #----------INTERFACE-PATH-------------------------------------------------------------#

    dalton_path = ''
    molpro_path = ''
    orca_path = ''

    if molecule_data['interface'] == 'DALTON':
        dalton_match = re.search(r'dalton_path\s+(\S+)', content, re.IGNORECASE)

        if not(dalton_match):
            error_msg(f'NO DALTON EXECUTABLE PATH SPECIFIED, add \n "dalton_path <path_address> \n in input file. \nExiting"')
            sys.exit(1)
        else:
            dalton_path = dalton_match.group(1)

    elif molecule_data['interface'] == 'MOLPRO':
        molpro_match = re.search(r'molpro_path\s+(\S+)', content, re.IGNORECASE)
        
        if not(molpro_match):
            error_msg('NO MOLPRO EXECUTABLE PATH SPECIFIED, add \n "molpro_path <path_address> \n in input file"')
            sys.exit(1)
        else:
            molpro_path = molpro_match.group(1)
            
    elif molecule_data['interface'] == 'ORCA':
        orca_match = re.search(r'orca_path\s+(\S+)', content, re.IGNORECASE)
        
        if not(orca_match):
            error_msg('NO ORCA EXECUTABLE PATH SPECIFIED, add \n "orca_path <path_address> \n in input file"')
            sys.exit(1)
        else:
            orca_path = orca_match.group(1)

    basis_path = f"{gammcor_path}/bazy-od-marcina/{basis}.txt"

    project_paths = {
        'basis_path': basis_path,
        'interface': interface,
        'gammcor_path': gammcor_path,
        'dalton_path': dalton_path,
        'molpro_path': molpro_path,
        'orca_path': orca_path,
        'input_file': input_file,
        'my_template_path': my_template_path
        }
    
    return molecule_data, project_paths

def create_output_dir_debug(project_paths):

    output_dir = os.path.join(project_paths['input_dir'], project_paths['input_name'])

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    project_paths['output_dir'] = output_dir

def create_output_dir(project_paths):

    input_name = project_paths['input_name']

    original_dir_path = os.path.join(project_paths['input_dir'], input_name)

    output_dir = original_dir_path

    suffix = 0

    while os.path.exists(output_dir):
        output_dir = os.path.join(project_paths['input_dir'], f"{input_name}_{suffix}")
        suffix += 1
        if suffix == 69:
            print('nice')
        elif suffix == 420:
            print('nice')
        elif suffix > 999:
            raise Exception("Achievement unlocked, over 999 directories with the same name")


    if os.path.exists(original_dir_path):
        shutil.move(original_dir_path, output_dir)
        output_dir = os.path.join(project_paths['input_dir'], input_name)
        os.mkdir(output_dir)
    else:
        os.mkdir(output_dir)

    project_paths['output_dir'] = output_dir



def run_gammcor(molecule_data, basis_set_data, project_paths):

    input_file = project_paths['input_file']
    input_name = project_paths['input_name']
    input_dir = project_paths['input_dir']
    output_dir = project_paths['output_dir']
    gammcor_path = project_paths['gammcor_path']
    
    os.chdir(input_dir)
    
    gammcor_input_name = "input.inp"
    gammcor_input = os.path.join(output_dir, gammcor_input_name)

    shutil.copy(input_file, gammcor_input)

    #----sort geometry

    update_gammcor_input(gammcor_input, molecule_data, basis_set_data, project_paths)


    os.chdir(output_dir)

    command = f'{gammcor_path}/gammcor > {output_dir}/{input_name}.out'

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    stdout, stderr = process.communicate()
    if process.returncode != 0:
        error_msg(f"Gammcor eroor check output, Exiting: {stderr.decode()}")
    else:
        inf_msg("Gammcor executed succesfully.", color = 'blue')

    os.chdir(input_dir)
    sys.exit(0)

def update_gammcor_input(gammcor_input, molecule_data, basis_set_data, project_paths):

    new_xyz = "\n".join([
        f"{len(molecule_data['geometry'])}",
        *[f"{entry['atom']} {' '.join(f'{x:.6f}' for x in entry['coordinates'])}"
          for entry in molecule_data['geometry']]
    ])
    nbasis = sum(
        (2 * angular_momentum_dict[angmom] + 1) * basis_set_data[atom['atom']][angmom]['num_contr']
        for atom in molecule_data['geometry']
        for angmom in basis_set_data[atom['atom']]['angmoms']
    )

    with open(gammcor_input, 'r') as file:
        content = file.read()

    content = re.sub(r'xyz\n.*?end', f'xyz\n{new_xyz}\nend', content, flags=re.DOTALL)

    if re.search(r'NBasis\s+\d+', content, re.IGNORECASE):
        content = re.sub(r'NBasis\s+\d+', f'NBasis {nbasis}', content, flags=re.IGNORECASE)
    else:
        content = re.sub(
            r'(calculation\s*\n)',
            f'\\1 NBasis {nbasis}\n',
            content,
            flags=re.IGNORECASE
        )

    basis_path = project_paths['basis_path']
    content = re.sub(
        r'(calculation\s*\n.*?)\b(basis\s+\S+)',
        f'\\1BasisPath {basis_path}',
        content,
        flags=re.DOTALL | re.IGNORECASE
    )

    with open(gammcor_input, 'w') as file:
        file.write(content)    

def run_dalton(project_paths):

    dalton_dir = project_paths['dalton_dir']
    input_name = project_paths['input_name']
    input_dir = project_paths['input_dir']
    dalton_path = project_paths['dalton_path']
    
    dalton_files_list = ['SIRIFC', 'AMFI_SYMINFO.TXT', 'AOONEINT','AOTWOINT', 'DALTON.MOPUN', \
                         'SIRIUS.RST', 'rdm2.dat', 'rdms2.dat', 'rdms1.dat']

    os.chdir(dalton_dir)

    if not debug:
        arguments = " ".join(dalton_files_list)
        command = f'{dalton_path} -get "{arguments}" {input_name}'
    
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        stdout, stderr = process.communicate()
        if process.returncode != 0:
            error_msg(f"DALTON ERROR, check output, Exiting: {stderr.decode()}")
        else:
            inf_msg("DALTON executed succesfully.")


    destination_dir = os.path.join(dalton_dir, '..')

    for dalton_file in dalton_files_list:

        this_name = f"{input_name}.{dalton_file}"
        new_name = f"{dalton_file}"
        
        if os.path.exists(this_name):
            os.rename(this_name, new_name)

            source_path = os.path.join(dalton_dir, new_name)
            destination_path = os.path.join(destination_dir, new_name)
    
            shutil.copy(source_path, destination_path)

        else:
            inf_msg(f"Warning: File '{this_name}' does not exist and was not copied.")

    os.chdir(input_dir)
        
def create_dalton_dir(project_paths):

    input_name = project_paths['input_name']
    dalton_dir = os.path.join(project_paths['output_dir'], 'dal_files')
    project_paths['dalton_dir'] = dalton_dir

    if not os.path.exists(dalton_dir):
        os.makedirs(dalton_dir)

    if project_paths['my_template_path'] is not False:
        template_dal_path = project_paths['my_template_path']
        #= os.path.join(project_paths['input_dir'], f'{input_name}.dal')
        if not os.path.exists(template_dal_path):
            error_msg(f"Template file {template_dal_path} not found.Please provide your own .dal template or set my_template_path to False to use the default template.\nExiting")
            sys.exit(1)
    else:
        template_dal_path = os.path.join(project_paths['gammcor_path'], 'templates', 'dalton-temp.dal')

    output_dal = os.path.join(dalton_dir, f"{input_name}.dal")
    output_mol = os.path.join(dalton_dir, f"{input_name}.mol")
    shutil.copy(template_dal_path, output_dal)

    project_paths['.dal'] = output_dal
    project_paths['.mol'] = output_mol

    if debug:
        dalton_files_list = ['SIRIFC', 'AMFI_SYMINFO.TXT', 'AOONEINT','AOTWOINT', 'DALTON.MOPUN', \
                             'SIRIUS.RST', 'rdm2.dat', 'rdms2.dat', 'rdms1.dat']
        for file_ext in dalton_files_list:
            open(os.path.join(dalton_dir, f"{input_name}.{file_ext}"), 'w').close()

def create_molpro_dir(project_paths):

    input_name = project_paths['input_name']
    molpro_directory = os.path.join(project_paths['output_dir'], 'molpro_files')

    if not os.path.exists(molpro_directory):
        os.makedirs(molpro_directory)

    if project_paths['my_templates'] is True:
        template_mol_path = os.path.join(project_paths['input_dir'], f'{input_name}-molpro.inp')
        if not os.path.exists(template_dal_path):
            error_msg(f"Error: Template file {template_mol_path} not found.Please provide your own molpro template or set my_templates to False to use the default template.\nExiting")
            sys.exit(1)
    else:
        template_mol_path = os.path.join(project_paths['gammcor_path'], 'templates', 'molpro.inp')

    output_molpro = os.path.join(molpro_directory, f"{input_name}.dal")
    shutil.copy(template_mol_path, output_molpro)

    project_paths['molpro_input'] = output_molpro
            
            
def create_dalton_input(molecule_data, dalton_data, basis_set_data, project_paths):

    dal_file = project_paths['.dal']
    mol_file = project_paths['.mol']
    
    with open(dal_file, 'r') as file:
        content = file.read()
        
    content = content.format(
        symmetry=dalton_data['symmetry'],
        spin_mul=dalton_data['mult'],
        inactive=' '.join(map(str, dalton_data['inactive'])), 
        electrons=dalton_data['electrons'],
        cas=dalton_data['cas'],
        state=dalton_data['state']
    )

    with open(dal_file, 'w') as file:
        file.write(content)

    basis_path = project_paths['basis_path']
    basis = molecule_data['basis']

    dalton_basis = format_dalton(molecule_data, basis_set_data)
    
    atomtypes = len(basis_set_data)

    atom_basis = ""

    with open(mol_file, 'w') as f:
        f.write("INTGRL\n")
        f.write(f"{basis} basis set\n")
        f.write("Generated by Python script\n")
        f.write(f"Atomtypes={atomtypes} Generators=0 0 {molecule_data['units']} Integrals=1.00D-15\n")
        f.write(dalton_basis)
        f.write("\n")

    inf_msg('Dalton files created succesfully.')

def test_print(basis_set_data):

    for elem in basis_set_data:
        print('element:', elem)
        print('blocks:', basis_set_data[elem]['angmoms'])
        for angular_momentum in basis_set_data[elem]['angmoms']:
            block = basis_set_data[elem][angular_momentum]
            print(angular_momentum)
            for key in block:
                print(key, block[key])
            for l in range(0, len(block['contractions'])):
                print('kontrakcja', l+1, block['contractions'][l], block['idx'][l])
            print('')

def main():

    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        if not input_file.endswith('.inp'):
            input_file += '.inp' 
    else:
        error_msg(f'PROVIDE INPUT FILE AS ARGUMENT')
        sys.exit(0)
    
    molecule_data, project_paths = parse_input_file(input_file)
    initialize_project_structure(project_paths)
    
    basis_set_data = parse_basis_content(molecule_data, project_paths)

    # test_print(basis_set_data)



    if molecule_data['interface'] == 'DALTON':

        dalton_data = parse_dalton(project_paths)

        create_dalton_input(molecule_data, dalton_data, basis_set_data, project_paths)

        run_dalton(project_paths)

        run_gammcor(molecule_data, basis_set_data, project_paths)

        
    elif molecule_data['interface'] == 'MOLPRO':
        error_msg(f'making molpro input - not implemented')
        sys.exit(1)
    elif molecule_data['interface'] == 'ORCA':
        error_msg('making ORCA input - not implemented')
        sys.exit(1)
        



if __name__ == "__main__":
    main()
