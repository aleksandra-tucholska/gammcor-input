#!/usr/bin/env python3
import sys
import os
import re
from atomic_number import elements_dict
from atomic_number import elements_short
import shutil
from copy import deepcopy
from collections import defaultdict

dalton_path = "/dalton"
molpro_path = "/molpro"
orca_path = "/orca"



def parse_basis_content(data):

    file_path = data['basis_path']
    atom_list = data['atom_list']

    with open(file_path, 'r') as file:
        lines = file.readlines()


    
    data_section_started = False
    element_name = None
    elements_data = {}  
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
                elements_data[element_name] = {'limits':{}, 'blocks':[]}  
                elements_data[element_name]['limits']['startline'] = i+1
                n_of_elems += 1
                
                endline = None


                for j in range(i + 1, len(lines)):
                    
                    if re.match(r'^[A-Z]+', lines[j].strip()) and not re.search(r'\d', lines[j]):
                        endline = j - 1
                        break
                    elif "$END" in lines[j]:
                        endline = j - 1
                        break
                
                elements_data[element_name]['limits']['endline'] = endline


        
    keys = list(elements_data.keys())

    for l in range(len(keys)):
        elem = keys[l]

        st = elements_data[elem]['limits']['startline']
        nd = elements_data[elem]['limits']['endline']


        selected_lines = lines[st:nd]
        for i, line in enumerate(selected_lines):


            if line.startswith(('S', 'P', 'D', 'F','G','H','I','J','K','L','M','N')):

                coefficients_list = []
                block_type = line[0]
                parts = line.split()
            
                if block_type not in elements_data[elem]:
                    elements_data[elem]['blocks'].append(block_type)
                    elements_data[elem][block_type] = {
                        'num_prim': int(parts[1]),     # Największa liczba w bloku
                        'num_contr': 0,           # Liczba podbloków
                        'exponents': [],            # Lista eksponentów 
                        'coefficients': []          # Lista współczynników 
                    }
                else:
                    if int(parts[1]) > elements_data[elem][block_type]['num_prim']:
                        elements_data[elem][block_type]['num_prim'] =  int(parts[1])

                elements_data[elem][block_type]['num_contr'] += 1  # Liczymy podblok

                coefficients = []
                exponents = []
                k = 0
                for next_line in selected_lines[i+1:]:
                    if next_line.startswith(('S', 'P', 'D', 'F','G','H')):
                        elements_data[elem][block_type]['exponents'].append(exponents)
                        elements_data[elem][block_type]['coefficients'].append(coefficients)
                        break
                    elif i+st+k+2 == nd:
                        parsed_values = re.findall(r'[\d\.\-Ee\+]+', next_line)
                        if len(parsed_values) >= 3:
                            exponents.append(parsed_values[1])
                            coefficients.append(parsed_values[2])
                            k+=1
                        elements_data[elem][block_type]['exponents'].append(exponents)
                        elements_data[elem][block_type]['coefficients'].append(coefficients)
                        break
                    else:
                        parsed_values = re.findall(r'[\d\.\-Ee\+]+', next_line)
                        if len(parsed_values) >= 3:
                            exponents.append(parsed_values[1])
                            coefficients.append(parsed_values[2])
                            k+=1
                            


    for elem in elements_data:
        del elements_data[elem]['limits']
        for block_name in elements_data[elem]['blocks']:
            sorted_data= sort_and_index(elements_data[elem][block_name]['exponents'], elements_data[elem][block_name]['coefficients'])
            elements_data[elem][block_name]['exp'] = sorted_data['exp']
            elements_data[elem][block_name]['num_prim'] = len(elements_data[elem][block_name]['exp'])
            elements_data[elem][block_name]['contractions'] = sorted_data['contractions']
            elements_data[elem][block_name]['idx'] = sorted_data['idx']
            del elements_data[elem][block_name]['coefficients']
            del elements_data[elem][block_name]['exponents']
            block = elements_data[elem][block_name]




        

    
    return elements_data



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
    

def format_dalton(data, elements_data):
    # Formats the parsed BASIS GAMESS US data into the Dalton format

    dalton_formatted = ""

    geometry = data['geometry']

    atomdict = {}
    
    for i, elem in enumerate(geometry):
        if elem['atom'] not in atomdict.keys():
            atomdict[elem['atom']] = [i]
        else:
            atomdict[elem['atom']].append(i)

    new_geom = []
    for elem in atomdict:
        for i in range(0, len(atomdict[elem])):
            new_geom.append(geometry[atomdict[elem][i]])

    geometry = new_geom
    basis_for_this_atom = []

    for atom in atomdict:
        instances = len(atomdict[atom])

    k = 0
    for dct in geometry:
        elem = dct['atom']

        instances = len(atomdict[elem])


        bl = len(elements_data[elem]['blocks'])
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
            for block_name in elements_data[elem]['blocks']:
                mx = elements_data[elem][block_name]['num_prim']

                bl = elements_data[elem][block_name]['num_contr']
                dalton_formatted += f"H  {mx}    {bl}\n"
                for ln in range(0, mx):

                    exp = elements_data[elem][block_name]['exp'][ln]
                    first_offset = " " *5
                    current_line = ""
                    line = ""
                    for x in range(0, bl):
                        coef = elements_data[elem][block_name]['contractions'][x]
                        coef_idx = elements_data[elem][block_name]['idx'][x]


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
    return f"{standarized}.dat"

def standarize_interface(interface):

    standarized = None
    if re.match(r'(?i)dal.*', interface):
        standarized = 'DALTON'        
    elif re.match(r'(?i)mol.*', interface):
        standarized = 'MOLPRO'
    elif re.match(r'(?i)orc.*', interface):
        standarized = 'ORCA'
    else:
        print('NO SUCH INTERFACE IMPLEMENTED', interface)
        sys.exit(1)
        
    return standarized

def standarize_units(units):

    standarized = None
    if re.match(r'(?i)a.*', units):
        standarized = 'ANGS'
    elif re.match(r'(?i)b.*', units):
        standarized = 'BOHR'
    else:
        print('DEFAULT UNITS USED [BOHR]', units)
        sys.exit(1)
        
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
            print('INCORRECT GEOMETRY in LINE:', line)
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
                print(f'ERROR: Atom {atom_symbol} missing coordinates')
                sys.exit(1)

        coord_tuple = tuple(coords)
        if coord_tuple in coords_set:
            print(f'ERROR: Atom {atom_symbol} have duplicate coordinates {coords}')
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

def parse_dalton(input_file):
    
    with open(input_file, 'r') as f:
        content = f.read()



    dalton_match = re.search(r'^\s*dalton\s*\n(.*?)\nend', content, re.DOTALL | re.IGNORECASE | re.MULTILINE)
    if not dalton_match:
        print('INTERFACE DALTON BUT NO DALTON BLOCK SPECIFIED, EXITING\n')
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
            print('ERROR: Symmetry value exceeds number of symmetries, exiting\n')
            sys.exit(1)
    else:
        symmetry = 1


    
    if multip_match:
        mult = int(multip_match.group(1))
    else:
        print('MULTIPLICITY NOT SPECIFIED, ASSUMING SINGLET, MULT = 1\n')
        mult = 1

    if state_match:
        state = int(state_match.group(1))
    else:
        state = 1  


    if inactive_match:
        inactive_list = inactive_match.group(1).split()

        if len(inactive_list) != nsym:
            print(f'ERROR: Expected {nsym} values for inactive, got {len(inactive_list)}')
            sys.exit(1)
        else:
            inactive = inactive_match.group(1).strip()            
            
    if cas_match:
        cas_list = cas_match.group(1).split()

        if len(cas_list) != nsym:
            print(f'ERROR: Expected {nsym} values for cas, got {len(cas_list)}')
            sys.exit(1)
        else:
            cas = cas_match.group(1).strip()            

    if electrons_match:
        electrons = int(electrons_match.group(1))
    else:
        print('ERROR: ELECRTONS NOT SPECIFIED, EXITING\n')
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
        content = f.read()

    #----------BASIS--------------------------------------------------------------------#
    basis_match = re.search(r'^\s*basis\s*(\S+)', content, re.IGNORECASE | re.MULTILINE)
    if not(basis_match):
        print('NO BASIS SET SPECIFIED, EXITING \n')
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
        print('NO UNITS SPECIFIED, USING DEFAULT [BOHR]\n')
        units = 'BOHR'
    else:
        units = units_match.group(1)
        units = standarize_units(units)

    #----------INTERFACE----------------------------------------------------------------#    
    interface_match = re.search(r'interface\s*(\S+)', content, re.IGNORECASE)
    if not(interface_match):
        print('NO INTERFACE SPECIFIED, EXITING\n')
        sys.exit(1)
    else:
        interface = interface_match.group(1)
        interface = standarize_interface(interface)

    #----------GAMMCOR-PATH-------------------------------------------------------------#    
    gammcor_match = re.search(r'gammcor_path\s+(\S+)', content, re.IGNORECASE)

    if not(gammcor_match):
        print('NO GAMMCOR PATH SPECIFIED, EXITING \n')
        print('add some default path?')
        sys.exit(1)
    else:
        gammcor_path = gammcor_match.group(1)
        print('Gammcor Path:', gammcor_path)

    #----------GEOMETRY---------------------------------------------------------------#    
    geometry_match =  re.search(r'xyz\n(.*?)\nend', content, re.DOTALL | re.IGNORECASE)
    if not(geometry_match):
        print('NO GEOMETRY SPECIFIED, EXITING\n')
        sys.exit(1)
    else:
        geometry = geometry_match.group(1).strip()    
        geometry_data, atom_list = read_geometry(geometry, units, interface)

    #------------------------------------------------------------------------------#

    basis_path = f"{gammcor_path}/basis/{basis}"

    data = {
        'basis': basis,
        'basis_path': basis_path,
        'gammcor_path': gammcor_path,
        'interface': interface,
        'units': units,
        'geometry': geometry_data,
        'charge': charge,
        'atom_list': atom_list
    }

    
    return data


def create_dalton_input(data, dalton_data, elements_data):
    current_path = os.getcwd()
    dal_directory = os.path.join(current_path, 'dal_files')
    print('dal_directory', dal_directory)
    print('')

    if not os.path.exists(dal_directory):
        os.makedirs(dal_directory)

    gp = data['gammcor_path']
    template_dal_path = os.path.join(gp, 'templates', 'dalton-temp.dal')
    print('template_dal_path', template_dal_path)
    print('')

    output_dal_name = 'nazwa.dal'
    output_mol_name = 'nazwa.mol'
    output_dal = os.path.join(dal_directory, output_dal_name)
    output_mol = os.path.join(dal_directory, output_mol_name)
    print(' output_dal', output_dal)
    print(' output_moll', output_mol)
    print('')

    shutil.copy(template_dal_path, output_dal)

    with open(output_dal, 'r') as file:
        content = file.read()
        
    content = content.format(
        symmetry=dalton_data['symmetry'],
        spin_mul=dalton_data['mult'],
        inactive=' '.join(map(str, dalton_data['inactive'])), 
        electrons=dalton_data['electrons'],
        cas=dalton_data['cas'],
        state=dalton_data['state']
    )
    with open(output_dal, 'w') as file:
        file.write(content)

    print(f'Plik {output_dal_name} został utworzony w {dal_directory}.')

    
    basis_path = data['basis_path']
    basis = data['basis']
    dalton_basis = format_dalton(data, elements_data)

    
    atomtypes = len(elements_data)

    atom_basis = ""


    

    with open(output_mol, 'w') as f:
        f.write("INTGRL\n")
        f.write(f"{basis} basis set\n")
        f.write("Generated by Python script\n")
        f.write(f"Atomtypes={atomtypes} Generators=0 0 {data['units']} Integrals=1.00D-15\n")
        f.write(dalton_basis)
        f.write("\n")


def main():

    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        print(input_file)
    else:
        print('PROVIDE INPUT FILE AS ARGUMENT')
        sys.exit(0)


    data = parse_input_file(input_file)
    elements_data = parse_basis_content(data)

    for elem in elements_data:
        print('element:', elem)
        print('blocks:', elements_data[elem]['blocks'])
        for block_name in elements_data[elem]['blocks']:
            block = elements_data[elem][block_name]
            print(block_name)
            for key in block:
                print(key, block[key])
            for l in range(0, len(block['contractions'])):
                print('kontrakcja', l+1, block['contractions'][l], block['idx'][l])
            print('')

    

    if data['interface'] == 'DALTON':

        dalton_data = parse_dalton(input_file)
        
        print('making DALTON file')
        create_dalton_input(data, dalton_data, elements_data)
        
    elif data['interface'] == 'MOLPRO':
        print('making molpro input - not implemented')
        sys.exit(1)
    elif data['interface'] == 'ORCA':
        print('making ORCA input - not implemented')
        sys.exit(1)
        



if __name__ == "__main__":
    main()
