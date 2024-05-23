import os

class OutcarParser:
    """Class to parse a OUTCAR file"""

    def __init__(self, filename):
        """parse OUTCAR and find positions of atoms and energy at each geometry"""
        self.filename = filename
        self.data = []
        self.energies = []
        self.positions = []
        self.magnetizations = []
        if not os.path.exists('POSCAR'):
            poscar = 'CONTCAR'
        else:
            poscar= 'POSCAR'
        self.poscar = PoscarParser(poscar)
        self.atom_count = self.poscar.number_of_atoms()

        with open(self.filename, 'r') as file:
            lines = file.readlines()
            lenght = len(lines)
            
            for i in range(lenght):
                if i % 10000 == 0:
                    print('reading OUTCAR file; line: ', i, f' out of {lenght}', end='\r')
                line = lines[i].strip()
                if line.startswith('POSITION'):
                    section_position = []
                    i += 2
                    current_i = i
                    for i in range(current_i, current_i + self.atom_count):
                        line = lines[i]
                        position_data = line.split()
                        atom_data = [float(x) for x in position_data[:3]]
                        section_position.append(atom_data)
                    self.positions.append(section_position)
                elif line.startswith('FREE ENERGIE'):
                    i += 2
                    energy_data = lines[i].strip().split()[4]
                    self.energies.append(float(energy_data))
                elif line.startswith('magnetization (x)'):
                    section_magnetization = []
                    i += 4
                    current_i = i
                    for i in range(current_i, current_i + self.atom_count):
                        line = lines[i]
                        magnetization_data = line.split()
                        atom_mag = float(magnetization_data[-1])
                        section_magnetization.append(atom_mag)
                    self.magnetizations.append(section_magnetization)
                if line.startswith('Voluntary context switches:'):
                    self.magnetizations.pop()
            if self.positions==[]:
                for i in range(lenght):
                    line = lines[i].strip()
                    if line.startswith('position of ions in cartesian'):
                        i += 1
                        current_i = i
                        for i in range(current_i, current_i + self.atom_count):
                            line = lines[i]
                            position_data = line.split()
                            atom_data = [float(x) for x in position_data[:3]]
                            section_position.append(atom_data)
                        self.positions.append(section_position)
                self.magnetizations = [["N/A" for atom in self.atom_count]]
                

        print('\n')

    def find_coordinates(self):
        """returns coordinates of each electronically converged calculation step"""
        return self.positions

    def find_energy(self):
        """returns converged energy in eV"""
        return self.energies
    
    def print_magnetization(self):
        return self.magnetizations
    
    def find_magnetization(self):
        search_string = 'magnetization'
        with open(self.filename, 'r') as file:
            file.seek(0, 2)
            file.seek(file.tell() - 6500, 0)
            file.readline()
            lines = file.readlines()
            file.seek(file.tell() - 6500, 0)
            file.readline()
            i=0
            lines_mag = []
            while i < len(lines):
                i+=1
                line = file.readline()
                if search_string in line:
                    j = 0
                    while j < 10:
                        j+=1
                        line=file.readline()
                        line=line.strip()
                        if line.startswith('1'):
                            lines_mag = file.readlines()
                            lines_mag=[line.strip() for line in lines_mag]
                            break
                    break

        lines_mag.insert(0,line)
        lines_mag=[el.split() for el in lines_mag]
        lines_mag=lines_mag[:self.poscar.number_of_atoms()]
        mag_values=[lst[-1] for lst in lines_mag]
        return mag_values


class PoscarParser:
    """class to parse POSCAR / CONTCAR files"""

    def __init__(self, filename):
        self.filename = filename
        self.atom_symbols_exists = None
        self.dynamic_exists = None
        self.total_atoms = 0
        with open(self.filename, 'r') as file:
            self.lines = file.readlines()
            self.scale = self.scale_factor()
            self.unit_cell_vectors()
            self.atomic_symbols()
            self.dynamics()

    def title(self):
        title = self.lines[0].strip()
        return title

    def scale_factor(self):
        scale_factor = float(self.lines[1].strip())
        return scale_factor

    def unit_cell_vectors(self):
        unit_cell_vectors = []
        for line in self.lines[2:5]:
            unit_cell_vectors.append([float(value)*self.scale_factor() for value in line.split()])
        return unit_cell_vectors

    @staticmethod
    def is_integer(string):
        try:
            int(string)
            return True
        except ValueError:
            return False

    def atomic_symbols(self):
        atom_symbols = []
        if PoscarParser.is_integer(self.lines[5].split()[0]):
            self.atom_symbols_exists = False
            if os.path.exists('POTCAR'):
                potcar = 'POTCAR'
            elif os.path.exists('../POTCAR'):
                potcar = '../POTCAR'
            else:
                raise FileNotFoundError('Error! No POTCAR file found!')
            with open(potcar, 'r') as file:
                lines = file.readlines()
                lenght = len(lines)
                first = lines[0].split()
                atom_symbols.append(first[1])

                for i in range(1, lenght - 10):
                    if lines[i].strip().startswith('End of Dataset'):
                        i += 1
                        method_line = lines[i].split()
                        atom_symbols.append(method_line[1])
        else:
            self.atom_symbols_exists = True
            atom_symbols = self.lines[5].split()
        return atom_symbols

    def list_atomic_symbols(self):
        symbol_list = [s for s, c in zip(self.atomic_symbols(), self.atom_counts()) for _ in range(c)]
        return symbol_list

    def atom_counts(self):
        if self.atom_symbols_exists:
            counts = [int(value) for value in self.lines[6].split()]
        else:
            counts = [int(value) for value in self.lines[5].split()]
        return counts

    def number_of_atoms(self):
        self.total_atoms = sum(self.atom_counts())
        return sum(self.atom_counts())

    def symbol_and_number(self):
        sym_num_list = []
        for symbol, number in zip(self.list_atomic_symbols(), range(1, self.number_of_atoms() + 1)):
            sym_num_list.append(str(symbol) + str(number))
        return sym_num_list

    def dynamics(self):
        if self.atom_symbols_exists and self.lines[7].strip()[0].lower() == "s":
            self.dynamic_exists = True
            return self.lines[7].strip()
        elif self.atom_symbols_exists is False and self.lines[6].strip()[0].lower() == "s":
            self.dynamic_exists = True
            return self.lines[6].strip()
        else:
            self.dynamic_exists = False
            return 'no dynamics'

    def coordinate_type(self):
        if self.atom_symbols_exists and self.dynamic_exists:
            return self.lines[8].strip()
        elif self.atom_symbols_exists and self.dynamic_exists is False:
            return self.lines[7].strip()
        elif self.atom_symbols_exists is False and self.dynamic_exists:
            return self.lines[7].strip()
        elif self.atom_symbols_exists is False and self.dynamic_exists is False:
            return self.lines[6].strip()

    def parse_coordinates(self):
        start_line = None
        if self.atom_symbols_exists and self.dynamic_exists:
            start_line = 9
        elif self.atom_symbols_exists and self.dynamic_exists is False:
            start_line = 8
        elif self.atom_symbols_exists is False and self.dynamic_exists:
            start_line = 8
        elif self.atom_symbols_exists is False and self.dynamic_exists is False:
            start_line = 7
        coordinates = []
        constrain = []
        for line in self.lines[start_line:start_line + self.number_of_atoms()]:
            values = line.split()
            coordinates.append([float(value) for value in values[:3]])
            if len(values) == 6:
                constrain.append(values[3])
            else: 
                constrain.append('n/a')
        if self.coordinate_type() == "Direct":  # convert from direct to cartesian
            coords_cart = []
            for coor in coordinates:
                x = coor[0] * self.unit_cell_vectors()[0][0]
                y = coor[1] * self.unit_cell_vectors()[1][1]
                z = coor[2] * self.unit_cell_vectors()[2][2]
                coords_cart.append([x, y, z])
            return coords_cart, constrain
        else:
            return coordinates, constrain

    def coordinates(self):
        return self.parse_coordinates()[0]

    def constrains(self):
        constrain = self.parse_coordinates()[1]
        if len(constrain[0]) == 0:
            return ['n/a']*self.total_atoms
        else:
            return self.parse_coordinates()[1]

if __name__ == "__main__":
    outcar = OutcarParser('OUTCAR')
    print(outcar.print_magnetization()[-1])
    print('done')
