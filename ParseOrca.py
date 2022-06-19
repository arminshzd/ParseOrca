import json
import numpy as np
from pathlib import Path

class ParseOrca:
    def __init__(self, rep_dir) -> None:
        """Instantiate a ParseOrca object with report at rep_dir. Immediately checks whether the calculations where successful.
        Developed for ORCA 5.0.1.
        Args:
            rep_dir (Path): path to report file
        """
        if type(rep_dir) is not type(Path()):
            rep_dir = Path(rep_dir)

        with open(rep_dir, 'r') as f:
            self.rep_str = f.read()

        self.rep_cnt = self.rep_str.split('\n')

        self._get_n_atoms()
        self.success = self.successful_termination()

        self.ptable = {1 : "H", 2  : "He", 3  : "Li", 4  : "Be", 5  : "Be", 6  : "C", 7  : "N", 8  : "O", 9  : "F", 10 : "Ne", 11 : "Na",\
                       12 : "Mg", 13 : "Al", 14 : "Si", 15 : "P", 16 : "Si", 17 : "Cl", 18 : "Ar", 19 : "K", 20 : "Ca", 21 : "Sc", 22 : "Ti",\
                       23 : "V", 24 : "Cr", 25 : "Mn", 26 : "Fe", 27 : "Co", 28 : "Ni", 29 : "Cu", 30 : "Zn", 31 : "Ga", 32 : "Ge", 33 : "As",\
                       34 : "Se", 35 : "Br", 36 : "Kr", 37 : "Rb", 38 : "Sr", 39 : "Y", 40 : "Zr", 41 : "Nb", 42 : "Mo", 43 : "Tc", 44 : "Ru",\
                       45 : "Rh", 46 : "Pd", 47 : "Ag", 48 : "Cd", 49 : "In", 50 : "Sn", 51 : "Sb", 52 : "Te", 53 : "In", 54 : "Xe",  55 : "Cs",\
                       56 : "Ba", 57 : "La", 58 : "Ce", 59 : "Pr", 60 : "Nd", 61 : "Pm", 62 : "Sm", 63 : "Eu", 64 : "Gd", 65 : "Tb", 66 : "Dy",\
                       67 : "Ho", 68 : "Er", 69 : "Tm", 70 : "Yb", 71 : "Lu", 72 : "Hf", 73 : "Ta", 74 : "W ", 75 : "Re", 76 : "Os", 77 : "Ir",\
                       78 : "Pt", 79 : "Au", 80 : "Hg", 81 : "Tl", 82 : "Pb", 83 : "Bi", 84 : "Po", 85 : "At", 86 : "Rn"}
        self.ptable_r = {}
        for i, val in enumerate(list(self.ptable.values())):
            self.ptable_r[val] = list(self.ptable.keys())[i]
            
        self.final_coords = None
        self.final_struct = None
        self.optimized_coords = None
        self.optimized_struct = None
        self.atom_list_sym = None
        self.atom_list_num = None
    
    def _get_n_atoms(self):
        for line in self.rep_cnt:
            if "Number of atoms" in line:
                line_list = line.split()
                self.natoms = int(line_list[-1])
    
    def successful_termination(self):
        """Checks whether the optimization was successful. Searches whithin the whole report as a string.
        Returns:
            bool: True if successful False otherwise.
        """
        if ("HURRAY" in self.rep_str) and ("ERROR !!!" not in self.rep_str):
            return True
        else:
            return False

    def get_optimized_coords(self):
        """Extract the optimized coords from log file. ORCA's output is in au (Bohr)
        but the output is converted to Angstroms.
        Returns:
            np.array: optimized (or final if the optimization failed) coordinates.
        """
        bohr2ang = 0.529177
        
        match_term = "  #           XYZ [au]              \
        r0(AA) [Ang.]  CN      C6(AA)     C8(AA)    C10(AA) [au]"
        match_list = match_term.split()
        match_list.sort()
        for i, line in enumerate(self.rep_cnt[::-1]):
            line_list = line.split()
            line_list.sort()
            if line_list == match_list:
                start_line = i
                break
        
        structure = self.rep_cnt[-start_line:-start_line + self.natoms]
        structure = [i.split()[:-5] for i in structure]

        for atom in structure:
            atom[0] = int(atom[0])
            atom[1] = float(atom[1])*bohr2ang
            atom[2] = float(atom[2])*bohr2ang
            atom[3] = float(atom[3])*bohr2ang
            atom[4] = atom[4].title()
        
        coords = [atom[1:-1] for atom in structure]
        for i, atom in enumerate(structure):
            structure[i] = [atom[0], coords[i], atom[-1]]
        
        self.final_struct = structure
        self.final_coords = np.array(coords)
        self.atom_list_sym = [atom[-1] for atom in structure]
        self.atom_list_num = [self.ptable_r[atom[-1]] for atom in structure]

        if not self.success:
            print('WARNING! Optimization failed or did not converge. Latest coordinates will be returned instead.')
        else:
            self.optimized_struct = self.final_struct
            self.optimized_coords = self.final_coords

        return self.optimized_coords

    def get_frequencies(self):
        """Returns a list of frequencies. The 6 zero frequencies are included in the list.

        Returns:
            list: list of frequencies (cm**-1)
        """
        freq_key = "Scaling factor for frequencies"
        for i, line in enumerate(self.rep_cnt[::-1]):
            if freq_key in line:
                startline = i - 1
                break

        nmodes = self.natoms * 3
        endline = -startline + nmodes

        freqs_txt = self.rep_cnt[-startline: endline]
        freqs = []
        for line in freqs_txt:
            freqs.append(float(line.split()[1]))

        return freqs

    def get_energies(self):
        """Returns electronic energy, enthalpy, entropy*temperature, and Gibbs free energy from report

        Returns:
            Float: Electronic energy (Ha)
            Float: Enthalpy (Ha)
            Float: Entropy*Temperature (Ha)
            Float: Gibbs free energy (Ha)
        """
        
        G_key = "Final Gibbs free energy"
        TS_key = "Final entropy term"
        H_key = "Total Enthalpy"
        E_key = "Electronic energy"

        keys_found = 0
        for line in self.rep_cnt[::-1]:
            if E_key in line:
                keys_found += 1
                E = float(line.split()[-2])

            elif H_key in line:
                keys_found += 1
                H = float(line.split()[-2])
            
            elif TS_key in line:
                keys_found += 1
                TS = float(line.split()[-4])
            
            elif G_key in line:
                keys_found += 1
                G = float(line.split()[-2])
            
            if keys_found == 4:
                break
        
        return E, H, TS, G

    def generate_report(self, out_dir, dict_key):
        if type(out_dir) is not type(Path()):
            out_dir = Path(out_dir)

        opt_coords = self.get_optimized_coords().tolist()
        freqs = self.get_frequencies()
        E, H, TS, G = self.get_energies()
        results = {"Successful Job completion" :  self.successful_termination(),
        "Stationary Point Coordinates" :  opt_coords,
        "# Imaginary Frequencies" :  len([freq for freq in freqs if freq<0]),
        "Imaginary Frequencies (cm**-1)": [freq for freq in freqs if freq<0],
        "Electronic Energy (Ha)" :  E,
        "G (Ha)" :  G,
        "H (Ha)" : H,
        "TS (Ha)" : TS}
        
        if out_dir.exists():
            with open(out_dir, 'r') as f:
                species_results = json.load(f)
            species_results[dict_key] = results
        else:
            species_results = {dict_key: results}

        with open(out_dir, 'w') as f:
            json.dump(species_results, f, indent=2)
        
        return
