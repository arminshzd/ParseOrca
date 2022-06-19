# ParseOrca: Python tool to parse ORCA output files

ParseOrca is a small Python script that reads ORCA DFT optimization input files and returns useful information for furthur calculations in Python scripts.

## Usage

Create a ParseOrca object and pass the path to the ORCA log file as an input. The ParseOrca class has the following methods:

- `successful_termination()`: Returns `True` if calculation was successful and `False` otherwise.
- `get_optimized_coords()`: Returns a numpy array of the optimized structure in Angestroms.
- `get_frequcies()`: Returns a list of frequencies. If no frequency caculation is found, returns an empty list.
- `get_energies()`: Returns a tuple of electronic energy, enthalpy, entropy*temperature, and Gibbs free energy in Hartrees (in this order).
-`generate_report(out_dir, dict_key)`: Generates a report from the log file in the form of an dictionary with `dict_key` as the key and saves it as a JSON at `out_dir` directory.

`Example` directory contains a test case for a single water molecule calculation. 