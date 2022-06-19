from ParseOrca import ParseOrca

if __name__ == "__main__":
    POproc = ParseOrca("/Users/arminsh/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Github_Repo/Parse_ORCA/Example/H2O.log")
    POproc.generate_report("/Users/arminsh/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Github_Repo/Parse_ORCA/Example/H2O_report.json", "H2O")