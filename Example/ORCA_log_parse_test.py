from ParseOrca import ParseOrca

if __name__ == "__main__":
    POproc = ParseOrca("./H2O.log")
    POproc.generate_report("./H2O_report.json", "H2O")