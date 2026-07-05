import tkinter as tk
import sys, os

import solver, interface

def _resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


def on_solve():
    inputs = ui.get_inputs_on_solve()
    
    results, solver_errors = solver.generate_aerospike_contour(inputs)

    ui.results, ui.errors = results, solver_errors
    if solver_errors:
        ui.show_errors()
        return
    
    ui.show_results()
    

def on_save():
    ui.save_file()


if __name__ == "__main__":
    root = tk.Tk()
    root.title('Aerospike Nozzle Geometry')
    root.geometry('320x370')
    root.resizable(False, True)

    icon_path = _resource_path(os.path.join("assets", "icon.ico"))
    root.iconbitmap(icon_path)

    ui = interface.Interface(root=root, on_solve=on_solve, on_save=on_save, results=None, errors=None)

    root.mainloop()