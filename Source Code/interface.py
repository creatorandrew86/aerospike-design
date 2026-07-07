import matplotlib
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, font


class Interface:
    def __init__(self, root: tk, results: dict, errors: list[str], on_solve: callable, on_save: callable):
        self.root = root
        self.on_solve = on_solve
        self.on_save = on_save
        self.results = results
        self.errors = errors
        self.build_interface()


    def get_inputs_on_solve(self) -> dict:
        return {
            "oxidizer": self.input_oxidizer.get()    if self.input_oxidizer.get() else None,
            "fuel":     self.input_fuel.get()        if self.input_fuel.get() else None,
            "Pc":       float(self.input_Pc.get())   if self.input_Pc.get() else None,
            "unit_Pc":  self.unit_Pc.get()           if self.unit_Pc.get() else None,  
            "MR":       float(self.input_MR.get())   if self.input_MR.get() else None,
            "eps":      float(self.input_eps.get())  if self.input_eps.get() else None,

            "effective_throat_eps": float(self.input_effective_throat_eps.get())  if self.input_effective_throat_eps.get() else None,
            "truncate_percent":     float(self.input_truncate_percentage.get())   if self.input_truncate_percentage.get() else None,
            "aerospike_resolution": int(self.input_aerospike_resolution.get())    if self.input_aerospike_resolution.get() else None,
        }


    def build_interface(self):
        UNDERLINE_FONT = font.Font(family="Arial", size=9, underline=1)
        OXIDIZER_ITEMS = ["LOX", "GOX", "N2O4", "N2O", "IRFNA", "H2O2", "Peroxide90", "Peroxide98", "MON3", "MON15", "MON25"]
        FUEL_ITEMS = ["RP1", "LH2", "CH4", "MMH", "N2H4", "UDMH", "A50", "Ethanol", "Methanol", "GH2", "GCH4", "JetA", "JP10"]

        ttk.Label(self.root, text="Engine definition", font=UNDERLINE_FONT).grid(column=0, row=0, padx=7, pady=4, sticky='W')

        # Oxidizer
        ttk.Label(self.root, text='Oxidizer:').grid(column=0, row=1, padx=2, sticky='W')
        self.input_oxidizer = ttk.Combobox(self.root, values=OXIDIZER_ITEMS, width=15)
        self.input_oxidizer.grid(column=1, row=1, padx=5, pady=1, sticky='W')
        # Fuel
        ttk.Label(self.root, text='Fuel:').grid(column=0, row=2, padx=2, sticky='W')
        self.input_fuel = ttk.Combobox(self.root, values=FUEL_ITEMS, width=15)
        self.input_fuel.grid(column=1, row=2, padx=5, pady=1, sticky='W')

        # Chamber Pressure
        ttk.Label(self.root, text='Chamber Pressure:').grid(column=0, row=3, padx=2, sticky='W')
        self.input_Pc = ttk.Entry(self.root, width=8)
        self.input_Pc.grid(column=1, row=3, padx=5, pady=1, sticky='W')

        self.unit_Pc = ttk.Combobox(self.root, values=["bar", "Pa", "atm", "kPa", "psi", "MPa"], width=4)
        self.unit_Pc.grid(column=1, row=3, padx=70, pady=1, sticky='W')
        self.unit_Pc.set("bar")

        # Mixture Ratio
        ttk.Label(self.root, text='Mixture Ratio:').grid(column=0, row=4, padx=2, sticky='W')
        self.input_MR = ttk.Entry(self.root, width=18)
        self.input_MR.grid(column=1, row=4, padx=5, pady=1, sticky='W')
        # Area Ratio
        ttk.Label(self.root, text='Design Area Ratio:').grid(column=0, row=5, padx=2, pady=5, sticky='W')
        self.input_eps = ttk.Entry(self.root, width=18)
        self.input_eps.grid(column=1, row=5, padx=5, pady=1, sticky='W')




        ttk.Label(self.root, text='Aerospike Definition', font=UNDERLINE_FONT).grid(column=0, row=6, padx=7, pady=4, sticky='W')

        # Effective Throat Area Ratio
        ttk.Label(self.root, text='Effective Throat Area Ratio:').grid(column=0, row=7, padx=2, sticky='W')
        self.input_effective_throat_eps = ttk.Entry(self.root, width=18)
        self.input_effective_throat_eps.grid(column=1, row=7, padx=5, pady=1, sticky='W')

        # Truncation
        ttk.Label(self.root, text='Truncate at (%):').grid(column=0, row=8, padx=2, sticky='W')
        self.input_truncate_percentage = ttk.Entry(self.root, width=18)
        self.input_truncate_percentage.grid(column=1, row=8, padx=5, pady=1, sticky='W')

        # Resolution
        ttk.Label(self.root, text='Aerospike Resolution:').grid(column=0, row=9, padx=2, sticky='W')
        self.input_aerospike_resolution = ttk.Entry(self.root, width=18)
        self.input_aerospike_resolution.grid(column=1, row=9, padx=5, pady=1, sticky='W')

        # Solver
        ttk.Button(self.root, text='Run', command=self.on_solve).grid(column=0, row=11, pady=7)
        ttk.Button(self.root, text='Save Contour', command=self.on_save).grid(column=1, row=11, pady=7)


    def show_results(self):
        import matplotlib
        matplotlib.use("TkAgg")
        import matplotlib.pyplot as plt
        
        x, y = self.results["x"], self.results["R_x"]

        plt.figure()
        plt.plot(y, x, color='black')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('x')
        plt.grid()
        plt.show()


    def show_errors(self):
        for error in self.errors:
            messagebox.showerror("Error", error)


    def save_file(self):
        try:
            x, R, eps = self.results["x"], self.results["R_x"], self.results["eps_x"]
        except Exception:
            messagebox.showerror("Error", "You must run the program before attempting to save the results.")
            return

        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".dat", filetypes=[("DAT files", "*.dat")])
            if not file_path:
                return
            
            # Open the file and write the points
            with open(file_path, 'w') as file:
                for x, y, eps in zip(x, R, eps):
                    file.write(f"{x} {y} {eps}\n")
            
            messagebox.showinfo("File Saved", "File Successfully Saved")
        
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while saving the file: {e}")
            return