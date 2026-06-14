import math
import tkinter as tk
from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt
from tkinter import ttk, messagebox, filedialog

global R_show, x_show, eps_dat

'''

    **********************************************************************************************************************************
                                                                    MAIN FUNCTION
    **********************************************************************************************************************************
    
'''

def aerospike():
    global R_show, x_show, eps_dat

    
    def P_M_angle(mach, gam):  # P-M angle (Prandtl-Meyer angle)
        term1 = math.sqrt((gam - 1) / (gam + 1) * (mach ** 2 - 1))  # First term
        term2 = math.sqrt((gam + 1) / (gam - 1)) * math.atan(term1)  # Second Term
        term3 = math.atan(math.sqrt(mach ** 2 - 1))  # Third term
        nu = term2 - term3  # P-M angle in radians
        return math.degrees(nu)

    #Function to get Mach number from Prandtl-Meyer angle
    def get_mach(gam, nu):
        
        mach_guess = 1
        while True:
            nu_guess = P_M_angle(mach_guess, gam)
            if math.fabs(nu_guess - nu) < 0.7:
                break
            mach_guess += 0.01
        return mach_guess


    #Take inputs from the interface
    Pc = float(chamber_pressure.get())              #Chamber pressure
    MR = float(mixture_ratio.get())                 #Mixture Ratio
    theta_throat = float(throat_angle.get())        #Throat angle
    eps_exit = float(eps_in.get())                  #Exit Area Ratio
    mfr = float(mfri.get())                         #Nozzle Design Mass flux
    
    if float(truncate_val.get()) == 0:
        truncate = 1
    else:
        truncate = float(truncate_val.get())
    
    Ox = ox_box.get()                               #Oxidizer
    Fuel = fuel_box.get()                           #Fuel


    #Setup the CEA Case
    c = CEA_Obj(oxName = Ox, fuelName = Fuel, temperature_units = 'K',
                pressure_units = 'bar', sonic_velocity_units='m/s', density_units = 'kg/m^3')

    #Function to get epsilon(Area Ratio) from Mach number
    def get_eps(M):
        
        eps_guess = 1
        while True:
            mach_guess = c.get_MachNumber(Pc = Pc, MR = MR, eps = eps_guess)
            if math.fabs(mach_guess - M) < 0.1:
                break
            eps_guess += 0.005
        return mach_guess

    
    #Process the output from CEA
    rho_th = c.get_Densities(Pc = Pc, MR = MR, eps = 2)[1]          #Gas density at the throat (kg/m^3)
    a_th = c.get_SonicVelocities(Pc = Pc, MR = MR, eps = 2)[1]      #Gas sound velocity at the throat (m/s)
    M_exit = c.get_MachNumber(Pc = Pc, MR = MR, eps = eps_exit)     #Gas Mach number at the nozzle exit
    gam = c.get_Chamber_MolWt_gamma(Pc = Pc, MR = MR)[1]            #Gas gamma value at combustion

    
    #Process geometric parameters
    area_throat = (mfr/(a_th*rho_th))*10000          #!!! Throat is different(M = 1) || cm^2
    throat_area = round(area_throat, 2)
    

    #Note: The "Aerospike throat" is the point at which the aerospike begins the expansion, the minimum Mach number at the
    #aerospike throat is M=1, but it can be higher.
    
    #P-M angles at exit and aerospike throat
    nu_exit = P_M_angle(M_exit, gam)
    nu_throat = nu_exit - theta_throat

    M_throat = get_mach(gam, nu_throat)                             #Aerospike throat Mach number
    eps_throat = get_eps(M_throat)                                  #Aerospike throat area ratio
    Re = round(math.sqrt(eps_exit*area_throat/math.pi), 2)          #Aerospike radius

    #Initialize lists, "R" means radius, "x" means length, and "eps" means epsilon(Area Ratio)
    R_list, R_show, R_dot = [], [], []
    x_list, x_show, x_dot = [], [], []
    eps_dat, eps_list = [], []




    max_it = 100

    #Iteration is done over epsilon
    d_eps = (eps_exit - eps_throat)/max_it
    
    for i in range(max_it):

        #Get epsilon and Mach number values at the current analysis point
        eps = eps_throat + i*d_eps
        M = c.get_MachNumber(Pc = Pc, MR = MR, eps = eps)

        #Get the P-M angle and Mach angle
        nu_x = P_M_angle(M, gam)
        mu_x = math.degrees(math.asin(1/M))
        term_0 = eps*area_throat/math.pi

        #Calculate the radius and length of the current analysis point
        R_x = math.sqrt(pow(Re, 2) - term_0)
        x = (Re - R_x)/math.tan(math.radians(nu_exit - nu_x + mu_x))

        #Append the radius, epsilon and length to lists
        R_list.append(R_x)
        x_list.append(x)
        eps_list.append(eps)

    #Define the values for the truncation
    for i in range(len(R_list)):
        if x_list[i] < x_list[-1]*truncate/100:
            x_show.append(x_list[i])
            R_show.append(R_list[i])
            eps_dat.append(eps_list[i])

        elif x_list[i] >= x_list[-1]*truncate/100:
            x_dot.append(x_list[i])
            R_dot.append(R_list[i])

    #Calculate general engine parameters
    Pressure_ratio = c.get_PcOvPe(Pc = Pc, MR = MR, eps = eps_exit)
    Isp = round(c.estimate_Ambient_Isp(Pc = Pc, MR = MR, eps = eps_exit, Pamb = Pc/Pressure_ratio)[0], 2)
    length = round(float(x_list[-1]), 2)
    truncated_length = round(float(x_list[-1]) * truncate/100, 2)

    #Calculate Output Values
    v_eff = Isp * 9.81
    T = round(v_eff*mfr/1000, 2)

    '''

    **********************************************************************************************************************************
                                                              MAIN FUNCTION OUTPUT
    **********************************************************************************************************************************
    
    '''

    win = tk.Tk()
    win.title('Results')
    win.geometry('250x250')

    ttk.Label(win, text='Geometric Parameters').grid(column=0, row=0, padx=5, pady=5)

    ttk.Label(win, text=f'Radius: {Re} cm').grid(column=0, row=1, padx=2, pady=2, sticky='W')
    ttk.Label(win, text=f'Length: {length} cm').grid(column=0, row=2, padx=2, pady=2, sticky='W')
    ttk.Label(win, text=f'Truncated Length: {truncated_length} cm').grid(column=0, row=3, padx=2, pady=2, sticky='W')
    ttk.Label(win, text=f'Throat Area: {throat_area} cm\u00B2').grid(column=0, row=4, padx=2, pady=2, sticky='W')

    ttk.Label(win, text='Engine Parameters').grid(column=0, row=5, padx=5, pady=5)
    
    ttk.Label(win, text=f'Isp: {Isp} s').grid(column=0, row=6, padx=2, pady=2, sticky='W')
    ttk.Label(win, text=f'Thrust: {T} kN').grid(column=0, row=7, padx=2, pady=2, sticky='W')
    ttk.Label(win, text=f'Aen/At: {eps_throat}').grid(column=0, row=8, padx=2, pady=2, sticky='W')
    

    plt.figure(figsize=(10, 5))
    plt.plot(x_show, R_show, linestyle='-', color='black', label='Aerospike')
    plt.plot(x_dot, R_dot, linestyle='--', color='black', label='Full')
    plt.title('Aerospike contour')
    plt.xlabel('x(cm)')
    plt.ylabel('R(cm)')
    plt.grid(False)
    plt.show()





'''

    **********************************************************************************************************************************
                                                                    PROGRAM OUTPUT
    **********************************************************************************************************************************
    
'''    

def save_file():
    try:
        file_path = filedialog.asksaveasfilename(defaultextension=".dat", filetypes=[("DAT files", "*.dat")])
        if not file_path:
            return  # User cancelled the save dialog
        
        # Open the file and write x_points and y_points
        with open(file_path, 'w') as file:
            for x, y, eps in zip(x_show, R_show, eps_dat):
                file.write(f"{x} {y} {eps}\n")
        
        messagebox.showinfo("Success", "Points saved successfully to the .dat file.")
    
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while saving the file: {e}")    





'''

    **********************************************************************************************************************************
                                                                    USER INTERFACE
    **********************************************************************************************************************************
    
'''
ox_list = ['LOX', 'H2O2']
fuel_list = ['H2', 'C2H6', 'CH4']

#Initalize the window
root = tk.Tk()
root.title('Aerospike Nozzle Geometry')
root.geometry('350x350')


''' ********************************************************* CEA INPUTS ********************************************************* '''

ttk.Label(root, text='CEA Inputs').grid(column=0, row=0, padx=5, pady=2)


#Oxidizer Input
ttk.Label(root, text='Oxidizer:').grid(column=0, row=1, padx=5, sticky='W')
ox_box = ttk.Combobox(root, values=ox_list)
ox_box.grid(column=1, row=1, padx=5, pady=5)
ox_box.set('LOX')

#Fuel Input
ttk.Label(root, text='Fuel:').grid(column=0, row=2, padx=5, sticky='W')
fuel_box = ttk.Combobox(root, values=fuel_list)
fuel_box.grid(column=1, row=2, padx=5, pady=5)
fuel_box.set('H2')


#Chamber Pressure input
ttk.Label(root, text='Chamber Pressure (bar):').grid(column=0, row=3, padx=5, sticky='W')
chamber_pressure = ttk.Entry(root)
chamber_pressure.grid(column=1, row=3, padx=5, pady=5 , sticky='W')
chamber_pressure.insert(0, '120')

#Mixture Ratio input
ttk.Label(root, text='Mixture Ratio:').grid(column=0, row=4, padx=5, sticky='W')
mixture_ratio = ttk.Entry(root)
mixture_ratio.grid(column=1, row=4, padx=5, pady=5, sticky='W')
mixture_ratio.insert(0, '6')


''' **************************************************** AEROSPIKE PARAMETERS ****************************************************'''

ttk.Label(root, text='Aerospike Parameters').grid(column=0, columnspan=2, row=5, padx=10, pady=7, sticky='W')

#Nozzle Mass Flux Input
ttk.Label(root, text='Nozzle Mass Flow:').grid(column=0, row=6, padx=5, pady=3, sticky='W')
mfri = ttk.Entry(root)
mfri.insert(0, '4')
mfri.grid(column=1, row=6, padx=5, pady=5, sticky='W')

#Design epsilon input
ttk.Label(root, text='Design area ratio:').grid(column=0, row=7, padx=5, pady=3, sticky='W')
eps_in = ttk.Entry(root)
eps_in.insert(0, '270')
eps_in.grid(column=1, row=7, padx=5, pady=5, sticky='W')

#Truncation percentage across length
ttk.Label(root, text='Truncate at(%):').grid(column=0, row=8, padx=5, pady=3, sticky='W')
truncate_val = ttk.Entry(root)
truncate_val.insert(0, '40')
truncate_val.grid(column=1, row=8, padx=5, pady=5, sticky='W')

#Throat angle input
ttk.Label(root, text='Throat Angle (deg):').grid(column=0, row=9, padx=5, sticky='W')
throat_angle = ttk.Entry(root)
throat_angle.grid(column=1, row=9, padx=5, pady=5, sticky='W')
throat_angle.insert(0, '70')

''' ************************************************** SOLVER BUTTONS ************************************************** '''

#Main Function button
calculate = ttk.Button(root, text='Generate', command=aerospike)
calculate.grid(column=0, row=10, padx=5, pady=5, sticky='W')

#File generation button
save_file = ttk.Button(root, text='Get Aerospike .dat file', command = save_file)
save_file.grid(column=1, row=10, padx=5, pady=5, sticky='W')


#Run the Tkinter event loop
root.mainloop()

