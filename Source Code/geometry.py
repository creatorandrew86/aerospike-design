from rocketcea.cea_obj_w_units import CEA_Obj
from scipy.optimize import brentq
import numpy as np

# Helper functions
def _prandtl_meyer_angle(mach: float, gamma: float) -> float:
    import numpy as np

    mach = max(mach, 1)
    pm_angle = np.sqrt((gamma + 1)/(gamma - 1)) * np.arctan(np.sqrt((gamma - 1)/(gamma + 1) * (mach**2 - 1))) - np.arctan(np.sqrt(mach**2 - 1))
    return pm_angle

def _eps_from_mach(mach: float, Pc: float, MR: float, c: CEA_Obj) -> float:
    mach = max(mach, 1 + 1e-6)

    def residual(eps):
        guess_mach = c.get_MachNumber(Pc=Pc, MR=MR, eps=eps)
        return mach - guess_mach
    
    # brentq with lower boud eps = 1, upper bound eps = 1000
    return brentq(residual, 1, 1000)

def _mach_angle(mach: float) -> float:
    import numpy as np

    mach = max(mach, 1)
    return np.arcsin(1 / mach)



def generate_aerospike_contour(inputs: dict) -> tuple[dict, list[str]]:
    errors = []

    oxidizer = inputs["oxidizer"]
    fuel     = inputs["fuel"]
    MR       = inputs["MR"]
    exit_eps = inputs["eps"]
    Pc       = inputs["Pc"]

    effective_throat_eps = inputs["effective_throat_eps"]
    truncate_percent     = inputs["truncate_percent"] / 100
    N                    = inputs["aerospike_resolution"]

    sizing         = inputs["sizing"]
    aerospike_type = inputs["aerospike_type"]
    mfr            = inputs["mfr"]
    radius         = inputs["radius"]
    width          = inputs["width"]
    length         = inputs["length"]


    results = {
        "aerospike_type": aerospike_type,
        "mfr": None,
        "radius": None,
        "length": None,
        "width": None,
        "x": [None] * N,
        "R_x": [None] * N,
        "eps_x": [None] * N,
    }

    
    # CEA Object
    c = CEA_Obj(oxName=oxidizer, fuelName=fuel, pressure_units="Pa", sonic_velocity_units="m/s", density_units="kg/m^3", temperature_units="K")
    

    # Aerospike Throat Mach
    try:
        entry_mach = max(c.get_MachNumber(Pc=Pc, MR=MR, eps=effective_throat_eps), 1 + 1e-6)
    except Exception as e:
        errors.append(f"Entry Mach number calculation failed with {e}")
        return {}, errors
    

    # Aerospike Exit Conditions
    try:
        exit_mach  = max(c.get_MachNumber(Pc=Pc, MR=MR, eps=exit_eps), 1 + 1e-6)
        exit_gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=exit_eps)[1]
        exit_pm_angle = _prandtl_meyer_angle(exit_mach, exit_gamma)

    except Exception as e:
        errors.append(f"Exit conditions calculation failed with {e} ")
        return {}, errors


    # Loop - generates contour assuming throat_length = 1
    for i in range(N):
        mach = entry_mach + (i / N) * (exit_mach - entry_mach)

        # Station Values
        try:
            # Angles
            eps = _eps_from_mach(mach, Pc, MR, c)
            gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=eps)[1]
            mach_angle = _mach_angle(mach)

            pm_angle   = _prandtl_meyer_angle(mach, gamma)
            flow_angle = mach_angle + (exit_pm_angle - pm_angle)

            # Length and coordinates
            characteristic_length = mach * eps

            x = characteristic_length * -np.cos(flow_angle)
            R = characteristic_length * np.sin(flow_angle)

        except Exception as e:
            errors.append(f"Station {i+1} calculation failed with {e}")
            return {}, errors


        results["x"][i]     = x
        results["R_x"][i]   = R
        results["eps_x"][i] = eps


    # Non-dimensionalization of coordinates
    if not sizing:
        max_length = max(np.max(np.abs(results["x"])), np.max(np.abs(results["R_x"])))
        
        results["x"]   = [v / max_length for v in results["x"]]
        results["R_x"] = [v / max_length for v in results["R_x"]]


    # Sizing
    if sizing:
        if aerospike_type == "toroidal":
            if mfr is None:
                # Radius + eps -> throat area -> mass flow
                throat_length = radius / np.max(np.abs(results["R_x"]))

                results["x"]   = [v * throat_length for v in results["x"]]
                results["R_x"] = [v * throat_length for v in results["R_x"]]

                try:
                    mfr = _get_mass_flow_toroidal(radius, c, Pc, MR, exit_eps, effective_throat_eps)
                except Exception as e:
                    errors.append(f"Mass flow calculation failed with {e}")
                    return {}, errors
                

            if radius is None:
                # Mass flow -> throat area + area ratio -> exit area -> radius
                try:
                    radius = _get_radius_toroidal(mfr, c, Pc, MR, exit_eps, effective_throat_eps)
                except Exception as e:
                    errors.append(f"Cowl radius calculation failed with {e}")
                    return {}, errors

                throat_length = radius / np.max(np.abs(results["R_x"]))

                results["x"]   = [v * throat_length for v in results["x"]]
                results["R_x"] = [v * throat_length for v in results["R_x"]]

            results["mfr"]    = mfr
            results["radius"] = radius

        
        if aerospike_type == "linear":
            if mfr is None:
                # Width (radius equivalent) + area ratio -> throat area -> mass flow
                throat_length = width / np.max(np.abs(results["R_x"]))

                results["x"]   = [v * throat_length for v in results["x"]]
                results["R_x"] = [v * throat_length for v in results["R_x"]]

                try:
                    mfr = _get_mass_flow_linear(length, width, c, Pc, MR, exit_eps, effective_throat_eps)
                except Exception as e:
                    errors.append(f"Mass flow calculation failed with {e}")
                    return {}, errors
                

            if length is None:
                # Width + mfr -> throat area + area ratio -> exit area -> length
                throat_length = width / np.max(np.abs(results["R_x"]))

                results["x"]   = [v * throat_length for v in results["x"]]
                results["R_x"] = [v * throat_length for v in results["R_x"]]

                try:
                    length = _get_length_linear(mfr, width, c, Pc, MR, exit_eps, effective_throat_eps)
                except Exception as e:
                    errors.append(f"Aerospike length calculation failed with {e}")
                    return {}, errors
                
            
            if width is None:
                # Mass flow + length -> throat area + area ratio -> exit area -> width
                try:
                    width = _get_width_linear(mfr, length, c, Pc, MR, exit_eps, effective_throat_eps)
                except Exception as e:
                    errors.append(f"Width calculation failed with {e}")
                    return {}, errors

                throat_length = width / np.max(np.abs(results["R_x"]))

                results["x"]   = [v * throat_length for v in results["x"]]
                results["R_x"] = [v * throat_length for v in results["R_x"]]

            results["mfr"]    = mfr
            results["length"] = length
            results["width"]  = width

    
    # Truncating
    results["x"] = results["x"][:int(len(results["x"]) * truncate_percent)]
    results["R_x"] = results["R_x"][:int(len(results["R_x"]) * truncate_percent)]
    results["eps_x"] = results["eps_x"][:int(len(results["eps_x"]) * truncate_percent)]

    return results, errors




# Torroidal aerospike values
def _get_mass_flow_toroidal(radius: float, c: CEA_Obj, Pc: float, MR: float, exit_eps: float, throat_eps: float) -> float:
    if throat_eps == 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=exit_eps)[1]
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=exit_eps)[1]

    if throat_eps > 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=throat_eps)[2] * c.get_MachNumber(Pc=Pc, MR=MR, eps=throat_eps)
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=throat_eps)[2]
                    

    throat_area = np.pi * pow(radius, 2) / (exit_eps / throat_eps)
    mfr = throat_area * effective_throat_velocity * effective_throat_rho

    return mfr

def _get_radius_toroidal(mfr: float, c: CEA_Obj, Pc: float, MR: float, exit_eps: float, throat_eps: float) -> float:
    if throat_eps == 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=exit_eps)[1]
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=exit_eps)[1]

    if throat_eps > 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=throat_eps)[2] * c.get_MachNumber(Pc=Pc, MR=MR, eps=throat_eps)
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=throat_eps)[2]

    throat_area = mfr / (effective_throat_velocity * effective_throat_rho)
    exit_area = throat_area * (exit_eps / throat_eps)
    radius = np.sqrt(exit_area / np.pi)

    return radius




# Linear aerospike values
def _get_mass_flow_linear(length: float, width: float, c: CEA_Obj, Pc: float, MR: float, exit_eps: float, throat_eps: float) -> float:
    if throat_eps == 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=exit_eps)[1]
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=exit_eps)[1]

    if throat_eps > 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=throat_eps)[2] * c.get_MachNumber(Pc=Pc, MR=MR, eps=throat_eps)
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=throat_eps)[2]

    throat_area = length * width / (exit_eps / throat_eps)
    mfr = throat_area * effective_throat_velocity * effective_throat_rho

    return mfr

def _get_length_linear(mfr: float, width: float, c: CEA_Obj, Pc: float, MR: float, exit_eps: float, throat_eps: float) -> float:
    if throat_eps == 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=exit_eps)[1]
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=exit_eps)[1]

    if throat_eps > 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=throat_eps)[2] * c.get_MachNumber(Pc=Pc, MR=MR, eps=throat_eps)
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=throat_eps)[2]

    throat_area = mfr / (effective_throat_velocity * effective_throat_rho)
    exit_area = throat_area * (exit_eps / throat_eps)
    length = exit_area / width

    return length

def _get_width_linear(mfr: float, length: float, c: CEA_Obj, Pc: float, MR: float, exit_eps: float, throat_eps: float) -> float:
    if throat_eps == 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=exit_eps)[1]
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=exit_eps)[1]

    if throat_eps > 1:
        effective_throat_velocity = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=throat_eps)[2] * c.get_MachNumber(Pc=Pc, MR=MR, eps=throat_eps)
        effective_throat_rho      = c.get_Densities(Pc=Pc, MR=MR, eps=throat_eps)[2]

    throat_area = mfr / (effective_throat_velocity * effective_throat_rho)
    exit_area = throat_area * (exit_eps / throat_eps)
    width = exit_area / length

    return width