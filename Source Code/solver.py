from rocketcea.cea_obj_w_units import CEA_Obj
from scipy.optimize import brentq


# Unit Conversion
def _pressure_conversion(value: float, unit: str) -> float:
    match unit:
        case "Pa" :  return value
        case "kPa" : return value * 1e2
        case "MPa" : return value * 1e6
        case "bar" : return value * 1e5
        case "atm" : return value * 101325
        case "psi" : return value * 6894.76


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
    import numpy as np
    
    errors = []

    if any(value is None for _, value in inputs.items()):
        errors.append("Input Error. Review your inputs before proceeding.")
        return {}, errors

    oxidizer = inputs["oxidizer"]
    fuel     = inputs["fuel"]
    MR       = inputs["MR"]
    exit_eps = inputs["eps"]
    Pc       = _pressure_conversion(inputs["Pc"], inputs["unit_Pc"])

    effective_throat_eps = inputs["effective_throat_eps"]
    truncate_percent     = inputs["truncate_percent"] / 100
    N                    = inputs["aerospike_resolution"]



    results = {
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
        return errors
    

    # Aerospike Exit Conditions
    try:
        exit_mach  = max(c.get_MachNumber(Pc=Pc, MR=MR, eps=exit_eps), 1 + 1e-6)

        exit_gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=exit_eps)[1]
        exit_pm_angle = _prandtl_meyer_angle(exit_mach, exit_gamma)
        exit_mach_angle = _mach_angle(exit_mach)

    except Exception as e:
        errors.append(f"Exit conditions calculation failed with {e} ")
        return errors


    # Loop - generates contour assuming throat_length = 1
    for i in range(N):
        mach = entry_mach + (i / N) * (exit_mach - entry_mach)

        # Station Values
        try:
            # Angles
            eps = _eps_from_mach(mach, Pc, MR, c)
            gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=eps)[1]

            pm_angle   = _prandtl_meyer_angle(mach, gamma)
            flow_angle = exit_mach_angle + (exit_pm_angle - pm_angle)


            # Length of the current characteristic
            characteristic_length = mach * eps

            # Coordinates
            x = characteristic_length * -np.cos(flow_angle)
            R = characteristic_length * np.sin(flow_angle)

        except Exception as e:
            errors.append(f"Station {i+1} calculation failed with {e}")
            return errors

        # Append to results
        results["x"][i]     = x
        results["R_x"][i]   = R
        results["eps_x"][i] = eps


    # Non-dimensionalization of coordinates
    max_length = max(np.max(np.abs(results["x"])), np.max(np.abs(results["R_x"])))
    
    results["x"]   = [v / max_length for v in results["x"]]
    results["R_x"] = [v / max_length for v in results["R_x"]] 
    
    # Truncating
    results["x"] = results["x"][:int(len(results["x"]) * truncate_percent)]
    results["R_x"] = results["R_x"][:int(len(results["R_x"]) * truncate_percent)]
    results["eps_x"] = results["eps_x"][:int(len(results["eps_x"]) * truncate_percent)]

    return results, errors