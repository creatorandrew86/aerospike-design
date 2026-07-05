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

def _length_conversion(value: float, unit: str) -> float:
    match unit:
        case "m":   return value
        case "cm":  return value / 100
        case "in":  return value * 0.0254
        case "ft":  return value / 3.28

def _mass_flow_conversion(value, unit):
    match unit:
        case "kg/s":   return value
        case "lb/s":   return value * 0.453592



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
    mfr      = _mass_flow_conversion(inputs["mfr"], inputs["unit_mfr"])
    
    radius               = _length_conversion(inputs["radius"], inputs["unit_radius"])
    throat_angle         = np.radians(inputs["throat_angle"])

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
    
    # Throat area calculation
    try:
        throat_c = c.get_SonicVelocities(Pc=Pc, MR=MR, eps=exit_eps)[1]
        throat_rho = c.get_Densities(Pc=Pc, MR=MR, eps=exit_eps)[1]

        throat_area = mfr / (throat_c * throat_rho)

    except ZeroDivisionError as e:
        errors.append(f"Throat Area calculation failed with {e}")
        return errors


    try:
        if throat_angle == np.radians(90):
            throat_length = throat_area / (2 * np.pi * radius)

        elif (throat_angle >= np.radians(30) and throat_angle < np.radians(90)):
            throat_length = (-2 * np.pi * radius + np.sqrt(pow(-2 * np.pi * radius, 2) + 4 * np.pi * np.sin(throat_angle) * throat_area)) \
                            / (2 * np.pi * np.sin(throat_angle))
            
        elif throat_angle < np.radians(90):
            errors.append("Throat angle cannot be smaller than 30 degrees")

        elif throat_angle > np.radians(90):
            errors.append("Throat angle cannot be greater than 90 degrees")
    
    except Exception as e:
        errors.append(f"Throat Length Calculation failed with {e}")
        return errors



    # Aerospike Throat Conditions
    try:
        entry_mach = max(c.get_MachNumber(Pc=Pc, MR=MR, eps=effective_throat_eps), 1 + 1e-6)

        entry_gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=effective_throat_eps)[1]
        entry_pm_angle = _prandtl_meyer_angle(entry_mach, entry_gamma)
        entry_mach_angle = _mach_angle(entry_mach)

    except Exception as e:
        errors.append(f"Entry conditions calculations failed with {e}")
        return errors
    

    try:
        exit_mach  = max(c.get_MachNumber(Pc=Pc, MR=MR, eps=exit_eps), 1 + 1e-6)

        exit_gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=exit_eps)[1]
        exit_pm_angle = _prandtl_meyer_angle(exit_mach, exit_gamma)

    except Exception as e:
        errors.append(f"Exit conditions calculation failed with {e} ")
        return errors

    # Loop
    for i in range(N):
        mach = entry_mach + (i / N) * (exit_mach - entry_mach)

        # Station Values
        try:
            eps = _eps_from_mach(mach, Pc, MR, c)
            gamma = c.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=eps)[1]

            # Angles
            pm_angle   = _prandtl_meyer_angle(mach, gamma)
            mach_angle = _mach_angle(mach)
            flow_angle = (exit_pm_angle - entry_pm_angle) + (mach_angle - entry_mach_angle) - pm_angle - (np.radians(90) - throat_angle)

            # Lengths
            characteristic_length = mach * eps * throat_length

            # Coordinates
            x = characteristic_length * np.cos(flow_angle)
            R = characteristic_length * np.sin(flow_angle)

        except Exception as e:
            errors.append(f"Station {i+1} calculation failed with {e}")
            return errors

        # Append to results
        results["x"][i]     = x
        results["R_x"][i]   = R
        results["eps_x"][i] = eps


    results["x"] = results["x"][:int(len(results["x"]) * truncate_percent)]
    results["R_x"] = results["R_x"][:int(len(results["R_x"]) * truncate_percent)]
    results["eps_x"] = results["eps_x"][:int(len(results["eps_x"]) * truncate_percent)]

    return results, errors