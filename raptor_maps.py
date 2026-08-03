"""Generate reusable maps for the Raptor-like FullFlow example.

``raptor.py`` is intentionally a fast runtime network: expensive real-fluid and
chemical-equilibrium calculations are performed here once and stored in
``raptor.h5`` as rectangular interpolation maps.  This script does not solve the
engine, impose a target flow, resize hardware, or write a cycle operating point.
It generates only the component data consumed by the steady network.

Map families
------------
* methane and LOX pump head/torque maps;
* fuel-rich and oxygen-rich turbine torque/flow-parameter maps;
* real-fluid methane pressure-enthalpy properties for regenerative cooling;
* GRCop-42 thermal conductivity versus temperature;
* reactant enthalpy maps on the CEA/ThermoProp reference basis;
* pressure-enthalpy equilibrium maps for both preburners, both turbine exhausts,
  and the main chamber;
* entropy-pressure equilibrium maps for ideal turbine expansion and the nozzle.

Why generation can be slow
--------------------------
The pressure-enthalpy equilibrium maps require an outer temperature root solve
at each grid cell.  Each root-function evaluation performs a complete TP
chemical-equilibrium solve, so rebuilding all maps can take much longer than
running the engine.  A complete ``raptor.h5`` is included with the repository;
normal users do not need to regenerate it.

Safe restart behavior
---------------------
``REBUILD_MAPS`` defaults to ``False``.  Completed map cells are resumed and
reused, which makes accidental reruns inexpensive.  Set it to ``True`` only
after intentionally changing map axes, outputs, reference chemistry, or
component data.  A rebuild overwrites the corresponding HDF5 map groups.

The script requests only the outputs read by ``raptor.py``.  Dynamic-only
properties such as internal energy, wall density, and wall heat capacity are
not generated because the public example is steady-state only.

All map inputs and outputs use SI units.
"""

import numpy as np
from scipy.optimize import root_scalar

from fullplot import Axis, generate_map
from thermoprop import (
    CombustionGas,
    Equilibrium,
    Fluid,
    Material,
    Propellant,
    Reactants,
)

from raptor_data import (
    BAR,
    CONVERGING_PRESSURE_GUESS,
    EXIT_PRESSURE_GUESS,
    FILE_NAME,
    GRCOP_MATERIAL_NAME,
    FUEL_PREBURNER_MAP_MIXTURE_RATIO,
    FUEL_TURBINE_EXIT_MAP_PRESSURE,
    FUEL_TURBINE_EXIT_MAP_TEMPERATURE,
    FUEL_TURBINE_MAP_FLOW_PARAMETER,
    FUEL_TURBINE_MAP_PRESSURE_RATIO,
    FUEL_TURBINE_MAP_TORQUE,
    LOX_DISCHARGE_PRESSURE_GUESS,
    LOX_MANIFOLD_TEMPERATURE_GUESS,
    LOX_PUMP_MAP_TORQUE,
    LOX_PUMP_MAP_HEAD,
    LOX_PUMP_MAP_SPEED,
    LOX_PUMP_MAP_VOLUMETRIC_FLOW,
    MAIN_CHAMBER_MAP_MIXTURE_RATIO,
    MAIN_CHAMBER_MAP_PRESSURE,
    MAIN_CHAMBER_MAP_TEMPERATURE,
    METHANE_PUMP_MAP_TORQUE,
    METHANE_PUMP_MAP_HEAD,
    METHANE_PUMP_MAP_SPEED,
    METHANE_PUMP_MAP_VOLUMETRIC_FLOW,
    OXYGEN_PREBURNER_MAP_MIXTURE_RATIO,
    OXYGEN_PREBURNER_METHANE_MAP_TEMPERATURE,
    OXYGEN_TURBINE_EXIT_MAP_PRESSURE,
    OXYGEN_TURBINE_EXIT_MAP_TEMPERATURE,
    OXYGEN_TURBINE_MAP_FLOW_PARAMETER,
    OXYGEN_TURBINE_MAP_PRESSURE_RATIO,
    OXYGEN_TURBINE_MAP_TORQUE,
    PREBURNER_MAP_PRESSURE,
    HOT_METHANE_MAP_TEMPERATURE,
    THROAT_PRESSURE_GUESS,
    UPPER_NOZZLE_PRESSURE_GUESS,
)


# The repository ships with completed maps.  Leave this False to resume or
# reuse those groups.  Set it True only when deliberately regenerating every
# map after changing the data or map definitions below.
REBUILD_MAPS = False
MAP_RESUME = not REBUILD_MAPS
MAP_OVERWRITE = REBUILD_MAPS


def value(obj, *names):
    """Return the first available numeric attribute from ``obj``.

    ThermoProp objects use consistent public properties, but a few releases
    expose common quantities under alternate names such as ``cp`` versus
    ``specific_heat_cp``.  This small compatibility helper keeps the actual map
    callbacks readable while still failing loudly when a required property is
    unavailable.
    """
    for name in names:
        result = getattr(obj, name, None)
        if result is not None:
            return float(result)
    raise AttributeError(f"{type(obj).__name__} does not provide {names}.")


GAS_OUTPUT_ALIASES = {
    "temperature": ("temperature",),
    "enthalpy": ("enthalpy",),
    "entropy": ("entropy",),
    "density": ("density",),
    "specific_heat_cp": ("specific_heat_cp", "cp"),
    "gamma": ("gamma", "specific_heat_ratio"),
    "gas_constant": ("gas_constant",),
    "speed_of_sound": ("speed_of_sound",),
    "dynamic_viscosity": ("dynamic_viscosity", "viscosity"),
    "conductivity": ("conductivity", "thermal_conductivity"),
    "prandtl": ("prandtl", "prandtl_number"),
}


def gas_outputs(gas, output_names=None):
    """Extract only the requested scalar properties from a gas state.

    Reading transport and derivative properties can be appreciably more
    expensive than reading temperature or enthalpy.  Limiting extraction to the
    outputs stored by each map avoids calculating properties that FullPlot would
    immediately discard.
    """
    if output_names is None:
        output_names = GAS_OUTPUT_ALIASES

    return {
        name: value(gas, *GAS_OUTPUT_ALIASES[name])
        for name in output_names
    }


def equilibrium_outputs_at_enthalpy(
    reactants,
    pressure,
    enthalpy,
    temperature_bounds,
    guess_temperature,
    output_names,
):
    """Return an equilibrium state at a requested pressure and enthalpy.

    ThermoProp supplies TP equilibrium directly.  This helper performs the
    pressure-enthalpy inversion offline by solving for the temperature whose TP
    equilibrium enthalpy matches ``enthalpy``.  The resulting map lets
    ``raptor.py`` use pressure and enthalpy as normal network inputs without an
    extra temperature iteration inside every FullFlow residual evaluation.

    ``guess_temperature`` is used to split the broad safety bracket before
    Brent's method is called.  The method remains bracketed and robust, while a
    physically reasonable guess often reduces the number of expensive
    equilibrium evaluations.
    """
    minimum_temperature, maximum_temperature = temperature_bounds
    guess_temperature = min(
        max(float(guess_temperature), minimum_temperature),
        maximum_temperature,
    )

    def enthalpy_error(temperature):
        """Return equilibrium enthalpy minus the requested target enthalpy."""
        state = Equilibrium(
            reactants=reactants,
            mode="tp",
            pressure=pressure,
            temperature=temperature,
            guess_temperature=temperature,
        )
        return value(state, "enthalpy") - enthalpy

    minimum_error = enthalpy_error(minimum_temperature)
    maximum_error = enthalpy_error(maximum_temperature)

    if minimum_error == 0.0:
        temperature = minimum_temperature
    elif maximum_error == 0.0:
        temperature = maximum_temperature
    elif minimum_error * maximum_error > 0.0:
        raise ValueError(
            f"Requested enthalpy {enthalpy:.6g} J/kg is outside the equilibrium "
            f"map range at {pressure:.6g} Pa."
        )
    else:
        guess_error = enthalpy_error(guess_temperature)

        if guess_error == 0.0:
            temperature = guess_temperature
        else:
            if minimum_error * guess_error < 0.0:
                bracket = (minimum_temperature, guess_temperature)
            elif guess_error * maximum_error < 0.0:
                bracket = (guess_temperature, maximum_temperature)
            else:
                bracket = (minimum_temperature, maximum_temperature)

            temperature = root_scalar(
                enthalpy_error,
                bracket=bracket,
                method="brentq",
                xtol=1.0e-5,
            ).root

    state = Equilibrium(
        reactants=reactants,
        mode="tp",
        pressure=pressure,
        temperature=temperature,
        guess_temperature=temperature,
    )
    return gas_outputs(state, output_names)


FUEL_PREBURNER_ENTHALPIES = np.array(
    [-4.10e6, -3.85e6, -3.60e6, -3.5168e6, -3.30e6, -2.95e6, -2.60e6]
)
OXYGEN_PREBURNER_ENTHALPIES = np.array(
    [-5.60e5, -5.00e5, -4.40e5, -3.95e5, -3.60e5, -2.80e5, -2.00e5, -1.50e5]
)
FUEL_TURBINE_EXHAUST_ENTHALPIES = np.array(
    [-4.25e6, -4.00e6, -3.75815e6, -3.50e6, -3.25e6, -3.00e6, -2.75e6, -2.55e6, -2.35e6]
)
OXYGEN_TURBINE_EXHAUST_ENTHALPIES = np.array(
    [-6.00e5, -5.50e5, -5.00e5, -4.892e5, -4.50e5, -4.00e5, -3.25e5, -2.50e5]
)
MAIN_CHAMBER_ENTHALPIES = np.array(
    [-3.00e6, -2.50e6, -2.00e6, -1.50e6, -1.22297e6, -1.00e6, -5.00e5, 0.0, 1.00e6]
)

# Slightly denser pressure and mixture-ratio axes preserve the original TP-map
# interpolation closely after changing the runtime parameterization to PH.
FUEL_PREBURNER_PRESSURES = np.array([549.9, 580.0, 611.405, 642.0, 672.1]) * BAR
OXYGEN_PREBURNER_PRESSURES = np.array([549.9, 585.0, 621.476, 650.0, 672.1]) * BAR
FUEL_PREBURNER_MIXTURE_RATIOS = np.array([0.07, 0.085, 0.101922, 0.125, 0.15])
OXYGEN_PREBURNER_MIXTURE_RATIOS = np.array([48.0, 54.0, 59.3116, 64.0, 70.0])

FUEL_TURBINE_PH_PRESSURES = np.array([256.5, 300.0, 336.0, 380.0, 427.5]) * BAR
OXYGEN_TURBINE_PH_PRESSURES = np.array([285.0, 335.0, 384.576, 430.0, 475.0]) * BAR
FUEL_TURBINE_PH_MIXTURE_RATIOS = FUEL_PREBURNER_MIXTURE_RATIOS
OXYGEN_TURBINE_PH_MIXTURE_RATIOS = OXYGEN_PREBURNER_MIXTURE_RATIOS

CHAMBER_MIXTURE_RATIOS = [
    0.85 * MAIN_CHAMBER_MAP_MIXTURE_RATIO,
    MAIN_CHAMBER_MAP_MIXTURE_RATIO,
    1.15 * MAIN_CHAMBER_MAP_MIXTURE_RATIO,
]

MAIN_CHAMBER_PH_PRESSURES = np.array([260.0, 280.0, 299.873, 320.0, 340.0]) * BAR
MAIN_CHAMBER_PH_MIXTURE_RATIOS = np.array(
    [CHAMBER_MIXTURE_RATIOS[0], 3.20, 3.45499, 3.70, CHAMBER_MIXTURE_RATIOS[-1]]
)


PREBURNER_OUTPUTS = [
    "temperature",
    "entropy",
    "gamma",
    "gas_constant",
]

MAIN_CHAMBER_OUTPUTS = [
    "temperature",
    "entropy",
    "density",
    "specific_heat_cp",
    "gas_constant",
    "speed_of_sound",
    "dynamic_viscosity",
    "prandtl",
]

NOZZLE_OUTPUTS = [
    "temperature",
    "enthalpy",
    "density",
    "gas_constant",
    "speed_of_sound",
    "dynamic_viscosity",
    "prandtl",
]


# Only generate the turbine properties that raptor.py actually reads.
TURBINE_EXHAUST_OUTPUTS = [
    "temperature",
    "gamma",
    "gas_constant",
]
TURBINE_IDEAL_OUTPUTS = [
    "temperature",
    "enthalpy",
]



def reference_preburner(mixture_ratio, methane_temperature):
    """Build one reference preburner equilibrium state for map composition.

    The returned gas composition fixes the elemental/product basis used when
    constructing the main-chamber stream maps.  It does not prescribe the
    solved preburner pressure, temperature, or mixture ratio in ``raptor.py``.
    """
    return Equilibrium(
        reactants=Reactants(
            fuels=Propellant("ch4", temperature=methane_temperature),
            oxidizers=Propellant("lox", temperature=LOX_MANIFOLD_TEMPERATURE_GUESS),
            mixture_ratio=mixture_ratio,
        ),
        mode="hp",
        pressure=PREBURNER_MAP_PRESSURE,
    )


FuelPreburnerReference = reference_preburner(
    FUEL_PREBURNER_MAP_MIXTURE_RATIO,
    HOT_METHANE_MAP_TEMPERATURE,
)
OxygenPreburnerReference = reference_preburner(
    OXYGEN_PREBURNER_MAP_MIXTURE_RATIO,
    OXYGEN_PREBURNER_METHANE_MAP_TEMPERATURE,
)

fuel_preburner_composition = FuelPreburnerReference.gas_composition
oxygen_preburner_composition = OxygenPreburnerReference.gas_composition


def chamber_reactants(stream_mixture_ratio):
    """Return the two turbine-product streams used by chamber equilibrium.

    ``stream_mixture_ratio`` is oxygen-rich-stream mass flow divided by
    fuel-rich-stream mass flow.  The reference stream compositions come from
    the corresponding preburner equilibrium states above.
    """
    fuel_stream = CombustionGas(
        fuel_preburner_composition,
        basis="mass",
        pressure=FUEL_TURBINE_EXIT_MAP_PRESSURE,
        temperature=FUEL_TURBINE_EXIT_MAP_TEMPERATURE,
    )
    oxygen_stream = CombustionGas(
        oxygen_preburner_composition,
        basis="mass",
        pressure=OXYGEN_TURBINE_EXIT_MAP_PRESSURE,
        temperature=OXYGEN_TURBINE_EXIT_MAP_TEMPERATURE,
    )
    return Reactants(
        fuels=fuel_stream,
        oxidizers=oxygen_stream,
        mixture_ratio=stream_mixture_ratio,
    )


# -----------------------------------------------------------------------------
# Regenerative-cooling property maps
# -----------------------------------------------------------------------------

REGEN_METHANE_OUTPUTS = [
    "temperature",
    "density",
    "specific_heat_cp",
    "dynamic_viscosity",
    "conductivity",
    "prandtl",
]


def regen_methane_state(pressure, target_enthalpy):
    """Return real-fluid methane properties at one pressure-enthalpy point.

    The map spans the cryogenic-to-hot supercritical cooling path.  Any invalid
    pressure/enthalpy pair is re-raised with engineering units so a user can
    correct map coverage rather than silently accepting extrapolated nonsense.
    """
    try:
        methane = Fluid(
            "Methane",
            pressure=pressure,
            enthalpy=target_enthalpy,
        )
    except Exception as error:
        raise ValueError(
            "Invalid regenerative-cooling methane state: "
            f"pressure={pressure / BAR:.3f} bar, "
            f"enthalpy={target_enthalpy / 1.0e6:.6f} MJ/kg."
        ) from error

    return {
        "temperature": value(methane, "temperature"),
        "enthalpy": value(methane, "enthalpy"),
        "density": value(methane, "density"),
        "specific_heat_cp": value(methane, "specific_heat_cp", "cp"),
        "dynamic_viscosity": value(
            methane,
            "dynamic_viscosity",
            "viscosity",
        ),
        "conductivity": value(
            methane,
            "conductivity",
            "thermal_conductivity",
        ),
        "prandtl": value(methane, "prandtl", "prandtl_number"),
    }


def grcop42_state(temperature):
    """Return GRCop-42 thermal conductivity at ``temperature``."""
    material = Material(
        GRCOP_MATERIAL_NAME,
        temperature=temperature,
        allow_extrapolation=True,
    )
    return {
        "thermal_conductivity": value(material, "thermal_conductivity", "conductivity"),
    }


# Map generation writes directly into raptor.h5.  With REBUILD_MAPS=False,
# FullPlot keeps completed cells and evaluates only missing work.
print("Generating regenerative-cooling property maps...")
generate_map(
    FILE_NAME,
    group="regen_methane_ph",
    axes=[
        # The solved cooling circuit operates from roughly 690 to 880 bar.
        # Extending the rectangular PH grid to 950-1000 bar paired those high
        # pressures with the lowest enthalpy point, which is outside CoolProp's
        # valid methane domain. 600-900 bar covers the actual network with
        # margin while keeping every generated pressure-enthalpy pair valid.
        Axis.linear("pressure", 600.0 * BAR, 900.0 * BAR, 7, units="Pa"),
        Axis.linear(
            "target_enthalpy",
            0.15e6,
            2.25e6,
            17,
            units="J/kg",
        ),
    ],
    evaluate=regen_methane_state,
    outputs=REGEN_METHANE_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="grcop42_t",
    axes=[
        Axis.linear("temperature", 100.0, 2600.0, 31, units="K"),
    ],
    evaluate=grcop42_state,
    outputs=["thermal_conductivity"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)


# -----------------------------------------------------------------------------
# Pump maps
# -----------------------------------------------------------------------------


def pump_map(
    rotor_speed,
    volumetric_flow,
    design_speed,
    design_flow,
    design_head,
    design_torque,
):
    """Evaluate the smooth representative pump head and torque surface.

    This is a compact public map shape centered on one documented design point,
    not a proprietary Raptor pump map.  FullFlow later interpolates the generated
    grid and calculates pump efficiency from hydraulic and shaft power.
    """
    speed_ratio = rotor_speed / design_speed
    flow_ratio = volumetric_flow / design_flow
    flow_coefficient_ratio = flow_ratio / speed_ratio

    head_denominator = 0.80 + 0.20 * flow_coefficient_ratio**2
    torque_shape = (
        flow_coefficient_ratio
        * (1.0 + 0.70 * (flow_coefficient_ratio - 1.0) ** 2)
        / head_denominator
    )

    head_rise = design_head * speed_ratio**2 / head_denominator
    torque = design_torque * speed_ratio**2 * torque_shape

    return {
        "head_rise": head_rise,
        "torque": torque,
    }


def methane_pump_map(rotor_speed, volumetric_flow):
    """Evaluate the methane-pump map at one speed and volumetric flow."""
    return pump_map(
        rotor_speed,
        volumetric_flow,
        METHANE_PUMP_MAP_SPEED,
        METHANE_PUMP_MAP_VOLUMETRIC_FLOW,
        METHANE_PUMP_MAP_HEAD,
        METHANE_PUMP_MAP_TORQUE,
    )


def lox_pump_map(rotor_speed, volumetric_flow):
    """Evaluate the LOX-pump map at one speed and volumetric flow."""
    return pump_map(
        rotor_speed,
        volumetric_flow,
        LOX_PUMP_MAP_SPEED,
        LOX_PUMP_MAP_VOLUMETRIC_FLOW,
        LOX_PUMP_MAP_HEAD,
        LOX_PUMP_MAP_TORQUE,
    )


print("Generating pump maps...")
generate_map(
    FILE_NAME,
    group="methane_pump",
    axes=[
        Axis.linear(
            "rotor_speed",
            0.75 * METHANE_PUMP_MAP_SPEED,
            1.25 * METHANE_PUMP_MAP_SPEED,
            9,
            units="rpm",
        ),
        Axis.linear(
            "volumetric_flow",
            0.75 * METHANE_PUMP_MAP_VOLUMETRIC_FLOW,
            1.25 * METHANE_PUMP_MAP_VOLUMETRIC_FLOW,
            9,
            units="m3/s",
        ),
    ],
    evaluate=methane_pump_map,
    outputs=["head_rise", "torque"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="lox_pump",
    axes=[
        Axis.linear(
            "rotor_speed",
            0.75 * LOX_PUMP_MAP_SPEED,
            1.25 * LOX_PUMP_MAP_SPEED,
            9,
            units="rpm",
        ),
        Axis.linear(
            "volumetric_flow",
            0.75 * LOX_PUMP_MAP_VOLUMETRIC_FLOW,
            1.25 * LOX_PUMP_MAP_VOLUMETRIC_FLOW,
            9,
            units="m3/s",
        ),
    ],
    evaluate=lox_pump_map,
    outputs=["head_rise", "torque"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)


# -----------------------------------------------------------------------------
# Turbine maps
# -----------------------------------------------------------------------------


def turbine_map(
    pressure_ratio,
    speed_ratio,
    design_pressure_ratio,
    design_torque,
    design_flow_parameter,
):
    """Evaluate the representative turbine torque and flow-parameter surface.

    The map supplies aerodynamic performance data only.  ``GasTurbine`` uses
    the solved shaft speed and thermodynamic states to calculate shaft power,
    actual enthalpy drop, and efficiency.
    """
    pressure_drive = max(
        (1.0 - 1.0 / pressure_ratio)
        / (1.0 - 1.0 / design_pressure_ratio),
        1.0e-12,
    )
    return {
        "torque": (
            design_torque
            * pressure_drive
            * (1.05 - 0.05 * speed_ratio)
        ),
        "flow_parameter": (
            design_flow_parameter
            * pressure_drive**0.25
            * (1.02 - 0.02 * speed_ratio)
        ),
    }


def fuel_turbine_map(pressure_ratio, speed_ratio):
    """Evaluate the fuel-rich turbine map at one operating point."""
    return turbine_map(
        pressure_ratio,
        speed_ratio,
        FUEL_TURBINE_MAP_PRESSURE_RATIO,
        FUEL_TURBINE_MAP_TORQUE,
        FUEL_TURBINE_MAP_FLOW_PARAMETER,
    )


def oxygen_turbine_map(pressure_ratio, speed_ratio):
    """Evaluate the oxygen-rich turbine map at one operating point."""
    return turbine_map(
        pressure_ratio,
        speed_ratio,
        OXYGEN_TURBINE_MAP_PRESSURE_RATIO,
        OXYGEN_TURBINE_MAP_TORQUE,
        OXYGEN_TURBINE_MAP_FLOW_PARAMETER,
    )


print("Generating turbine maps...")
generate_map(
    FILE_NAME,
    group="fuel_turbine",
    axes=[
        Axis.linear(
            "pressure_ratio",
            0.80 * FUEL_TURBINE_MAP_PRESSURE_RATIO,
            1.20 * FUEL_TURBINE_MAP_PRESSURE_RATIO,
            9,
        ),
        Axis.linear("speed_ratio", 0.75, 1.25, 9),
    ],
    evaluate=fuel_turbine_map,
    outputs=["torque", "flow_parameter"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="oxygen_turbine",
    axes=[
        Axis.linear(
            "pressure_ratio",
            0.80 * OXYGEN_TURBINE_MAP_PRESSURE_RATIO,
            1.20 * OXYGEN_TURBINE_MAP_PRESSURE_RATIO,
            9,
        ),
        Axis.linear("speed_ratio", 0.75, 1.25, 9),
    ],
    evaluate=oxygen_turbine_map,
    outputs=["torque", "flow_parameter"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)


# -----------------------------------------------------------------------------
# Reactant enthalpy maps
# -----------------------------------------------------------------------------

# Fluid supplies real-fluid hydraulic states in raptor.py.  The combustion
# balances use Propellant and Equilibrium on the CEA reference basis.
#
# LOX remains a compressed liquid, so pressure and temperature are retained
# for its reactant enthalpy. Both methane reactant streams now come from the
# regenerative circuit: the mixed return feeds the fuel-rich preburner and the
# nozzle-return tap feeds the oxygen-rich preburner. These methane streams are
# above the critical temperature, so both use the CEA gas-species enthalpy map,
# which depends on temperature only.


def hot_methane_reactant_enthalpy(temperature):
    """Return CEA-reference methane enthalpy for a hot reactant stream."""
    methane = Propellant("ch4", temperature=temperature)
    return {"enthalpy": value(methane, "enthalpy")}


def lox_reactant_enthalpy(pressure, temperature):
    """Return CEA-reference LOX enthalpy at the post-valve state."""
    lox = Propellant(
        "lox",
        pressure=pressure,
        temperature=temperature,
    )
    return {"enthalpy": value(lox, "enthalpy")}


print("Generating reactant enthalpy maps...")
generate_map(
    FILE_NAME,
    group="hot_methane_reactant_cea_t",
    axes=[
        Axis.linear("temperature", 250.0, 850.0, 25, units="K"),
    ],
    evaluate=hot_methane_reactant_enthalpy,
    outputs=["enthalpy"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="lox_reactant_pt",
    axes=[
        Axis.linear(
            "pressure",
            0.80 * LOX_DISCHARGE_PRESSURE_GUESS,
            1.20 * LOX_DISCHARGE_PRESSURE_GUESS,
            5,
            units="Pa",
        ),
        Axis.linear("temperature", 90.0, 135.0, 9, units="K"),
    ],
    evaluate=lox_reactant_enthalpy,
    outputs=["enthalpy"],
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)


# -----------------------------------------------------------------------------
# Preburner equilibrium maps
# -----------------------------------------------------------------------------


def preburner_ph(pressure, mixture_ratio, target_enthalpy):
    """Return a preburner HP-equilibrium state through offline TP inversion."""
    reactants = Reactants(
        fuels=Propellant("ch4", temperature=300.0),
        oxidizers=Propellant("lox", temperature=100.0),
        mixture_ratio=mixture_ratio,
    )
    return equilibrium_outputs_at_enthalpy(
        reactants,
        pressure,
        target_enthalpy,
        temperature_bounds=(400.0, 1400.0),
        guess_temperature=850.0,
        output_names=PREBURNER_OUTPUTS,
    )


print("Generating pressure-enthalpy preburner equilibrium maps...")
generate_map(
    FILE_NAME,
    group="fuel_preburner_equilibrium_ph",
    axes=[
        Axis.values(
            "pressure",
            FUEL_PREBURNER_PRESSURES,
            units="Pa",
        ),
        Axis.values(
            "mixture_ratio",
            FUEL_PREBURNER_MIXTURE_RATIOS,
        ),
        Axis.values(
            "target_enthalpy",
            FUEL_PREBURNER_ENTHALPIES,
            units="J/kg",
        ),
    ],
    evaluate=preburner_ph,
    outputs=PREBURNER_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="oxygen_preburner_equilibrium_ph",
    axes=[
        Axis.values(
            "pressure",
            OXYGEN_PREBURNER_PRESSURES,
            units="Pa",
        ),
        Axis.values(
            "mixture_ratio",
            OXYGEN_PREBURNER_MIXTURE_RATIOS,
        ),
        Axis.values(
            "target_enthalpy",
            OXYGEN_PREBURNER_ENTHALPIES,
            units="J/kg",
        ),
    ],
    evaluate=preburner_ph,
    outputs=PREBURNER_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)


# -----------------------------------------------------------------------------
# Turbine thermodynamic maps
# -----------------------------------------------------------------------------

# The preburner and turbine maps use the same CEA equilibrium thermodynamic
# basis.  No CombustionGas reference conversion or fixed-composition bridge is
# used between them.  The actual outlet maps are pressure-enthalpy equilibrium states.  The
# ideal outlet maps are SP equilibrium states at the turbine discharge pressure
# and the inlet entropy.  This makes h_in, h_out, and h_out,s directly
# comparable in GasTurbine's efficiency calculation.

FUEL_TURBINE_MIXTURE_RATIOS = [
    0.07,
    FUEL_PREBURNER_MAP_MIXTURE_RATIO,
    0.15,
]
OXYGEN_TURBINE_MIXTURE_RATIOS = [
    48.0,
    OXYGEN_PREBURNER_MAP_MIXTURE_RATIO,
    70.0,
]

FUEL_TURBINE_TEMPERATURES = [
    600.0,
    700.0,
    800.0,
    900.0,
    1000.0,
    1125.0,
]
OXYGEN_TURBINE_TEMPERATURES = [
    525.0,
    600.0,
    675.0,
    750.0,
    825.0,
    925.0,
]


def turbine_reactants(mixture_ratio):
    """Return the CH4/O2 elemental inventory for turbine equilibrium maps."""
    # Reactant temperatures do not affect a TP or SP equilibrium composition;
    # they only provide the CH4/O2 elemental inventory and mixture ratio.
    return Reactants(
        fuels=Propellant("ch4", temperature=300.0),
        oxidizers=Propellant("lox", temperature=100.0),
        mixture_ratio=mixture_ratio,
    )


def turbine_ph(pressure, mixture_ratio, target_enthalpy):
    """Return the actual turbine-exhaust equilibrium state at ``P`` and ``h``."""
    return equilibrium_outputs_at_enthalpy(
        turbine_reactants(mixture_ratio),
        pressure,
        target_enthalpy,
        temperature_bounds=(400.0, 1350.0),
        guess_temperature=775.0,
        output_names=TURBINE_EXHAUST_OUTPUTS,
    )


def turbine_sp(pressure, mixture_ratio, inlet_entropy):
    """Return the ideal isentropic turbine-exhaust equilibrium state."""
    equilibrium = Equilibrium(
        reactants=turbine_reactants(mixture_ratio),
        mode="sp",
        pressure=pressure,
        entropy=inlet_entropy,
        guess_temperature=750.0,
    )
    return gas_outputs(equilibrium, TURBINE_IDEAL_OUTPUTS)


def turbine_entropy_axis(mixture_ratios, temperatures):
    """Build an entropy axis that covers expected preburner inlet states."""
    values = []
    for mixture_ratio in mixture_ratios:
        for temperature in temperatures:
            state = Equilibrium(
                reactants=turbine_reactants(mixture_ratio),
                mode="tp",
                pressure=PREBURNER_MAP_PRESSURE,
                temperature=temperature,
                guess_temperature=temperature,
            )
            values.append(value(state, "entropy"))

    margin = 100.0
    return np.linspace(min(values) - margin, max(values) + margin, 7)


fuel_turbine_entropy_axis = turbine_entropy_axis(
    FUEL_TURBINE_MIXTURE_RATIOS,
    [650.0, 760.0, 860.0, 980.0, 1100.0],
)
oxygen_turbine_entropy_axis = turbine_entropy_axis(
    OXYGEN_TURBINE_MIXTURE_RATIOS,
    [500.0, 650.0, 760.0, 900.0, 1050.0, 1200.0],
)

print("Generating equilibrium turbine pressure-enthalpy and ideal-SP maps...")
generate_map(
    FILE_NAME,
    group="fuel_turbine_exhaust_ph",
    axes=[
        Axis.values(
            "pressure",
            FUEL_TURBINE_PH_PRESSURES,
            units="Pa",
        ),
        Axis.values("mixture_ratio", FUEL_TURBINE_PH_MIXTURE_RATIOS),
        Axis.values(
            "target_enthalpy",
            FUEL_TURBINE_EXHAUST_ENTHALPIES,
            units="J/kg",
        ),
    ],
    evaluate=turbine_ph,
    outputs=TURBINE_EXHAUST_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="oxygen_turbine_exhaust_ph",
    axes=[
        Axis.values(
            "pressure",
            OXYGEN_TURBINE_PH_PRESSURES,
            units="Pa",
        ),
        Axis.values("mixture_ratio", OXYGEN_TURBINE_PH_MIXTURE_RATIOS),
        Axis.values(
            "target_enthalpy",
            OXYGEN_TURBINE_EXHAUST_ENTHALPIES,
            units="J/kg",
        ),
    ],
    evaluate=turbine_ph,
    outputs=TURBINE_EXHAUST_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)

generate_map(
    FILE_NAME,
    group="fuel_turbine_ideal_sp",
    axes=[
        Axis.values(
            "pressure",
            np.array([0.75, 1.00, 1.25])
            * FUEL_TURBINE_EXIT_MAP_PRESSURE,
            units="Pa",
        ),
        Axis.values("mixture_ratio", FUEL_TURBINE_MIXTURE_RATIOS),
        Axis.values(
            "inlet_entropy",
            fuel_turbine_entropy_axis,
            units="J/kg-K",
        ),
    ],
    evaluate=turbine_sp,
    outputs=TURBINE_IDEAL_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)
generate_map(
    FILE_NAME,
    group="oxygen_turbine_ideal_sp",
    axes=[
        Axis.values(
            "pressure",
            np.array([0.75, 1.00, 1.25])
            * OXYGEN_TURBINE_EXIT_MAP_PRESSURE,
            units="Pa",
        ),
        Axis.values("mixture_ratio", OXYGEN_TURBINE_MIXTURE_RATIOS),
        Axis.values(
            "inlet_entropy",
            oxygen_turbine_entropy_axis,
            units="J/kg-K",
        ),
    ],
    evaluate=turbine_sp,
    outputs=TURBINE_IDEAL_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)


# -----------------------------------------------------------------------------
# Main chamber and equilibrium nozzle maps
# -----------------------------------------------------------------------------


def main_chamber_ph(pressure, stream_mixture_ratio, target_enthalpy):
    """Return the adiabatic mixed-stream chamber equilibrium state."""
    return equilibrium_outputs_at_enthalpy(
        chamber_reactants(stream_mixture_ratio),
        pressure,
        target_enthalpy,
        temperature_bounds=(2500.0, 4500.0),
        guess_temperature=3800.0,
        output_names=MAIN_CHAMBER_OUTPUTS,
    )


print("Generating main chamber pressure-enthalpy map...")
generate_map(
    FILE_NAME,
    group="main_chamber_ph",
    axes=[
        Axis.values(
            "pressure",
            MAIN_CHAMBER_PH_PRESSURES,
            units="Pa",
        ),
        Axis.values("stream_mixture_ratio", MAIN_CHAMBER_PH_MIXTURE_RATIOS),
        Axis.values(
            "target_enthalpy",
            MAIN_CHAMBER_ENTHALPIES,
            units="J/kg",
        ),
    ],
    evaluate=main_chamber_ph,
    outputs=MAIN_CHAMBER_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)

def nozzle_equilibrium(pressure, chamber_entropy, stream_mixture_ratio):
    """Return one equilibrium, isentropic chamber/nozzle station state."""
    equilibrium = Equilibrium(
        reactants=chamber_reactants(stream_mixture_ratio),
        mode="sp",
        pressure=pressure,
        entropy=chamber_entropy,
    )
    return gas_outputs(equilibrium, NOZZLE_OUTPUTS)


station_pressures = np.array(
    [
        MAIN_CHAMBER_MAP_PRESSURE,
        CONVERGING_PRESSURE_GUESS,
        THROAT_PRESSURE_GUESS,
        UPPER_NOZZLE_PRESSURE_GUESS,
        EXIT_PRESSURE_GUESS,
    ]
)
pressure_axis = np.unique(
    np.concatenate(
        [
            np.geomspace(0.25 * BAR, 1.05 * MAIN_CHAMBER_MAP_PRESSURE, 9),
            station_pressures,
        ]
    )
)
pressure_axis.sort()

chamber_entropies = []
for mixture_ratio in CHAMBER_MIXTURE_RATIOS:
    chamber = Equilibrium(
        reactants=chamber_reactants(mixture_ratio),
        mode="tp",
        pressure=MAIN_CHAMBER_MAP_PRESSURE,
        temperature=MAIN_CHAMBER_MAP_TEMPERATURE,
    )
    chamber_entropies.append(value(chamber, "entropy"))

entropy_axis = np.linspace(
    min(chamber_entropies) - 150.0,
    max(chamber_entropies) + 150.0,
    3,
)

print("Generating equilibrium SP nozzle map...")
generate_map(
    FILE_NAME,
    group="nozzle_sp",
    axes=[
        Axis(
            name="pressure",
            values=pressure_axis,
            units="Pa",
            spacing="log",
        ),
        Axis.values(
            "chamber_entropy",
            entropy_axis,
            units="J/kg-K",
        ),
        Axis.values(
            "stream_mixture_ratio",
            CHAMBER_MIXTURE_RATIOS,
        ),
    ],
    evaluate=nozzle_equilibrium,
    outputs=NOZZLE_OUTPUTS,
    resume=MAP_RESUME,
    overwrite=MAP_OVERWRITE,
    raise_errors=True,
)

print(f"Finished. Maps are stored in {FILE_NAME}.h5")
print("Run `python3 raptor.py` to solve the steady engine network.")
