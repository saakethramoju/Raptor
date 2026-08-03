"""Input data for the public Raptor-like FullFlow example.

The file is deliberately separated from ``raptor.py`` so that geometry,
reference map data, boundary conditions, and solver guesses can be reviewed or
changed without searching through the network assembly.  Nothing in this file
solves the engine.

Data categories
---------------
Boundary conditions
    Methane-tank pressure and temperature, LOX-tank pressure and temperature,
    and ambient pressure.  These are the only operating conditions imposed on
    the network.

Component data
    Pump/turbine map reference points, valve and injector geometry, cooling
    passages, wall thickness, and chamber/nozzle geometry.  These quantities
    define the simplified hardware represented by the network.

Map reference states
    Values used to choose the coverage and reference compositions of the
    generated HDF5 maps.  They are not operating-point constraints.

Solver starting guesses
    Initial values supplied to the nonlinear solver.  They are centered on the
    saved trim solution so the example runs quickly, but they are not targets,
    bounds, balances, or prescribed results.

Units and conventions
---------------------
All values are SI unless a variable name or comment states otherwise:
pressure in Pa, temperature in K, mass flow in kg/s, area in m^2, length in m,
head in m, torque in N*m, enthalpy in J/kg, and shaft speed in rpm.  ``BAR`` is
provided only as a readable conversion factor.

Model provenance
----------------
This is a Raptor-like educational cycle model, not a claim of proprietary
hardware geometry or map data.  Publicly unavailable quantities are documented
as explicit effective assumptions so users can replace them with their own
validated data.
"""

import math

# FullPlot/FullFlow append the ``.h5`` suffix where required.  The same HDF5
# file stores reusable maps and the latest saved FullFlow solution.
FILE_NAME = "raptor"

BAR = 100_000.0
G0 = 9.80665


# -----------------------------------------------------------------------------
# Boundary conditions
# -----------------------------------------------------------------------------

METHANE_TANK_PRESSURE = 4.0 * BAR
METHANE_TANK_TEMPERATURE = 111.0

LOX_TANK_PRESSURE = 4.0 * BAR
LOX_TANK_TEMPERATURE = 90.0

AMBIENT_PRESSURE = 101_325.0


# -----------------------------------------------------------------------------
# Pump-map data
#
# Head and torque are the map outputs. Pump efficiency is calculated by the
# ConstantDensityPump component from hydraulic power and shaft power.
# -----------------------------------------------------------------------------

METHANE_PUMP_MAP_SPEED = 36_000.0
METHANE_PUMP_MAP_VOLUMETRIC_FLOW = 0.3515021537856995
METHANE_PUMP_MAP_HEAD = 21_236.340032085456
METHANE_PUMP_MAP_TORQUE = 9_792.287043656543

LOX_PUMP_MAP_SPEED = 30_000.0
LOX_PUMP_MAP_VOLUMETRIC_FLOW = 0.46866740062102363
LOX_PUMP_MAP_HEAD = 6_173.94917818599
LOX_PUMP_MAP_TORQUE = 12_587.61865981648


# -----------------------------------------------------------------------------
# Turbine-map data
# -----------------------------------------------------------------------------

FUEL_TURBINE_MAP_PRESSURE_RATIO = 1.82
FUEL_TURBINE_MAP_TORQUE = METHANE_PUMP_MAP_TORQUE
FUEL_TURBINE_MAP_FLOW_PARAMETER = 0.001677793915715393

# Each turbine design torque matches the pump design torque on the same shaft.
# No shaft-efficiency multiplier is used. Pump and turbine efficiencies are
# calculated from the map torque, flow, head, pressure ratio, and thermodynamic
# states at the solved operating point.
#
# The oxygen turbine is centered at a larger expansion ratio than the early
# open-source 430-bar outlet estimate. At the calculated LOX-pump power, that
# smaller ratio did not provide enough ideal enthalpy drop.
OXYGEN_TURBINE_MAP_PRESSURE_RATIO = 1.60
OXYGEN_TURBINE_MAP_TORQUE = LOX_PUMP_MAP_TORQUE
OXYGEN_TURBINE_MAP_FLOW_PARAMETER = 0.0038897144649489342



# -----------------------------------------------------------------------------
# Property-map reference states
# -----------------------------------------------------------------------------

# These values choose map coverage and the fixed reference compositions used
# by the turbine-exhaust property maps. The actual preburner mixture ratios,
# temperatures, pressures, and flows are calculated by the network.
FUEL_PREBURNER_MAP_MIXTURE_RATIO = 0.1012
OXYGEN_PREBURNER_MAP_MIXTURE_RATIO = 58.1563
PREBURNER_MAP_PRESSURE = 611.0 * BAR

FUEL_TURBINE_EXIT_MAP_PRESSURE = 342.0 * BAR
FUEL_TURBINE_EXIT_MAP_TEMPERATURE = 808.0
OXYGEN_TURBINE_EXIT_MAP_PRESSURE = 380.0 * BAR
OXYGEN_TURBINE_EXIT_MAP_TEMPERATURE = 725.0

MAIN_CHAMBER_MAP_PRESSURE = 298.6 * BAR
MAIN_CHAMBER_MAP_TEMPERATURE = 3748.92
MAIN_CHAMBER_MAP_MIXTURE_RATIO = 3.44192007

# Reference temperature used only to construct the fixed fuel-rich composition
# axes in raptor_maps.py. The solved merged hot-methane temperature remains a
# network output.
HOT_METHANE_MAP_TEMPERATURE = 620.0

# The oxygen-rich preburner now receives methane tapped from the hot nozzle
# cooling return rather than cold pump-discharge methane. This value is used
# only to construct the fixed reference composition for generated maps.
OXYGEN_PREBURNER_METHANE_MAP_TEMPERATURE = 550.0


# -----------------------------------------------------------------------------
# Fixed valves, restrictions, and injectors
# -----------------------------------------------------------------------------

LIQUID_INJECTOR_CD = 0.90
GAS_INJECTOR_CD = 0.90

# Fully open tank isolation valves. These are fixed hardware CdA values, not
# commands or solve variables. Their areas are deliberately large so they add
# only a small inlet loss at the expected main-stage flow rates.
MAIN_METHANE_VALVE_CD = 0.90
MAIN_METHANE_VALVE_AREA = 4.00e-2
MAIN_LOX_VALVE_CD = 0.90
MAIN_LOX_VALVE_AREA = 1.20e-1

# Engine control valves shown in the cycle schematic.
#
# The main fuel valve is upstream of the entire regenerative circuit. The
# oxygen-preburner fuel valve (OPFV/OBFV) is supplied from the hot nozzle-return
# tap before the chamber and nozzle returns mix. The main oxidizer valve (MOV)
# is one equivalent valve representing the parallel integrated shutoff elements
# visible in public Raptor reconstructions. The fuel-preburner oxidizer valve
# remains a small branch from the LOX-pump discharge manifold.
#
# All areas below are fixed hardware CdA data. Their flows and pressure drops
# are solved by the network; none is a prescribed operating condition.
MAIN_FUEL_VALVE_CD = 0.85
MAIN_FUEL_VALVE_AREA = 5.18e-3

OXYGEN_PREBURNER_FUEL_VALVE_CD = 0.85
OXYGEN_PREBURNER_FUEL_VALVE_AREA = 3.50e-4

FUEL_PREBURNER_OXIDIZER_VALVE_CD = 0.85
FUEL_PREBURNER_OXIDIZER_VALVE_AREA = 3.11e-4

MAIN_OXIDIZER_VALVE_CD = 0.85
MAIN_OXIDIZER_VALVE_AREA = 1.80e-2

# Injector areas are fixed hardware geometry. The minority-stream injector
# areas are sized together with the upstream control valves so the calculated
# branch flows remain near the reference cycle while preserving the actual
# valve -> injector ordering shown in the schematic.
FUEL_PREBURNER_METHANE_AREA = 2.827e-3
FUEL_PREBURNER_LOX_AREA = 1.127e-4
OXYGEN_PREBURNER_METHANE_AREA = 2.20e-4
OXYGEN_PREBURNER_LOX_AREA = 4.22e-3

MAIN_FUEL_INJECTOR_AREA = 7.736e-3
# Sized so the oxygen-rich turbine can expand to roughly 380 bar while the
# main chamber remains near 300 bar at the design flow.
MAIN_OXYGEN_INJECTOR_AREA = 1.14110e-2


# -----------------------------------------------------------------------------
# Chamber and nozzle geometry
# -----------------------------------------------------------------------------

THROAT_AREA = 0.0430
EXPANSION_RATIO = 34.0
EXIT_AREA = THROAT_AREA * EXPANSION_RATIO

CHAMBER_DIAMETER = 0.30
CHAMBER_RADIUS = 0.5 * CHAMBER_DIAMETER
CHAMBER_AREA = math.pi * CHAMBER_RADIUS**2
CHAMBER_AREA_RATIO = CHAMBER_AREA / THROAT_AREA

THROAT_RADIUS = math.sqrt(THROAT_AREA / math.pi)
EXIT_RADIUS = math.sqrt(EXIT_AREA / math.pi)

CONVERGING_AREA_RATIO = 1.25
UPPER_NOZZLE_AREA_RATIO = 4.0
UPPER_NOZZLE_SPLIT_AREA_RATIO = 8.0
UPPER_NOZZLE_SPLIT_RADIUS = math.sqrt(
    THROAT_AREA * UPPER_NOZZLE_SPLIT_AREA_RATIO / math.pi
)

# Representative station used for heat transfer over the final nozzle section.
# The section spans A/At = 8 to 34, so using the geometric exit state for the
# entire section underpredicts its average gas-side heating.
EXIT_NOZZLE_THERMAL_RADIUS = 0.5 * (
    UPPER_NOZZLE_SPLIT_RADIUS + EXIT_RADIUS
)
EXIT_NOZZLE_THERMAL_AREA_RATIO = (
    math.pi * EXIT_NOZZLE_THERMAL_RADIUS**2 / THROAT_AREA
)


# -----------------------------------------------------------------------------
# Regenerative cooling and materials
# -----------------------------------------------------------------------------

# The methane cooling supply manifold is near the throat. It feeds two
# parallel branches:
#
#   1. a chamber branch running from the throat region toward the injector end;
#   2. a nozzle branch using alternating downflow/upflow channels. Half of the
#      nozzle channels carry methane toward the exit and the neighboring half
#      return it toward the throat-side nozzle return manifold.
#
# Public Raptor channel dimensions are not available. The values below are
# explicit effective geometry assumptions for this open-source cycle model.
# Pressure is never interpolated or imposed inside the cooling circuit: every
# listed passage is represented by a FlowTube and every intermediate pressure
# is solved from the local momentum equation. The two effective channel widths
# were selected together with the documented roughness so the detailed channel
# model remains near the public 886-bar methane-pump / 300-bar chamber cycle.
COOLING_CHANNEL_ROUGHNESS = 2.5e-7
CHAMBER_COOLING_ROUGHNESS = COOLING_CHANNEL_ROUGHNESS
NOZZLE_COOLING_ROUGHNESS = COOLING_CHANNEL_ROUGHNESS

CHAMBER_COOLING_CHANNEL_COUNT = 306
CHAMBER_COOLING_CHANNEL_WIDTH = 0.912e-3
CHAMBER_COOLING_CHANNEL_HEIGHT = 6.25e-3
CHAMBER_COOLING_FLOW_AREA = (
    CHAMBER_COOLING_CHANNEL_COUNT
    * CHAMBER_COOLING_CHANNEL_WIDTH
    * CHAMBER_COOLING_CHANNEL_HEIGHT
)
CHAMBER_COOLING_HYDRAULIC_DIAMETER = (
    2.0
    * CHAMBER_COOLING_CHANNEL_WIDTH
    * CHAMBER_COOLING_CHANNEL_HEIGHT
    / (
        CHAMBER_COOLING_CHANNEL_WIDTH
        + CHAMBER_COOLING_CHANNEL_HEIGHT
    )
)

NOZZLE_COOLING_CHANNEL_COUNT = 306
NOZZLE_COOLING_PASS_CHANNEL_COUNT = NOZZLE_COOLING_CHANNEL_COUNT // 2
NOZZLE_COOLING_CHANNEL_WIDTH = 2.477e-3
NOZZLE_COOLING_CHANNEL_HEIGHT = 6.25e-3
NOZZLE_COOLING_PASS_FLOW_AREA = (
    NOZZLE_COOLING_PASS_CHANNEL_COUNT
    * NOZZLE_COOLING_CHANNEL_WIDTH
    * NOZZLE_COOLING_CHANNEL_HEIGHT
)
NOZZLE_COOLING_HYDRAULIC_DIAMETER = (
    2.0
    * NOZZLE_COOLING_CHANNEL_WIDTH
    * NOZZLE_COOLING_CHANNEL_HEIGHT
    / (
        NOZZLE_COOLING_CHANNEL_WIDTH
        + NOZZLE_COOLING_CHANNEL_HEIGHT
    )
)

# The two branch outlets connect to the hot-methane merge header through short,
# large-area return manifolds. These are explicit geometric assumptions rather
# than guessed pressure drops. Their solved losses are expected to be small
# relative to the cooling-channel losses.
REGEN_RETURN_PIPE_DIAMETER = 0.120
REGEN_RETURN_PIPE_FLOW_AREA = math.pi * REGEN_RETURN_PIPE_DIAMETER**2 / 4.0
REGEN_RETURN_PIPE_ROUGHNESS = 2.5e-6
CHAMBER_RETURN_PIPE_LENGTH = 0.20
NOZZLE_RETURN_PIPE_LENGTH = 0.25

# Effective copper conduction thickness. The earlier 1.2 mm value caused the
# explicit Solid/Convection/Conduction network to absorb substantially more heat
# than the original post-MOV steady model. A 2.0 mm wall is still a direct
# geometry assumption—not a heat-transfer efficiency or correction factor.
GRCOP_LINER_THICKNESS = 2.0e-3
GRCOP_MATERIAL_NAME = "GRCop-42"


def frustum_area(radius_1, radius_2, length):
    """Return the lateral area of a cylindrical or conical-frustum section.

    Parameters
    ----------
    radius_1, radius_2:
        Inlet and outlet radii in metres.
    length:
        Axial section length in metres.

    Returns
    -------
    float
        Gas-side lateral area in square metres.  Equal radii reduce naturally
        to the cylindrical-wall area ``2*pi*r*length``.
    """
    slant_length = math.hypot(length, radius_2 - radius_1)
    return math.pi * (radius_1 + radius_2) * slant_length



THERMAL_SECTIONS = (
    {
        "name": "Chamber Barrel",
        "area_ratio": CHAMBER_AREA_RATIO,
        "thermal_area_ratio": CHAMBER_AREA_RATIO,
        "length": 0.44,
        "radius_in": CHAMBER_RADIUS,
        "radius_out": CHAMBER_RADIUS,
    },
    {
        "name": "Converging",
        "area_ratio": CONVERGING_AREA_RATIO,
        "thermal_area_ratio": CONVERGING_AREA_RATIO,
        "length": 0.16,
        "radius_in": CHAMBER_RADIUS,
        "radius_out": THROAT_RADIUS,
    },
    {
        "name": "Throat",
        "area_ratio": 1.0,
        "thermal_area_ratio": 1.0,
        "length": 0.08,
        "radius_in": THROAT_RADIUS,
        "radius_out": THROAT_RADIUS,
    },
    {
        "name": "Upper Nozzle",
        "area_ratio": UPPER_NOZZLE_AREA_RATIO,
        "thermal_area_ratio": UPPER_NOZZLE_AREA_RATIO,
        "length": 0.48,
        "radius_in": THROAT_RADIUS,
        "radius_out": UPPER_NOZZLE_SPLIT_RADIUS,
    },
    {
        "name": "Exit Nozzle",
        "area_ratio": EXPANSION_RATIO,
        "thermal_area_ratio": EXIT_NOZZLE_THERMAL_AREA_RATIO,
        "length": 0.62,
        "radius_in": UPPER_NOZZLE_SPLIT_RADIUS,
        "radius_out": EXIT_RADIUS,
    },
)

for section in THERMAL_SECTIONS:
    section["hot_gas_area"] = frustum_area(
        section["radius_in"],
        section["radius_out"],
        section["length"],
    )


CHAMBER_THERMAL_SECTIONS = tuple(
    section
    for section in THERMAL_SECTIONS
    if section["name"] in {"Chamber Barrel", "Converging", "Throat"}
)
NOZZLE_THERMAL_SECTIONS = tuple(
    section
    for section in THERMAL_SECTIONS
    if section["name"] in {"Upper Nozzle", "Exit Nozzle"}
)

CHAMBER_REGEN_LENGTH = sum(
    section["length"] for section in CHAMBER_THERMAL_SECTIONS
)
NOZZLE_ONE_WAY_LENGTH = sum(
    section["length"] for section in NOZZLE_THERMAL_SECTIONS
)
NOZZLE_REGEN_LENGTH = 2.0 * NOZZLE_ONE_WAY_LENGTH

# Physical hot-wall length. The nozzle coolant travels this nozzle length twice,
# but the wall itself is counted only once.
TOTAL_COOLED_LENGTH = sum(section["length"] for section in THERMAL_SECTIONS)


# -----------------------------------------------------------------------------
# Solver starting guesses only
# -----------------------------------------------------------------------------
#
# Public examples benefit from a reliable near-trim initial state, especially
# when a dense finite-difference Jacobian is used.  Updating these values after
# changing hardware data is encouraged, but doing so does not prescribe the
# result: every listed State remains a nonlinear iteration variable.

# These values initialize the nonlinear solve. They are not targets or
# prescribed operating conditions. They are centered on the converged Raptor
# trim point printed by raptor.py so the solver starts close to the expected
# operating state while still solving every pressure, flow, temperature, and
# shaft speed from the connected network.
METHANE_PUMP_FLOW_GUESS = 147.646
LOX_PUMP_FLOW_GUESS = 532.610

METHANE_PUMP_INLET_PRESSURE_GUESS = 380_100.0
LOX_PUMP_INLET_PRESSURE_GUESS = 389_400.0

METHANE_ROTOR_SPEED_GUESS = 36_165.7
LOX_ROTOR_SPEED_GUESS = 30_067.7

METHANE_DISCHARGE_PRESSURE_GUESS = 89_900_000.0
LOX_DISCHARGE_PRESSURE_GUESS = 70_120_000.0
METHANE_MANIFOLD_TEMPERATURE_GUESS = 145.6
LOX_MANIFOLD_TEMPERATURE_GUESS = 111.9

# Regenerative-cooling guesses are ordered in the physical flow direction used
# in raptor.py. Outlet enthalpies were obtained from the converged outlet
# pressure and temperature through the same regen_methane_ph property map used
# by the network.
CHAMBER_REGEN_FLOW_GUESS = 67.583
CHAMBER_COOLANT_PRESSURE_GUESSES = (
    86_981_000.0,  # throat-channel outlet
    83_179_000.0,  # converging-channel outlet
    69_814_000.0,  # chamber-barrel return manifold
)
CHAMBER_COOLANT_ENTHALPY_GUESSES = (
    426_158.0,
    792_846.0,
    1_739_397.0,
)

NOZZLE_REGEN_FLOW_GUESS = 80.063
NOZZLE_COOLANT_PRESSURE_GUESSES = (
    85_657_000.0,  # upper-nozzle downflow outlet
    81_140_000.0,  # physical nozzle turnaround
    75_518_000.0,  # exit-nozzle upflow outlet
    69_814_000.0,  # upper-nozzle upflow / OPFV tap
)
NOZZLE_COOLANT_ENTHALPY_GUESSES = (
    665_328.0,
    1_019_125.0,
    1_357_303.0,
    1_730_533.0,
)

# Copper-wall guesses follow the same physical flow order. The nozzle report
# prints one averaged wall temperature for each gas-side station, so the same
# converged station value is used for its downflow and upflow passes.
CHAMBER_HOT_WALL_TEMPERATURE_GUESSES = (
    2_237.0,  # throat
    2_161.7,  # converging
    2_060.9,  # chamber barrel
)
CHAMBER_COLD_WALL_TEMPERATURE_GUESSES = (
    988.2,
    1_040.6,
    1_100.3,
)
NOZZLE_HOT_WALL_TEMPERATURE_GUESSES = (
    1_382.4,  # upper-nozzle downflow
    762.7,    # exit-nozzle downflow
    762.7,    # exit-nozzle upflow
    1_382.4,  # upper-nozzle upflow
)
NOZZLE_COLD_WALL_TEMPERATURE_GUESSES = (
    862.7,
    621.6,
    621.6,
    862.7,
)

CHAMBER_RETURN_FLOW_GUESS = 67.583
NOZZLE_RETURN_FLOW_GUESS = 71.532
REGEN_OUTLET_PRESSURE_GUESS = 69_814_000.0

# A representative converged value is used for each repeated branch family.
# Every section still solves its own friction factor through Churchill.
CHAMBER_REGEN_FRICTION_FACTOR_GUESS = 0.01416
NOZZLE_REGEN_FRICTION_FACTOR_GUESS = 0.01205
REGEN_RETURN_FRICTION_FACTOR_GUESS = 0.0100

MAIN_FUEL_VALVE_OUTLET_PRESSURE_GUESS = 88_640_000.0
OXYGEN_PREBURNER_FUEL_VALVE_OUTLET_PRESSURE_GUESS = 67_370_000.0
FUEL_PREBURNER_OXIDIZER_VALVE_OUTLET_PRESSURE_GUESS = 68_940_000.0
MAIN_OXIDIZER_VALVE_OUTLET_PRESSURE_GUESS = 69_630_000.0

FUEL_PREBURNER_PRESSURE_GUESS = 60_920_000.0
FUEL_PREBURNER_TEMPERATURE_GUESS = 876.4
OXYGEN_PREBURNER_PRESSURE_GUESS = 61_720_000.0
OXYGEN_PREBURNER_TEMPERATURE_GUESS = 785.4

FUEL_TURBINE_EXIT_PRESSURE_GUESS = 33_420_000.0
FUEL_TURBINE_EXIT_TEMPERATURE_GUESS = 818.8
OXYGEN_TURBINE_EXIT_PRESSURE_GUESS = 38_570_000.0
OXYGEN_TURBINE_EXIT_TEMPERATURE_GUESS = 716.2

CHAMBER_PRESSURE_GUESS = 29_948_000.0
CHAMBER_TEMPERATURE_GUESS = 3_810.9

# These internal nozzle pressure guesses were already centered on the same
# converged nozzle solution and are retained because the compact terminal report
# does not print their solved values individually.
CONVERGING_PRESSURE_GUESS = 24_955_700.0
THROAT_PRESSURE_GUESS = 17_250_500.0
UPPER_NOZZLE_PRESSURE_GUESS = 1_492_260.0
EXIT_NOZZLE_THERMAL_PRESSURE_GUESS = 187_250.0
EXIT_PRESSURE_GUESS = 95_500.0
