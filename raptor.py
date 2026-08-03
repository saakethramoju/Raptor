"""Solve a simplified Raptor-like full-flow staged-combustion engine.

This file is intentionally written as a complete FullFlow user example.  The
model is assembled in the same order that an engineer would follow the actual
engine flow path: tanks, prevalves, pumps, control valves, regenerative cooling,
preburners, turbines, main injectors, chamber, nozzle, and wall heat transfer.
Shared :class:`fullflow.State` objects connect component outlets to downstream
component inlets, so the code itself is the network schematic.

Modeling objective
------------------
The goal is a transparent, steady, one-dimensional cycle model—not an exact
reconstruction of proprietary flight hardware.  Only the tank states and
ambient pressure are operating boundary conditions.  Pump and turbine maps,
valve and injector geometry, cooling-channel geometry, chamber geometry, and
nozzle area ratios are component data.  The network solves the engine mass
flows, pressures, shaft speeds, regenerative split, preburner mixture ratios,
turbine discharge states, chamber pressure, nozzle flow, and wall temperatures.

Full-flow topology represented here
-----------------------------------
* Separate methane and oxygen turbopumps are driven by separate turbines.
* Regeneratively heated methane feeds the fuel-rich preburner, with a hot
  nozzle-return tap feeding the oxygen-rich preburner through the OPFV.
* Most LOX feeds the oxygen-rich preburner; a smaller branch feeds the
  fuel-rich preburner through the FPOV.
* Every preburner product stream passes through its turbine and then enters the
  main chamber.  No propellant is dumped overboard.
* The methane cooling system is an explicit chamber branch plus a nozzle
  downflow/turnaround/upflow branch with solved return lines and a solved merge
  pressure.

Steady-state closure
--------------------
The model avoids redundant solve variables.  Exact steady adiabatic enthalpy
relations are propagated directly through pump-discharge manifolds, valves,
preburners, turbine exhaust plenums, the hot-methane merge, and the main
chamber.  Coolant enthalpies remain independent unknowns because each cooling
section receives a solved wall heat rate.

Only three general-purpose ``Balance`` components are required:

1. fuel-rich turbine torque equals methane-pump torque;
2. oxygen-rich turbine torque equals LOX-pump torque;
3. throat Mach number equals one.

All other equations come from ordinary FullFlow components: valve/orifice
relations, pump and turbine maps, FlowTube momentum equations, algebraic node
continuity, coolant energy balances, wall heat balances, and equilibrium nozzle
stations.  No state bounds, prescribed engine flow rates, imposed chamber or
preburner pressures, mechanical-efficiency multipliers, or hidden tuning
balances are used.

Regenerative-cooling treatment
------------------------------
Each cooling section contains a real-fluid pressure-enthalpy methane map, a
Churchill friction model, a ``FlowTube`` momentum relation, and an algebraic
mass/energy node.  A two-temperature GRCop-42 wall is assembled from gas-side
convection, copper conduction, and coolant-side convection.  The thermal model
is one-way coupled at cycle level: heat is added to the wall/coolant network,
but it is not subtracted from the equilibrium chamber/nozzle gas expansion.

Running the example
-------------------
The repository includes ``raptor.h5`` with all pre-generated property maps, so
normal use is simply::

    python3 raptor.py

Map generation is intentionally separate because the equilibrium maps are
expensive to build.  See ``raptor_maps.py`` and ``README.md`` before rebuilding
them.  ``diagram.py`` contains a detailed text flow diagram of the saved trim
solution.

All calculations use SI units internally.  The terminal report converts
pressures to bar and powers to MW only for readability.
"""

import math
import time

from fullflow import (
    AdiabaticFlow,
    AdiabaticWallTemperature,
    Balance,
    Bartz,
    Churchill,
    CompressibleOrifice,
    ConstantDensityPump,
    Conduction,
    Convection,
    FlowTube,
    DischargeCoefficient,
    GasTurbine,
    Gnielinski,
    Lookup,
    Map,
    Network,
    Solid,
    State,
    SteadyState,
    TemperatureRecoveryFactor,
    EckertReferenceTemperature,
    Volume,
)
from thermoprop import Fluid

from raptor_data import (
    AMBIENT_PRESSURE,
    BAR,
    CHAMBER_AREA,
    CHAMBER_DIAMETER,
    CHAMBER_PRESSURE_GUESS,
    CONVERGING_AREA_RATIO,
    CONVERGING_PRESSURE_GUESS,
    CHAMBER_COOLING_CHANNEL_COUNT,
    CHAMBER_COOLING_CHANNEL_HEIGHT,
    CHAMBER_COOLING_CHANNEL_WIDTH,
    CHAMBER_COOLING_FLOW_AREA,
    CHAMBER_COOLING_HYDRAULIC_DIAMETER,
    CHAMBER_REGEN_FLOW_GUESS,
    CHAMBER_REGEN_FRICTION_FACTOR_GUESS,
    CHAMBER_REGEN_LENGTH,
    CHAMBER_COOLANT_PRESSURE_GUESSES,
    CHAMBER_COOLANT_ENTHALPY_GUESSES,
    CHAMBER_HOT_WALL_TEMPERATURE_GUESSES,
    CHAMBER_COLD_WALL_TEMPERATURE_GUESSES,
    CHAMBER_COOLING_ROUGHNESS,
    CHAMBER_RETURN_FLOW_GUESS,
    CHAMBER_RETURN_PIPE_LENGTH,
    CHAMBER_THERMAL_SECTIONS,
    EXIT_AREA,
    EXIT_NOZZLE_THERMAL_AREA_RATIO,
    EXIT_NOZZLE_THERMAL_PRESSURE_GUESS,
    EXPANSION_RATIO,
    EXIT_PRESSURE_GUESS,
    FILE_NAME,
    FUEL_PREBURNER_LOX_AREA,
    FUEL_PREBURNER_METHANE_AREA,
    FUEL_PREBURNER_OXIDIZER_VALVE_AREA,
    FUEL_PREBURNER_OXIDIZER_VALVE_CD,
    FUEL_PREBURNER_OXIDIZER_VALVE_OUTLET_PRESSURE_GUESS,
    FUEL_PREBURNER_PRESSURE_GUESS,
    FUEL_TURBINE_EXIT_PRESSURE_GUESS,
    G0,
    GAS_INJECTOR_CD,
    GRCOP_LINER_THICKNESS,
    LIQUID_INJECTOR_CD,
    LOX_DISCHARGE_PRESSURE_GUESS,
    LOX_PUMP_FLOW_GUESS,
    LOX_PUMP_INLET_PRESSURE_GUESS,
    LOX_PUMP_MAP_SPEED,
    LOX_ROTOR_SPEED_GUESS,
    LOX_TANK_PRESSURE,
    LOX_TANK_TEMPERATURE,
    MAIN_LOX_VALVE_AREA,
    MAIN_LOX_VALVE_CD,
    MAIN_OXIDIZER_VALVE_AREA,
    MAIN_OXIDIZER_VALVE_CD,
    MAIN_OXIDIZER_VALVE_OUTLET_PRESSURE_GUESS,
    MAIN_METHANE_VALVE_AREA,
    MAIN_METHANE_VALVE_CD,
    MAIN_FUEL_INJECTOR_AREA,
    MAIN_FUEL_VALVE_AREA,
    MAIN_FUEL_VALVE_CD,
    MAIN_FUEL_VALVE_OUTLET_PRESSURE_GUESS,
    MAIN_OXYGEN_INJECTOR_AREA,
    METHANE_DISCHARGE_PRESSURE_GUESS,
    METHANE_PUMP_FLOW_GUESS,
    METHANE_PUMP_INLET_PRESSURE_GUESS,
    METHANE_PUMP_MAP_SPEED,
    METHANE_ROTOR_SPEED_GUESS,
    METHANE_TANK_PRESSURE,
    METHANE_TANK_TEMPERATURE,
    OXYGEN_PREBURNER_FUEL_VALVE_AREA,
    OXYGEN_PREBURNER_FUEL_VALVE_CD,
    OXYGEN_PREBURNER_FUEL_VALVE_OUTLET_PRESSURE_GUESS,
    OXYGEN_PREBURNER_LOX_AREA,
    OXYGEN_PREBURNER_METHANE_AREA,
    OXYGEN_PREBURNER_PRESSURE_GUESS,
    OXYGEN_TURBINE_EXIT_PRESSURE_GUESS,
    NOZZLE_COOLING_CHANNEL_COUNT,
    NOZZLE_COOLING_CHANNEL_HEIGHT,
    NOZZLE_COOLING_CHANNEL_WIDTH,
    NOZZLE_COOLING_HYDRAULIC_DIAMETER,
    NOZZLE_COOLING_PASS_CHANNEL_COUNT,
    NOZZLE_COOLING_PASS_FLOW_AREA,
    NOZZLE_ONE_WAY_LENGTH,
    NOZZLE_REGEN_FLOW_GUESS,
    NOZZLE_REGEN_FRICTION_FACTOR_GUESS,
    NOZZLE_REGEN_LENGTH,
    NOZZLE_COOLANT_PRESSURE_GUESSES,
    NOZZLE_COOLANT_ENTHALPY_GUESSES,
    NOZZLE_HOT_WALL_TEMPERATURE_GUESSES,
    NOZZLE_COLD_WALL_TEMPERATURE_GUESSES,
    NOZZLE_COOLING_ROUGHNESS,
    NOZZLE_RETURN_FLOW_GUESS,
    NOZZLE_RETURN_PIPE_LENGTH,
    NOZZLE_THERMAL_SECTIONS,
    REGEN_OUTLET_PRESSURE_GUESS,
    REGEN_RETURN_FRICTION_FACTOR_GUESS,
    REGEN_RETURN_PIPE_DIAMETER,
    REGEN_RETURN_PIPE_FLOW_AREA,
    REGEN_RETURN_PIPE_ROUGHNESS,
    THERMAL_SECTIONS,
    THROAT_AREA,
    THROAT_PRESSURE_GUESS,
    TOTAL_COOLED_LENGTH,
    UPPER_NOZZLE_AREA_RATIO,
    UPPER_NOZZLE_SPLIT_AREA_RATIO,
    UPPER_NOZZLE_PRESSURE_GUESS,
)


program_start_time = time.perf_counter()

# A FullFlow network is a container for components and shared State objects.
# Components are registered with the network when they are constructed.  The
# model therefore remains a normal Python script rather than a separate graph
# definition or configuration file.
Raptor = Network("Raptor")


# -----------------------------------------------------------------------------
# Tank boundary conditions
# -----------------------------------------------------------------------------

methane_tank_pressure = State(METHANE_TANK_PRESSURE)
lox_tank_pressure = State(LOX_TANK_PRESSURE)
ambient_pressure = State(AMBIENT_PRESSURE)

MethaneTank = Lookup(
    "Methane Tank",
    Raptor,
    Fluid,
    "Methane",
    pressure=methane_tank_pressure,
    temperature=METHANE_TANK_TEMPERATURE,
)
LOXTank = Lookup(
    "LOX Tank",
    Raptor,
    Fluid,
    "Oxygen",
    pressure=lox_tank_pressure,
    temperature=LOX_TANK_TEMPERATURE,
)


# -----------------------------------------------------------------------------
# Main propellant valves and pump inlet nodes
# -----------------------------------------------------------------------------

# These are numerical starting values. The pump components solve the flows from
# their maps and connected pressures; the inlet nodes enforce valve/pump
# continuity. No mass flow is prescribed.
methane_pump_mass_flow = State(METHANE_PUMP_FLOW_GUESS)
lox_pump_mass_flow = State(LOX_PUMP_FLOW_GUESS)

methane_pump_inlet_pressure = State(METHANE_PUMP_INLET_PRESSURE_GUESS)
lox_pump_inlet_pressure = State(LOX_PUMP_INLET_PRESSURE_GUESS)

MethanePumpInletFluid = Lookup(
    "Methane Pump Inlet Fluid",
    Raptor,
    Fluid,
    "Methane",
    pressure=methane_pump_inlet_pressure,
    temperature=METHANE_TANK_TEMPERATURE,
)
LOXPumpInletFluid = Lookup(
    "LOX Pump Inlet Fluid",
    Raptor,
    Fluid,
    "Oxygen",
    pressure=lox_pump_inlet_pressure,
    temperature=LOX_TANK_TEMPERATURE,
)

MethaneMainValve = DischargeCoefficient(
    "Methane Tank Prevalve",
    Raptor,
    upstream_pressure=methane_tank_pressure,
    downstream_pressure=methane_pump_inlet_pressure,
    density=MethaneTank.density,
    discharge_coefficient=MAIN_METHANE_VALVE_CD,
    cross_sectional_area=MAIN_METHANE_VALVE_AREA,
)
LOXMainValve = DischargeCoefficient(
    "LOX Tank Prevalve",
    Raptor,
    upstream_pressure=lox_tank_pressure,
    downstream_pressure=lox_pump_inlet_pressure,
    density=LOXTank.density,
    discharge_coefficient=MAIN_LOX_VALVE_CD,
    cross_sectional_area=MAIN_LOX_VALVE_AREA,
)

MethanePumpInletNode = Volume(
    "Methane Pump Inlet Node",
    Raptor,
    pressure=methane_pump_inlet_pressure,
    mass_flow_in=MethaneMainValve.mass_flow,
    mass_flow_out=methane_pump_mass_flow,
)
LOXPumpInletNode = Volume(
    "LOX Pump Inlet Node",
    Raptor,
    pressure=lox_pump_inlet_pressure,
    mass_flow_in=LOXMainValve.mass_flow,
    mass_flow_out=lox_pump_mass_flow,
)


# -----------------------------------------------------------------------------
# Pump states, maps, and components
# -----------------------------------------------------------------------------

methane_rotor_speed = State(METHANE_ROTOR_SPEED_GUESS)
lox_rotor_speed = State(LOX_ROTOR_SPEED_GUESS)

methane_discharge_pressure = State(METHANE_DISCHARGE_PRESSURE_GUESS)
lox_discharge_pressure = State(LOX_DISCHARGE_PRESSURE_GUESS)

MethanePumpMap = Map.from_hdf5(
    "Methane Pump Map",
    Raptor,
    filename=FILE_NAME,
    group="methane_pump",
    inputs={
        "rotor_speed": methane_rotor_speed,
        "volumetric_flow": methane_pump_mass_flow / MethanePumpInletFluid.density,
    },
    extrapolate=True,
)
LOXPumpMap = Map.from_hdf5(
    "LOX Pump Map",
    Raptor,
    filename=FILE_NAME,
    group="lox_pump",
    inputs={
        "rotor_speed": lox_rotor_speed,
        "volumetric_flow": lox_pump_mass_flow / LOXPumpInletFluid.density,
    },
    extrapolate=True,
)

MethanePump = ConstantDensityPump(
    "Methane Pump",
    Raptor,
    mass_flow=methane_pump_mass_flow,
    rotor_speed=methane_rotor_speed,
    head_rise=MethanePumpMap.head_rise,
    density=MethanePumpInletFluid.density,
    torque=MethanePumpMap.torque,
    upstream_pressure=methane_pump_inlet_pressure,
    discharge_pressure=methane_discharge_pressure,
    upstream_total_enthalpy=MethanePumpInletFluid.enthalpy,
)
LOXPump = ConstantDensityPump(
    "LOX Pump",
    Raptor,
    mass_flow=lox_pump_mass_flow,
    rotor_speed=lox_rotor_speed,
    head_rise=LOXPumpMap.head_rise,
    density=LOXPumpInletFluid.density,
    torque=LOXPumpMap.torque,
    upstream_pressure=lox_pump_inlet_pressure,
    discharge_pressure=lox_discharge_pressure,
    upstream_total_enthalpy=LOXPumpInletFluid.enthalpy,
)

# Pump discharge enthalpy follows directly from the pump component. The
# discharge manifolds solve only the pressure required by flow continuity.
MethaneManifoldFluid = Lookup(
    "Methane Discharge Manifold Fluid",
    Raptor,
    Fluid,
    "Methane",
    pressure=methane_discharge_pressure,
    enthalpy=MethanePump.discharge_total_enthalpy,
)
LOXManifoldFluid = Lookup(
    "LOX Discharge Manifold Fluid",
    Raptor,
    Fluid,
    "Oxygen",
    pressure=lox_discharge_pressure,
    enthalpy=LOXPump.discharge_total_enthalpy,
)


# -----------------------------------------------------------------------------
# Engine control valves and regenerative-cooling network
# -----------------------------------------------------------------------------

# These pressure states are ordinary network unknowns. Valve areas and channel
# geometry are fixed component properties; no valve flow or cooling-branch flow
# is prescribed.
main_fuel_valve_outlet_pressure = State(
    MAIN_FUEL_VALVE_OUTLET_PRESSURE_GUESS
)
oxygen_preburner_fuel_valve_outlet_pressure = State(
    OXYGEN_PREBURNER_FUEL_VALVE_OUTLET_PRESSURE_GUESS
)
fuel_preburner_oxidizer_valve_outlet_pressure = State(
    FUEL_PREBURNER_OXIDIZER_VALVE_OUTLET_PRESSURE_GUESS
)
main_oxidizer_valve_outlet_pressure = State(
    MAIN_OXIDIZER_VALVE_OUTLET_PRESSURE_GUESS
)

MainFuelValve = DischargeCoefficient(
    "Main Fuel Valve",
    Raptor,
    upstream_pressure=methane_discharge_pressure,
    downstream_pressure=main_fuel_valve_outlet_pressure,
    density=MethaneManifoldFluid.density,
    discharge_coefficient=MAIN_FUEL_VALVE_CD,
    cross_sectional_area=MAIN_FUEL_VALVE_AREA,
)
FuelPreburnerOxidizerValve = DischargeCoefficient(
    "Fuel Preburner Oxidizer Valve",
    Raptor,
    upstream_pressure=lox_discharge_pressure,
    downstream_pressure=fuel_preburner_oxidizer_valve_outlet_pressure,
    density=LOXManifoldFluid.density,
    discharge_coefficient=FUEL_PREBURNER_OXIDIZER_VALVE_CD,
    cross_sectional_area=FUEL_PREBURNER_OXIDIZER_VALVE_AREA,
)
MainOxidizerValve = DischargeCoefficient(
    "Main Oxidizer Valve",
    Raptor,
    upstream_pressure=lox_discharge_pressure,
    downstream_pressure=main_oxidizer_valve_outlet_pressure,
    density=LOXManifoldFluid.density,
    discharge_coefficient=MAIN_OXIDIZER_VALVE_CD,
    cross_sectional_area=MAIN_OXIDIZER_VALVE_AREA,
)

# The liquid valves are adiabatic throttles. Their outlet temperatures are
# obtained from pressure and unchanged specific enthalpy.
MainFuelValveOutletFluid = Lookup(
    "Main Fuel Valve Outlet Fluid",
    Raptor,
    Fluid,
    "Methane",
    pressure=main_fuel_valve_outlet_pressure,
    enthalpy=MethaneManifoldFluid.enthalpy,
)
FuelPreburnerOxidizerValveOutletFluid = Lookup(
    "Fuel Preburner Oxidizer Valve Outlet Fluid",
    Raptor,
    Fluid,
    "Oxygen",
    pressure=fuel_preburner_oxidizer_valve_outlet_pressure,
    enthalpy=LOXManifoldFluid.enthalpy,
)
MainOxidizerValveOutletFluid = Lookup(
    "Main Oxidizer Valve Outlet Fluid",
    Raptor,
    Fluid,
    "Oxygen",
    pressure=main_oxidizer_valve_outlet_pressure,
    enthalpy=LOXManifoldFluid.enthalpy,
)

# The chamber branch flows from the throat toward the injector end. The nozzle
# branch flows down the upper and exit nozzle, turns at the exit, then returns
# through neighboring channels. Each listed coolant node is the outlet of one
# physical channel section.
chamber_path_sections = tuple(reversed(CHAMBER_THERMAL_SECTIONS))
nozzle_path_sections = (
    [("Nozzle Downflow", section) for section in NOZZLE_THERMAL_SECTIONS]
    + [
        ("Nozzle Upflow", section)
        for section in reversed(NOZZLE_THERMAL_SECTIONS)
    ]
)

chamber_coolant_nodes = []
for section, pressure_guess, enthalpy_guess in zip(
    chamber_path_sections,
    CHAMBER_COOLANT_PRESSURE_GUESSES,
    CHAMBER_COOLANT_ENTHALPY_GUESSES,
):
    coolant_pressure = State(pressure_guess)
    coolant_enthalpy = State(enthalpy_guess)
    Coolant = Map.from_hdf5(
        f"Chamber {section['name']} Outlet Methane",
        Raptor,
        filename=FILE_NAME,
        group="regen_methane_ph",
        inputs={
            "pressure": coolant_pressure,
            "target_enthalpy": coolant_enthalpy,
        },
        outputs=["density", "temperature"],
        extrapolate=True,
    )
    chamber_coolant_nodes.append(
        {
            "pass_name": "Chamber",
            "section": section,
            "pressure": coolant_pressure,
            "enthalpy": coolant_enthalpy,
            "fluid": Coolant,
        }
    )

nozzle_coolant_nodes = []
for (pass_name, section), pressure_guess, enthalpy_guess in zip(
    nozzle_path_sections,
    NOZZLE_COOLANT_PRESSURE_GUESSES,
    NOZZLE_COOLANT_ENTHALPY_GUESSES,
):
    coolant_pressure = State(pressure_guess)
    coolant_enthalpy = State(enthalpy_guess)
    Coolant = Map.from_hdf5(
        f"{pass_name} {section['name']} Outlet Methane",
        Raptor,
        filename=FILE_NAME,
        group="regen_methane_ph",
        inputs={
            "pressure": coolant_pressure,
            "target_enthalpy": coolant_enthalpy,
        },
        outputs=["density", "temperature"],
        extrapolate=True,
    )
    nozzle_coolant_nodes.append(
        {
            "pass_name": pass_name,
            "section": section,
            "pressure": coolant_pressure,
            "enthalpy": coolant_enthalpy,
            "fluid": Coolant,
        }
    )

# One independent mass-flow state is used for every hydraulic branch. The
# intervening algebraic nodes enforce continuity instead of forcing the same
# State into several momentum equations.
chamber_segment_flows = [
    State(CHAMBER_REGEN_FLOW_GUESS)
    for _ in chamber_coolant_nodes
]
nozzle_segment_flows = [
    State(NOZZLE_REGEN_FLOW_GUESS)
    for _ in nozzle_coolant_nodes
]
chamber_return_mass_flow = State(CHAMBER_RETURN_FLOW_GUESS)
nozzle_return_to_merge_mass_flow = State(NOZZLE_RETURN_FLOW_GUESS)

chamber_regen_mass_flow = chamber_segment_flows[0]
nozzle_regen_mass_flow = nozzle_segment_flows[0]
regen_supply_mass_flow = chamber_regen_mass_flow + nozzle_regen_mass_flow

# Build each chamber cooling section from its actual inlet and outlet states.
# The midpoint property lookup is a control-volume average derived from solved
# endpoint states; it is not an imposed pressure or enthalpy profile.
#
# FlowTube is retained because it handles changing real-fluid density through
# the supercritical cooling circuit. In a SteadyState solve, FullFlow drives its
# mass-flow derivative to zero, which is exactly the steady momentum equation;
# no line state is time-integrated in this model.
chamber_regen_pipes = []
for index, node in enumerate(chamber_coolant_nodes):
    upstream_pressure = (
        main_fuel_valve_outlet_pressure
        if index == 0
        else chamber_coolant_nodes[index - 1]["pressure"]
    )
    upstream_enthalpy = (
        MainFuelValveOutletFluid.enthalpy
        if index == 0
        else chamber_coolant_nodes[index - 1]["enthalpy"]
    )
    UpstreamFluid = (
        MainFuelValveOutletFluid
        if index == 0
        else chamber_coolant_nodes[index - 1]["fluid"]
    )

    MeanFluid = Map.from_hdf5(
        f"Chamber {node['section']['name']} Mean Methane",
        Raptor,
        filename=FILE_NAME,
        group="regen_methane_ph",
        inputs={
            "pressure": 0.5 * (upstream_pressure + node["pressure"]),
            "target_enthalpy": 0.5 * (
                upstream_enthalpy + node["enthalpy"]
            ),
        },
        outputs=[
            "temperature",
            "conductivity",
            "specific_heat_cp",
            "dynamic_viscosity",
        ],
        extrapolate=True,
    )
    friction_factor = State(CHAMBER_REGEN_FRICTION_FACTOR_GUESS)
    Friction = Churchill(
        f"Chamber {node['section']['name']} Friction",
        Raptor,
        mass_flow=chamber_segment_flows[index],
        friction_factor=friction_factor,
        hydraulic_diameter=CHAMBER_COOLING_HYDRAULIC_DIAMETER,
        dynamic_viscosity=MeanFluid.dynamic_viscosity,
        cross_sectional_area=CHAMBER_COOLING_FLOW_AREA,
        roughness=CHAMBER_COOLING_ROUGHNESS,
    )
    Pipe = FlowTube(
        f"Chamber {node['section']['name']} Channels",
        Raptor,
        mass_flow=chamber_segment_flows[index],
        upstream_static_pressure=upstream_pressure,
        downstream_static_pressure=node["pressure"],
        length=node["section"]["length"],
        hydraulic_diameter=CHAMBER_COOLING_HYDRAULIC_DIAMETER,
        cross_sectional_area=CHAMBER_COOLING_FLOW_AREA,
        upstream_density=UpstreamFluid.density,
        downstream_density=node["fluid"].density,
        friction_factor=friction_factor,
        upstream_static_enthalpy=upstream_enthalpy,
    )
    node["mean_fluid"] = MeanFluid
    node["friction"] = Friction
    node["pipe"] = Pipe
    chamber_regen_pipes.append(Pipe)

# Build the four nozzle channel sections, including a real outlet node at the
# physical downflow-to-upflow turnaround.
nozzle_regen_pipes = []
for index, node in enumerate(nozzle_coolant_nodes):
    upstream_pressure = (
        main_fuel_valve_outlet_pressure
        if index == 0
        else nozzle_coolant_nodes[index - 1]["pressure"]
    )
    upstream_enthalpy = (
        MainFuelValveOutletFluid.enthalpy
        if index == 0
        else nozzle_coolant_nodes[index - 1]["enthalpy"]
    )
    UpstreamFluid = (
        MainFuelValveOutletFluid
        if index == 0
        else nozzle_coolant_nodes[index - 1]["fluid"]
    )

    MeanFluid = Map.from_hdf5(
        f"{node['pass_name']} {node['section']['name']} Mean Methane",
        Raptor,
        filename=FILE_NAME,
        group="regen_methane_ph",
        inputs={
            "pressure": 0.5 * (upstream_pressure + node["pressure"]),
            "target_enthalpy": 0.5 * (
                upstream_enthalpy + node["enthalpy"]
            ),
        },
        outputs=[
            "temperature",
            "conductivity",
            "specific_heat_cp",
            "dynamic_viscosity",
        ],
        extrapolate=True,
    )
    friction_factor = State(NOZZLE_REGEN_FRICTION_FACTOR_GUESS)
    Friction = Churchill(
        f"{node['pass_name']} {node['section']['name']} Friction",
        Raptor,
        mass_flow=nozzle_segment_flows[index],
        friction_factor=friction_factor,
        hydraulic_diameter=NOZZLE_COOLING_HYDRAULIC_DIAMETER,
        dynamic_viscosity=MeanFluid.dynamic_viscosity,
        cross_sectional_area=NOZZLE_COOLING_PASS_FLOW_AREA,
        roughness=NOZZLE_COOLING_ROUGHNESS,
    )
    Pipe = FlowTube(
        f"{node['pass_name']} {node['section']['name']} Channels",
        Raptor,
        mass_flow=nozzle_segment_flows[index],
        upstream_static_pressure=upstream_pressure,
        downstream_static_pressure=node["pressure"],
        length=node["section"]["length"],
        hydraulic_diameter=NOZZLE_COOLING_HYDRAULIC_DIAMETER,
        cross_sectional_area=NOZZLE_COOLING_PASS_FLOW_AREA,
        upstream_density=UpstreamFluid.density,
        downstream_density=node["fluid"].density,
        friction_factor=friction_factor,
        upstream_static_enthalpy=upstream_enthalpy,
    )
    node["mean_fluid"] = MeanFluid
    node["friction"] = Friction
    node["pipe"] = Pipe
    nozzle_regen_pipes.append(Pipe)

ChamberReturnManifoldFluid = chamber_coolant_nodes[-1]["fluid"]
NozzleReturnTapFluid = nozzle_coolant_nodes[-1]["fluid"]
chamber_return_enthalpy = chamber_coolant_nodes[-1]["enthalpy"]
nozzle_return_enthalpy = nozzle_coolant_nodes[-1]["enthalpy"]
nozzle_turnaround_pressure = nozzle_coolant_nodes[1]["pressure"]

# The OPFV branches from the solved nozzle-return tap before the remainder of
# the nozzle flow reaches the hot-methane merge header.
OxygenPreburnerFuelValve = DischargeCoefficient(
    "Oxygen Preburner Fuel Valve",
    Raptor,
    upstream_pressure=nozzle_coolant_nodes[-1]["pressure"],
    downstream_pressure=oxygen_preburner_fuel_valve_outlet_pressure,
    density=NozzleReturnTapFluid.density,
    discharge_coefficient=OXYGEN_PREBURNER_FUEL_VALVE_CD,
    cross_sectional_area=OXYGEN_PREBURNER_FUEL_VALVE_AREA,
)
OxygenPreburnerFuelValveOutletFluid = Lookup(
    "Oxygen Preburner Fuel Valve Outlet Fluid",
    Raptor,
    Fluid,
    "Methane",
    pressure=oxygen_preburner_fuel_valve_outlet_pressure,
    enthalpy=nozzle_return_enthalpy,
)

# The two branch returns are explicit FlowTubes. Separate branch property
# lookups at the merge pressure preserve each stream's enthalpy before mixing.
regen_outlet_pressure = State(REGEN_OUTLET_PRESSURE_GUESS)

ChamberReturnAtMergeFluid = Map.from_hdf5(
    "Chamber Return Methane At Merge Pressure",
    Raptor,
    filename=FILE_NAME,
    group="regen_methane_ph",
    inputs={
        "pressure": regen_outlet_pressure,
        "target_enthalpy": chamber_return_enthalpy,
    },
    outputs=["density"],
    extrapolate=True,
)
NozzleReturnAtMergeFluid = Map.from_hdf5(
    "Nozzle Return Methane At Merge Pressure",
    Raptor,
    filename=FILE_NAME,
    group="regen_methane_ph",
    inputs={
        "pressure": regen_outlet_pressure,
        "target_enthalpy": nozzle_return_enthalpy,
    },
    outputs=["density"],
    extrapolate=True,
)

ChamberReturnMeanFluid = Map.from_hdf5(
    "Chamber Return Pipe Mean Methane",
    Raptor,
    filename=FILE_NAME,
    group="regen_methane_ph",
    inputs={
        "pressure": 0.5 * (
            chamber_coolant_nodes[-1]["pressure"]
            + regen_outlet_pressure
        ),
        "target_enthalpy": chamber_return_enthalpy,
    },
    outputs=["dynamic_viscosity"],
    extrapolate=True,
)
NozzleReturnMeanFluid = Map.from_hdf5(
    "Nozzle Return Pipe Mean Methane",
    Raptor,
    filename=FILE_NAME,
    group="regen_methane_ph",
    inputs={
        "pressure": 0.5 * (
            nozzle_coolant_nodes[-1]["pressure"]
            + regen_outlet_pressure
        ),
        "target_enthalpy": nozzle_return_enthalpy,
    },
    outputs=["dynamic_viscosity"],
    extrapolate=True,
)

chamber_return_friction_factor = State(
    REGEN_RETURN_FRICTION_FACTOR_GUESS
)
nozzle_return_friction_factor = State(
    REGEN_RETURN_FRICTION_FACTOR_GUESS
)
ChamberReturnFriction = Churchill(
    "Chamber Return Pipe Friction",
    Raptor,
    mass_flow=chamber_return_mass_flow,
    friction_factor=chamber_return_friction_factor,
    hydraulic_diameter=REGEN_RETURN_PIPE_DIAMETER,
    dynamic_viscosity=ChamberReturnMeanFluid.dynamic_viscosity,
    cross_sectional_area=REGEN_RETURN_PIPE_FLOW_AREA,
    roughness=REGEN_RETURN_PIPE_ROUGHNESS,
)
NozzleReturnFriction = Churchill(
    "Nozzle Return Pipe Friction",
    Raptor,
    mass_flow=nozzle_return_to_merge_mass_flow,
    friction_factor=nozzle_return_friction_factor,
    hydraulic_diameter=REGEN_RETURN_PIPE_DIAMETER,
    dynamic_viscosity=NozzleReturnMeanFluid.dynamic_viscosity,
    cross_sectional_area=REGEN_RETURN_PIPE_FLOW_AREA,
    roughness=REGEN_RETURN_PIPE_ROUGHNESS,
)
ChamberReturnPipe = FlowTube(
    "Chamber Return Pipe",
    Raptor,
    mass_flow=chamber_return_mass_flow,
    upstream_static_pressure=chamber_coolant_nodes[-1]["pressure"],
    downstream_static_pressure=regen_outlet_pressure,
    length=CHAMBER_RETURN_PIPE_LENGTH,
    hydraulic_diameter=REGEN_RETURN_PIPE_DIAMETER,
    cross_sectional_area=REGEN_RETURN_PIPE_FLOW_AREA,
    upstream_density=ChamberReturnManifoldFluid.density,
    downstream_density=ChamberReturnAtMergeFluid.density,
    friction_factor=chamber_return_friction_factor,
    upstream_static_enthalpy=chamber_return_enthalpy,
)
NozzleReturnPipe = FlowTube(
    "Nozzle Return Pipe",
    Raptor,
    mass_flow=nozzle_return_to_merge_mass_flow,
    upstream_static_pressure=nozzle_coolant_nodes[-1]["pressure"],
    downstream_static_pressure=regen_outlet_pressure,
    length=NOZZLE_RETURN_PIPE_LENGTH,
    hydraulic_diameter=REGEN_RETURN_PIPE_DIAMETER,
    cross_sectional_area=REGEN_RETURN_PIPE_FLOW_AREA,
    upstream_density=NozzleReturnTapFluid.density,
    downstream_density=NozzleReturnAtMergeFluid.density,
    friction_factor=nozzle_return_friction_factor,
    upstream_static_enthalpy=nozzle_return_enthalpy,
)

regen_merge_mass_flow = (
    chamber_return_mass_flow + nozzle_return_to_merge_mass_flow
)
regen_mixed_inlet_enthalpy = (
    chamber_return_mass_flow * chamber_return_enthalpy
    + nozzle_return_to_merge_mass_flow * nozzle_return_enthalpy
) / regen_merge_mass_flow

RegenOutletFluid = Map.from_hdf5(
    "Merged Hot Methane",
    Raptor,
    filename=FILE_NAME,
    group="regen_methane_ph",
    inputs={
        "pressure": regen_outlet_pressure,
        "target_enthalpy": regen_mixed_inlet_enthalpy,
    },
    outputs=["density", "temperature"],
    extrapolate=True,
)

# Each cooling section is a solved real-fluid pressure-enthalpy node. Heat rates
# are attached after the gas-side nozzle states and local heat-transfer
# coefficients have been created below.
for index, node in enumerate(chamber_coolant_nodes):
    inlet_enthalpy = (
        MainFuelValveOutletFluid.enthalpy
        if index == 0
        else chamber_coolant_nodes[index - 1]["enthalpy"]
    )
    outlet_flow = (
        chamber_segment_flows[index + 1]
        if index + 1 < len(chamber_segment_flows)
        else chamber_return_mass_flow
    )
    node["volume"] = Volume(
        f"Chamber {node['section']['name']} Coolant Node",
        Raptor,
        pressure=node["pressure"],
        enthalpy=node["enthalpy"],
        total_enthalpy_in=inlet_enthalpy,
        mass_flow_in=chamber_segment_flows[index],
        mass_flow_out=outlet_flow,
        energy_variable="enthalpy",
    )

for index, node in enumerate(nozzle_coolant_nodes):
    inlet_enthalpy = (
        MainFuelValveOutletFluid.enthalpy
        if index == 0
        else nozzle_coolant_nodes[index - 1]["enthalpy"]
    )
    if index + 1 < len(nozzle_segment_flows):
        outlet_flow = nozzle_segment_flows[index + 1]
    else:
        outlet_flow = (
            OxygenPreburnerFuelValve.mass_flow
            + nozzle_return_to_merge_mass_flow
        )
    node["volume"] = Volume(
        f"{node['pass_name']} {node['section']['name']} Coolant Node",
        Raptor,
        pressure=node["pressure"],
        enthalpy=node["enthalpy"],
        total_enthalpy_in=inlet_enthalpy,
        mass_flow_in=nozzle_segment_flows[index],
        mass_flow_out=outlet_flow,
        energy_variable="enthalpy",
    )


# -----------------------------------------------------------------------------
# Reactant enthalpy lookups
# -----------------------------------------------------------------------------

# Fluid supplies the real-fluid hydraulic states. Combustion energy balances
# use the CEA reference carried by Propellant and Equilibrium. Both methane
# reactant streams are hot, post-regen methane, so their CEA enthalpies use the
# temperature-only gas-species map. The LOX reactants retain their actual
# post-valve pressure and temperature.
HotMethaneReactant = Map.from_hdf5(
    "Hot Methane Reactant Enthalpy",
    Raptor,
    filename=FILE_NAME,
    group="hot_methane_reactant_cea_t",
    inputs={
        "temperature": RegenOutletFluid.temperature,
    },
    extrapolate=True,
)
OxygenPreburnerMethaneReactant = Map.from_hdf5(
    "Oxygen Preburner Hot Methane Reactant Enthalpy",
    Raptor,
    filename=FILE_NAME,
    group="hot_methane_reactant_cea_t",
    inputs={
        "temperature": OxygenPreburnerFuelValveOutletFluid.temperature,
    },
    extrapolate=True,
)
FuelPreburnerLOXReactant = Map.from_hdf5(
    "Fuel Preburner LOX Reactant Enthalpy",
    Raptor,
    filename=FILE_NAME,
    group="lox_reactant_pt",
    inputs={
        "pressure": fuel_preburner_oxidizer_valve_outlet_pressure,
        "temperature": FuelPreburnerOxidizerValveOutletFluid.temperature,
    },
    extrapolate=True,
)
OxygenPreburnerLOXReactant = Map.from_hdf5(
    "Oxygen Preburner LOX Reactant Enthalpy",
    Raptor,
    filename=FILE_NAME,
    group="lox_reactant_pt",
    inputs={
        "pressure": main_oxidizer_valve_outlet_pressure,
        "temperature": MainOxidizerValveOutletFluid.temperature,
    },
    extrapolate=True,
)


# -----------------------------------------------------------------------------
# Preburner states and injector components
# -----------------------------------------------------------------------------

fuel_preburner_pressure = State(FUEL_PREBURNER_PRESSURE_GUESS)
oxygen_preburner_pressure = State(OXYGEN_PREBURNER_PRESSURE_GUESS)

FuelPreburnerMethaneInjector = DischargeCoefficient(
    "Fuel Preburner Methane Injector",
    Raptor,
    upstream_pressure=regen_outlet_pressure,
    downstream_pressure=fuel_preburner_pressure,
    density=RegenOutletFluid.density,
    discharge_coefficient=LIQUID_INJECTOR_CD,
    cross_sectional_area=FUEL_PREBURNER_METHANE_AREA,
)
FuelPreburnerLOXInjector = DischargeCoefficient(
    "Fuel Preburner LOX Injector",
    Raptor,
    upstream_pressure=fuel_preburner_oxidizer_valve_outlet_pressure,
    downstream_pressure=fuel_preburner_pressure,
    density=FuelPreburnerOxidizerValveOutletFluid.density,
    discharge_coefficient=LIQUID_INJECTOR_CD,
    cross_sectional_area=FUEL_PREBURNER_LOX_AREA,
)
OxygenPreburnerMethaneInjector = DischargeCoefficient(
    "Oxygen Preburner Methane Injector",
    Raptor,
    upstream_pressure=oxygen_preburner_fuel_valve_outlet_pressure,
    downstream_pressure=oxygen_preburner_pressure,
    density=OxygenPreburnerFuelValveOutletFluid.density,
    discharge_coefficient=LIQUID_INJECTOR_CD,
    cross_sectional_area=OXYGEN_PREBURNER_METHANE_AREA,
)
OxygenPreburnerLOXInjector = DischargeCoefficient(
    "Oxygen Preburner LOX Injector",
    Raptor,
    upstream_pressure=main_oxidizer_valve_outlet_pressure,
    downstream_pressure=oxygen_preburner_pressure,
    density=MainOxidizerValveOutletFluid.density,
    discharge_coefficient=LIQUID_INJECTOR_CD,
    cross_sectional_area=OXYGEN_PREBURNER_LOX_AREA,
)

fpb_methane_flow = FuelPreburnerMethaneInjector.mass_flow
fpb_lox_flow = FuelPreburnerLOXInjector.mass_flow
opb_methane_flow = OxygenPreburnerMethaneInjector.mass_flow
opb_lox_flow = OxygenPreburnerLOXInjector.mass_flow

# The chamber and nozzle returns mix adiabatically at one solved pressure.
# Their mixed enthalpy is determined directly by the two inlet streams.
RegenMergeHeader = Volume(
    "Hot Methane Merge Header",
    Raptor,
    pressure=regen_outlet_pressure,
    mass_flow_in=regen_merge_mass_flow,
    mass_flow_out=fpb_methane_flow,
)

RegenSupplyManifold = Volume(
    "Throat Regen Supply Manifold",
    Raptor,
    pressure=main_fuel_valve_outlet_pressure,
    mass_flow_in=MainFuelValve.mass_flow,
    mass_flow_out=regen_supply_mass_flow,
)
OxygenPreburnerFuelValveOutletNode = Volume(
    "Oxygen Preburner Fuel Valve Outlet Node",
    Raptor,
    pressure=oxygen_preburner_fuel_valve_outlet_pressure,
    mass_flow_in=OxygenPreburnerFuelValve.mass_flow,
    mass_flow_out=opb_methane_flow,
)
FuelPreburnerOxidizerValveOutletNode = Volume(
    "Fuel Preburner Oxidizer Valve Outlet Node",
    Raptor,
    pressure=fuel_preburner_oxidizer_valve_outlet_pressure,
    mass_flow_in=FuelPreburnerOxidizerValve.mass_flow,
    mass_flow_out=fpb_lox_flow,
)
MainOxidizerValveOutletNode = Volume(
    "Main Oxidizer Valve Outlet Node",
    Raptor,
    pressure=main_oxidizer_valve_outlet_pressure,
    mass_flow_in=MainOxidizerValve.mass_flow,
    mass_flow_out=opb_lox_flow,
)


# -----------------------------------------------------------------------------
# Pump discharge manifolds
# -----------------------------------------------------------------------------

# These adiabatic manifolds close only mass continuity at their solved
# pressures.  Pump discharge enthalpy is already known from each pump component,
# so adding separate manifold enthalpy variables would duplicate the exact
# steady energy relation.
MethaneDischargeManifold = Volume(
    "Methane Discharge Manifold",
    Raptor,
    pressure=methane_discharge_pressure,
    mass_flow_in=methane_pump_mass_flow,
    mass_flow_out=MainFuelValve.mass_flow,
)
LOXDischargeManifold = Volume(
    "LOX Discharge Manifold",
    Raptor,
    pressure=lox_discharge_pressure,
    mass_flow_in=lox_pump_mass_flow,
    mass_flow_out=(
        FuelPreburnerOxidizerValve.mass_flow + MainOxidizerValve.mass_flow
    ),
)


# -----------------------------------------------------------------------------
# Preburner equilibrium and continuity nodes
# -----------------------------------------------------------------------------

fuel_preburner_inlet_flow = fpb_methane_flow + fpb_lox_flow
oxygen_preburner_inlet_flow = opb_methane_flow + opb_lox_flow

# Mixture ratios are direct derived states from the injector component outputs.
# They are not independent solver variables or operating-point targets.
fuel_preburner_mixture_ratio = fpb_lox_flow / fpb_methane_flow
oxygen_preburner_mixture_ratio = opb_lox_flow / opb_methane_flow

fuel_preburner_inlet_enthalpy = (
    fpb_methane_flow * HotMethaneReactant.enthalpy
    + fpb_lox_flow * FuelPreburnerLOXReactant.enthalpy
) / fuel_preburner_inlet_flow
oxygen_preburner_inlet_enthalpy = (
    opb_methane_flow * OxygenPreburnerMethaneReactant.enthalpy
    + opb_lox_flow * OxygenPreburnerLOXReactant.enthalpy
) / oxygen_preburner_inlet_flow

# Each preburner is adiabatic at steady state, so its equilibrium enthalpy is
# the mass-weighted reactant enthalpy. The pressure remains a network unknown.
FuelPreburnerMap = Map.from_hdf5(
    "Fuel-Rich Preburner Equilibrium",
    Raptor,
    filename=FILE_NAME,
    group="fuel_preburner_equilibrium_ph",
    inputs={
        "pressure": fuel_preburner_pressure,
        "mixture_ratio": fuel_preburner_mixture_ratio,
        "target_enthalpy": fuel_preburner_inlet_enthalpy,
    },
    extrapolate=True,
)
OxygenPreburnerMap = Map.from_hdf5(
    "Oxygen-Rich Preburner Equilibrium",
    Raptor,
    filename=FILE_NAME,
    group="oxygen_preburner_equilibrium_ph",
    inputs={
        "pressure": oxygen_preburner_pressure,
        "mixture_ratio": oxygen_preburner_mixture_ratio,
        "target_enthalpy": oxygen_preburner_inlet_enthalpy,
    },
    extrapolate=True,
)

FuelPreburner = Volume(
    "Fuel-Rich Preburner",
    Raptor,
    pressure=fuel_preburner_pressure,
    mass_flow_in=fuel_preburner_inlet_flow,
)
OxygenPreburner = Volume(
    "Oxygen-Rich Preburner",
    Raptor,
    pressure=oxygen_preburner_pressure,
    mass_flow_in=oxygen_preburner_inlet_flow,
)


# -----------------------------------------------------------------------------
# Turbines and shaft torque balances
# -----------------------------------------------------------------------------

fuel_turbine_exit_pressure = State(FUEL_TURBINE_EXIT_PRESSURE_GUESS)
oxygen_turbine_exit_pressure = State(OXYGEN_TURBINE_EXIT_PRESSURE_GUESS)

FuelTurbineMap = Map.from_hdf5(
    "Fuel Turbine Map",
    Raptor,
    filename=FILE_NAME,
    group="fuel_turbine",
    inputs={
        "pressure_ratio": fuel_preburner_pressure / fuel_turbine_exit_pressure,
        "speed_ratio": methane_rotor_speed / METHANE_PUMP_MAP_SPEED,
    },
    extrapolate=True,
)
OxygenTurbineMap = Map.from_hdf5(
    "Oxygen Turbine Map",
    Raptor,
    filename=FILE_NAME,
    group="oxygen_turbine",
    inputs={
        "pressure_ratio": (
            oxygen_preburner_pressure / oxygen_turbine_exit_pressure
        ),
        "speed_ratio": lox_rotor_speed / LOX_PUMP_MAP_SPEED,
    },
    extrapolate=True,
)

# GasTurbine calculates efficiency from the isentropic enthalpy drop supplied
# here.  The ideal maps are native ThermoProp SP-equilibrium maps generated
# from the same CH4/O2 mixture-ratio basis as the preburner PH maps.  The inlet
# enthalpy, ideal outlet enthalpy, and actual outlet enthalpy therefore share one
# CEA thermodynamic reference basis.  No extra temperature states or
# user-written turbine energy balances are needed.
FuelTurbineIdealExhaust = Map.from_hdf5(
    "Fuel Turbine Ideal SP State",
    Raptor,
    filename=FILE_NAME,
    group="fuel_turbine_ideal_sp",
    inputs={
        "pressure": fuel_turbine_exit_pressure,
        "mixture_ratio": fuel_preburner_mixture_ratio,
        "inlet_entropy": FuelPreburnerMap.entropy,
    },
    extrapolate=True,
)
OxygenTurbineIdealExhaust = Map.from_hdf5(
    "Oxygen Turbine Ideal SP State",
    Raptor,
    filename=FILE_NAME,
    group="oxygen_turbine_ideal_sp",
    inputs={
        "pressure": oxygen_turbine_exit_pressure,
        "mixture_ratio": oxygen_preburner_mixture_ratio,
        "inlet_entropy": OxygenPreburnerMap.entropy,
    },
    extrapolate=True,
)

fuel_turbine_ideal_enthalpy_change = (
    fuel_preburner_inlet_enthalpy - FuelTurbineIdealExhaust.enthalpy
)
oxygen_turbine_ideal_enthalpy_change = (
    oxygen_preburner_inlet_enthalpy - OxygenTurbineIdealExhaust.enthalpy
)

FuelTurbine = GasTurbine(
    "Fuel-Rich Turbine",
    Raptor,
    rotor_speed=methane_rotor_speed,
    torque=FuelTurbineMap.torque,
    flow_parameter=FuelTurbineMap.flow_parameter,
    upstream_total_pressure=fuel_preburner_pressure,
    upstream_total_temperature=FuelPreburnerMap.temperature,
    downstream_pressure=fuel_turbine_exit_pressure,
    gas_constant=FuelPreburnerMap.gas_constant,
    specific_heat_ratio=FuelPreburnerMap.gamma,
    upstream_total_enthalpy=fuel_preburner_inlet_enthalpy,
    ideal_total_enthalpy_change=fuel_turbine_ideal_enthalpy_change,
    mass_flow=FuelPreburner.mass_flow_out,
)
OxygenTurbine = GasTurbine(
    "Oxygen-Rich Turbine",
    Raptor,
    rotor_speed=lox_rotor_speed,
    torque=OxygenTurbineMap.torque,
    flow_parameter=OxygenTurbineMap.flow_parameter,
    upstream_total_pressure=oxygen_preburner_pressure,
    upstream_total_temperature=OxygenPreburnerMap.temperature,
    downstream_pressure=oxygen_turbine_exit_pressure,
    gas_constant=OxygenPreburnerMap.gas_constant,
    specific_heat_ratio=OxygenPreburnerMap.gamma,
    upstream_total_enthalpy=oxygen_preburner_inlet_enthalpy,
    ideal_total_enthalpy_change=oxygen_turbine_ideal_enthalpy_change,
    mass_flow=OxygenPreburner.mass_flow_out,
)


# The turbine and pump maps act directly on the same shaft. Shaft speed is a
# solved network variable, and each necessary balance simply enforces equal
# turbine and pump torque. No rotor inertia is needed in a steady-state model.
MethaneShaftBalance = Balance(
    "Methane Turbopump Shaft Torque Balance",
    Raptor,
    variable=methane_rotor_speed,
    function=FuelTurbineMap.torque / MethanePumpMap.torque - 1.0,
)
LOXShaftBalance = Balance(
    "LOX Turbopump Shaft Torque Balance",
    Raptor,
    variable=lox_rotor_speed,
    function=OxygenTurbineMap.torque / LOXPumpMap.torque - 1.0,
)


# -----------------------------------------------------------------------------
# Turbine exhaust properties and main injectors
# -----------------------------------------------------------------------------

# Each turbine discharge enthalpy follows directly from the shaft work
# calculated by GasTurbine. The exhaust-plenum pressures remain solved nodes.
FuelTurbineExhaust = Map.from_hdf5(
    "Fuel Turbine Exhaust Properties",
    Raptor,
    filename=FILE_NAME,
    group="fuel_turbine_exhaust_ph",
    inputs={
        "pressure": fuel_turbine_exit_pressure,
        "mixture_ratio": fuel_preburner_mixture_ratio,
        "target_enthalpy": FuelTurbine.discharge_total_enthalpy,
    },
    extrapolate=True,
)
OxygenTurbineExhaust = Map.from_hdf5(
    "Oxygen Turbine Exhaust Properties",
    Raptor,
    filename=FILE_NAME,
    group="oxygen_turbine_exhaust_ph",
    inputs={
        "pressure": oxygen_turbine_exit_pressure,
        "mixture_ratio": oxygen_preburner_mixture_ratio,
        "target_enthalpy": OxygenTurbine.discharge_total_enthalpy,
    },
    extrapolate=True,
)

chamber_pressure = State(CHAMBER_PRESSURE_GUESS)

MainFuelInjector = CompressibleOrifice(
    "Main Fuel Gas Injector",
    Raptor,
    upstream_total_pressure=fuel_turbine_exit_pressure,
    upstream_total_temperature=FuelTurbineExhaust.temperature,
    downstream_pressure=chamber_pressure,
    discharge_coefficient=GAS_INJECTOR_CD,
    cross_sectional_area=MAIN_FUEL_INJECTOR_AREA,
    gas_constant=FuelTurbineExhaust.gas_constant,
    specific_heat_ratio=FuelTurbineExhaust.gamma,
    upstream_static_enthalpy=FuelTurbine.discharge_total_enthalpy,
    upstream_static_temperature=FuelTurbineExhaust.temperature,
)
MainOxygenInjector = CompressibleOrifice(
    "Main Oxygen Gas Injector",
    Raptor,
    upstream_total_pressure=oxygen_turbine_exit_pressure,
    upstream_total_temperature=OxygenTurbineExhaust.temperature,
    downstream_pressure=chamber_pressure,
    discharge_coefficient=GAS_INJECTOR_CD,
    cross_sectional_area=MAIN_OXYGEN_INJECTOR_AREA,
    gas_constant=OxygenTurbineExhaust.gas_constant,
    specific_heat_ratio=OxygenTurbineExhaust.gamma,
    upstream_static_enthalpy=OxygenTurbine.discharge_total_enthalpy,
    upstream_static_temperature=OxygenTurbineExhaust.temperature,
)

fuel_main_mass_flow = MainFuelInjector.mass_flow
oxygen_main_mass_flow = MainOxygenInjector.mass_flow

FuelTurbineExhaustPlenum = Volume(
    "Fuel Turbine Exhaust Plenum",
    Raptor,
    pressure=fuel_turbine_exit_pressure,
    mass_flow_in=FuelTurbine.mass_flow,
    mass_flow_out=fuel_main_mass_flow,
)
OxygenTurbineExhaustPlenum = Volume(
    "Oxygen Turbine Exhaust Plenum",
    Raptor,
    pressure=oxygen_turbine_exit_pressure,
    mass_flow_in=OxygenTurbine.mass_flow,
    mass_flow_out=oxygen_main_mass_flow,
)


# -----------------------------------------------------------------------------
# Main chamber
# -----------------------------------------------------------------------------

main_inlet_flow = fuel_main_mass_flow + oxygen_main_mass_flow
main_inlet_enthalpy = (
    fuel_main_mass_flow * MainFuelInjector.total_enthalpy
    + oxygen_main_mass_flow * MainOxygenInjector.total_enthalpy
) / main_inlet_flow

# The chamber map uses oxygen-rich-stream flow divided by fuel-rich-stream
# flow. Both flows are calculated by the main injector components.
main_stream_mixture_ratio = oxygen_main_mass_flow / fuel_main_mass_flow

# The two turbine exhaust streams mix adiabatically in the main chamber. The
# chamber pressure is solved from injector inflow and choked nozzle outflow.
MainChamberMap = Map.from_hdf5(
    "Main Chamber Equilibrium",
    Raptor,
    filename=FILE_NAME,
    group="main_chamber_ph",
    inputs={
        "pressure": chamber_pressure,
        "stream_mixture_ratio": main_stream_mixture_ratio,
        "target_enthalpy": main_inlet_enthalpy,
    },
    extrapolate=True,
)

MainChamber = Volume(
    "Main Chamber",
    Raptor,
    pressure=chamber_pressure,
    mass_flow_in=main_inlet_flow,
)

# -----------------------------------------------------------------------------
# Equilibrium nozzle
# -----------------------------------------------------------------------------

converging_pressure = State(CONVERGING_PRESSURE_GUESS)
throat_pressure = State(THROAT_PRESSURE_GUESS)
upper_nozzle_pressure = State(UPPER_NOZZLE_PRESSURE_GUESS)
exit_nozzle_thermal_pressure = State(EXIT_NOZZLE_THERMAL_PRESSURE_GUESS)
exit_pressure = State(EXIT_PRESSURE_GUESS)

ConvergingMap = Map.from_hdf5(
    "Converging Equilibrium Nozzle State",
    Raptor,
    filename=FILE_NAME,
    group="nozzle_sp",
    inputs={
        "pressure": converging_pressure,
        "chamber_entropy": MainChamberMap.entropy,
        "stream_mixture_ratio": main_stream_mixture_ratio,
    },
    extrapolate=True,
)
ThroatMap = Map.from_hdf5(
    "Throat Equilibrium Nozzle State",
    Raptor,
    filename=FILE_NAME,
    group="nozzle_sp",
    inputs={
        "pressure": throat_pressure,
        "chamber_entropy": MainChamberMap.entropy,
        "stream_mixture_ratio": main_stream_mixture_ratio,
    },
    extrapolate=True,
)
UpperNozzleMap = Map.from_hdf5(
    "Upper Nozzle Equilibrium Nozzle State",
    Raptor,
    filename=FILE_NAME,
    group="nozzle_sp",
    inputs={
        "pressure": upper_nozzle_pressure,
        "chamber_entropy": MainChamberMap.entropy,
        "stream_mixture_ratio": main_stream_mixture_ratio,
    },
    extrapolate=True,
)
ExitNozzleThermalMap = Map.from_hdf5(
    "Exit Nozzle Thermal State",
    Raptor,
    filename=FILE_NAME,
    group="nozzle_sp",
    inputs={
        "pressure": exit_nozzle_thermal_pressure,
        "chamber_entropy": MainChamberMap.entropy,
        "stream_mixture_ratio": main_stream_mixture_ratio,
    },
    extrapolate=True,
)
ExitMap = Map.from_hdf5(
    "Exit Nozzle Equilibrium Nozzle State",
    Raptor,
    filename=FILE_NAME,
    group="nozzle_sp",
    inputs={
        "pressure": exit_pressure,
        "chamber_entropy": MainChamberMap.entropy,
        "stream_mixture_ratio": main_stream_mixture_ratio,
    },
    extrapolate=True,
)

converging_area = THROAT_AREA * CONVERGING_AREA_RATIO
upper_nozzle_area = THROAT_AREA * UPPER_NOZZLE_AREA_RATIO
exit_nozzle_thermal_area = THROAT_AREA * EXIT_NOZZLE_THERMAL_AREA_RATIO

# The throat component calculates the engine nozzle flow. Its output is wired
# directly to the chamber outlet; no nozzle mass flow is prescribed.
ThroatFlow = AdiabaticFlow(
    "Equilibrium Nozzle Throat Flow",
    Raptor,
    upstream_static_enthalpy=main_inlet_enthalpy,
    downstream_static_enthalpy=ThroatMap.enthalpy,
    downstream_density=ThroatMap.density,
    downstream_cross_sectional_area=THROAT_AREA,
    mass_flow=MainChamber.mass_flow_out,
)

throat_velocity = (
    ThroatFlow.mass_flow / (ThroatMap.density * THROAT_AREA)
)
throat_mach = abs(throat_velocity) / ThroatMap.speed_of_sound

# This is the only explicit nozzle balance. Throat pressure is otherwise free,
# so the balance does not compete with a Volume equation.
NozzleChokedBoundary = Balance(
    "Equilibrium Nozzle Choked Boundary",
    Raptor,
    variable=throat_pressure,
    function=throat_mach - 1.0,
)


# The remaining station components calculate the mass flow supported by each
# local area and SP state. A one-equation Volume varies only that station's
# pressure until its flow equals the throat flow.
ConvergingFlow = AdiabaticFlow(
    "Equilibrium Nozzle Converging Flow",
    Raptor,
    upstream_static_enthalpy=main_inlet_enthalpy,
    downstream_static_enthalpy=ConvergingMap.enthalpy,
    downstream_density=ConvergingMap.density,
    downstream_cross_sectional_area=converging_area,
)
ConvergingNode = Volume(
    "Equilibrium Nozzle Converging Node",
    Raptor,
    pressure=converging_pressure,
    mass_flow_in=ThroatFlow.mass_flow,
    mass_flow_out=ConvergingFlow.mass_flow,
)

UpperNozzleFlow = AdiabaticFlow(
    "Equilibrium Nozzle Upper Flow",
    Raptor,
    upstream_static_enthalpy=main_inlet_enthalpy,
    downstream_static_enthalpy=UpperNozzleMap.enthalpy,
    downstream_density=UpperNozzleMap.density,
    downstream_cross_sectional_area=upper_nozzle_area,
)
UpperNozzleNode = Volume(
    "Equilibrium Nozzle Upper Node",
    Raptor,
    pressure=upper_nozzle_pressure,
    mass_flow_in=ThroatFlow.mass_flow,
    mass_flow_out=UpperNozzleFlow.mass_flow,
)

# This station represents the average gas state over the final cooled nozzle
# section. The actual exit state below is still used for thrust and Isp.
ExitNozzleThermalFlow = AdiabaticFlow(
    "Equilibrium Nozzle Thermal Exit Flow",
    Raptor,
    upstream_static_enthalpy=main_inlet_enthalpy,
    downstream_static_enthalpy=ExitNozzleThermalMap.enthalpy,
    downstream_density=ExitNozzleThermalMap.density,
    downstream_cross_sectional_area=exit_nozzle_thermal_area,
)
ExitNozzleThermalNode = Volume(
    "Equilibrium Nozzle Thermal Exit Node",
    Raptor,
    pressure=exit_nozzle_thermal_pressure,
    mass_flow_in=ThroatFlow.mass_flow,
    mass_flow_out=ExitNozzleThermalFlow.mass_flow,
)

ExitFlow = AdiabaticFlow(
    "Equilibrium Nozzle Exit Flow",
    Raptor,
    upstream_static_enthalpy=main_inlet_enthalpy,
    downstream_static_enthalpy=ExitMap.enthalpy,
    downstream_density=ExitMap.density,
    downstream_cross_sectional_area=EXIT_AREA,
)
ExitNode = Volume(
    "Equilibrium Nozzle Exit Node",
    Raptor,
    pressure=exit_pressure,
    mass_flow_in=ThroatFlow.mass_flow,
    mass_flow_out=ExitFlow.mass_flow,
)

chamber_velocity = (
    ThroatFlow.mass_flow / (MainChamberMap.density * CHAMBER_AREA)
)
converging_velocity = (
    ConvergingFlow.mass_flow / (ConvergingMap.density * converging_area)
)
upper_nozzle_velocity = (
    UpperNozzleFlow.mass_flow
    / (UpperNozzleMap.density * upper_nozzle_area)
)
exit_nozzle_thermal_velocity = (
    ExitNozzleThermalFlow.mass_flow
    / (ExitNozzleThermalMap.density * exit_nozzle_thermal_area)
)
exit_velocity = ExitFlow.mass_flow / (ExitMap.density * EXIT_AREA)

chamber_mach = abs(chamber_velocity) / MainChamberMap.speed_of_sound
converging_mach = (
    abs(converging_velocity) / ConvergingMap.speed_of_sound
)
upper_nozzle_mach = (
    abs(upper_nozzle_velocity) / UpperNozzleMap.speed_of_sound
)
exit_nozzle_thermal_mach = (
    abs(exit_nozzle_thermal_velocity)
    / ExitNozzleThermalMap.speed_of_sound
)
exit_mach = abs(exit_velocity) / ExitMap.speed_of_sound

nozzle_stations = {
    "Chamber Barrel": {
        "pressure": chamber_pressure,
        "temperature": MainChamberMap.temperature,
        "density": MainChamberMap.density,
        "gas_constant": MainChamberMap.gas_constant,
        "dynamic_viscosity": MainChamberMap.dynamic_viscosity,
        "prandtl": MainChamberMap.prandtl,
        "velocity": chamber_velocity,
        "mach": chamber_mach,
    },
    "Converging": {
        "pressure": converging_pressure,
        "temperature": ConvergingMap.temperature,
        "density": ConvergingMap.density,
        "gas_constant": ConvergingMap.gas_constant,
        "dynamic_viscosity": ConvergingMap.dynamic_viscosity,
        "prandtl": ConvergingMap.prandtl,
        "velocity": converging_velocity,
        "mach": converging_mach,
    },
    "Throat": {
        "pressure": throat_pressure,
        "temperature": ThroatMap.temperature,
        "density": ThroatMap.density,
        "gas_constant": ThroatMap.gas_constant,
        "dynamic_viscosity": ThroatMap.dynamic_viscosity,
        "prandtl": ThroatMap.prandtl,
        "velocity": throat_velocity,
        "mach": throat_mach,
    },
    "Upper Nozzle": {
        "pressure": upper_nozzle_pressure,
        "temperature": UpperNozzleMap.temperature,
        "density": UpperNozzleMap.density,
        "gas_constant": UpperNozzleMap.gas_constant,
        "dynamic_viscosity": UpperNozzleMap.dynamic_viscosity,
        "prandtl": UpperNozzleMap.prandtl,
        "velocity": upper_nozzle_velocity,
        "mach": upper_nozzle_mach,
    },
    "Exit Nozzle": {
        "pressure": exit_nozzle_thermal_pressure,
        "temperature": ExitNozzleThermalMap.temperature,
        "density": ExitNozzleThermalMap.density,
        "gas_constant": ExitNozzleThermalMap.gas_constant,
        "dynamic_viscosity": ExitNozzleThermalMap.dynamic_viscosity,
        "prandtl": ExitNozzleThermalMap.prandtl,
        "velocity": exit_nozzle_thermal_velocity,
        "mach": exit_nozzle_thermal_mach,
    },
}


# -----------------------------------------------------------------------------
# Regenerative heat-transfer model
# -----------------------------------------------------------------------------

# Each hydraulic coolant node also closes a steady energy balance. Every hot-gas
# zone has two solved GRCop-42 Solid temperatures: the gas-side copper surface
# and the coolant-side copper surface. Gas convection, copper conduction, and
# coolant convection are separate FullFlow components. Heat entering the
# coolant node is exactly the heat leaving the cold copper node.
def cooling_pass(
    node,
    channel_count,
    channel_width,
    channel_height,
    flow_area,
    hydraulic_diameter,
    hot_wall_temperature_guess,
    cold_wall_temperature_guess,
    hot_area_fraction=1.0,
):
    """Assemble one hydraulic/thermal regenerative-cooling section.

    Parameters
    ----------
    node:
        Dictionary describing the section outlet pressure/enthalpy state, the
        local property map, friction model, FlowTube, pass name, and geometric
        section data.  These dictionaries are created explicitly above while
        following the physical coolant path.
    channel_count, channel_width, channel_height:
        Effective rectangular-channel geometry for this pass.
    flow_area, hydraulic_diameter:
        Total flow area and hydraulic diameter used by the friction and
        convection correlations.
    hot_wall_temperature_guess, cold_wall_temperature_guess:
        Initial guesses only.  FullFlow solves both wall temperatures from the
        two algebraic heat balances.
    hot_area_fraction:
        Fraction of the gas-side surface assigned to this coolant pass.  Each
        nozzle downflow/upflow pass receives half of the physical nozzle wall.

    Returns
    -------
    dict
        References to the created FullFlow components and derived quantities.
        The returned dictionary is used only for reporting, tracking, and
        assembling section totals; it does not replace any network equation.

    Notes
    -----
    The local heat path is

    ``hot gas -> gas convection -> hot copper -> conduction -> cold copper``
    ``-> coolant convection -> methane control volume``.

    The two ``Solid`` objects have no thermal mass, so in this steady model each
    one contributes a zero-net-heat algebraic temperature equation.
    """
    pass_name = node["pass_name"]
    section = node["section"]
    section_name = section["name"]
    Station = nozzle_stations[section_name]
    LocalCoolant = node["mean_fluid"]

    CoolantHTC = Gnielinski(
        f"{pass_name} {section_name} Coolant Gnielinski",
        Raptor,
        mass_flow=node["pipe"].mass_flow,
        hydraulic_diameter=hydraulic_diameter,
        friction_factor=node["friction"].friction_factor,
        fluid_conductivity=LocalCoolant.conductivity,
        fluid_specific_heat=LocalCoolant.specific_heat_cp,
        fluid_dynamic_viscosity=LocalCoolant.dynamic_viscosity,
        cross_sectional_area=flow_area,
    )

    RecoveryFactor = TemperatureRecoveryFactor(
        f"{pass_name} {section_name} Recovery Factor",
        Raptor,
        prandtl_number=Station["prandtl"],
        turbulent=True,
    )
    AdiabaticWall = AdiabaticWallTemperature(
        f"{pass_name} {section_name} Adiabatic Wall Temperature",
        Raptor,
        total_temperature=MainChamberMap.temperature,
        static_temperature=Station["temperature"],
        recovery_factor=RecoveryFactor.recovery_factor,
    )

    hot_wall_temperature = State(hot_wall_temperature_guess)
    cold_wall_temperature = State(cold_wall_temperature_guess)

    HotCopperMaterial = Map.from_hdf5(
        f"{pass_name} {section_name} Hot GRCop-42",
        Raptor,
        filename=FILE_NAME,
        group="grcop42_t",
        inputs={"temperature": hot_wall_temperature},
        extrapolate=True,
    )
    ColdCopperMaterial = Map.from_hdf5(
        f"{pass_name} {section_name} Cold GRCop-42",
        Raptor,
        filename=FILE_NAME,
        group="grcop42_t",
        inputs={"temperature": cold_wall_temperature},
        extrapolate=True,
    )
    local_flow_area = THROAT_AREA * section["thermal_area_ratio"]
    local_hot_gas_diameter = math.sqrt(
        4.0 * local_flow_area / math.pi
    )
    hot_area = hot_area_fraction * section["hot_gas_area"]

    # Channel ribs are treated as straight copper fins. Nozzle downflow and
    # upflow use alternating channels, so each pass receives half of the hot
    # circumference and half of the physical channel count.
    channel_pitch = (
        hot_area_fraction
        * math.pi
        * local_hot_gas_diameter
        / channel_count
    )
    fin_thickness = channel_pitch - channel_width
    if fin_thickness <= 0.0:
        raise ValueError(
            f"{pass_name} {section_name}: channel width exceeds channel pitch."
        )

    HotCopper = Solid(
        f"{pass_name} {section_name} Hot Copper Wall",
        Raptor,
        temperature=hot_wall_temperature,
    )
    ColdCopper = Solid(
        f"{pass_name} {section_name} Cold Copper Wall",
        Raptor,
        temperature=cold_wall_temperature,
    )

    # The Bartz mean-temperature correction now uses the solved gas-side wall
    # temperature through FullFlow's Eckert reference-temperature component.
    ReferenceTemperature = EckertReferenceTemperature(
        f"{pass_name} {section_name} Eckert Reference Temperature",
        Raptor,
        wall_temperature=HotCopper.temperature,
        static_temperature=Station["temperature"],
        adiabatic_wall_temperature=AdiabaticWall.adiabatic_wall_temperature,
    )
    reference_density = (
        Station["pressure"]
        / (Station["gas_constant"] * ReferenceTemperature.reference_temperature)
    )
    reference_viscosity = (
        Station["dynamic_viscosity"]
        * (
            ReferenceTemperature.reference_temperature
            / Station["temperature"]
        )
        ** 0.70
    )

    HotGasHTC = Bartz(
        f"{pass_name} {section_name} Bartz",
        Raptor,
        mass_flow=ThroatFlow.mass_flow,
        hydraulic_diameter=local_hot_gas_diameter,
        chamber_specific_heat_cp=MainChamberMap.specific_heat_cp,
        chamber_prandtl_number=MainChamberMap.prandtl,
        chamber_dynamic_viscosity=MainChamberMap.dynamic_viscosity,
        local_freestream_density=Station["density"],
        mean_temperature_density=reference_density,
        mean_temperature_dynamic_viscosity=reference_viscosity,
    )

    fin_parameter = (
        (
            2.0 * CoolantHTC.convection_coefficient
            / (
                ColdCopperMaterial.thermal_conductivity
                * fin_thickness
            )
        )
        ** 0.5
        * channel_height
    )
    fin_efficiency = fin_parameter.tanh() / fin_parameter
    effective_coolant_perimeter = channel_count * (
        channel_width
        + 2.0 * fin_efficiency * channel_height
    )
    effective_coolant_area = (
        effective_coolant_perimeter * section["length"]
    )

    GasToHotWall = Convection(
        f"{pass_name} {section_name} Gas to Hot Copper Convection",
        Raptor,
        surface_temperature=HotCopper.temperature,
        fluid_temperature=AdiabaticWall.adiabatic_wall_temperature,
        convective_area=hot_area,
        convection_coefficient=HotGasHTC.convection_coefficient,
    )

    copper_conductivity = (
        2.0
        * HotCopperMaterial.thermal_conductivity
        * ColdCopperMaterial.thermal_conductivity
        / (
            HotCopperMaterial.thermal_conductivity
            + ColdCopperMaterial.thermal_conductivity
        )
    )
    HotToColdCopper = Conduction(
        f"{pass_name} {section_name} Copper Wall Conduction",
        Raptor,
        temperature1=HotCopper.temperature,
        temperature2=ColdCopper.temperature,
        thermal_conductivity=copper_conductivity,
        length=GRCOP_LINER_THICKNESS,
        conductive_area=hot_area,
    )

    ColdWallToCoolant = Convection(
        f"{pass_name} {section_name} Cold Copper to Coolant Convection",
        Raptor,
        surface_temperature=ColdCopper.temperature,
        fluid_temperature=LocalCoolant.temperature,
        convective_area=effective_coolant_area,
        convection_coefficient=CoolantHTC.convection_coefficient,
    )

    # FullFlow sign conventions:
    #   Convection is positive into the surface when the fluid is hotter.
    #   Conduction is positive into temperature1 when temperature2 is hotter.
    # These connections make both copper nodes close their own heat balances.
    HotCopper.heat_rate = (
        GasToHotWall.heat_rate
        + HotToColdCopper.heat_rate
    )
    ColdCopper.heat_rate = (
        -HotToColdCopper.heat_rate
        + ColdWallToCoolant.heat_rate
    )
    node["volume"].heat_rate = -ColdWallToCoolant.heat_rate

    gas_to_wall_heat = GasToHotWall.heat_rate
    wall_conduction_heat = -HotToColdCopper.heat_rate
    wall_to_coolant_heat = -ColdWallToCoolant.heat_rate

    return {
        "pass_name": pass_name,
        "name": section_name,
        "node": node,
        "coolant_pressure": 0.5 * (
            node["pipe"].upstream_static_pressure
            + node["pressure"]
        ),
        "coolant_temperature": LocalCoolant.temperature,
        "hot_gas": Station,
        "heat_rate": wall_to_coolant_heat,
        "gas_to_wall_heat": gas_to_wall_heat,
        "wall_conduction_heat": wall_conduction_heat,
        "wall_to_coolant_heat": wall_to_coolant_heat,
        "hot_wall_energy_error": HotCopper.heat_rate,
        "cold_wall_energy_error": ColdCopper.heat_rate,
        "hot_wall_temperature": HotCopper.temperature,
        "cold_wall_temperature": ColdCopper.temperature,
        "hot_copper": HotCopper,
        "cold_copper": ColdCopper,
        "gas_convection": GasToHotWall,
        "wall_conduction": HotToColdCopper,
        "coolant_convection": ColdWallToCoolant,
        "fin_efficiency": fin_efficiency,
        "effective_coolant_area": effective_coolant_area,
        "coolant_reynolds_number": node["friction"].reynolds_number,
        "coolant_friction_factor": node["friction"].friction_factor,
    }


cooling_passes = []
cooling_by_name = {}

for node, hot_guess, cold_guess in zip(
    chamber_coolant_nodes,
    CHAMBER_HOT_WALL_TEMPERATURE_GUESSES,
    CHAMBER_COLD_WALL_TEMPERATURE_GUESSES,
):
    Pass = cooling_pass(
        node,
        CHAMBER_COOLING_CHANNEL_COUNT,
        CHAMBER_COOLING_CHANNEL_WIDTH,
        CHAMBER_COOLING_CHANNEL_HEIGHT,
        CHAMBER_COOLING_FLOW_AREA,
        CHAMBER_COOLING_HYDRAULIC_DIAMETER,
        hot_guess,
        cold_guess,
    )
    cooling_passes.append(Pass)

    section = node["section"]
    cooling_by_name[section["name"]] = {
        **Pass,
        "branch": "Chamber",
        "heat_flux": Pass["heat_rate"] / section["hot_gas_area"],
        "thermal_area_ratio": section["thermal_area_ratio"],
    }

nozzle_passes_by_name = {
    section["name"]: [] for section in NOZZLE_THERMAL_SECTIONS
}
for node, hot_guess, cold_guess in zip(
    nozzle_coolant_nodes,
    NOZZLE_HOT_WALL_TEMPERATURE_GUESSES,
    NOZZLE_COLD_WALL_TEMPERATURE_GUESSES,
):
    Pass = cooling_pass(
        node,
        NOZZLE_COOLING_PASS_CHANNEL_COUNT,
        NOZZLE_COOLING_CHANNEL_WIDTH,
        NOZZLE_COOLING_CHANNEL_HEIGHT,
        NOZZLE_COOLING_PASS_FLOW_AREA,
        NOZZLE_COOLING_HYDRAULIC_DIAMETER,
        hot_guess,
        cold_guess,
        hot_area_fraction=0.5,
    )
    cooling_passes.append(Pass)
    nozzle_passes_by_name[node["section"]["name"]].append(Pass)

for section in NOZZLE_THERMAL_SECTIONS:
    Downflow, Upflow = nozzle_passes_by_name[section["name"]]
    heat_rate = Downflow["heat_rate"] + Upflow["heat_rate"]

    cooling_by_name[section["name"]] = {
        "name": section["name"],
        "branch": "Nozzle round trip",
        "coolant_temperature": 0.5 * (
            Downflow["coolant_temperature"]
            + Upflow["coolant_temperature"]
        ),
        "hot_gas": Downflow["hot_gas"],
        "heat_rate": heat_rate,
        "heat_flux": heat_rate / section["hot_gas_area"],
        "hot_wall_temperature": 0.5 * (
            Downflow["hot_wall_temperature"]
            + Upflow["hot_wall_temperature"]
        ),
        "cold_wall_temperature": 0.5 * (
            Downflow["cold_wall_temperature"]
            + Upflow["cold_wall_temperature"]
        ),
        "fin_efficiency": 0.5 * (
            Downflow["fin_efficiency"]
            + Upflow["fin_efficiency"]
        ),
        "effective_coolant_area": (
            Downflow["effective_coolant_area"]
            + Upflow["effective_coolant_area"]
        ),
        "thermal_area_ratio": section["thermal_area_ratio"],
    }

cooling_sections = [
    cooling_by_name[section["name"]] for section in THERMAL_SECTIONS
]

chamber_regen_heat = State(0.0)
nozzle_regen_heat = State(0.0)
for section in cooling_sections:
    if section["branch"] == "Chamber":
        chamber_regen_heat = chamber_regen_heat + section["heat_rate"]
    else:
        nozzle_regen_heat = nozzle_regen_heat + section["heat_rate"]

total_regen_heat = chamber_regen_heat + nozzle_regen_heat

# These are diagnostics only. The actual fluid energy equations are the seven
# local coolant-node balances, and the wall energy equations are the fourteen
# GRCop-42 Solid balances above.
chamber_regen_energy_error = (
    chamber_return_mass_flow
    * (
        chamber_return_enthalpy
        - MainFuelValveOutletFluid.enthalpy
    )
    - chamber_regen_heat
)
nozzle_regen_energy_error = (
    nozzle_regen_pipes[-1].mass_flow
    * (
        nozzle_return_enthalpy
        - MainFuelValveOutletFluid.enthalpy
    )
    - nozzle_regen_heat
)


# -----------------------------------------------------------------------------
# Performance, saved outputs, and solve
# -----------------------------------------------------------------------------

# Momentum thrust plus the exit-pressure correction.  The nozzle mass flow and
# exit state are solved by the connected chamber/nozzle network; neither thrust
# nor specific impulse participates in closing the engine operating point.
thrust = (
    ExitFlow.mass_flow * exit_velocity
    + (exit_pressure - ambient_pressure) * EXIT_AREA
)
specific_impulse = thrust / (ExitFlow.mass_flow * G0)

# ``track`` stores selected engineering outputs in the steady-state HDF5 run
# group.  It does not add a solve variable or equation.  The component data are
# also exported automatically, while these tracks provide stable, descriptive
# names for post-processing and for the accompanying solution diagram.
Raptor.track("Methane Tank Prevalve Mass Flow [kg/s]", MethaneMainValve.mass_flow)
Raptor.track("Methane Tank Prevalve Outlet Pressure [Pa]", methane_pump_inlet_pressure)
Raptor.track("LOX Tank Prevalve Mass Flow [kg/s]", LOXMainValve.mass_flow)
Raptor.track("LOX Tank Prevalve Outlet Pressure [Pa]", lox_pump_inlet_pressure)
Raptor.track("Main Fuel Valve Mass Flow [kg/s]", MainFuelValve.mass_flow)
Raptor.track(
    "Main Fuel Valve Outlet Pressure [Pa]",
    main_fuel_valve_outlet_pressure,
)
Raptor.track(
    "Oxygen Preburner Fuel Valve Mass Flow [kg/s]",
    OxygenPreburnerFuelValve.mass_flow,
)
Raptor.track(
    "Oxygen Preburner Fuel Valve Outlet Pressure [Pa]",
    oxygen_preburner_fuel_valve_outlet_pressure,
)
Raptor.track(
    "Fuel Preburner Oxidizer Valve Mass Flow [kg/s]",
    FuelPreburnerOxidizerValve.mass_flow,
)
Raptor.track(
    "Fuel Preburner Oxidizer Valve Outlet Pressure [Pa]",
    fuel_preburner_oxidizer_valve_outlet_pressure,
)
Raptor.track(
    "Main Oxidizer Valve Mass Flow [kg/s]",
    MainOxidizerValve.mass_flow,
)
Raptor.track(
    "Main Oxidizer Valve Outlet Pressure [Pa]",
    main_oxidizer_valve_outlet_pressure,
)
Raptor.track("Methane Pump Mass Flow [kg/s]", methane_pump_mass_flow)
Raptor.track("Methane Pump Efficiency", MethanePump.efficiency)
Raptor.track("LOX Pump Mass Flow [kg/s]", lox_pump_mass_flow)
Raptor.track("LOX Pump Efficiency", LOXPump.efficiency)
Raptor.track("Methane Rotor Speed [rpm]", methane_rotor_speed)
Raptor.track("LOX Rotor Speed [rpm]", lox_rotor_speed)
Raptor.track("Methane Pump Discharge Pressure [Pa]", methane_discharge_pressure)
Raptor.track("LOX Pump Discharge Pressure [Pa]", lox_discharge_pressure)
Raptor.track("Regen Supply Pressure [Pa]", main_fuel_valve_outlet_pressure)
Raptor.track("Regen Supply Mass Flow [kg/s]", regen_supply_mass_flow)
Raptor.track("Chamber Regen Mass Flow [kg/s]", chamber_regen_mass_flow)
Raptor.track("Nozzle Regen Mass Flow Before OPFV [kg/s]", nozzle_regen_mass_flow)
Raptor.track(
    "Nozzle Return To Merge Mass Flow [kg/s]",
    nozzle_return_to_merge_mass_flow,
)
Raptor.track("Regen Merge Mass Flow [kg/s]", regen_merge_mass_flow)
Raptor.track("Nozzle Turnaround Pressure [Pa]", nozzle_turnaround_pressure)
Raptor.track("Hot Methane Merge Pressure [Pa]", regen_outlet_pressure)
Raptor.track(
    "Chamber Return Manifold Temperature [K]",
    ChamberReturnManifoldFluid.temperature,
)
Raptor.track(
    "Nozzle Return Tap Temperature [K]",
    NozzleReturnTapFluid.temperature,
)
Raptor.track("Hot Methane Merge Temperature [K]", RegenOutletFluid.temperature)
Raptor.track(
    "Chamber Return Manifold Enthalpy [J/kg]",
    chamber_return_enthalpy,
)
Raptor.track(
    "Nozzle Return Tap Enthalpy [J/kg]",
    nozzle_return_enthalpy,
)
Raptor.track("Hot Methane Merge Enthalpy [J/kg]", regen_mixed_inlet_enthalpy)
Raptor.track("Chamber Regen Heat Load [W]", chamber_regen_heat)
Raptor.track("Nozzle Regen Heat Load [W]", nozzle_regen_heat)
Raptor.track(
    "Chamber Regen Energy Error [W]",
    chamber_regen_energy_error,
)
Raptor.track(
    "Nozzle Regen Energy Error [W]",
    nozzle_regen_energy_error,
)
Raptor.track("Fuel Preburner Pressure [Pa]", fuel_preburner_pressure)
Raptor.track("Fuel Preburner Mixture Ratio", fuel_preburner_mixture_ratio)
Raptor.track("Fuel Preburner Temperature [K]", FuelPreburnerMap.temperature)
Raptor.track("Oxygen Preburner Pressure [Pa]", oxygen_preburner_pressure)
Raptor.track("Oxygen Preburner Mixture Ratio", oxygen_preburner_mixture_ratio)
Raptor.track("Oxygen Preburner Temperature [K]", OxygenPreburnerMap.temperature)
Raptor.track("Fuel Turbine Efficiency", FuelTurbine.efficiency)
Raptor.track("Fuel Turbine Actual Enthalpy Drop [J/kg]", FuelTurbine.shaft_power / FuelTurbine.mass_flow)
Raptor.track("Fuel Turbine Ideal Enthalpy Drop [J/kg]", fuel_turbine_ideal_enthalpy_change)
Raptor.track("Fuel Turbine Ideal Exit Temperature [K]", FuelTurbineIdealExhaust.temperature)
Raptor.track("Fuel Turbine Exit Pressure [Pa]", fuel_turbine_exit_pressure)
Raptor.track("Fuel Turbine Exit Temperature [K]", FuelTurbineExhaust.temperature)
Raptor.track("Oxygen Turbine Efficiency", OxygenTurbine.efficiency)
Raptor.track("Oxygen Turbine Actual Enthalpy Drop [J/kg]", OxygenTurbine.shaft_power / OxygenTurbine.mass_flow)
Raptor.track("Oxygen Turbine Ideal Enthalpy Drop [J/kg]", oxygen_turbine_ideal_enthalpy_change)
Raptor.track("Oxygen Turbine Ideal Exit Temperature [K]", OxygenTurbineIdealExhaust.temperature)
Raptor.track("Oxygen Turbine Exit Pressure [Pa]", oxygen_turbine_exit_pressure)
Raptor.track(
    "Oxygen Turbine Exit Temperature [K]",
    OxygenTurbineExhaust.temperature,
)
Raptor.track("Chamber Pressure [Pa]", chamber_pressure)
Raptor.track("Chamber Temperature [K]", MainChamberMap.temperature)
Raptor.track("Chamber Stream Mixture Ratio", main_stream_mixture_ratio)
Raptor.track("Nozzle Mass Flow [kg/s]", ThroatFlow.mass_flow)
Raptor.track("Throat Mach", throat_mach)
Raptor.track("Exit Pressure [Pa]", exit_pressure)
Raptor.track("Exit Mach", exit_mach)
Raptor.track("Regenerative Heat Load [W]", total_regen_heat)
Raptor.track("Thrust [N]", thrust)
Raptor.track("Specific Impulse [s]", specific_impulse)

for section in cooling_sections:
    name = section["name"]
    Raptor.track(f"Cooling Hot Gas {name} Pressure [Pa]", section["hot_gas"]["pressure"])
    Raptor.track(f"Cooling Hot Gas {name} Temperature [K]", section["hot_gas"]["temperature"])
    Raptor.track(f"Cooling Hot Gas {name} Mach", section["hot_gas"]["mach"])
    Raptor.track(f"{name} Heat Rate [W]", section["heat_rate"])
    Raptor.track(f"{name} Heat Flux [W/m^2]", section["heat_flux"])
    Raptor.track(
        f"{name} Hot Wall Temperature [K]",
        section["hot_wall_temperature"],
    )
    Raptor.track(
        f"{name} Cold Wall Temperature [K]",
        section["cold_wall_temperature"],
    )
    Raptor.track(
        f"{name} Coolant Temperature [K]",
        section["coolant_temperature"],
    )
    Raptor.track(f"{name} Channel Fin Efficiency", section["fin_efficiency"])
    Raptor.track(
        f"{name} Effective Coolant Area [m^2]",
        section["effective_coolant_area"],
    )

for cooling_pass_data in cooling_passes:
    pass_label = (
        f"{cooling_pass_data['pass_name']} "
        f"{cooling_pass_data['name']}"
    )
    Raptor.track(
        f"{pass_label} Heat Rate [W]",
        cooling_pass_data["heat_rate"],
    )
    Raptor.track(
        f"{pass_label} Coolant Temperature [K]",
        cooling_pass_data["coolant_temperature"],
    )
    Raptor.track(
        f"{pass_label} Coolant Outlet Pressure [Pa]",
        cooling_pass_data["node"]["pressure"],
    )
    Raptor.track(
        f"{pass_label} Coolant Mass Flow [kg/s]",
        cooling_pass_data["node"]["pipe"].mass_flow,
    )
    Raptor.track(
        f"{pass_label} Coolant Reynolds Number",
        cooling_pass_data["coolant_reynolds_number"],
    )
    Raptor.track(
        f"{pass_label} Coolant Friction Factor",
        cooling_pass_data["coolant_friction_factor"],
    )
    Raptor.track(
        f"{pass_label} Gas To Wall Heat Rate [W]",
        cooling_pass_data["gas_to_wall_heat"],
    )
    Raptor.track(
        f"{pass_label} Copper Conduction Heat Rate [W]",
        cooling_pass_data["wall_conduction_heat"],
    )
    Raptor.track(
        f"{pass_label} Wall To Coolant Heat Rate [W]",
        cooling_pass_data["wall_to_coolant_heat"],
    )
    Raptor.track(
        f"{pass_label} Hot Copper Energy Error [W]",
        cooling_pass_data["hot_wall_energy_error"],
    )
    Raptor.track(
        f"{pass_label} Cold Copper Energy Error [W]",
        cooling_pass_data["cold_wall_energy_error"],
    )

Raptor.track("Chamber Return Pipe Mass Flow [kg/s]", ChamberReturnPipe.mass_flow)
Raptor.track("Nozzle Return Pipe Mass Flow [kg/s]", NozzleReturnPipe.mass_flow)
Raptor.track("Chamber Return Pipe Friction Factor", ChamberReturnFriction.friction_factor)
Raptor.track("Nozzle Return Pipe Friction Factor", NozzleReturnFriction.friction_factor)

# Start from the guesses in raptor_data.py and solve the complete square
# nonlinear system.  The guesses are numerical aids only; the final values are
# determined by the component equations.  Two derived-state passes make the
# script robust to component ordering while retaining a conventional user-style
# network definition.
solve_start_time = time.perf_counter()
SteadyState(Raptor).solve(
    filename=FILE_NAME,
    verbose=True,
    jacobian_method="2-point",
    state_max_passes=2,
)
solve_runtime = time.perf_counter() - solve_start_time
model_runtime = time.perf_counter() - program_start_time

print("\nRaptor steady solution")
print(f"  Methane flow: {methane_pump_mass_flow.value:.3f} kg/s")
print(f"  LOX flow: {lox_pump_mass_flow.value:.3f} kg/s")
print(f"  Total flow: {(methane_pump_mass_flow.value + lox_pump_mass_flow.value):.3f} kg/s")
print(f"  Overall O/F: {(lox_pump_mass_flow.value / methane_pump_mass_flow.value):.4f}")
print(f"  Methane rotor: {methane_rotor_speed.value:.1f} rpm")
print(f"  LOX rotor: {lox_rotor_speed.value:.1f} rpm")
print(
    f"  Methane tank prevalve: {methane_tank_pressure.value / BAR:.3f} -> "
    f"{methane_pump_inlet_pressure.value / BAR:.3f} bar, "
    f"loss={(methane_tank_pressure.value - methane_pump_inlet_pressure.value) / BAR:.3f} bar"
)
print(
    f"  LOX tank prevalve: {lox_tank_pressure.value / BAR:.3f} -> "
    f"{lox_pump_inlet_pressure.value / BAR:.3f} bar, "
    f"loss={(lox_tank_pressure.value - lox_pump_inlet_pressure.value) / BAR:.3f} bar"
)
print(
    f"  Methane pump: discharge={methane_discharge_pressure.value / BAR:.1f} bar, "
    f"power={MethanePump.shaft_power.value / 1.0e6:.2f} MW, "
    f"efficiency={MethanePump.efficiency.value:.4f}"
)
print(
    f"  LOX pump: discharge={lox_discharge_pressure.value / BAR:.1f} bar, "
    f"power={LOXPump.shaft_power.value / 1.0e6:.2f} MW, "
    f"efficiency={LOXPump.efficiency.value:.4f}"
)
print(
    f"  Fuel preburner: {fuel_preburner_pressure.value / BAR:.1f} bar, "
    f"MR={fuel_preburner_mixture_ratio.value:.4f}, "
    f"T={FuelPreburnerMap.temperature.value:.1f} K"
)
print(
    f"  Oxygen preburner: {oxygen_preburner_pressure.value / BAR:.1f} bar, "
    f"MR={oxygen_preburner_mixture_ratio.value:.2f}, "
    f"T={OxygenPreburnerMap.temperature.value:.1f} K"
)
print(
    f"  Fuel turbine exit: "
    f"{fuel_turbine_exit_pressure.value / BAR:.1f} bar, "
    f"{FuelTurbineExhaust.temperature.value:.1f} K, "
    f"power={FuelTurbine.shaft_power.value / 1.0e6:.2f} MW, "
    f"efficiency={FuelTurbine.efficiency.value:.4f}, "
    f"dh={FuelTurbine.shaft_power.value / FuelTurbine.mass_flow.value / 1000.0:.1f} kJ/kg, "
    f"dh_s={fuel_turbine_ideal_enthalpy_change.value / 1000.0:.1f} kJ/kg, "
    f"ideal T={FuelTurbineIdealExhaust.temperature.value:.1f} K"
)
print(
    f"  Oxygen turbine exit: "
    f"{oxygen_turbine_exit_pressure.value / BAR:.1f} bar, "
    f"{OxygenTurbineExhaust.temperature.value:.1f} K, "
    f"power={OxygenTurbine.shaft_power.value / 1.0e6:.2f} MW, "
    f"efficiency={OxygenTurbine.efficiency.value:.4f}, "
    f"dh={OxygenTurbine.shaft_power.value / OxygenTurbine.mass_flow.value / 1000.0:.1f} kJ/kg, "
    f"dh_s={oxygen_turbine_ideal_enthalpy_change.value / 1000.0:.1f} kJ/kg, "
    f"ideal T={OxygenTurbineIdealExhaust.temperature.value:.1f} K"
)
print(
    f"  Chamber: {chamber_pressure.value / BAR:.2f} bar, "
    f"{MainChamberMap.temperature.value:.1f} K, "
    f"stream MR={main_stream_mixture_ratio.value:.3f}"
)
print(
    f"  Main fuel valve: {methane_discharge_pressure.value / BAR:.1f} -> "
    f"{main_fuel_valve_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(methane_discharge_pressure.value - main_fuel_valve_outlet_pressure.value) / BAR:.1f} bar"
)
print(
    f"  Regen supply/merge: {main_fuel_valve_outlet_pressure.value / BAR:.1f} -> "
    f"{regen_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(main_fuel_valve_outlet_pressure.value - regen_outlet_pressure.value) / BAR:.1f} bar"
)
print(
    f"  Regen supply split: chamber={chamber_regen_mass_flow.value:.3f} kg/s, "
    f"nozzle before OPFV={nozzle_regen_mass_flow.value:.3f} kg/s, "
    f"total={regen_supply_mass_flow.value:.3f} kg/s"
)
print(
    f"  Nozzle return split: OPFV={OxygenPreburnerFuelValve.mass_flow.value:.3f} kg/s, "
    f"to merge={nozzle_return_to_merge_mass_flow.value:.3f} kg/s, "
    f"merge total={regen_merge_mass_flow.value:.3f} kg/s"
)
print(
    f"  Nozzle turnaround pressure: {nozzle_turnaround_pressure.value / BAR:.1f} bar"
)
print(
    f"  Chamber return pipe: {chamber_coolant_nodes[-1]['pressure'].value / BAR:.1f} -> "
    f"{regen_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(chamber_coolant_nodes[-1]['pressure'].value - regen_outlet_pressure.value) / BAR:.3f} bar"
)
print(
    f"  Nozzle return pipe: {nozzle_coolant_nodes[-1]['pressure'].value / BAR:.1f} -> "
    f"{regen_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(nozzle_coolant_nodes[-1]['pressure'].value - regen_outlet_pressure.value) / BAR:.3f} bar"
)
print(
    f"  Oxygen-preburner fuel valve: {nozzle_coolant_nodes[-1]['pressure'].value / BAR:.1f} -> "
    f"{oxygen_preburner_fuel_valve_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(nozzle_coolant_nodes[-1]['pressure'].value - oxygen_preburner_fuel_valve_outlet_pressure.value) / BAR:.1f} bar, "
    f"inlet T={NozzleReturnTapFluid.temperature.value:.1f} K"
)
print(
    f"  Main oxidizer valve: {lox_discharge_pressure.value / BAR:.1f} -> "
    f"{main_oxidizer_valve_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(lox_discharge_pressure.value - main_oxidizer_valve_outlet_pressure.value) / BAR:.1f} bar"
)
print(
    f"  Fuel-preburner oxidizer valve: {lox_discharge_pressure.value / BAR:.1f} -> "
    f"{fuel_preburner_oxidizer_valve_outlet_pressure.value / BAR:.1f} bar, "
    f"loss={(lox_discharge_pressure.value - fuel_preburner_oxidizer_valve_outlet_pressure.value) / BAR:.1f} bar"
)
print(
    f"  Regen temperatures: inlet={MainFuelValveOutletFluid.temperature.value:.1f} K, "
    f"chamber return={ChamberReturnManifoldFluid.temperature.value:.1f} K, "
    f"nozzle return tap={NozzleReturnTapFluid.temperature.value:.1f} K, "
    f"merged={RegenOutletFluid.temperature.value:.1f} K"
)
print(
    f"  Hot-methane pressure budget: pump={methane_discharge_pressure.value / BAR:.1f} bar, "
    f"MFV outlet={main_fuel_valve_outlet_pressure.value / BAR:.1f} bar, "
    f"regen outlet={regen_outlet_pressure.value / BAR:.1f} bar, "
    f"preburner={fuel_preburner_pressure.value / BAR:.1f} bar"
)
print(
    f"  Regenerative heat load: total={total_regen_heat.value / 1.0e6:.2f} MW, "
    f"chamber={chamber_regen_heat.value / 1.0e6:.2f} MW, "
    f"nozzle={nozzle_regen_heat.value / 1.0e6:.2f} MW"
)
print(
    f"  Regen branch energy errors: "
    f"chamber={chamber_regen_energy_error.value:.3e} W, "
    f"nozzle={nozzle_regen_energy_error.value:.3e} W"
)
print(
    f"  Nozzle: mdot={ThroatFlow.mass_flow.value:.3f} kg/s, "
    f"Mt={throat_mach.value:.6f}, "
    f"Pe={exit_pressure.value / BAR:.3f} bar, "
    f"Me={exit_mach.value:.3f}"
)
print(f"  Thrust: {thrust.value / 1.0e6:.3f} MN")
print(f"  Specific impulse: {specific_impulse.value:.2f} s")

print("\nValve and injector geometry")
print(
    f"  Methane tank prevalve: Cd={MAIN_METHANE_VALVE_CD:.3f}, "
    f"area={MAIN_METHANE_VALVE_AREA:.6f} m^2"
)
print(
    f"  LOX tank prevalve: Cd={MAIN_LOX_VALVE_CD:.3f}, "
    f"area={MAIN_LOX_VALVE_AREA:.6f} m^2"
)
print(
    f"  Main fuel valve: Cd={MAIN_FUEL_VALVE_CD:.3f}, "
    f"area={MAIN_FUEL_VALVE_AREA:.6f} m^2"
)
print(
    f"  Oxygen-preburner fuel valve: Cd={OXYGEN_PREBURNER_FUEL_VALVE_CD:.3f}, "
    f"area={OXYGEN_PREBURNER_FUEL_VALVE_AREA:.7f} m^2"
)
print(
    f"  Fuel-preburner oxidizer valve: Cd={FUEL_PREBURNER_OXIDIZER_VALVE_CD:.3f}, "
    f"area={FUEL_PREBURNER_OXIDIZER_VALVE_AREA:.7f} m^2"
)
print(
    f"  Main oxidizer valve: Cd={MAIN_OXIDIZER_VALVE_CD:.3f}, "
    f"area={MAIN_OXIDIZER_VALVE_AREA:.7f} m^2"
)
print(
    f"  Liquid injector areas [m^2]: FPB CH4={FUEL_PREBURNER_METHANE_AREA:.7f}, "
    f"FPB LOX={FUEL_PREBURNER_LOX_AREA:.7f}, "
    f"OPB CH4={OXYGEN_PREBURNER_METHANE_AREA:.7f}, "
    f"OPB LOX={OXYGEN_PREBURNER_LOX_AREA:.7f}"
)
print(
    f"  Main gas injector areas [m^2]: fuel={MAIN_FUEL_INJECTOR_AREA:.6f}, "
    f"oxygen={MAIN_OXYGEN_INJECTOR_AREA:.6f}"
)

print("\nNozzle and cooling geometry")
print(
    f"  Chamber: D={CHAMBER_DIAMETER:.4f} m, "
    f"A={CHAMBER_AREA:.6f} m^2, Ac/At={CHAMBER_AREA / THROAT_AREA:.4f}"
)
print(
    f"  Throat: A={THROAT_AREA:.6f} m^2, "
    f"D={2.0 * math.sqrt(THROAT_AREA / math.pi):.4f} m"
)
print(
    f"  Exit: A={EXIT_AREA:.6f} m^2, "
    f"D={2.0 * math.sqrt(EXIT_AREA / math.pi):.4f} m, "
    f"Ae/At={EXPANSION_RATIO:.2f}"
)
print(
    f"  Intermediate station area ratios: converging={CONVERGING_AREA_RATIO:.3f}, "
    f"upper nozzle={UPPER_NOZZLE_AREA_RATIO:.3f}, "
    f"upper split={UPPER_NOZZLE_SPLIT_AREA_RATIO:.3f}, "
    f"exit thermal={EXIT_NOZZLE_THERMAL_AREA_RATIO:.3f}"
)
print(f"  Physical cooled length: {TOTAL_COOLED_LENGTH:.3f} m")
print(
    f"  Chamber branch: L={CHAMBER_REGEN_LENGTH:.3f} m, "
    f"N={CHAMBER_COOLING_CHANNEL_COUNT}, "
    f"width={CHAMBER_COOLING_CHANNEL_WIDTH * 1.0e3:.3f} mm, "
    f"height={CHAMBER_COOLING_CHANNEL_HEIGHT * 1.0e3:.2f} mm, "
    f"area={CHAMBER_COOLING_FLOW_AREA:.6f} m^2, "
    f"Dh={CHAMBER_COOLING_HYDRAULIC_DIAMETER * 1.0e3:.3f} mm"
)
print(
    f"  Nozzle branch: one-way L={NOZZLE_ONE_WAY_LENGTH:.3f} m, "
    f"round-trip L={NOZZLE_REGEN_LENGTH:.3f} m, "
    f"N_total={NOZZLE_COOLING_CHANNEL_COUNT}, "
    f"N_per_pass={NOZZLE_COOLING_PASS_CHANNEL_COUNT}, "
    f"width={NOZZLE_COOLING_CHANNEL_WIDTH * 1.0e3:.3f} mm, "
    f"height={NOZZLE_COOLING_CHANNEL_HEIGHT * 1.0e3:.2f} mm, "
    f"pass area={NOZZLE_COOLING_PASS_FLOW_AREA:.6f} m^2, "
    f"Dh={NOZZLE_COOLING_HYDRAULIC_DIAMETER * 1.0e3:.3f} mm"
)
print(
    f"  Return manifolds: D={REGEN_RETURN_PIPE_DIAMETER:.3f} m, "
    f"area={REGEN_RETURN_PIPE_FLOW_AREA:.6f} m^2, "
    f"chamber L={CHAMBER_RETURN_PIPE_LENGTH:.3f} m, "
    f"nozzle L={NOZZLE_RETURN_PIPE_LENGTH:.3f} m"
)
for section in THERMAL_SECTIONS:
    print(
        f"  {section['name']}: L={section['length']:.3f} m, "
        f"end A/At={section['area_ratio']:.3f}, "
        f"thermal A/At={section['thermal_area_ratio']:.3f}, "
        f"Din={2.0 * section['radius_in']:.4f} m, "
        f"Dout={2.0 * section['radius_out']:.4f} m"
    )

print("\nRegenerative coolant section results")
print(
    "  section                              "
    "Pin [bar]   Pout [bar]   mdot [kg/s]   Tout [K]      Re        f"
)
for cooling_pass_data in cooling_passes:
    node = cooling_pass_data["node"]
    label = f"{cooling_pass_data['pass_name']} {cooling_pass_data['name']}"
    print(
        f"  {label:<36} "
        f"{node['pipe'].upstream_static_pressure.value / BAR:>9.2f} "
        f"{node['pressure'].value / BAR:>11.2f} "
        f"{node['pipe'].mass_flow.value:>13.3f} "
        f"{node['fluid'].temperature.value:>10.1f} "
        f"{node['friction'].reynolds_number.value:>9.3e} "
        f"{node['friction'].friction_factor.value:>8.5f}"
    )

print("\nWall temperatures by hot-gas station")
print("  station                 A/At   hot wall [K]   cold wall [K]   coolant [K]   heat flux [MW/m^2]   fin eta")
for section in cooling_sections:
    print(
        f"  {section['name']:<22} "
        f"{section['thermal_area_ratio']:>5.2f} "
        f"{section['hot_wall_temperature'].value:>12.1f} "
        f"{section['cold_wall_temperature'].value:>15.1f} "
        f"{section['coolant_temperature'].value:>13.1f} "
        f"{section['heat_flux'].value / 1.0e6:>18.2f} "
        f"{section['fin_efficiency'].value:>9.3f}"
    )

print("\nCopper-wall thermal network results")
print(
    "  section                              "
    "qgas [MW]   qcond [MW]   qcool [MW]   hot residual [W]   cold residual [W]"
)
for cooling_pass_data in cooling_passes:
    label = f"{cooling_pass_data['pass_name']} {cooling_pass_data['name']}"
    print(
        f"  {label:<36} "
        f"{cooling_pass_data['gas_to_wall_heat'].value / 1.0e6:>9.3f} "
        f"{cooling_pass_data['wall_conduction_heat'].value / 1.0e6:>12.3f} "
        f"{cooling_pass_data['wall_to_coolant_heat'].value / 1.0e6:>11.3f} "
        f"{cooling_pass_data['hot_wall_energy_error'].value:>18.3e} "
        f"{cooling_pass_data['cold_wall_energy_error'].value:>19.3e}"
    )

print("\nNozzle round-trip pass details")
for cooling_pass_data in cooling_passes:
    if not cooling_pass_data["pass_name"].startswith("Nozzle"):
        continue
    print(
        f"  {cooling_pass_data['pass_name']:<16} "
        f"{cooling_pass_data['name']:<14} "
        f"coolant={cooling_pass_data['coolant_temperature'].value:.1f} K, "
        f"heat={cooling_pass_data['heat_rate'].value / 1.0e6:.2f} MW"
    )

print(f"\nSolve runtime: {solve_runtime:.3f} s")
print(f"Total model runtime: {model_runtime:.3f} s")
print(f"Results and maps: {FILE_NAME}.h5")
