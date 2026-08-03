# Raptor 2 Full-Flow Staged-Combustion Model

This repository is a detailed steady-state **FullFlow** example for a simplified Raptor-like full-flow staged-combustion (FFSC) methane/oxygen engine.

The model is intentionally written as ordinary user code. The Python file follows the physical engine flow path, creates recognizable components and nodes, and lets the connected network solve its own pressures, flows, shaft speeds, mixture ratios, heat loads, chamber state, and nozzle performance.

It is not an exact or proprietary reconstruction of a specific Raptor block revision. Publicly unavailable pump maps, turbine maps, valve areas, injector areas, and cooling-channel geometry are explicit model inputs that can be replaced with better data.

## What the model solves

Only the following are operating boundary conditions:

- methane-tank pressure and temperature;
- LOX-tank pressure and temperature;
- ambient pressure.

The network solves, among other quantities:

- total methane and LOX flow;
- both pump-inlet and pump-discharge pressures;
- both shaft speeds;
- all valve and injector flows;
- the chamber/nozzle regenerative split;
- every coolant-section pressure and enthalpy;
- fuel-rich and oxygen-rich preburner pressures and mixture ratios;
- turbine discharge pressures and temperatures;
- main-chamber pressure and equilibrium temperature;
- nozzle station pressures and Mach numbers;
- hot- and cold-wall temperatures;
- thrust and specific impulse.

No engine flow rate, pump discharge pressure, preburner pressure, chamber pressure, cooling split, turbine exit pressure, or shaft speed is prescribed.

## Topology represented

The model preserves the defining FFSC organization:

1. A methane pump and LOX pump operate on separate shafts.
2. Regeneratively heated methane feeds the fuel-rich preburner.
3. A hot nozzle-return methane tap feeds the oxygen-rich preburner through the OPFV.
4. Most LOX feeds the oxygen-rich preburner through the MOV branch.
5. A smaller LOX branch feeds the fuel-rich preburner through the FPOV.
6. The fuel-rich turbine drives the methane pump.
7. The oxygen-rich turbine drives the LOX pump.
8. Both complete turbine exhaust streams enter the main chamber.
9. The main chamber discharges through a solved choked equilibrium nozzle.

The regenerative circuit contains an explicit chamber branch, nozzle downflow, nozzle turnaround, nozzle upflow, OPFV tap, return lines, and hot-methane merge header.

## Repository contents

| File | Purpose |
|---|---|
| `raptor.py` | Builds and solves the complete steady FullFlow network. |
| `raptor_data.py` | Boundary conditions, component data, geometry, map reference values, and initial guesses. |
| `raptor_maps.py` | Generates the HDF5 pump, turbine, real-fluid, material, and equilibrium maps. |
| `raptor.h5` | Included pre-generated maps and the latest saved steady solution. |
| `diagram.py` | Detailed text flow diagram and high-precision saved solution record. |
| `requirements.txt` | Python dependencies used directly by this example. |

## Installation

Create and activate a virtual environment, then install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

## Run the engine model

The repository already includes the required maps:

```bash
python3 raptor.py
```

The script prints a detailed operating-point report and saves the latest result into `raptor.h5` under the FullFlow steady-state run group.

A representative saved solution is approximately:

```text
Methane flow:       147.646 kg/s
LOX flow:           532.610 kg/s
Chamber pressure:   299.483 bar
Thrust:             2.355 MN
Specific impulse:   352.946 s
```

## View the detailed solution diagram

Open `diagram.py` as a text file, or print it in the terminal:

```bash
python3 diagram.py
```

The diagram records the saved flow path with high-precision pressures, flows, shaft powers, preburner states, turbine states, nozzle stations, regenerative hydraulics, wall temperatures, component counts, and solver diagnostics.

## Regenerating the maps

Normal users should **not** need to regenerate the maps. The pressure-enthalpy equilibrium maps are expensive because every grid cell requires a temperature root solve wrapped around repeated chemical-equilibrium calculations.

To intentionally rebuild all maps:

1. Open `raptor_maps.py`.
2. Change:

   ```python
   REBUILD_MAPS = False
   ```

   to:

   ```python
   REBUILD_MAPS = True
   ```

3. Run:

   ```bash
   python3 raptor_maps.py
   ```

4. Set `REBUILD_MAPS` back to `False` afterward.

With `REBUILD_MAPS = False`, FullPlot resumes and reuses completed map cells rather than overwriting them.

## Modeling approach

### Component-oriented closure

Only three general-purpose `Balance` components are used:

- fuel-rich turbine torque equals methane-pump torque;
- oxygen-rich turbine torque equals LOX-pump torque;
- throat Mach number equals one.

Everything else is closed by ordinary components:

- `DischargeCoefficient` for liquid valves and injectors;
- `CompressibleOrifice` for the two main gas injectors;
- `ConstantDensityPump` for map-based pump pressure rise and power;
- `GasTurbine` for map-based turbine flow and shaft work;
- `FlowTube` for regenerative and return-line momentum relations;
- algebraic `Volume` nodes for continuity and coolant energy balances;
- algebraic `Solid` nodes for wall heat balances;
- `Bartz`, `Gnielinski`, `Churchill`, `Convection`, and `Conduction` for heat transfer and friction;
- equilibrium property maps for the preburners, turbines, chamber, and nozzle.

### Direct steady enthalpy propagation

The model does not create redundant enthalpy unknowns where steady adiabatic conservation determines the result exactly. Pump discharge manifolds, throttling valves, preburners, turbine exhaust plenums, the hot-methane merge, and the main chamber use direct steady enthalpy relations.

Coolant enthalpies remain solved because every cooling section receives a solved wall heat rate.

### No imposed efficiencies

Pump and turbine efficiencies are calculated outputs. The pump maps provide head and torque; the turbine maps provide torque and flow parameter. Shaft power and ideal/actual enthalpy changes determine the reported efficiencies.

### No hidden bounds or operating targets

The public model does not apply `State` bounds or target equations to force a desired operating point. The values in the `Solver starting guesses only` section of `raptor_data.py` are numerical initial conditions, not constraints.

## Important limitations

This is a cycle-level engineering model. Major deliberate simplifications include:

- steady, one-dimensional, lumped flow;
- representative public pump and turbine surfaces;
- equilibrium preburner and chamber chemistry;
- equilibrium isentropic nozzle stations;
- no finite-rate chemistry or combustion-stability model;
- no seal leakage, bearing flow, auxiliary shaft loads, or turbine blade cooling;
- no inducer/cavitation model;
- no detailed injector-element or atomization model;
- no distributed manifolds or secondary cooling passages;
- one-way cycle thermal coupling: heat is added to the wall/coolant network but is not subtracted from the equilibrium nozzle gas path.

## Editing the model

For most studies, change values in `raptor_data.py` first. The comments distinguish:

- boundary conditions;
- fixed component data;
- map reference states;
- solver starting guesses.

After changing geometry or map definitions, regenerate the affected maps before solving.

## License

No license has been selected in this repository. Add the license you want before publishing if you intend to grant reuse, modification, or redistribution rights.
