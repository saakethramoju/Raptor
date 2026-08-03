# Raptor 2 Full-Flow Staged-Combustion Model

This repository is a detailed steady-state **FullFlow** example for a simplified Raptor-like full-flow staged-combustion (FFSC) methane/oxygen engine.

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

## Editing the model

For most studies, change values in `raptor_data.py` first. The comments distinguish:

- boundary conditions;
- fixed component data;
- map reference states;
- solver starting guesses.

After changing geometry or map definitions, regenerate the affected maps before solving.

## Codebases:

- FullFlow: [github.com/saakethramoju/FullFlow](https://github.com/saakethramoju/FullFlow)
- ThermoProp: [github.com/saakethramoju/ThermoProp](https://github.com/saakethramoju/ThermoProp)
- FullPlot: [https://github.com/saakethramoju/FullPlot](https://github.com/saakethramoju/FullPlot)
