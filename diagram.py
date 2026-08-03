r"""
RAPTOR-LIKE FULL-FLOW STAGED-COMBUSTION ENGINE
DETAILED SAVED SOLUTION FLOW DIAGRAM

This module is a static engineering record for the converged solution stored in
``raptor.h5``.  It contains no engine model equations.  Open it as a text file,
or run ``python3 diagram.py`` to print the complete diagram in a terminal.

The model is a transparent, steady, one-dimensional Raptor-like FFSC example.
It is not presented as a proprietary reconstruction of a specific Raptor block.
Pressures and flows are solved from tank boundary conditions and component data;
no engine mass flow, pump discharge pressure, preburner pressure, chamber
pressure, turbine exit pressure, cooling split, or shaft speed is prescribed.

Pressure convention
-------------------
Cycle/plenum values are the pressures used by the lumped cycle components.
Nozzle-station pressures are static pressures from the equilibrium SP map.
All numerical values below were read from the included saved solution.


OVERALL FULL-FLOW CYCLE
=======================

                                METHANE SHAFT
                         +--------------------------+
                         |                          |
CH4 TANK -> PREVALVE -> CH4 PUMP -> DISCHARGE MANIFOLD -> MFV
                         ^                          |
                         |                          v
                         |                 REGEN SUPPLY MANIFOLD
                         |                   /              \
                         |                  /                \
                         |       CHAMBER COOLING        NOZZLE ROUND TRIP
                         |                  \                /
                         |                   \-- RETURN ----/
                         |                          |
                         |              +-----------+-----------+
                         |              |                       |
                         |              v                       v
                         |        HOT-METHANE MERGE            OPFV
                         |              |                       |
                         |              v                       v
                         +------ FUEL-RICH TURBINE       OXYGEN-RICH PREBURNER
                                        ^                       |
                                        |                       v
                                  FUEL-RICH PREBURNER    OXYGEN-RICH TURBINE
                                        ^                       |
                                        |                       |
LOX TANK -> PREVALVE -> LOX PUMP -> DISCHARGE MANIFOLD          |
                         ^             /          \              |
                         |            /            \             |
                         |          FPOV           MOV            |
                         |            |             |             |
                         |            v             v             |
                         |      FPB LOX INJ.   OPB LOX INJ.       |
                         +------------|-------------|-------------+
                                      |             |
                                      v             v
                              FUEL-RICH GAS   OXYGEN-RICH GAS
                                  INJECTOR        INJECTOR
                                      \             /
                                       \           /
                                        MAIN CHAMBER
                                             |
                                             v
                                      CHOKED NOZZLE


MODEL CLOSURE
=============

Operating boundary conditions
    CH4 tank pressure:        4.0000000000 bar
    CH4 tank temperature:     111.0000000000 K
    LOX tank pressure:        4.0000000000 bar
    LOX tank temperature:     90.0000000000 K
    Ambient pressure:         1.0132500000 bar

Fixed component data
    pump head/torque maps
    turbine torque/flow-parameter maps
    valve and injector CdA geometry
    regenerative-channel dimensions and roughness
    chamber, throat, and nozzle geometry
    equilibrium and real-fluid property maps
    GRCop-42 thermal-conductivity map

Only three explicit general-purpose Balance components
    1. fuel-rich turbine torque / methane-pump torque - 1 = 0
    2. oxygen-rich turbine torque / LOX-pump torque - 1 = 0
    3. throat Mach number - 1 = 0

All other closure comes from ordinary components
    valve and injector pressure-flow relations
    pump-map pressure rise
    turbine-map flow and shaft work
    FlowTube momentum equations
    algebraic node mass continuity
    coolant-node mass and energy balances
    algebraic hot/cold copper heat balances
    equilibrium chamber and nozzle property maps

No State lower or upper bounds are applied.
No pump, turbine, shaft, or mechanical efficiency is imposed as a multiplier.
Reported efficiencies are calculated outputs from flow, head, torque, speed,
and thermodynamic enthalpy changes.

Saved steady-solver diagnostics
    Iteration variables:      60
    Residual equations:       60
    Function evaluations:     4
    Jacobian evaluations:     4
    Maximum |residual|:       1.087784767151e-06
    Cost:                     7.826899382750e-13
    Optimality:               2.585504676944e-01
    Solver time in HDF5:      7.311095750 s
    Termination:              `xtol` termination condition is satisfied.


LATEST STEADY OPERATING POINT
=============================

    CH4 flow:                 147.645677270243 kg/s
    LOX flow:                 532.610057664678 kg/s
    Total flow:               680.255734934921 kg/s
    Overall O/F:              3.607352870141
    Chamber pressure:         299.482927169875 bar
    Chamber temperature:      3810.939108886642 K
    Chamber stream O/F:       3.444676510620
    Thrust:                   2354512.324748931 N
    Thrust:                   2.354512324749 MN
    Equilibrium Isp:          352.945850788387 s


METHANE FLOW PATH
=================

CH4 TANK
    Pressure:                 4.0000000000 bar
    Temperature:              111.00000000 K
    Flow:                     147.645677270243 kg/s
    |
    v
METHANE TANK PREVALVE
    Cd:                       0.9000000000
    Area:                     0.040000000000 m^2
    Pressure:                 4.0000000000 -> 3.8014637712 bar
    Pressure loss:            0.1985362288 bar
    Flow:                     147.645677270243 kg/s
    |
    v
METHANE PUMP  <==============================  FUEL-RICH TURBINE / SHAFT
    Inlet pressure:           3.8014637712 bar
    Discharge pressure:       898.9975859905 bar
    Flow:                     147.645677270243 kg/s
    Volumetric flow:          0.348556289114 m^3/s
    Inlet density:            423.592062118146 kg/m^3
    Head rise:                21550.120534775851 m
    Shaft torque:             9809.005255495127 N*m
    Shaft speed:              36165.703415597825 rpm
    Shaft power absorbed:     37.149288609428 MW
    Calculated efficiency:    0.839925204681
    Discharge enthalpy:       249696.810371248692 J/kg
    |
    v
METHANE DISCHARGE MANIFOLD
    Algebraic pressure/continuity node
    No storage volume and no separate enthalpy solve variable
    Enthalpy is propagated from the methane-pump discharge
    |
    v
MAIN FUEL VALVE — MFV
    Cd:                       0.8500000000
    Area:                     0.005180000000 m^2
    Pressure:                 898.9975859905 -> 886.4050647251 bar
    Pressure loss:            12.5925212654 bar
    Flow:                     147.645677270260 kg/s
    |
    v
REGEN SUPPLY MANIFOLD
    Pressure:                 886.4050647251 bar
    Flow:                     147.645677270243 kg/s
    Inlet temperature:        145.6000000000 K
    |
    +-- CHAMBER BRANCH:       67.582513994591 kg/s
    |       throat -> converging -> chamber barrel -> return
    |
    +-- NOZZLE ROUND TRIP:    80.063163275652 kg/s
            upper downflow -> exit downflow -> turnaround
            -> exit upflow -> upper upflow -> OPFV tap

NOZZLE RETURN / OPFV TAP
    Pressure:                 698.1417699578 bar
    Temperature:              610.9455631416 K
    Total branch flow:        80.063163275652 kg/s
    |
    +-- OPFV -> OPB CH4:      8.531227661169 kg/s
    |       outlet pressure:  673.7344068394 bar
    |
    +-- TO HOT-CH4 MERGE:     71.531935614484 kg/s

CHAMBER RETURN
    Pressure before line:     698.1352010629 bar
    Merge pressure:           698.1185283324 bar
    Pressure loss:            0.0166727305 bar
    Flow:                     67.582513994591 kg/s
    Temperature:              613.2660355511 K

HOT-METHANE MERGE HEADER
    Pressure:                 698.1185283324 bar
    Temperature:              612.0731915306 K
    Mixed enthalpy:           1734864.5735085758 J/kg
    Flow:                     139.114449609075 kg/s
    Algebraic pressure/continuity node; mixed enthalpy is calculated directly
    |
    v
FUEL-PREBURNER METHANE INJECTOR
    Flow:                     139.114449609075 kg/s
    Area:                     0.002827000000 m^2
    |
    v
FUEL-RICH PREBURNER
    Pressure:                 609.2455870672 bar
    Temperature:              876.4213366146 K
    CH4 flow:                 139.114449609075 kg/s
    LOX flow:                 13.935099286003 kg/s
    Total flow:               153.049548895077 kg/s
    O/F mixture ratio:        0.100170034998
    |
    v
FUEL-RICH TURBINE
    Inlet pressure:           609.2455870672 bar
    Exit pressure:            334.1927190347 bar
    Exit temperature:         818.8047396487 K
    Flow:                     153.049548895077 kg/s
    Torque produced:          9809.005255495127 N*m
    Shaft power produced:     37.149288609428 MW
    Actual enthalpy drop:     242.727200946511 kJ/kg
    Ideal enthalpy drop:      257.973177496341 kJ/kg
    Ideal exit temperature:   815.2885645435 K
    Calculated efficiency:    0.940900923508
    |
    v
MAIN FUEL-GAS INJECTOR
    Upstream pressure:        334.1927190347 bar
    Chamber pressure:         299.4829271699 bar
    Area:                     0.007736000000 m^2
    Flow:                     153.049548895079 kg/s
    Choked:                   False


OXYGEN FLOW PATH
================

LOX TANK
    Pressure:                 4.0000000000 bar
    Temperature:              90.00000000 K
    Flow:                     532.610057664678 kg/s
    |
    v
LOX TANK PREVALVE
    Cd:                       0.9000000000
    Area:                     0.120000000000 m^2
    Pressure:                 4.0000000000 -> 3.8935902029 bar
    Pressure loss:            0.1064097971 bar
    Flow:                     532.610057664678 kg/s
    |
    v
LOX PUMP  <==================================  OXYGEN-RICH TURBINE / SHAFT
    Inlet pressure:           3.8935902029 bar
    Discharge pressure:       701.1773192796 bar
    Flow:                     532.610057664678 kg/s
    Volumetric flow:          0.466078127176 m^3/s
    Inlet density:            1142.748450547677 kg/m^3
    Head rise:                6222.117590660235 m
    Shaft torque:             12587.718415418760 N*m
    Shaft speed:              30067.667927215443 rpm
    Shaft power absorbed:     39.634682396425 MW
    Calculated efficiency:    0.819960385472
    Discharge enthalpy:       -59117.597006354321 J/kg
    |
    v
LOX DISCHARGE MANIFOLD
    Algebraic pressure/continuity node
    No storage volume and no separate enthalpy solve variable
    Enthalpy is propagated from the LOX-pump discharge
    |
    +-- FPOV -> FPB LOX
    |       pressure:         701.1773192796 -> 689.4073987074 bar
    |       valve flow:       13.935099286003 kg/s
    |       injector area:    0.000112700000 m^2
    |
    +-- MOV -> OPB LOX
            pressure:         701.1773192796 -> 696.3096649840 bar
            valve flow:       518.674958378673 kg/s
            injector area:    0.004220000000 m^2

OXYGEN-RICH PREBURNER
    Pressure:                 617.2282903742 bar
    Temperature:              785.4330726028 K
    CH4 flow:                 8.531227661169 kg/s
    LOX flow:                 518.674958378673 kg/s
    Total flow:               527.206186039843 kg/s
    O/F mixture ratio:        60.797223914152
    |
    v
OXYGEN-RICH TURBINE
    Inlet pressure:           617.2282903742 bar
    Exit pressure:            385.7383450848 bar
    Exit temperature:         716.2492900069 K
    Flow:                     527.206186039843 kg/s
    Torque produced:          12587.718415418760 N*m
    Shaft power produced:     39.634682396425 MW
    Actual enthalpy drop:     75.178712704691 kJ/kg
    Ideal enthalpy drop:      97.301982914887 kJ/kg
    Ideal exit temperature:   699.5579936478 K
    Calculated efficiency:    0.772632894547
    |
    v
MAIN OXYGEN-GAS INJECTOR
    Upstream pressure:        385.7383450848 bar
    Chamber pressure:         299.4829271699 bar
    Area:                     0.011411000000 m^2
    Flow:                     527.206186039842 kg/s
    Choked:                   False


MAIN CHAMBER AND NOZZLE
=======================

FUEL-RICH TURBINE EXHAUST               OXYGEN-RICH TURBINE EXHAUST
    153.049548895077 kg/s                            527.206186039843 kg/s
    334.1927190347 bar, 818.80473965 K              385.7383450848 bar, 716.24929001 K
                   \                                  /
                    \                                /
                     +------------------------------+
                                    |
                                    v
                             MAIN CHAMBER
    Algebraic pressure/continuity node
    Pressure:                 299.482927169875 bar
    Temperature:              3810.939108886642 K
    Flow:                     680.255734934917 kg/s
    Stream O/F:               3.444676510620
    Overall engine O/F:       3.607352870141
    Mixed inlet enthalpy is passed directly to the HP-equilibrium map
                                    |
                                    v
                             EQUILIBRIUM NOZZLE
    Throat mass flow is not prescribed.
    Throat pressure varies until Mach = 1.

Nozzle station results
    station                      A/At      pressure [bar]    temperature [K]           Mach
Chamber barrel           1.64385662      299.482927170      3810.93910889   0.3618700564
Converging               1.25000000      249.633569286      3743.84034545   0.5646941046
Throat                   1.00000000      172.575405664      3613.89765824   1.0000000000
Upper nozzle             4.00000000       14.949517148      2872.05435333   2.5148255440
Exit thermal            18.74621125        1.878168397      2270.77802308   3.4641681443
Physical exit           34.00000000        0.956004680      2060.86522910   3.7803780018

Geometry
    Chamber diameter:        0.300000000000 m
    Chamber area:            0.070685834706 m^2
    Chamber A/At:            1.643856621064
    Throat area:             0.043000000000 m^2
    Throat diameter:         0.233985684228 m
    Exit area:               1.462000000000 m^2
    Exit diameter:           1.364359268816 m
    Exit A/At:               34.000000000000

Performance at ambient pressure
    Nozzle flow:             680.255734934917 kg/s
    Exit pressure:           0.956004679983 bar
    Exit Mach:               3.780378001771
    Ambient pressure:        1.013250000000 bar
    Thrust:                  2354512.324748931 N
    Specific impulse:        352.945850788387 s


REGENERATIVE COOLING NETWORK
============================

Physical routing
    MFV outlet -> parallel split
        chamber branch: throat -> converging -> chamber barrel -> return
        nozzle branch: upper downflow -> exit downflow -> turnaround
                       -> exit upflow -> upper upflow -> OPFV tap/return

    The OPFV methane branch leaves the nozzle-return tap before the remaining
    nozzle flow joins the chamber return at the hot-methane merge header.

Each cooling section contains
    one independent branch mass-flow State
    one outlet-pressure State
    one outlet-enthalpy State
    one real-fluid methane pressure-enthalpy Map
    one Churchill friction component
    one FlowTube momentum relation
    one algebraic coolant mass balance
    one algebraic coolant energy balance
    one Gnielinski coolant-side coefficient
    one Bartz gas-side coefficient
    one hot-copper Solid temperature balance
    one cold-copper Solid temperature balance
    one copper Conduction component
    two Convection components

FlowTube is retained because the cooling methane is dense/supercritical and its
pressure loss depends on local inlet/outlet density, enthalpy, velocity,
friction, area, and geometry.  In this steady solve the line length participates
in the momentum equation; no transient line state is advanced.

Cooling hydraulic and energy results
    section                         Pin [bar] -> Pout [bar]    mdot [kg/s]     Tout [K]       hout [J/kg]   rho [kg/m3]             Re             f       Q [MW]
Chamber throat                 886.405065 ->   869.809080    67.58251399   203.760705    426038.314670    397.353009   625909.00898   0.014786792   11.9176022
Chamber converging             869.809080 ->   831.785649    67.58251399   325.917921    792900.148658    306.957259  1128834.50024   0.014163776   24.7934450
Chamber barrel                 831.785649 ->   698.135201    67.58251399   613.266036   1739270.981586    167.917218  1932571.35371   0.013779703   63.9581201
Nozzle downflow upper          886.405065 ->   856.570414    80.06316328   283.177405    665259.658402    337.755904  1666387.66679   0.012482539   33.2712762
Nozzle downflow exit           856.570414 ->   811.396017    80.06316328   398.951770   1018971.633148    263.284930  2887002.29126   0.012046690   28.3192996
Nozzle upflow exit             811.396017 ->   755.180904    80.06316328   503.575541   1357219.696838    209.785798  3635824.65717   0.011907154   27.0812100
Nozzle upflow upper            755.180904 ->   698.141770    80.06316328   610.945563   1730701.452058    168.460570  3998746.73228   0.011855989   29.9021307

Individual coolant-pass wall results
    section                        hot wall [K]     cold wall [K]        Q [MW]
Chamber throat                  2236.96597282      988.16083487   11.91760218
Chamber converging              2161.67826164     1040.58456770   24.79344503
Chamber barrel                  2060.93663105     1100.29461841   63.95812006
Nozzle downflow upper           1290.37642486      748.73437651   33.27127615
Nozzle downflow exit             720.68299423      577.15968847   28.31929959
Nozzle upflow exit               804.65993685      666.01119522   27.08120995
Nozzle upflow upper             1474.38789640      976.69213494   29.90213075

Branch totals
    Chamber regenerative heat:     100.669167272110 MW
    Nozzle regenerative heat:      118.573916439356 MW
    Total regenerative heat:       219.243083711466 MW
    Chamber energy closure error:  -2.980232238770e-08 W
    Nozzle energy closure error:   3.278255462646e-07 W

Station-average thermal results
    station                  hot wall [K]     cold wall [K]   mean coolant [K]   heat flux [MW/m2]  fin efficiency
Chamber barrel            2060.93663105     1100.29461841       476.16223435      154.2310751286    0.1997582560
Converging                2161.67826164     1040.58456770       265.03789581      180.9328357056    0.1940361089
Throat                    2236.96597282      988.16083487       174.61362487      202.6561265218    0.1828763866
Upper nozzle              1382.38216063      862.71325573       386.11038240       85.4329430182    0.2485733417
Exit nozzle                762.67146554      621.58544185       396.98438544       24.4272365589    0.4594034839

Thermal-coupling qualification
    The gas-to-wall, copper-conduction, and wall-to-coolant rates close locally.
    The resulting heat is added to the methane circuit.  The same heat is not
    subtracted from the equilibrium chamber/nozzle gas path, so the present
    cycle model is one-way thermally coupled.


COOLING AND NOZZLE GEOMETRY
===========================

Chamber cooling branch
    Total cooled length:      0.680000000000 m
    Channel count:            306
    Channel width:            0.912000000000 mm
    Channel height:           6.250000000000 mm
    Total flow area:          0.001744200000 m^2
    Hydraulic diameter:       1.591734152471 mm
    Roughness:                0.250000000000 micrometres

Nozzle round-trip branch
    One-way cooled length:    1.100000000000 m
    Round-trip length:        2.200000000000 m
    Total physical channels:  306
    Channels per pass:        153
    Channel width:            2.477000000000 mm
    Channel height:           6.250000000000 mm
    Per-pass flow area:       0.002368631250 m^2
    Hydraulic diameter:       3.547897330125 mm
    Roughness:                0.250000000000 micrometres

Return manifolds
    Diameter:                 0.120000000000 m
    Flow area:                0.011309733553 m^2
    Roughness:                2.500000000000 micrometres
    Chamber return length:    0.200000000000 m
    Nozzle return length:     0.250000000000 m

Wall and total length
    GRCop-42 conduction thickness: 2.000000000000 mm
    Physical cooled length:        1.780000000000 m


CURRENT STEADY COMPONENT ORGANIZATION
=====================================

Component counts in the saved network
    Map                        53
    Volume                     25
    FlowTube                    9
    Churchill                   9
    Solid                      14
    Convection                 14
    Conduction                  7
    Gnielinski                  7
    Bartz                       7
    AdiabaticFlow               5
    DischargeCoefficient       10
    CompressibleOrifice         2
    ConstantDensityPump         2
    GasTurbine                  2
    Balance                     3

Interpretation of the node components
    Volume without a physical volume
        algebraic continuity node; pressure is the network unknown

    Coolant Volume without a physical volume
        algebraic mass/energy node; pressure and enthalpy are unknowns

    Solid without mass or specific heat
        algebraic zero-net-heat wall-temperature node

    FlowTube in a steady solve
        one-dimensional branch momentum and pressure-loss relation

There are no assumed plenum volumes, stored fluid masses, internal-energy
storage states, wall thermal masses, rotor inertias, schedules, or actuator
dynamics in this published steady model.


WHAT THE MODEL DOES AND DOES NOT CLAIM
======================================

The represented topology is a recognizable simplified full-flow
staged-combustion cycle:
    both propellants are fully gasified through separate preburner/turbine paths
    each turbine drives the pump on its corresponding shaft
    methane follows an explicit regenerative split, turnaround, tap, and merge
    oxygen-rich and fuel-rich turbine exhausts recombine in the main chamber
    nozzle flow is determined by the solved choked throat

The network is not claimed to be a proprietary reconstruction of a specific
Raptor block revision.  Valve dimensions, maps, channel geometry,
chamber/nozzle geometry, and property maps are model data and must be validated
against the intended hardware.

Deliberate cycle-level simplifications
    steady, one-dimensional, lumped flow
    equilibrium preburner and chamber chemistry
    equilibrium isentropic nozzle stations
    no finite-rate chemistry or combustion-stability model
    no leakage, seals, bearing loads, auxiliary shaft loads, or gear losses
    no explicit inducer/cavitation model
    no detailed injector-element pressure distribution or atomization model
    no distributed manifolds or secondary cooling passages
    no turbine blade-cooling model
    one-way gas-to-coolant thermal coupling
    representative public pump/turbine surfaces rather than proprietary maps

Within those assumptions, the network is organized as a connected engineering
model rather than a tuned algebraic fit.
"""