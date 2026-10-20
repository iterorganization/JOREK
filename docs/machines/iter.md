---
title: "ITER"
nav_order: 1
parent: "Machines"
layout: default
render_with_liquid: false
---

# Information relevant for simulating ITER

To obtain ITER equilibria, machine descriptions, and other simulation data through IMAS, see the [JOREK-IMAS documentation](../diagnostics/jorek-imas.md).

## Reference scenarios

[Zipped folder](assets/iter/201905_iter_scenarios.zip) containing data from **CORSICA simulations of ITER reference scenarios**. Contains G-files (i.e. eqdsk files, which can be converted into JOREK input using [eqdsk2jorek.f90](../howto/eqdsk2jorek), P-files (profiles) and T-files (time traces) for the following scenarios:

- 15MA/5.3T H-mode (D-T plasma)
- 15MA/5.3T L-mode (H plasma)
- 7.5MA/2.65T H-mode (He plasma)
- 5MA/1.77T H-mode (H plasma)

**Note:** the data should also be accessible through **IMAS**, using the name given in the Excel spreadsheet contained in the above folder, but I (Eric) have not tested this. Would the first person who tries this please confirm (by email to eric.nardon@cea.fr) that this works?

[Slides](assets/iter/iter_reference_scenario_20190702_v3.pdf) from A. Matsuyama showing some figures of the above scenarios and indicating **which time slices are recommended** by the ITER DMS Task Force for disruption simulations.

The **contact person** on the ITER side ("owner" of the data): S.H. Kim [sunhee.kim@iter.org](mailto:sunhee.kim@iter.org). If you have a request or a problem with the data, please also inform Eric [eric.nardon@cea.fr](mailto:eric.nardon@cea.fr). 
