---
title: "Reads and plot ''macroscopic_vars.dat''"
nav_order: 15
render_with_liquid: false
parent: "Running JOREK"
---

#### Reads and plot ''macroscopic_vars.dat''
This is a Python analog of [plot_live_data.sh](plot_live_data.sh.md) in case `gnuplot` isn't stable, e.g. when using x-terminal.

Usage:
```
plot_live_data.py
```
Options:
```
  -h, --help  show this help message and exit
  -legend     Print legend
  -q Q        Plot given quantity
  -f F        File name
  -l          List plottable quantities
  -ps         Save plot into eps file
  -nology     Linear y axis
```
