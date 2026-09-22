---
title: "plot_grids.sh"
nav_order: 14
render_with_liquid: false
parent: "Running JOREK"
---

## plot_grids.sh
This utility plots the grid to the screen, to a png file or a postscript file.

  ```
  Usage: plot_grids.sh [options]

    -h         Print this usage information
    -o <grid>  Plot only the specified grid
    -ps        Plot to a .ps file (default: plot to screen)
    -png       Plot to a .png file (default: plot to screen)
    -r <resol> Resolution of a png plot (default: '800x800')

  Example: plot_grids.sh -o xpoint -png -r 1600x1600
           plot_grids.sh -o initial -ps -r 1600x1600
  ```
