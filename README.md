# CRaSZe-AntS
CRaSZe-AntS is a hybrid algorithm for <ins>**C**</ins>lose Enough Orienteering Problem (CEOP) that combines the advantages of the <ins>**Ra**</ins>ndomized <ins>**S**</ins>teiner <ins>**Z**</ins>on<ins>**e**</ins> Discretization scheme and the <ins>**Ant**</ins> Colony <ins>**S**</ins>ystem (ACS). For solving CEOP with different cost functions to collect prizes, CRaSZe-AntS involves different operators. This repository contains the explanation how CRaSZe-AntS works in CEOP and CEOP-$\mathcal{N}$, and experimental results achieved by CRaSZe-AntS and benchmark algorithms. 

## Work flow of CRaSZe-AntS
The work flow of CRaSZe-AntS for solving **CEOP** can be found [here](figures/CEOP.pdf).

The work flow of CRaSZe-AntS for solving **CEOP-$\mathcal{N}$** can be found [here](figures/CEOPN.pdf).

## Dataset
In this paper, we mainly use two types of instances with suffix `.ceop` and `.tddp`. Because these two types of instances comes from the original CETSP dataset, we follow a similar file structure as CETSP dataset. Specifically, for `.ceop` instances, we have the following file structure:

| Xs | Ys | Zs (not used) | radius | prize |
|----|----|---------------|--------|-------|

For `.tddp` instances, we have the following file structure:

| Xs | Ys | Zs (not used) | radius | prize | parcel weight impact factor |
|----|----|---------------|--------|-----------------------------|-------|

## Experimental results
For experiment results of solving SOP, i.e. examine the discretization performance, please refer to [here](results/SOP).

For experiment results of solving CEOP, i.e. examine the discretization performance, please refer to [here](results/CEOP).

For experiment results of solving TDDP, i.e. examine the discretization performance, please refer to [here](results/TDDP).

Note that in each `solutions/` directory, the solution files have the following format:

| timestamp | number of waypoints | cost | prize | waypoints (x1, y1, x2, y2, ..., xn, yn) |
|---|---|---|---|-----------------------------------------|

In `steiner-zone-vertices/` directory, SZ vertex files have the following format:

| vertex ID | SZ prize | Xs | Ys |
|-----------|----------|----|----|

Similarly, in `discrete-points/` directory, discrete point files have the following format:

| point ID | circle prize | Xs | Ys |
|----------|--------------|----|----|

and convex hull files have the following format:

| The target circle ID in original dataset that are convex hull vertices |
|------------------------------------------------------------------------|