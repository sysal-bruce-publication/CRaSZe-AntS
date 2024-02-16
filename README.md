# CRaSZe-AntS
Please see more details in our [paper](https://arxiv.org/abs/2310.04257). Citation:
```
@article{qian2023solving,
  title={On Solving Close Enough Orienteering Problem with Overlapped Neighborhoods},
  author={Qian, Qiuchen and Wang, Yanran and Boyle, David},
  journal={arXiv preprint arXiv:2310.04257},
  year={2023}
}
```

CRaSZe-AntS is a hybrid algorithm for <ins>**C**</ins>lose Enough Orienteering Problem (CEOP) that combines the advantages of the <ins>**Ra**</ins>ndomized <ins>**S**</ins>teiner <ins>**Z**</ins>on<ins>**e**</ins> Discretization scheme and the <ins>**Ant**</ins> Colony <ins>**S**</ins>ystem (ACS). For solving CEOP with different cost functions to collect prizes, CRaSZe-AntS involves different operators. This repository explains how CRaSZe-AntS works in CEOP and CEOP- $\mathcal{N}$ and experimental results achieved by CRaSZe-AntS and benchmark algorithms. 

## Overall design of CRaSZe-AntS
![GraphicalAbstract](figures/GraphicalAbstract.pdf).

## Dataset
This paper mainly uses two types of instances with suffixes `.ceop` and `.tddp`. Because these two types of instances come from the original CETSP dataset, we follow a similar file structure as the CETSP dataset. Specifically, for `.ceop` instances, we have the following file structure:

| Xs | Ys | Zs (not used) | radius | prize |
|----|----|---------------|--------|-------|

For `.tddp` instances, we have the following file structure:

| Xs | Ys | Zs (not used) | radius | prize | parcel weight impact factor |
|----|----|---------------|--------|-----------------------------|-------|

## Experimental results
For experiment results of solving SOP, i.e., examine the discretization performance, please refer to [here](results/SOP).

Please refer to [here](results/CEOP) for experiment results of solving CEOP.

For experiment results of solving TDDP, please refer to [here](results/TDDP).

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

## Evolution process visualization of CRaSZe-AntS and BE-PSO-IACS
As stated in the paper, we show the full evolution process of two algorithms in solving TDDP instance _bubbles1_:

### CRaSZe-AntS

![CRaSZe-AntS](figures/CRaSZe-AntS.gif)

### BE-PSO-IACS

![BE-PSO-IACS](figures/BE-PSO-IACS.gif)
