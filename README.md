# CRaSZe-AntS
CRaSZe-AntS is a hybrid algorithm for <ins>**C**</ins>lose Enough Orienteering Problem (CEOP) that combines the advantages of the <ins>**Ra**</ins>ndomized <ins>**S**</ins>teiner <ins>**Z**</ins>on<ins>**e**</ins> Discretization scheme and the <ins>**Ant**</ins> Colony <ins>**S**</ins>ystem (ACS). For solving CEOP with different cost functions to collect prizes, CRaSZe-AntS involves different operators. This repository explains how CRaSZe-AntS works in CEOP and CEOP- $\mathcal{N}$ and experimental results achieved by CRaSZe-AntS and benchmark algorithms. CRaSZe-AntS can solve three problems:
1. SOP: The [Set Orienteering Problem](https://www.sciencedirect.com/science/article/pii/S0377221717310202). Here our SOP allows at most **one** visit to the vertex in each group.
2. CEOP: The [Close Enough Orienteering Problem](https://ieeexplore.ieee.org/document/7844421).
3. TDDP: The Truck-and-Drone Delivery Problem. Please check our paper for the formal formulation.

Please see more details in our [paper](https://arxiv.org/abs/2310.04257). Citation:
```
@article{qian2023solving,
  title={On Solving Close Enough Orienteering Problem with Overlapped Neighborhoods},
  author={Qian, Qiuchen and Wang, Yanran and Boyle, David},
  journal={arXiv preprint arXiv:2310.04257},
  year={2023}
}
```

## Overall design of CRaSZe-AntS
![GraphicalAbstract](figures/GraphicalAbstract.png)

## Quick Start

### Prerequisites
The source code is written in ISO C++20. The binary code has been tested on Red Hat Enterprise Linux (version 8.5). The dependencies are listed as follows:
* g++ (GCC) 8.5.0 20210514 (Red Hat 8.5.0-4)
* [nlohmann/json](https://github.com/nlohmann/json)
* Python 3.11.3 (for scripts)

### Build
Because we use `fstream` library for data I/O, to compile, please use the following command:
```bash
g++ -std=c++2a -o crazyants <path-to-all-src>*.cpp -lstdc++fs
```

### Instructions
1. Create directory structure as follows:
```
workspace/
├── crazyants (binary file)
├── configs/
│   ├── bubbles.tddp.json
├── data/
│   ├── tddp/
│   │   ├── bubbles1.tddp
│   │   ├── ...
│   │   ├── sz-vtxs/
│   │   │   ├── bubbles1.szppt
│   │   │   ├── ...
│   │   ├── sols/
│   │   │   ├── bubbles1.tddp_sol
│   │   │   ├── ...
├── scripts/
│   ├── myscript.sh
│   ├── applyBudgetRatioCEOP.py
├── src/
│   ├── tddp/
│   │   ├── PSO.h
│   │   ├── PSO.cpp
│   │   ├── ...
```
2. Run CRaSZe-AntS with **script**, check [example script](scripts/example_script.sh) for more details:
```bash
cd scripts
./example_script.sh
```
3. Check results in `workspace/data/tddp/sols/`.

## Evolution process visualization
Below we show the full evolution process of two algorithms in solving TDDP instance _bubbles1_:

### CRaSZe-AntS

![CRaSZe-AntS](figures/CRaSZe-AntS.gif)

### BE-PSO-IACS

![BE-PSO-IACS](figures/BE-PSO-IACS.gif)