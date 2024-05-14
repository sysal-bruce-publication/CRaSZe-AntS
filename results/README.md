# Results
For experiment results of solving SOP, i.e., examine the discretization performance, please refer to [here](./SOP).

Please refer to [here](./CEOP) for experiment results of solving CEOP.

For experiment results of solving TDDP, please refer to [here](./TDDP).

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