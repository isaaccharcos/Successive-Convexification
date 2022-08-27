# Successive Convexification
<p align="center">
  <img src="https://github.com/isaaccharcos/Successive-Convexification/blob/main/SCvx%20Iterations.gif" />
</p>

Implementation of [Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time](https://arxiv.org/abs/1802.03827) by Michael Szmuk, Behcet Acikmese.
## Requirements
- [CVX](http://cvxr.com/cvx/download/) for MATLAB.
```
cd ~/MATLAB/cvx
cvx_setup
```
## How to run
1. Run `SCvx.m` to generate `results.mat` file which contains solved trajectory information.
2. Run `plot_trajectory.m` to generate plots and gifs of any of the trajectories in `results.mat`.
3. Run `animate_iterations.m` to generate gif of all computed SCvx iterations.
