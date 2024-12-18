# Stable Dynamical Systems on Riemannian Manifolds (SDS-RM)

An implementation of the SDS-RM approach described in [Saveriano, Abu-dakka, and Kyrki, 2023](https://www.sciencedirect.com/science/article/pii/S0921889023001495).

## Demos description
- `main_R_LASA_SPD.m`: a demo to run SDS-RM on the Riemannian LASA (R-LASA) Handwriting dataset using Symmetric and Positive Definite (SPD) matrices.
- `main_R_LASA_UQ.m`: a demo to run SDS-RM on the Riemannian LASA (R-LASA) Handwriting dataset using Unit Quaternions (UQ).

## Software Requirements
The code is developed and tested under `Matlab2023b`.

## References
Please acknowledge the authors in any academic publication that used parts of these codes.
```
@article{Saveriano2023Learning,
  title={Learning stable robotic skills on Riemannian manifolds},
  author = {Matteo Saveriano and Fares J. Abu-Dakka and Ville Kyrki},
journal = {Robotics and Autonomous Systems},
volume = {169},
pages = {104510},
year = {2023}
}
```

## Third-party material
Third-party code and dataset have been included in this repository for convenience.

- **LASA Handwriting dataset**: please acknowledge the authors in any academic publications that have made use of the LASA HandWriting dataset by citing: *S. M. Khansari-Zadeh and A. Billard, "Learning Stable Non-Linear Dynamical Systems with Gaussian Mixture Models", IEEE Transaction on Robotics, 2011*.

- **GMR library**: please acknowledge the authors in any academic publications that have made use of the GMR library by citing: *S. Calinon et al., "On Learning, Representing and Generalizing a Task in a Humanoid Robot", IEEE Transactions on Systems, Man and Cybernetics, Part B., 2006*.

## Note
This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
