## Program Overview

This program implements Chorin's projection method to solve the incompressible Navier-Stokes equations on an unstructured triangular mesh using the finite volume method. The convection fluxes uses hybrid differencing to switch between centred and upwind, depending on the local Peclet number. The high-order upwind differencing employs a TVD scheme with the Superbee limiter.

![](https://github.com/user-attachments/assets/d5e43c2a-450d-4ac3-b2fa-a5695acf615f)

## Program Files

```
root/
│
├── data/                   # saved data folder
│   └── speed.mp4           # speed animation
│
├── mesh/                   # mesh folder
│   ├── body.txt            # body coordinates (user input)
│   ├── body.py             # generates body.txt
│   ├── geo.py              # generates mesh.geo from body.txt
│   ├── su2.py              # generates mesh.su2 from mesh.geo
│   └── mesh.f90            # generates mesh.txt from mesh.su2
│
├── mods/                   # modules folder
│   ├── mod_mesh.f90        # mesh type and procedures
│   ├── mod_config.f90      # configuration type and procedures
│   ├── mod_solve.f90       # solver procedures
│   └── mod_utils.f90       # utility procedures
│
├── run.sh                  # script to run program
├── main.f90                # script to run solver
├── read.py                 # script to read and plot saved data
├── config.txt              # configuration for solver (user input)
└── requirements.txt        # dependencies for solver
```

Clone repository:

```bash
git clone https://github.com/obdwinston/Incompressible-Flow.git && cd Incompressible-Flow
```

Execute program (for macOS users):

```bash
chmod +x run.sh && ./run.sh
```

For Windows users, you need to modify `run.sh` accordingly before executing the program. For custom bodies, coordinates in `mesh/body.txt` should be `x y` space-delimited and in clockwise order, with no repeated points or intersecting lines. To visualise the mesh, download Gmsh [here](https://gmsh.info/#Download). Open the `.geo` file and select `Mesh > 2D` to generate the mesh.

## Solver Verification

https://github.com/user-attachments/assets/a05195a0-c8fa-4a11-b4ce-0961b1f41902

## Solver Theory

![](https://github.com/user-attachments/assets/a6281bda-5e6a-4647-8b80-30ce33694fb9)

## References

[1] Moukalled et al. (2016). _The Finite Volume Method in Computational Fluid Dynamics._  
[2] Mazumder (2016). _Numerical Methods for Partial Differential Equations._  
[3] Versteeg et al. (2015). _An Introduction to Computational Fluid Dynamics._  
[4] Saad (2020). _Computational Fluid Dynamics [Lecture Series](https://www.youtube.com/watch?v=sSqtgi0zqT8&list=PLEaLl6Sf-KIC7oet7zvNfW03aocrIq-s4) (Spring 2019)._
