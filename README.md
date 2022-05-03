# PRIME Algorithm for Circular Cylinder

![Figure 1](https://user-images.githubusercontent.com/104728656/166450612-c7372f9b-b101-48fe-b2d3-8f51e8ba4bc4.png)

![Figure 2](https://user-images.githubusercontent.com/104728656/166154801-38149220-ac17-4c17-b148-3ff7f406d703.png)

- PRessure Implicit Momentum Explicit (PRIME) algorithm by Maliska and Raithby (1984)
- 2D steady incompressible viscous laminar flow
- finite volume method on collocated grid
- Rhie-Chow interpolation to suppress chequerboard oscillations
- volume-based interpolation for cell face values
- Green-Gauss method for gradient computation
- high resolution scheme for convection term
- TVD framework with MUSCL limiter function for anti-diffusive convection term
- second-order central differencing scheme for diffusion term
- over-relaxed approach for cross-diffusion term
- deferred correction for anti-diffusive convection and cross-diffusion terms
- Gauss-Seidel method for pressure correction equation
- boundary conditions
  - u = U and v = 0 at inlet boundary
  - u = 0 and v = 0 at wall and object boundaries
  - du/dn = 0 and dv/dn = 0 at outlet boundary
  - p = 0 at outlet boundary
  - dp/dn = 0 at inlet, wall, and object boundaries
- GMSH program required for mesh generation
  - label 'Physical Curve' boundaries as 'WALL', 'OBJECT', 'INLET', and 'OUTLET'
  - Geometry → Physical Groups → Add → Curve
  - export and save as .SU2 mesh file
  - File → Export → Save → Save all elements
- OpenFOAM verification for speed distribution
  - SimpleFOAM (Laminar) → U<sub>min</sub> = 0.0039074 m/s, U<sub>max</sub> = 1.68126 m/s
  - PRIME Algorithm → U<sub>min</sub> = 0.0057518 m/s, U<sub>max</sub> = 1.68645 m/s

![Figure 3](https://user-images.githubusercontent.com/104728656/166154832-c7ff5ac7-df72-4c74-b11e-06c73f407f01.png)

![Figure 4](https://user-images.githubusercontent.com/104728656/166428131-4c008a11-1295-40a7-af53-d61df99c3171.png)
