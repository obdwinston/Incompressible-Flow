# PRIME Algorithm for Circular Cylinder

![Figure 1](https://user-images.githubusercontent.com/104728656/167181000-57b17068-c12f-4cba-989d-fbc538ff4e6c.png)

- PRessure Implicit Momentum Explicit (PRIME) algorithm by Maliska and Raithby (1984)
- 2D steady incompressible viscous laminar flow
- finite volume method on collocated unstructured triangle mesh
- Rhie-Chow interpolation to suppress chequerboard oscillations
- volume-based interpolation for cell face values
- Green-Gauss method for gradient computation
- high resolution scheme for convection term
- TVD framework with MUSCL limiter function for anti-diffusive convection term
- second-order central differencing scheme for diffusion term
- over-relaxed approach for cross-diffusion term
- deferred correction for anti-diffusive convection and cross-diffusion terms
- Gauss-Seidel method for pressure correction equation
- flow conditions
  - Re = 20, U = 1.0 m/s, ν = 0.01 m²/s
  - H = 1.0 m, L = 3.0 m, D = 0.2 m
- boundary conditions
  - u = U and v = 0 at inlet boundary
  - u = 0 and v = 0 at wall and object boundaries
  - du/dn = 0 and dv/dn = 0 at outlet boundary
  - p = 0 at outlet boundary
  - dp/dn = 0 at inlet, wall, and object boundaries
- GMSH program required for mesh generation
  - label 'Physical Curve' boundaries as 'WALL', 'OBJECT', 'INLET', and 'OUTLET'
  - Geometry → Physical Groups → Add → Curve
  - export and save as .su2 mesh file
  - File → Export → Save → Save all elements
- OpenFOAM verification for speed distribution
  - SimpleFOAM (Laminar) → U<sub>min</sub> = 0.0039074 m/s, U<sub>max</sub> = 1.68126 m/s
  - PRIME Algorithm → U<sub>min</sub> = 0.0057518 m/s, U<sub>max</sub> = 1.68645 m/s

![Figure 2](https://user-images.githubusercontent.com/104728656/166672492-2b16dc79-8307-4d7a-9168-df67b9aacc0d.png)

![Figure 3](https://user-images.githubusercontent.com/104728656/166154832-c7ff5ac7-df72-4c74-b11e-06c73f407f01.png)

![Figure 4](https://user-images.githubusercontent.com/104728656/166466869-2388b035-9774-4404-86e9-08338cc54db9.png)

# Projection Method for Square Cylinder

https://user-images.githubusercontent.com/104728656/167190203-dae0d4c2-c8b0-4d87-b88a-5a7cd427a54a.mp4

- fractional step projection method by Chorin (1967)
- 2D unsteady incompressible viscous laminar flow
- finite volume method on staggered mesh
- first-order upwind scheme for convection term
- second-order central differencing scheme for diffusion term
- Gauss-Seidel method for pressure correction equation
- flow conditions
  - Re = 100, U = 5.0 m/s, ν = 0.01 m²/s
  - H = 1.0 m, L = 3.0 m, D = 0.2 m
- boundary conditions
  - u = U and v = 0 at inlet boundary
  - u = 0 and v = 0 at wall and object boundaries
  - du/dn = 0 and dv/dn = 0 at outlet boundary
  - p = 0 at outlet boundary
  - dp/dn = 0 at inlet, wall, and object boundaries
- OpenFOAM verification for vortex shedding frequency
  - IcoFOAM → St ≈ 0.22
  - Projection Method → St ≈ 0.20

![Figure 1](https://user-images.githubusercontent.com/104728656/166877617-8a88d143-b129-4485-8883-5bfbba8ff17c.png)

![Figure 2](https://user-images.githubusercontent.com/104728656/167041265-ce5f4591-1622-456f-bef0-c0d16196355e.png)

![Figure 3](https://user-images.githubusercontent.com/104728656/167041276-be49a522-9739-4686-8f5c-9c220dc8989e.png)

![Figure 4](https://user-images.githubusercontent.com/104728656/167041325-4de01f47-c319-4e33-9d88-e7efdb2d312f.png)
