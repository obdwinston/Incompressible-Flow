# Projection Method

## Unstructured Grid

- Python code available
- Fortran code coming soon

https://user-images.githubusercontent.com/104728656/230766064-31efa065-57d8-4287-a444-6dd2405f1aba.mp4

<br>

![image](https://user-images.githubusercontent.com/104728656/230762199-63650e6e-1055-4c30-b822-fc7292dc9cde.png)

## Structured Grid

https://user-images.githubusercontent.com/104728656/167190203-dae0d4c2-c8b0-4d87-b88a-5a7cd427a54a.mp4

<br>

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

![image](https://user-images.githubusercontent.com/104728656/166877617-8a88d143-b129-4485-8883-5bfbba8ff17c.png)

![image](https://user-images.githubusercontent.com/104728656/167041265-ce5f4591-1622-456f-bef0-c0d16196355e.png)

![image](https://user-images.githubusercontent.com/104728656/167041276-be49a522-9739-4686-8f5c-9c220dc8989e.png)

![image](https://user-images.githubusercontent.com/104728656/167041325-4de01f47-c319-4e33-9d88-e7efdb2d312f.png)

# PRIME Algorithm

![image](https://user-images.githubusercontent.com/104728656/194483295-bceea236-fdaa-492e-8df8-692a40c5afd6.png)

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

![image](https://user-images.githubusercontent.com/104728656/166672492-2b16dc79-8307-4d7a-9168-df67b9aacc0d.png)

![image](https://user-images.githubusercontent.com/104728656/166154832-c7ff5ac7-df72-4c74-b11e-06c73f407f01.png)

![image](https://user-images.githubusercontent.com/104728656/230038008-f9b0f2f6-62cf-465d-a5d5-bf26551c865c.png)
