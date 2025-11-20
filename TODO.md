
## TODO list for the Advection equation solver
**TODO:**

- [ ] modularization
    - [ ] Mesh
    - [ ] SimulationInfo
    - [ ] InitialConditions
    - [ ] BoundaryConditions
    - [ ] NumericalSchemes
        - [ ] cfl conditions
- [ ] better mesh representation
    - [x] change the array of structer paradigm to structure of array
    - [x] copy constructor
    - [x] implement flattened arrays 
       - [ ] Z-order curves
- [x] VTK library for data storage and visualization
    - [ ] rewrite so that no library is needed
- [ ] error handling
- [ ] implicit grid - possible bc of regular 2D grid

***

**DONE:**

- [x] upwind method
- [x] lax-friedrichs method
- [x] lax-wendroff method
- [x] SimulationInfo struct
    - [x] constructor
- [x] implement the numerical solver
- [x] implement boundary/initial conditions
