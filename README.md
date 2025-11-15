# Project for the course Programmnig of Engineering Applicatioins
Problem: solve the specified advection equation using some finite difference scheme.

$$u_t+u_x -2x\cdot u_y= 0$$\
$$u(x,0) = 100\cdot exp(-\frac{(x+1)^2+y^2}{0.01}), \quad u(x,t)|_{\partial \Omega} = 0$$\
$$t \in <0,2>, \quad \Omega = \lbrace [x,y]:x\in \langle -2,2\rangle, y \in \langle -2,2\rangle \rbrace$$

Track the maximum value of $u$ in time and the value $\int_{\Omega} u\ dS$ in time.

# Documentation

## SimulationInfo

SimulationInfo is a struct holding information about the simulatioin. It comprises the following data:

- step_t : time step
- step_x : spatial step in the 'x' direction
- step_y : spatial step in the 'y' direction
- x_min : lower bound of the domain in the 'x' direction
- x_max : upper bound of the domain in the 'x' direction
- y_min : lower bound of the domain in the 'y' direction
- y_min : upper bound of the domain in the 'y' direction
- number_x : number of discretization points in the 'x' direction
- number_y : number of discretization points in the 'y' direction

All the data is constant so it has to be defined in the constructor. This struct is initialized in the following way:
```
SimulatinInfo sim_info( step_t,
                        step_x,
                        step_y,
                        x_min,
                        y_min,
                        x_max,
                        y_max,
                        number_x,
                        number_y);
```

## Mesh

The Mesh class holds the spatial coordinates and the data associated with them. Its attributes are:
- array of data points
- array of x coordinates
- array of y coordinates
- number of rows
- number of columns

All the arrays are flattened using the Column-major indexing. Methods of the class are:
- constructor using the SimulationInfo struct
```
SimulationInfo sim_info;
Mesh data( sim_info );
```
- copy constructor
```
Mesh data;
Mesh data_copy( data );
```
- `float value_ref( size_t i, size_t j )` : getter for the value referenc at index i,j. Indexes should be of type `size_t`. A value gotten in this way can be changed.
```
Mesh mesh( sim_info );
float u_ij;
u_ij = mesh.value_ref( i, j );
u_ij = 5;
```
- `float value( size_t i, size_t j )` : getter for the value (constant) at index i,j. Indices should be ot type `size_t`.
```
Mesh mesh( sim_info );
u_ij = mesh.value( i, j );
```
- `size_t getCols()` : returns the number of columns of the Mesh.
- `size_t getRows()` : returns the number of rows of the Mesh.
- `void construct_regular_grid( float dx, float dy )` : constructs the grid for constant spatial steps.
- `void write_data( string file_name )` : writes the data in the .vti format. The '.vti' should be included in the file_name variable.

## NumericalSolver
NumericalSolver< NumericalScheme >( Mesh &mesh, SimulationInfo sim_info ) is a function that moves the data to the next time step.
```
SimulationInfo sim_info( ... );
Mesh mesh( sim_info );
NumericalSolver< Upwind >( mesh, sim_info );
```
There are three numerical schemes available: `Lax_Wendroff`, `Lax_Friedrichs`, `Upwind`

### 1. Lax-Friedrichs scheme

### 2. Lax-Wendroff scheme
The equation can be discretized using the Lax-Wendroff scheme. The scheme will have the following form in 2D:

$$U_{i,j}^{n+1} = U_{i,j}^n - \frac{dt}{2dx}(U_{i+1,j} - U_{i-1,j})-\frac{x_{i,j}dt}{dy}(U_{i,j+1}-U_{i,j-1})+\frac{dt^2}{2dx^2}(U_{i+1,j}-2U_{i,j}+U_{i-1,j})+\frac{2(x_{i,j}dt)^2}{dy^2}(U_{i,j+1}-2U_{i,j}+U_{i,j-1})+\frac{x_{i,j}dt^2}{2dxdy}(U_{i+1,j+1}-U_{i+1,j-1} - U_{i-1,j+1}+U_{i-1,j-1})$$

The CFL condition for this scheme is:

$$dt \leq \frac{1}{\frac{1}{dx}+\frac{1}{2x_m dx}}$$,

where $x_m = \max |x|\quad x \in \Omega$.

### 3. Upwind scheme

## CFl
CFL is a templated function that computes maximal time step from the used numerical scheme.
```
float dt = CFL< Lax_Wendroff>( float dx, float dy );
```
Note that time step dt is needed for the SimulationStruct initialization.

## InitialConditions
InitialConditions is a templated function that applays specified initial conditions to the data stored in Mesh.
```
Mesh mesh;
InitialConditions< My_Initial_Conditions >( mesh );
```
There are currently only the `My_Initial_Conditions` available.

### 1. My_Initial_Conditions

## BoundaryConditions
BoundaryConditions is a templated function that applays specified boundary conditions to the data stored in Mesh.
```
Mesh mesh;
InitialConditions< Zeros >( mesh );
```
There are currently only the `Zeros` boundary conditioins available.

### 1. Zeros
This boundary condtion sets all the elements on the boundary to zero.

