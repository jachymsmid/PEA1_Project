# Project for the course Programmnig of Engineering Applicatioins
Problem: solve the specified advection equation using some finite difference scheme.

$$u_t+u_x -2x\cdot u_y= 0$$\
$$u(x,0) = 100\cdot exp(-\frac{(x+1)^2+y^2}{0.01}), \quad u(x,t)|_{\partial \Omega} = 0$$\
$$t \in <0,2>, \quad \Omega = \lbrace [x,y]:x\in \langle -2,2\rangle, y \in \langle -2,2\rangle \rbrace$$

Track the maximum value of $u$ in time and the value $\int_{\Omega} u\ dS$ in time.
### 1. Lax-Friedrichs scheme

### 2. Lax-Wendroff scheme
The equation can be discretized using the Lax-Wendroff scheme. The scheme will have the following form in 2D:

$$U_{i,j}^{n+1} = U_{i,j}^n - \frac{dt}{2dx}(U_{i+1,j} - U_{i-1,j})-\frac{x_{i,j}dt}{dy}(U_{i,j+1}-U_{i,j-1})+\frac{dt^2}{2dx^2}(U_{i+1,j}-2U_{i,j}+U_{i-1,j})+\frac{2(x_{i,j}dt)^2}{dy^2}(U_{i,j+1}-2U_{i,j}+U_{i,j-1})+\frac{x_{i,j}dt^2}{2dxdy}(U_{i+1,j+1}-U_{i+1,j-1} - U_{i-1,j+1}+U_{i-1,j-1})$$

The CFL condition for this scheme is:

$$dt \leq \frac{1}{\frac{1}{dx}+\frac{1}{2x_m dx}}$$,

where $x_m = \max |x|\quad x \in \Omega$.
### 3. Upwind scheme
