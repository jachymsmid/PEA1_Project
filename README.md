# Project for the course Programmnig of Engineering Applicatioins
Problem: solve the specified advection equation using some finite difference scheme.

$$u_t+u_x -2p\cdot xu_y= 0$$\
$$u(x,0) = 100\cdot exp(-\frac{(x+1)^2+y^2}{0.01}), \quad u(x,t)|_{\partial \Omega} = 0$$\
$$t \in <0,2>, \quad \Omega = \lbrace [x,y]:x\in <-2,2>, y \in <-2,2> \rbrace$$

Track the maximum value of $u$ in time and the value $\int_{\Omega} u dS$ in time.
### 1. Lax-Friedrichs scheme

### 2. Lax-Wendroff scheme

### 3. Upwind scheme
