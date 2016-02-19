Solving 1-D Heat Equation and Burger's Equation using Finite Difference Method
===================

## 1-D Heat equation
We solve the 1-D heat equation
```math
\begin{aligned}
u_t(x,t) &= v u_xx(x,t) + f(x,t) &\forall [x,t] \in [0,1]\times[0,1] \\
 u(0,t) &= u(1,t) = 0 &\forall t \in [0,1] \\
u(x,0) &= u_0(x)  &\forall x \in [0,1]
\end{aligned}
```
by the implicit Euler's method.

First, we partition the interval of $x$ using $m+1$ equally spaced nodes and the interval of $t$ using $n+1$ equally spaced nodes. At time $t_{j+1}$ we approximate $u_t$ and $u_{xx}$ using the 2-point backward finite difference formula and the 3-point centered finite difference formula, respectively. This gives

\begin{align}
\frac{u_i^{j+1} - u_i^{j}}{\Delta t} - \nu \frac{u_{i+1}^{j+1} - 2u_{i}^{j+1} + u_{i-1}^{j+1}}{(\Delta x)^2} &= f(x_i,t_{j+1}) \\
u_{i-1}^{j+1}\left(-\frac{\nu \Delta t}{(\Delta x)^2}\right) + u_i^{j+1} \left(1 + \frac{2\nu \Delta t}{(\Delta t)^2}\right) + u_{i+1}^{j+1} \left(-\frac{\nu \Delta t}{(\Delta x)^2}\right) &= u_i^j + \Delta t f(x_i,t_{j+1}) \\
\end{align}
for $i = 0,1,2,...,m$ and $j = 0,1,2,...,n$.

Let $\alpha = -\frac{\nu \Delta t}{(\Delta x)^2}$ and $\gamma = 1 - 2\alpha$, the above system of equations can be written as

\begin{equation*}
  \begin{pmatrix}
    \alpha & 1-2\alpha &\alpha     &     &     &\\
        &  \alpha & 1-2\alpha &\alpha     &     &\\
        &   & \ddots & \ddots &    \ddots &\\
        &     &      & \alpha & 1-2\alpha &\alpha  \\
        &     &      &     &\alpha & 1-2\alpha &\alpha
  \end{pmatrix}
  \begin{pmatrix}
    u_0^{j+1}\\
    u_1^{j+1}\\
    \vdots\\
     u_{m-1}^{j+1}\\
     u_m^{j+1}
  \end{pmatrix}=
  \begin{pmatrix}
     u_0^{j} + \Delta t f(x_0,t_{j+1})\\
    u_1^{j} + \Delta t f(x_1,t_{j+1})\\
    \vdots\\
     u_{m-1}^{j} + \Delta t f(x_{m-1},t_{j+1})\\
     u_m^{j} + \Delta t f(x_m,t_{j+1})
  \end{pmatrix}
\end{equation*}

To enforce Dirichlet boundary conditions at $x_0$ we replace the first row of the coefficient matrix with all zeros and one at the first position, and the first entry in the right hand side vector with the boundary condition. To enforce Dirichlet boundary conditions at $x_m$ we replace the last row of the coefficient matrix with all zeros and one at the last position, and the last entry in the right hand side vector with the boundary condition. We eliminate the boundary conditions by substituting the values of $u_0^{j+1}$ and $u_m^{j+1}$ into the second and penultimate equations, respectively. This gives an $(m-1) \times (m-1)$ system of equations

\begin{equation*}
  \begin{pmatrix}
    1-2\alpha &\alpha     & 0 &    &     &\\
     \alpha & 1-2\alpha &\alpha     &     &\\
         & \ddots & \ddots &    \ddots &\\
             &      & \alpha & 1-2\alpha &\alpha  \\
             &      &     &\alpha & 1-2\alpha
  \end{pmatrix}
  \begin{pmatrix}
    u_1^{j+1}\\
    u_2^{j+1}\\
    \vdots\\
     u_{m-2}^{j+1}\\
     u_{m-1}^{j+1}
  \end{pmatrix}=
  \begin{pmatrix}
    u_1^{j} + \Delta t f(x_1,t_{j+1}) - \alpha BC1(t_{j+1})\\
    u_2^{j} + \Delta t f(x_2,t_{j+1})\\
    \vdots\\
     u_{m-2}^{j} + \Delta t f(x_{m-2},t_{j+1})  - \alpha BC2(t_{j+1})\\
     u_{m-1}^{j} + \Delta t f(x_{m-1},t_{j+1})\\
  \end{pmatrix}
\end{equation*}
where $BC1$ and $BC2$ are the boundary conditions at $x_0$ and $x_m$, respectively.

At each time $t_{j+1}$ we solve the above system to get the solution vector at time $t_{j+1}$, for $j = 0,1,2,...,n$. This is implemented in $\tt heateq.m$.

The method is implicit and unconditionally stable. For $\Delta t = 0.25$ and $\Delta x = 0.2$, the approximation is within two accurate digits from the exact solution.

### Validation problem
Here we construct the validation problem. Let the solution $u(x,t) = e^{-t} \sin(\pi x)$. Note that this solution satisfies the specified boundary conditions. Then
\begin{align*}
u(x,0) &= \sin(\pi x)\\
u_t &= -\sin(\pi x) e^{-t} \text{  and}\\
u_xx &= -(\pi)^2 \sin(\pi x) \
\end{align*}
Substituting in the original heat equation gives

\begin{equation*}
f(x,t) = -\sin(\pi x) e^{-t} + \nu e^{-t} \sin(\pi x)
\end{equation*}

The constructed problem
\begin{align*}
u_t - \nu u_xx &= -\sin(\pi x) e^{-t} + \nu e^{-t} \sin(\pi x) &   &\forall x \in [0,1] \text{ and } \forall t \in [0,1] \\
u(0,t) &= u(1,t) = 0 &   & \forall t \in [0,1] \\
u(x,0) &= \sin(\pi x)  &   & \forall x \in [0,1]
\end{align*}
has exact solution $u(x,t) = e^{-t} \sin(\pi x)$.

### Convergence rate
We determine numerically the order of convergence for the local (truncation) error. The local (truncation) error is the error made at the discretization using $\Delta x$. The local (truncation) error at a fixed time $t^*$ in $x$ is defined as
\begin{equation*}
\epsilon = \max_{1\leq j \leq m} | u(x_j) - u_j| = C\Delta x^p = O(\Delta x^p)
\end{equation*}
where $p$ is the order of convergence. $p$ is determined as follows
\begin{equation*}
\log_{10} \epsilon = \log_{10} C\Delta x^p = \log_{10} C + p \log_{10} \Delta t \\
\end{equation*}
We compute $p$ over several grids with $\Delta x$ decreasing by a factor of $1/2$ at each grid. Then $p$ is
\begin{equation*}
p = \frac{log_{10}\epsilon_{\Delta x} -  \log_{10}\epsilon_{\Delta x / 2}}{\log_{10}\Delta x - \log_{10}(\Delta x/2)} = \frac{log_{10}\epsilon_{\Delta x} - \log_{10}\epsilon_{\Delta x /2}}{\log_{10}2}
\end{equation*}
We use $t^* = t_{final} = 0$ to eliminate the need of determining exact index $j+1$ for $t$ and simplify the code. For the error, we use the infinity norm of the error vector at time $t^* = 1$.

Convergence in $t$ is determined analogously.  

Convergence analysis is carried out in $\tt convHeateq.m$.

## Burger equation
Here we want to solve the 1-D Burger's equation
\begin{aligned}
u_t(x,t) &= \nu u_xx(x,t) - u(x,t)  u_x(x,t) + f(x,t) &\forall [x,t] \in [0,1]\times[0,1] \\
u(0,t) &= u(1,t) = 0 &\forall t \in [0,1] \\
u(x,0) &= u_0(x)  &\forall x \in [0,1]
\end{aligned}
using backward Euler's method.

We partition the interval of $x$ using $m+1$ equally spaced nodes and the interval of $t$ using $n+1$ equally spaced nodes.

At time $t_{j+1}$ we approximate $u_t$, $u_x$, and $u_{xx}$ using the 2-point backward finite difference formula, 2-point centered finite difference formula, and the 3-point centered finite difference formula, respectively. This gives
\begin{equation*}
\frac{u_i^{j+1} - u_i^{j}}{\Delta t} + u_i^{j+1} \frac{u_{i+1}^{j+1} - u_{i-1}^{j+1}}{2 \Delta x} - \nu \frac{u_{i+1}^{j+1} - 2u_{i}^{j+1} + u_{i-1}^{j+1}}{(\Delta x)^2} = f(x_i,t_{j+1})
\end{equation*}

Letting $\alpha = \frac{\Delta t}{2 \Delta x}$ and $\gamma = \frac{\nu \Delta t}{(\Delta x)^2}$, and grouping like terms, gives
\begin{equation*}
u_{i-1}^{j+1}\left(-\alpha u_i^{j+1} - \gamma \right) + u_i^{j+1} \left(1 + \gamma \right) + u_{i+1}^{j+1} \left(\alpha u_i^{j+1} - \gamma \right) - u_i^j - \Delta t f(x_i,t_{j+1}) = 0\\
\end{equation*}
for $i = 0,1,2,...,m$ and $j = 0,1,2,...,n$.


To enforce Dirichlet boundary conditions at $x_0$ we replace the first equation by $u_0^{j+1} - BC1 = 0$; to enforce Dirichlet boundary conditions at $x_m$ we replace the last equation by $u_m^{j+1} - BC2 = 0$, where $BC1$ and $BC2$ are the boundary conditions at $x_0$ and $x_m$, respectively.

The above system of equation is non-linear. For more robust approximations, we use Matlab's built-in $\tt fsolve$. At each time $t_{j+1}$ we solve the above system to get the solution vector at time $t_{j+1}$, for $j = 0,1,2,...,n$. This is implemented in $\tt iburger.m$.

The method used in the code is implicit and unconditionally stable. For $\Delta t = 0.25$ and $\Delta x = 0.2$, the approximation is within 2 accurate digits from the exact solution.

Convergence analysis is similar to that of the 1-D Heat equation and is carried out in $\tt convIburger.m$.

### Validation problem
Here we construct the validation problem. Let the solution $u(x,t) = e^{-t} \sin(\pi x)$. Note that this solution satisfies the specified boundary conditions. Then
\begin{align*}
u(x,0) &= \sin(\pi x)\\
u_t &= -\sin(\pi x) e^{-t} \\
u_x &= \pi e^{-t} \cos(\pi x) \\
u_{xx} &= -(\pi)^2 \sin(\pi x) \
\end{align*}
Substituting in the original Burger's equation gives
\begin{equation*}
f(x,t) = -e^{-t} \sin(\pi x) + \pi e^{-2t} \sin(\pi x) \cos(\pi x) + \nu \pi^2 e^{-t} \sin(\pi x) \end{equation*}

The constructed problem
\begin{align*}
u_t + u u_x - \nu u_{xx} &= -e^{-t} \sin(\pi x) + \pi e^{-2t} \sin(\pi x) \cos(\pi x) + \nu \pi^2 e^{-t} \sin(\pi x) &   &\forall x \in [0,1] \text{ and }  \forall t \in [0,1] \\
u(0,t) &= u(1,t) = 0 &   & \forall t \in [0,1] \\
u(x,0) &= \sin(\pi x)  &   & \forall x \in [0,1]
\end{align*}
has exact solution $u(x,t) = e^{-t} \sin(\pi x)$.
