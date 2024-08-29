# Parametric instability of the interface between two fluids

We briefly describe an implementation of Kumar \& Tuckerman [^1]. We consider
two layers of immiscible and incompressible fluids: fluid 1 (density $\rho_1$,
viscosity $\mu_1$) is defined between $z=z_0$ and $z=z_1$; fluid 2 (density
$\rho_2$, viscosity $\mu_2$) is defined between $z=z_1$ and $z=z_2$. The
interface between fluids 1 and 2 is $\zeta$ (surface tension coefficient
$\gamma$). The equations of motion are given by

$$
\rho_j \left[\partial_t + (\vec{U}_j\cdot\nabla) \right]\vec{U}_j = -\nabla(P_j) + mu_j\nabla^2\vec{U}_j - \rho_j G(t) \vec{e}_z
$$

$$
\label{2}
\nabla\cdot\vec{U}_j = 0
$$

where $\vec{U}=(u_j,v_j,w_j)$ and $G(t)=g (1 + F\cos(\omega t))$.

We will consider a horizontally infinite plane, whose normal modes are
trigonometric functions, e.g. $\sin (\vec{k}\cdot\vec{x} )$. The horizontal wave
number $\vec{k}=k_x\vec{i} + k_y\vec{j}$, where $k^2 = k_x^2 + k_y^2$, can take
any real value. We can expand the fields in terms of horizontal normal modes of
the Laplacian since the form of the equations is such that each mode is
decoupled from the others. This is the approach followed by Benjamin \& Ursell
[^2] for the ideal fluid case, and it remains valid for the viscous fluid
equations in the present case. We now simply replace 
$$
  \tag{1.2}
  \label{eq2}
  \begin{aligned}
    \label{eq2a}
    w_j^*(\vec{x},z,t) &= \sin(\vec{k}\cdot\vec{x}) w_j (z,t) 
    \\
    \label{eq2b}
    \zeta^*(\vec{x},t) &= \sin(\vec{k}\cdot\vec{x}) \zeta (t) 
  \end{aligned}
$$
and the differential operator $\nabla^2_H$ by the algebraic one $-k^2$

The complete linear stability problem reads
$$
  \tag{1.3}
  \label{eq3}
  \begin{aligned}
    [\partial_t - \nu_1 (\partial_{zz} - k^2)](\partial_zz - k^2) w_1 &= 0, && \text{ for } -h_1 \leq z < 0
    \\
    [\partial_t - \nu_2 (\partial_{zz} - k^2)](\partial_zz - k^2) w_2 &= 0, && \text{ for } 0 \leq z < h_2
  \end{aligned}
$$
The boundary conditions at the two plates are given by
$$
  \tag{1.4}
  \label{eq4}
  \begin{aligned}
    w_1 &= 0, && \text{ for } z = -h_1
    \\
    w_2 &= 0, && \text{ for } z = h_2
    \\
    \partial_z w_1 &= 0, && \text{ for } z = -h_1
    \\
    \partial_z w_2 &= 0, && \text{ for } z = h_2
  \end{aligned}
$$
and the conditions at the interface are
$$
  \tag{1.5}
  \label{eq5}
  \begin{aligned}
    w_1 - w_2 &= 0,
    \\
    \partial w_1 - \partial w_2 &= 0,
    \\
    \nu_1 (\partial_{zz} + k^2) w_1 - \nu_2 (\partial_{zz} + k^2) w_2 &= 0,
    \\
    [\rho_1 \partial_t - \nu_1 (\partial_{zz} - k^2) + 2 \nu_1 k^2] \partial_z w_1
    -[\rho_2 \partial_t - \nu_2 (\partial_{zz} - k^2) + 2 \nu_2 k^2] \partial_z w_2
    &= - [(\rho_1-\rho_2)g(t) - \sigma k^2] k^2 \zeta,
  \end{aligned}
$$
The kinematic condition at the interface reads
$$
  \tag{1.6}
  \label{eq6}
  \begin{aligned}
    \partial \zeta - w  \vert_{z=0} = 0
  \end{aligned}
$$
The above set of equations \eqref{eq3}-\eqref{eq6} constitute the full
hydrodynamic system, which we shall refer to as FHS.

## References

[^1]: Kumar K, Tuckerman LS. Parametric instability of the interface between two
    fluids. J. Fluid Mech. 1994;279:49-68. doi:10.1017/S0022112094003812 
[^2]: Benjamin TB and Ursell FJ. The stability of the plane free surface of a
    liquid in vertical periodic motion. Proc. R. Soc. Lond. 1954. A225505–515.
    doi:10.1098/rspa.1954.0218
