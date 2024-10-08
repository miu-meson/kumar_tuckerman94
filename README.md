# Parametric instability of the interface between two fluids

We briefly describe an implementation of Kumar \& Tuckerman [^1]. We consider
two layers of immiscible and incompressible fluids: fluid 1 (density $\rho_1$,
viscosity $\mu_1$) is defined between $z=z_0$ and $z=z_1$; fluid 2 (density
$\rho_2$, viscosity $\mu_2$) is defined between $z=z_1$ and $z=z_2$. The
interface between fluids 1 and 2 is $\zeta$ (surface tension coefficient
$\gamma$). The equations of motion are given by the incompressible Navier-Stokes
equations. We will consider a horizontally infinite plane, whose normal modes
are trigonometric functions, e.g. $\sin (\vec{k}\cdot\vec{x} )$. The horizontal
wave number $\vec{k}=k_x\vec{i} + k_y\vec{j}$, where $k^2 = k_x^2 + k_y^2$, can
take any real value. We can expand the fields in terms of horizontal normal
modes of the Laplacian since the form of the equations is such that each mode is
decoupled from the others. This is the approach followed by Benjamin \& Ursell
[^2] for the ideal fluid case, and it remains valid for the viscous fluid
equations in the present case. 

An example of figures 1a and 1b from Kumar \& Tuckerman [^1] generated by this 
code can be seen here.

![Stability boundary for ideal fluids for the parameters in [^1].. The tongues
correspond alternately to subharmonic (red) and harmonic (black)
responses.](KT94_fig1a.png)

![Stability boundary for FHS for the parameters in [^1].](KT94_fig1b.png)

## References

[^1]: Kumar K, Tuckerman LS. Parametric instability of the interface between two
    fluids. J. Fluid Mech. 1994;279:49-68. doi:10.1017/S0022112094003812 
[^2]: Benjamin TB and Ursell FJ. The stability of the plane free surface of a
    liquid in vertical periodic motion. Proc. R. Soc. Lond. 1954. A225505–515.
    doi:10.1098/rspa.1954.0218
