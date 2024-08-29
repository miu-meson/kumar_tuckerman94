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
equations in the present case. We now simply expand  

$$
\begin{aligned}
w_j(\vec{x},z,t) &= \sin(\vec{k}\cdot\vec{x}) w_j (z,t) 
\\
\zeta(\vec{x},t) &= \sin(\vec{k}\cdot\vec{x}) \zeta (t)
\end{aligned}
$$

and the differential operator $\nabla^2_H$ by the algebraic one $-k^2$.

The complete linear stability problem reads

$$
\begin{aligned}
[\partial_t - \nu_1 (\partial_{zz} - k^2)](\partial_zz - k^2) w_1 = 0, &&&\text{ for } &-h_1 \leq z < 0
\\
[\partial_t - \nu_2 (\partial_{zz} - k^2)](\partial_zz - k^2) w_2 = 0, &&&\text{ for } &0 \leq z < h_2
\end{aligned}
$$

The boundary conditions at the two plates are given by

$$
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
\begin{aligned}
w_1 - w_2 &= 0,
\\
\partial w_1 - \partial w_2 &= 0,
\\
\nu_1 (\partial_{zz} + k^2) w_1 - \nu_2 (\partial_{zz} + k^2) w_2 &= 0,
\\
[\rho_1 \partial_t - \nu_1 (\partial_{zz} - k^2) + 2 \nu_1 k^2] \partial_z w_1 - [\rho_2 \partial_t - \nu_2 (\partial_{zz} - k^2) + 2 \nu_2 k^2] \partial_z w_2 &= - [(\rho_1-\rho_2)g(t) - \sigma k^2] k^2 \zeta,
\end{aligned}
$$

The kinematic condition at the interface reads

$$
\partial \zeta - w  \vert_{z=0} = 0
$$

The above set of equations constitute the full hydrodynamic system, which we
shall refer to as FHS.

## Solutions of Floquet form

In standard fashion, we search for solutions of Floquet form, i.e.,
$w_j(t)=w_j(t+nT)$ with integer $n$ and $T=2\pi/\omega$, 

$$
w_j (z,t) = e^{(i\alpha+\lambda)t} \tilde{w}_j (z, t \mod T)
$$

where $i\alpha+\lambda$ is the Floquet exponent and
$e^{(i\alpha + \lambda)T}$ is the Floquet multiplier. The function is
periodic in time with period $T$, and may therefore be expanded in the
Fourier series

$$
  \tilde{w}_j = \sum_{n} w_{jn}(z) e^{in\omega t}
$$

The Floquet multipliers are eigenvalues of a real mapping: this implies that
they are either real or complex-conjugate pairs. In addition, $\alpha$ is
defined only modulo $\omega$, since integer multiples of $\omega$ may be
absorbed into $\tilde{w}_j$. Hence, we restrict consideration to the range $0
\leq \alpha < \omega$. The two cases $\alpha = 0$ and $\alpha = \omega/2$ are
called harmonic and subharmonic, respectively, and correspond to positive or
negative real Floquet multipliers, whereas $0 < \alpha < \omega/2$ corresponds
to a complex Floquet multiplier.

The relationship between Fourier modes with positive and negative $n$ depends on
the value of $\alpha$. In the harmonic and subharmonic cases, $\tilde{w}_j$ must
obey reality conditions $w_{j,-n} = w_{j,n}^*$ (harmonic) or $w_{j,-n} =
w_{j,n-1}^*$ (harmonic) (subharmonic), so that the series may be rewritten in
terms only of non-negative Fourier indices. Only the harmonic and subharmonic
cases are relevant to this linear stability analysis : complex Floquet
multipliers are always of magnitude less than or equal to one, and hence do not
correspond to growing solutions. 

## References

[^1]: Kumar K, Tuckerman LS. Parametric instability of the interface between two
    fluids. J. Fluid Mech. 1994;279:49-68. doi:10.1017/S0022112094003812 
[^2]: Benjamin TB and Ursell FJ. The stability of the plane free surface of a
    liquid in vertical periodic motion. Proc. R. Soc. Lond. 1954. A225505–515.
    doi:10.1098/rspa.1954.0218
