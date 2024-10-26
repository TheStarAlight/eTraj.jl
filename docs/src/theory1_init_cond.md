# [Theory: Initial Conditions](@id theory_init_cond)

Several theories on strong-field ionization can be utilized to provide the initial conditions of the classical electrons in the trajectory simulation scheme.
The initial condition consists of three properties:

- Initial position ``\rr_0``, i.e., the tunneling exit position;
- Initial momentum ``\kk_0`` [^note1];
- The corresponding ionization probability ``W`` carried by each electron sample, which depends on the time-dependent laser field and properties of the target atom/molecule.

[^note1]: Here we use ``\kk`` to denote the momentum/velocity of the electron when the laser is on, to distinguish it from the canonical momentum ``\pp=\kk-\AA``. After the laser turned off, we have ``\pp=\kk``.

In this section we briefly revisit the available theories implemented in `eTraj`.
Atomic units (a.u.) are used throughout unless stated otherwise.

```@contents
Pages = ["theory1_init_cond.md"]
Depth = 2
```

---------------------------

## [Strong-Field Approximation with Saddle-Point Approximation (SFA-SPA)](@id SFA)

The *Strong-Field Approximation (SFA)* is originated from the Keldysh theory of strong-field ionization [^Keldysh_1965] [^Faisal_1973] [^Reiss_1980].
Compared with the perturbative methods and adiabatic tunneling theories, the SFA is able to predict both the multi-photon and the tunneling process during the laser-atom interaction, as well as high-order non-perturbative phenomena such as the ATI because it fully includes the non-adiabatic effect of the laser-atom interaction.
The broad scope of SFA has contributed to its widespread application in theoretical investigations of strong-field ionization.

[^Keldysh_1965]: L. V. Keldysh, Ionization in the field of a strong electromagnetic wave, *Sov. Phys. JETP* **20**, 1307 (1965).

[^Faisal_1973]: F. H. M. Faisal, Multiple absorption of laser photons by atoms, *J. Phys. B: At. Mol. Phys.* **6**, L89 (1973). DOI: [10.1088/0022-3700/6/4/011](https://doi.org/10.1088/0022-3700/6/4/011)

[^Reiss_1980]:  H. R. Reiss, Effect of an intense electromagnetic field on a weakly bound system, *Phys. Rev. A* **22**, 1786 (1980). DOI: [10.1103/PhysRevA.22.1786](https://doi.org/10.1103/PhysRevA.22.1786)

Consider the electron evolving under a combined field of the Coulomb field ``V(\rr)`` of the parent ion and the laser field ``\FF(t)=-\pd_t \AA(t)``.
Under the length gauge (LG), its Hamiltonian reads
```math
\begin{equation}
    H^{\rm{LG}} = \frac12 \pp^2 + V(\rr) + \FF(t) \cdot \rr.
\end{equation}
```

Denoting ``\ket{\Psi_0}=\ket{\psi_0} \ee^{\ii \Ip t}`` as the unperturbed initial state with ionization potential of ``\Ip``,
``\ket{\Psi_{\pp}}`` as the continuum state of momentum ``\pp``, and
```math
\begin{equation}
    U(\tf,t_0) = \exp \left[ -\ii\int_{t_0}^{\tf} H^{\rm{LG}}(\tau) \dd\tau \right]
\end{equation}
```

the time-evolution operator, the transition amplitude between the initial state (at ``t_0``) and the final state of momentum ``\pp`` (at ``\tf``) is written as
```math
\begin{equation}
    M_{\pp} = \mel{\Psi_{\pp}}{U(\tf,t_0)}{\Psi_0}.
\end{equation}
```

Here lies the key idea of SFA: when the influence of the Coulomb field to the ionized electrons is weak compared with that of the external laser field, we may neglect the influence of the Coulomb field in the expression of ``M_{\pp}`` by replacing the time-evolution operator with a Coulomb-free one ``U_{\rm{f}}``, and meanwhile replacing the continuum state with the Volkov state ``\ket{\Psi_{\pp}^{\rm{V}}}`` which represents a free electron evolving under the same laser field:
```math
\begin{equation}
    M_{\pp} \approx \mel{\Psi_{\pp}^{\rm{V}}}{U_{\rm{f}}(\tf,t_0)}{\Psi_0},
\end{equation}
```
where the Volkov state under the LG is the product of a plane wave and a phase factor:
```math
\begin{equation}
    \ket{\Psi_{\pp}^{\rm{V}}} = \ket{\pp+\AA(t)} \exp \left[-\ii \int^t \frac12 [\pp+\AA(\tau)]^2 \dd\tau \right].
\end{equation}
```
In this way the ``M_{\pp}`` is expressed as
```math
\begin{equation}
    M_{\pp} = -\ii \int_{t_0}^{\tf} \mel{\pp+\AA(\tau)}{\FF(\tau)\cdot\rr}{\psi_0} \ee^{-\ii\Sp(\tau)} \dd\tau,
\end{equation}
```
and we note that here we have extracted the phase factor of ``\ket{\Psi_0}`` and combined it with that of the Volkov state ``\ket{\Psi_{\pp}^{\rm{V}}}``, giving the phase
```math
\begin{equation}
    \Sp(t) = -\int^t \left[ \frac12 [\pp+\AA(\tau)]^2 + \Ip \right] \dd\tau.
\end{equation}
```
Inserting
```math
\begin{equation}
    \frac{\pd}{\pd t} \bra{\pp+\AA(t)} = \ii \bra{\pp+\AA(t)} [\FF(t)\cdot\rr]
\end{equation}
```
into the above expression of ``M_{\pp}`` [Eq. (6)], after integration by parts, one obtains
```math
\begin{equation}
\begin{aligned}
    M_{\pp}
    &= -\int_{t_0}^{\tf} \dd\tau \frac{\pd}{\pd\tau}[\bk{\pp+\AA(\tau)}{\psi_0}] \ee^{-\ii\Sp(\tau)} \\
    &= -{\bk{\pp+\AA(\tau)}{\psi_0} \ee^{-\ii\Sp(\tau)}} \rvert_{t_0}^{\tf} + \int_{t_0}^{\tf} \dd\tau \bk{\pp+\AA(\tau)}{\psi_0} \cdot (-\ii\Sp'(\tau)) \ee^{-\ii\Sp(\tau)}.
\end{aligned}
\end{equation}
```

An additional saddle-point approximation (SPA) facilitates preparation of initial conditions of the electron trajectories.
The variation of phase factor ``\ee^{\ii\Sp(t)}`` is much more sensitive than that of the prefactor as ``t`` varies, which leads to a fact that the whole integrand in Eq. (9) oscillates in its complex phase and its values cancel out in most cases, except when the variation of the phase ``\Sp(t)`` becomes stable, i.e., at the saddle points.
The saddle points ``\ts=\tr+\ii\ti`` are the zeroes of the derivative of the complex function ``\Sp(t)``, which satisfy
```math
\begin{equation}
    -\Sp'(\ts) = \frac12 [\pp+\AA(\ts)]^2 + \Ip = 0.
\end{equation}
```
The second term of the r.h.s. of Eq. (9), i.e., the integral, has significant contribution only in the vicinities of the two end points ``t_0, \tf`` and the saddle points ``\ts``,
while the contribution near the two end points cancels out the first term.
Therefore, the ``M_{\pp}`` is now approximated with the integration around the saddle points:
```math
\begin{equation}
    M_{\pp} \approx \sum_{\ts} \int_{C_{\ts}} \dd\tau \bk{\pp+\AA(\tau)}{\psi_0} \cdot (-\ii\Sp'(\tau)) \ee^{-\ii\Sp(\tau)},
\end{equation}
```
with ``C_{\ts}`` the integration contour following the steepest-descent path related to ``\ts``.

Further evaluation of the prefactor ``\tilde{\psi}_0(\kk) \rvert_{\kk=\pp+\AA(t)} = \braket{\pp+\AA(t)}{\psi_0}`` (i.e., the momentum-space wavefunction) in the vicinity of the saddle points in Eq. (11) is essential before applying the SPA.
We assume the field points towards the ``+ z`` axis, for an atom target at ``(l,m)`` state with ionization potential ``\Ip``, its wavefunction behaves asymptotically as
```math
\begin{equation}
    \psi_0(\rr) \sim 2 C_{\kappa l} \kappa^{3/2} (\kappa r)^{n^*-1} \ee^{-\kappa r} Y_{lm}(\hat{\rr})
\end{equation}
```
[^Perelomov_1966] for ``\kappa r \gg 1``, with ``\kappa=\sqrt{2\Ip}``, ``n^*=Z/\kappa`` the effective principal quantum number, ``Z`` the charge of the residual ion, ``Y_{lm}`` the spherical harmonics, and ``C_{\kappa l}`` the asymptotic coefficient for atoms, which can be approximated using the Hartree approximation formula [^Hartree_1928]
```math
\begin{equation}
    C_{\kappa l}^2 = \frac{2^{2n^*-2}}{n^* (n^*+l)! (n^*-l-1)!}.
\end{equation}
```
For atomic hydrogen at ground state we have ``C_{\kappa l} = 1``.
Moreover, for non-integer ``n^*``, the formula can be naturally extended by replacing the factorials ``x!`` with Gamma functions ``\Gamma(x+1)``, i.e.,
```math
\begin{equation}
    C_{\kappa l}^2 = \frac{2^{2n^*-2}}{n^* \Gamma(n^*+l+1) \Gamma(n^*-l)}.
\end{equation}
```
In the vicinity of the saddle points, which corresponds to the case when ``k^2 \rightarrow -\kappa^2``, the expression of ``\tilde{\psi}_0(\kk)`` is determined by the asymptotic behavior of the wavefunction:
```math
\begin{equation}
    \tilde{\psi}_0(\kk) = \frac{C_{\kappa l}}{\sqrt{\pi}} \frac{2^{n^*+3/2}\kappa^{2n^*+1/2}\Gamma(n^*+1)}{(k^2+\kappa^2)^{n^*+1}} Y_{lm}(\hat{\kk}),
\end{equation}
```
where ``\Gamma`` is the gamma function [^note2].
Substituting the above expression into Eq. (11), making use of the definition of ``\Sp(t)`` [Eq. (7)], we obtain
```math
\begin{equation}
    M_{\pp} = \ii \frac{C_{\kappa l}}{\sqrt{\pi}} 2^{1/2}\kappa^{2n^*+1/2}\Gamma(n^*+1) \sum_{\ts} \int_{C_{\ts}} \frac{Y_{lm}(\hat{\kk}(\tau))}{[\Sp'(\tau)]^{n^*}} \ee^{\ii\Sp(\tau)} \dd\tau,
\end{equation}
```
where ``\hat{\kk}(\tau)`` is the complex unit vector along ``\kk(\tau)=\pp+\AA(\tau)``,
and the evaluation method of spherical harmonics with complex arguments is based on Appendix B of Ref. [^Pisanty_2017], see also note [^note3].
A modified version of SPA can be carried out to handle the case when the integrand has a singularity at ``\ts`` (see Appendix B of Ref. [^Gribakin_1997]):
```math
\begin{equation}
\begin{aligned}
    \int_{C_{\ts}} \frac{Y_{lm}(\hat{\kk}(\tau))}{[\Sp'(\tau)]^{n^*}} \ee^{\ii\Sp(\tau)} \dd\tau
    &\approx \frac{Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{n^*}} \int_{C_{\ts}} \frac{\ee^{\ii\Sp(\tau)}}{(\tau-\ts)^{n^*}} \dd\tau \\
    &\approx \frac{Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{n^*}} \cdot \ii^{n^*} \frac{\Gamma(n^*/2)}{2\Gamma(n^*)} \sqrt{\frac{2\pi}{-\ii\Sp''(\ts)}} [-2\ii\Sp''(\ts)]^{n^*/2} \ee^{\ii\Sp(\ts)}.
\end{aligned}
\end{equation}
```
In this way we find the expression of transition amplitude:
```math
\begin{equation}
    M_{\pp} = c_{n^*} C_{\kappa l} \sum_{\ts} \frac{Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{(n^*+1)/2}} \ee^{-\ii\Sp(\ts)},
\end{equation}
```
with ``c_{n^*} = \ii^{(n^*-5)/2} 2^{n^*/2+1} \kappa^{2n^*+1/2} \Gamma(n^*/2+1)`` the constant coefficient.

[^note2]: Here ``k`` is actually ``\sqrt{\kk\cdot\kk}`` and is not the conventional "norm" of the complex vector, which is ``\abs{\abs{\kk}}=\sqrt{\kk^*\cdot\kk}``. Here we normalize the complex vector ``\kk`` through ``\hat{\kk}=\kk/\sqrt{\kk\cdot\kk}``, according to Ref. [^Perelomov_1966].

[^note3]: The evaluation of spherical harmonic in the complex domain makes use of the fact that the solid harmonics ``S_{lm}(\rr) := r^l Y_{lm}(\hat{\rr})`` can be expressed in a polynomial of ``x,y,z``, see Ref. [^Caola_1978].

[^Perelomov_1966]: A. Perelomov, V. Popov, and M. Terent’ev, Ionization of atoms in an alternating electric field, *Sov. Phys. JETP* **23**, 924 (1966).

[^Hartree_1928]: D. R. Hartree, The wave mechanics of an atom with a non-coulomb central field. Part I. Theory and methods, *Math. Proc. Cambridge Philos. Soc.* **24**, 89 (1928). DOI: [10.1017/S0305004100011919](https://doi.org/10.1017/S0305004100011919)

[^Pisanty_2017]: E. Pisanty and Á. Jiménez-Galán, Strong-field approximation in a rotating frame: High-order harmonic emission from p states in bicircular fields, *Phys. Rev. A* **96**, 063401 (2017). DOI: [10.1103/PhysRevA.96.063401](https://doi.org/10.1103/PhysRevA.96.063401)

[^Caola_1978]:  M. J. Caola, Solid harmonics and their addition theorems, *J. Phys. A: Math. Gen.* **11**, L23 (1978). DOI: [10.1088/0305-4470/11/2/001](https://doi.org/10.1088/0305-4470/11/2/001)

[^Gribakin_1997]: G. F. Gribakin and M. Yu. Kuchiev, Multiphoton detachment of electrons from negative ions, *Phys. Rev. A* **55**, 3760 (1997). DOI: [10.1103/PhysRevA.55.3760](https://doi.org/10.1103/PhysRevA.55.3760)

The SFA phase ``\Sp(\ts)`` is obtained by solving the integral
```math
\begin{equation}
\begin{aligned}
    \Sp(\ts)
    &= -\int_{\ts}^{\infty} \dd\tau \lbrack \frac12 [\pp+\AA(\tau)]^2 + \Ip \rbrack \\
    &= \left( - \int_{\ts}^{\tr} - \int_{\tr}^{\infty} \right) \dd\tau \left[ \frac12 [\pp+\AA(\tau)]^2 + \Ip \right] \\
    &= S_{\pp,\rm{tun}} + S_{\pp,\rm{traj}},
\end{aligned}
\end{equation}
```
where the ``S_{\pp,\rm{tun}}`` and ``S_{\pp,\rm{traj}}`` represent the complex phase accumulated during the tunneling process and the trajectory motion in the continuum, respectively.
The ``S_{\pp,\rm{tun}}`` is accumulated during an imaginary period of time (from time ``\ts`` to ``\tr``), in which the electron passes through the potential barrier with an "imaginary" momentum,
its real part denotes the quantum phase, while its imaginary part is related to the ionization probability.

To utilize the SFA to give initial condition of the photoelectrons, we suppose that the classical electron is ejected at time ``\tr`` at tunneling exit ``\rr_0`` with momentum ``\kk_0=\kk(\tr)``.
The initial momentum ``\kk_0``, neglecting the Coulomb interaction with the nucleus, is related to the final momentum ``\pp`` through
```math
\begin{equation}
    \pp = \kk_0 - \int_{\tr}^{\infty} \FF(\tau) \dd\tau = \kk_0 - \AA(\tr).
\end{equation}
```

The initial position ``\rr_0``, i.e., the tunneling exit, is found by constructing a quantum tunneling trajectory.
The beginning of the trajectory, i.e., the tunneling entrance, has a real part of zero; the electron tunnels through the barrier during the time interval ``\ts`` to ``\tr`` and emerges as a classical electron at the tunneling exit ``\rr_0`` with real position and momentum.
In this way we obtain the expression of the initial position:
```math
\begin{equation}
    \rr_0^{\rm{SFA-SPA}} = \Re \int_{\ts}^{\tr} [\pp+\AA(\tau)] \dd\tau = \Im \int_0^{\ti} [\pp+\AA(\tr+\ii\tau)] \dd\tau = \Im \int_0^{\ti} \AA(\tr+\ii\tau) \dd\tau,
\end{equation}
```
here ``\Re`` and ``\Im`` are the real and imaginary part notation, respectively,
and the real ``\pp`` does not come into play in the integral.

The probability density (in the final momentum space ``\pp``) carried by the electron sample is
```math
\begin{equation}
    \dd W^{\rm{SFA-SPA}}/\dd \pp = \sum_{\ts} \abs{\mathcal{P}^{\rm{SFA-SPA}}_{\pp}(\ts)}^2 \exp(-2\Im S_{\pp,\rm{tun}}(\ts)),
\end{equation}
```
where we have gathered the coefficients to the prefactor
```math
\begin{equation}
\begin{aligned}
    \mathcal{P}^{\rm{SFA-SPA}}_{\pp}(\ts)
    = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{(n^*+1)/2}}
    = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{\left\{[\pp+\AA(\ts)]\cdot\FF(\ts)\right\}^{(n^*+1)/2}}.
\end{aligned}
\end{equation}
```

We note that the ionization probability in Eq. (22) is expressed in the coordinate of final momentum ``\pp=(p_x,p_y,p_z)``,
however, in the trajectory simulation, the initial electrons are sampled in the ``(\kkt,\tr)`` coordinate, with ``\kkt`` the initial transversal momentum.
Thus, adding a Jacobian in the prefix of the ionization probability is required if we sample the initial electrons within such coordinate.
Suppose the laser propagates in the ``z`` axis and polarizes in the ``xy`` plane, the transformed expression reads
```math
\begin{equation}
    \dd W^{\rm{SFA-SPA}}/\dd\kkt\dd\tr = \sum_{\ts} J(\kp,\tr) \abs{\mathcal{P}_{\pp}(\ts)}^2 \exp(-2\Im S_{\pp,\rm{tun}}(\ts)),
\end{equation}
```
where ``\kp`` is the projection of ``\kkt`` on the polarization plane (i.e., the ``xy`` plane), and the Jacobian is
```math
\begin{equation}
    J(\kp,\tr) = \abs{\frac{\pd(p_x,p_y)}{\pd(\kp,\tr)}} =
    \begin{vmatrix}
        \pd p_x/\pd\kp & \pd p_x/\pd\tr \\
        \pd p_y/\pd\kp & \pd p_y/\pd\tr \\
    \end{vmatrix}.
\end{equation}
```

---------------------------

## [SFA-SPA with Non-adiabatic Expansion (SFA-SPANE)](@id SFAAE)

For small Keldysh parameter ``\gamma=\omega\kappa/F_0`` (``\omega`` is the laser angular frequency and ``F_0`` is the peak field strength), the non-adiabatic effect is not significant, thus a non-adiabatic expansion scheme can be carried out to develop a modified theory based on the SFA-SPA, which is named after the *SFA-SPA with Non-adiabatic Expansion (SFA-SPANE)* [^Ni_2018] [^Mao_2022] [^Ma_2021] [^Ma_2024].
It includes the non-adiabatic effect to a large extent and is capable of giving similar results compared with that given by the SFA-SPA under small Keldysh parameters.

[^Ni_2018]: H. Ni, N. Eicke, C. Ruiz, J. Cai, F. Oppermann, N. I. Shvetsov-Shilovski, and L.-W. Pi, Tunneling criteria and a nonadiabatic term for strong-field ionization, *Phys. Rev. A* **98**, 013411 (2018). DOI: [10.1103/PhysRevA.98.013411](https://doi.org/10.1103/PhysRevA.98.013411)

[^Mao_2022]: X. Mao, H. Ni, X. Gong, J. Burgdörfer, and J. Wu, Subcycle-resolved strong-field tunneling ionization: Identification of magnetic dipole and electric quadrupole effects, *Phys. Rev. A* **106**, 063105 (2022). DOI: [10.1103/PhysRevA.106.063105](https://doi.org/10.1103/PhysRevA.106.063105)

[^Ma_2021]: Y. Ma, J. Zhou, P. Lu, H. Ni, and J. Wu, Influence of nonadiabatic, nondipole and quantum effects on the attoclock signal, J. Phys. B: At. Mol. Opt. Phys. 54, 144001 (2021). DOI: [10.1088/1361-6455/ac0d3e](https://doi.org/10.1088/1361-6455/ac0d3e)

[^Ma_2024]: Y. Ma, H. Ni, and J. Wu, Attosecond ionization time delays in strong-field physics, Chin. Phys. B 33, 13201 (2024). DOI: [10.1088/1674-1056/ad0e5d](https://doi.org/10.1088/1674-1056/ad0e5d)


The SFA-SPANE is applicable when the Keldysh parameter is small, and the non-adiabatic effect is insignificant, which corresponds to the small-``\ti`` case.
We expand the vector potential ``\AA(\ts)=\AA(\tr+\ii\ti)`` in the SFA-SPA around ``\ti=0``, up to the second order of ``\ti``:
```math
\begin{equation}
    \AA(\tr+\ii\ti) = \AA(\tr) - \ii\ti\FF(\tr) + \frac12 \ti^2 \FF'(\tr) + o(\ti^2).
\end{equation}
```
Inserting Eq. (26) into the saddle-point equation in the SFA-SPA [Eq. (10)] gives
```math
\begin{equation}
    \kk(\tr)\cdot\FF(\tr) \approx 0
\end{equation}
```
and
```math
\begin{equation}
    \ti \approx \sqrt{\frac{k^2(\tr)+\kappa^2}{F^2(\tr)-\kk(\tr)\cdot\FF'(\tr)}},
\end{equation}
```
which allow for more specific expression of the ionization probability and other quantities.

The initial position ``\rr_0`` in SFA-SPANE, is given by
```math
\begin{equation}
    \rr_0^{\rm{SFA-SPANE}} = \Im \int_0^{\ti} \AA(\tr+\ii\tau) \dd\tau = - \frac{\FF}{2} \frac{\kt^2+\kappa^2}{F^2-\kk_0\cdot\FF'}.
\end{equation}
```

The ``\Im S_{\pp,\rm{tun}}`` which is related to the ionization probability, in the SFA-SPANE, is
```math
\begin{equation}
\begin{aligned}
    \Im S_{\pp,\rm{tun}}
    &\approx \Im \!\! \int_{\tr}^{\ts} \!\!\!\! \dd\tau \left\{ \frac12 \left[ \pp + \AA(\tr) - \ii\ti\FF(\tr) + \frac12 \ti^2 \FF'(\tr) \right]^2 \!\!\! + \Ip \right\} \\
    &\approx \left[ \Ip+\frac12 k^2(\tr) \right]\ti - [F^2(\tr)-\kk(\tr)\cdot\FF'(\tr)]\frac{\ti^3}{6} \\
    &= \frac13 \frac{(k^2+\kappa^2)^{3/2}}{\sqrt{F^2-\kk_0\cdot\FF'}}.
\end{aligned}
\end{equation}
```
Then follows the ionization probability
```math
\begin{equation}
    \dd W^{\rm{SFA-SPANE}}/\dd\pp = \abs{\mathcal{P}^{\rm{SFA-SPANE}}_{\pp}(\ts)}^2 \exp \left[ -\frac23 \frac{(\kt^2+\kappa^2)^{3/2}}{\sqrt{F^2-\kk_0\cdot\FF'}} \right],
\end{equation}
```
where ``\kkt`` is actually equivalent to ``\kk(\tr)`` in the SFA-SPANE because of the zero longitudinal initial momentum condition in Eq. (27), which is derived from the saddle-point equation under adiabatic expansion.
The prefactor, reads
```math
\begin{equation}
    \mathcal{P}^{\rm{SFA-SPANE}}_{\pp}(\ts) = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{\left[ (\kt^2+\kappa^2)(F^2-\kk_0\cdot\FF') \right]^{(n^*+1)/4}}.
\end{equation}
```

---------------------------

## [Ammosov-Delone-Krainov (ADK)](@id ADK)

The *Ammosov-Delone-Krainov (ADK)* theory [^Ammosov_1986] [^Delone_1998] is used to study the adiabatic tunneling in the strong-field ionization, and is, in a sense, the adiabatic limit of the SFA.

[^Ammosov_1986]: M. Ammosov, N. Delone, and V. Krainov, Tunnel ionization of complex atoms and of atomic ions in an alternating electromagnetic field, *Sov. Phys. JETP* **64**, 1191 (1986).

[^Delone_1998]: N. B. Delone and V. P. Krainov, Tunneling and barrier-suppression ionization of atoms and ions in a laser radiation field, *Phys. Usp.* **41**, 469 (1998). DOI: [10.1070/PU1998v041n05ABEH000393](https://doi.org/10.1070/PU1998v041n05ABEH000393)


In the adiabatic limit, the laser field can be treated as static, thus we have ``\FF'(t)=0`` (higher order derivatives of ``\FF(t)`` remains zero as well).
Substituting it into the ionization probability of SFA-SPANE [Eq. (31)] gives
```math
\begin{equation}
    \dd W^{\rm{ADK}}/\dd\pp = \abs{\mathcal{P}^{\rm{ADK}}_{\pp}(\ts)}^2 \exp \left[ -\frac23 \frac{(\kt^2+\kappa^2)^{3/2}}{F} \right],
\end{equation}
```
where the prefactor reads
```math
\begin{equation}
    \mathcal{P}^{\rm{ADK}}_{\pp}(\ts) = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{\left[ (\kt^2+\kappa^2)F^2 \right]^{(n^*+1)/4}},
\end{equation}
```
with ``\ti=\sqrt{\kt^2+\kappa^2}/F``.
If we expand Eq. (33) under the small-``\kt`` limit, we obtain
```math
\begin{equation}
    \dd W^{\rm{ADK}}/\dd\pp \propto \exp \left( -\frac{2\kappa^3}{3F} \right) \exp \left( -\frac{\kappa\kt^2}{F} \right),
\end{equation}
```
which is actually the exponential term of the well-known ADK rate.
However, we note that the result of our approach, i.e., applying the adiabatic limit of the SFA-SPA, is slightly different from the actual ADK rate in the prefactor.
As a remedy, adding an additional Coulomb-correction (CC) factor to Eq. (33) fills the gap between them:
```math
\begin{equation}
    C^{\rm{CC}} = \left(\frac{2\kappa^3}{F}\right)^{n^*} \!\!\!\! \left(1+2\gamma/e\right)^{-2n^*} \left[\Gamma\left(\frac{n^*}{2}+1\right)\right]^{-2}.
\end{equation}
```

The tunneling exit is found with the same approach:
```math
\begin{equation}
    \rr_0^{\rm{ADK}} = \Im \int_0^{\ti} \AA(\tr+\ii\tau) \dd\tau = - \frac{\FF}{2} \frac{\kt^2+\kappa^2}{F^2},
\end{equation}
```
which we refer to as the "``\Ip/F``" model, but is slightly different from the original version because we replaced the ionization potential ``\Ip=\kappa^2/2`` with the effective one ``\Ip'=(\kappa^2+\kt^2)/2`` to account for the initial kinetic energy, which suits the adiabatic tunneling scenario better.


---------------------------

## [Molecular SFA-SPA/SFA-SPANE/ADK](@id MOSFA)

The atomic SFA theory and its adiabatic versions mentioned in the previous sections can be generalized naturally to molecular cases [^Muth_2000] [^Tong_2002] [^Kjeldsen_2004] [^Kjeldsen_2005].
Under the Born-Oppenheimer and the single-active-electron (SAE) approximation, the strong-field ionization of the molecules can be modeled as the interaction of the ionizing orbital (usually the highest occupied molecular orbital (HOMO)) ``\psi_0(\rr)`` with the effective potential of the parent ion and the laser field, which simplifies the problem.

[^Muth_2000]: J. Muth-Böhm, A. Becker, and F. H. M. Faisal, Suppressed molecular ionization for a class of diatomics in intense femtosecond laser fields, *Phys. Rev. Lett.* **85**, 2280 (2000). DOI: [10.1103/PhysRevLett.85.2280](https://doi.org/10.1103/PhysRevLett.85.2280)

[^Tong_2002]: X. M. Tong, Z. X. Zhao, and C. D. Lin, Theory of molecular tunneling ionization, *Phys. Rev. A* **66**, 33402 (2002).DOI: [10.1103/PhysRevA.66.033402](https://doi.org/10.1103/PhysRevA.66.033402)

[^Kjeldsen_2004]: T. K. Kjeldsen and L. B. Madsen, Strong-field ionization of N₂ : Length and velocity gauge strong-field approximation and tunnelling theory, *J. Phys. B: At. Mol. Opt. Phys.* **37**, 2033 (2004). DOI: [10.1088/0953-4075/37/10/003](https://doi.org/10.1088/0953-4075/37/10/003)

[^Kjeldsen_2005]: T. K. Kjeldsen, C. Z. Bisgaard, L. B. Madsen, and H. Stapelfeldt, Influence of molecular symmetry on strong-field ionization: Studies on ethylene, benzene, fluorobenzene, and chlorofluorobenzene, *Phys. Rev. A* **71**, 13418 (2005). DOI: [10.1103/PhysRevA.71.042508](https://doi.org/10.1103/PhysRevA.71.042508)

To generalize the atomic SFA to the molecular SFA, we start from the transition amplitude in Eq. (11).
In the molecular frame (MF), the asymptotic wavefunction can be expanded into spherical harmonics:
```math
\begin{equation}
    \psi_0^{\rm{MF}}(\rr) \sim \sum_{l,m} 2 C_{lm} \kappa^{3/2} (\kappa r)^{n^*-1} \ee^{-\kappa r} Y_{lm}(\hat{\rr}),
\end{equation}
```
where the ``C_{lm}`` are asymptotic coefficients, and we continue to adopt the ``n^*=Z/\kappa`` notation for simplicity, although it does not represent the effective principal quantum number anymore.
We assume that in the field frame (FF) the field ``\FF`` points towards the ``z`` axis, and the rotation ``\RRh`` from the FF to the MF can be defined via a set of Euler angles ``(\phi,\theta,\chi)`` within the ``z-y'-z''`` convention, which satisfies
```math
\begin{equation}
    \psi_0^{\rm{MF}}(\RRh\rr) = \psi_0^{\rm{FF}}(\rr).
\end{equation}
```
Utilizing the Wigner-``D`` matrix, the rotated spherical harmonic function can be expressed as a linear combination of spherical harmonics of the same order ``l``:
```math
\begin{equation}
    \RRh(\phi,\theta,\chi) Y_{lm} = \sum_{m'} D_{m'm}^l(\phi,\theta,\chi) Y_{lm'}
\end{equation}
```
and the asymptotic behavior of the wavefunction in the FF is found by inserting Eq. (40) into Eq. (38), which gives
```math
\begin{equation}
     \psi_0^{\rm{FF}}(\rr) \sim \sum_{l,m,m'} 2 C_{lm} D_{m'm}^l(\phi,\theta,\chi) \kappa^{3/2} (\kappa r)^{n^*-1} \ee^{-\kappa r} Y_{lm'}(\hat{\rr}).
\end{equation}
```

It is obvious that the molecular version of the theory differs from the atomic one only in the expression of prefactor ``\mathcal{P}_{\pp}(\ts)``, while the expression of the tunneling exit and the initial momentum are identical.
Following the same procedure in the previous sections, we obtain the prefactor ``\mathcal{P}_{\pp}`` that is applicable for molecules:
```math
\begin{equation}
    \mathcal{P}^{\rm{SFA-SPA}}_{\pp}(\ts) = c_{n^*} \frac{\sum_{l,m,m'} C_{lm} D_{m'm}^l(\phi,\theta,\chi) Y_{lm'}(\hat{\kk}(\ts))}{\left\{[\pp+\AA(\ts)]\cdot\FF(\ts)\right\}^{(n^*+1)/2}},
\end{equation}
```
```math
\begin{equation}
    \mathcal{P}^{\rm{SFA-SPANE}}_{\pp}(\ts) = c_{n^*} \frac{\sum_{l,m,m'} C_{lm} D_{m'm}^l(\phi,\theta,\chi) Y_{lm'}(\hat{\kk}(\ts))}{\left[ (\kt^2+\kappa^2)(F^2-\kk_0\cdot\FF') \right]^{(n^*+1)/4}},
\end{equation}
```
```math
\begin{equation}
    \mathcal{P}^{\rm{ADK}}_{\pp}(\ts) = c_{n^*} \frac{\sum_{l,m,m'} C_{lm} D_{m'm}^l(\phi,\theta,\chi) Y_{lm'}(\hat{\kk}(\ts))}{\left[ (\kt^2+\kappa^2)F^2 \right]^{(n^*+1)/4}}.
\end{equation}
```
We also note that after applying an additional Coulomb-correction factor [Eq. (36)], the ionization rate is aligned with the original MO-ADK theory [^Tong_2002] under the adiabatic and small-``\kt`` limit.

---------------------------

## [Weak-Field Asymptotic Theory (WFAT)](@id WFAT)

The *Weak-Field Asymptotic Theory (WFAT)* generalizes the tunneling ionization from isotropic atomic potentials to arbitrary molecular potentials [^tolstikhin_2011] [^madsen_2013] [^madsen_2017] [^dnestryan_2018].
Compared with the MO-ADK theory, the WFAT naturally accounts for the influence of the molecule's permanent dipole moment, and calculates the structure factors (a similar concept to the asymptotic coefficients ``C_{lm}`` in Eq. (38)) based on the wavefunction close to the core, rather than using the wavefunction in the asymptotic region, allowing for enhanced accuracy in numerical simulations.

[^tolstikhin_2011]: O. I. Tolstikhin, T. Morishita, and L. B. Madsen, Theory of tunneling ionization of molecules: Weak-field asymptotics including dipole effects, *Phys. Rev. A* **84**, 053423 (2011). DOI: [10.1103/PhysRevA.84.053423](https://doi.org/10.1103/PhysRevA.84.053423)

[^madsen_2013]: L. B. Madsen, F. Jensen, O. I. Tolstikhin, and T. Morishita, Structure factors for tunneling ionization rates of molecules, *Phys. Rev. A* **87**, 013406 (2013). DOI: [10.1103/PhysRevA.87.014062](https://doi.org/10.1103/PhysRevA.87.014062)

[^madsen_2017]: L. B. Madsen, F. Jensen, A. I. Dnestryan, and O. I. Tolstikhin, Structure factors for tunneling ionization rates of molecules: General hartree-fock-based integral representation, *Phys. Rev. A* **96**, 013423 (2017). DOI: [10.1103/PhysRevA.96.013423](https://doi.org/10.1103/PhysRevA.96.013423)

[^dnestryan_2018]: A. I. Dnestryan, O. I. Tolstikhin, L. B. Madsen, and F. Jensen, Structure factors for tunneling ionization rates of molecules: General grid-based methodology and convergence studies, *J. Chem. Phys.* **149**, 164107 (2018). DOI: [10.1063/1.5046902](https://doi.org/10.1063/1.5046902)


The formulation of the WFAT is based on the expansion in the parabolic coordinates.
The total ionization rate ``w``, is split into different parabolic channels:
```math
\begin{equation}
    w^{\rm{WFAT}} = \sum_\nu w_\nu,
\end{equation}
```
where ``w_\nu`` are partial rates of parabolic quantum number indices ``\nu=(n_\xi,m)`` with ``n_\xi=0,1,2,\cdots`` and ``m=0,\pm 1,\pm 2,\cdots``.
In the leading-order approximation of the WFAT, the partial rates can be separated into two factors, namely the structural part ``\abs{G_\nu(\theta,\chi)}^2`` and the field part ``\mathcal{W}_\nu(F)``:
```math
\begin{equation}
    w_\nu = \abs{G_\nu(\theta,\chi)}^2 \mathcal{W}_\nu (F).
\end{equation}
```
The field factor is expressed as
```math
\begin{equation}
    \mathcal{W}_\nu (F) = \frac{\kappa}{2} \left(\frac{4\kappa^2}{F}\right)^{2n^*-2n_\xi-\abs{m}-1} \ee^{-2\kappa^3/3F}.
\end{equation}
```
The structure factor ``G_\nu(\theta,\chi)`` is found by an integral related to the ionizing orbital and a reference function, which has significant contribution only in the vicinities of the nuclei and is insensitive to the wavefunction's asymptotic behavior:
```math
\begin{equation}
    G_\nu(\theta,\chi) = \ee^{-\kappa\mu_F} \int \dd\rr \ \Omega_\nu^*(\RRh^{-1}\rr) \hat{V}_{\rm{c}}(\rr) \psi_0(\rr),
\end{equation}
```
which is evaluated in the MF, with ``\psi_0(\rr)`` the wavefunction of the ionizing orbital;
```math
\begin{equation}
    \bm{\mu} = - \int\dd\rr \ \psi_0^*(\rr) \rr \psi_0(\rr)
\end{equation}
```
is the orbital dipole moment in the MF, with ``\mu_F`` being its component along the field direction;
```math
\begin{equation}
    \Omega_\nu(\rr) = \sum_{l=\abs{m}}^{\infty} \Omega_{lm}^{\nu}(\rr) = \sum_{l=\abs{m}}^{\infty} R_l^\nu(r) Y_{lm}(\hat{\rr})
\end{equation}
```
is a reference function which can be expanded into spherical harmonics, with its radial part expressed as
```math
\begin{equation}
    R_l^\nu(r) = \omega_l^\nu\ (\kappa r)^l\ \ee^{-\kappa r}\ \rm{M}(l+1-n^*,2l+2,2\kappa r),
\end{equation}
```
where ``\rm{M}(a,b,x)`` is the confluent hyper-geometric function and
```math
\begin{equation}
\begin{aligned}
    \omega_l^\nu = & \quad (-1)^{l+(\abs{m}-m)/2+1}\ 2^{l+3/2}\ \kappa^{n^*-(\abs{m}+1)/2-n_\xi}\\
                   & \times \sqrt{(2l+1)(l+m)!(l-m)!(\abs{m}+n_\xi)!n_\xi!}\ \frac{l!}{(2l+1)!}\\
                   & \times \!\!\!\!\!\! \sum_{k=0}^{\min{(n_\xi,l-\abs{m})}} \!\!\!\!\!\!\!\!\!\! \frac{\Gamma(l+1-n^*+n_\xi-k)}{k!(l-k)!(\abs{m}+k)!(l-\abs{m}-k)!(n_\xi-k)!}
\end{aligned}
\end{equation}
```
is the normalization coefficient;
``\hat{V}_\rm{c}(\rr)=\hat{V}(\rr)+Z/r`` is the core potential with the Coulomb tail removed, where ``Z`` is the asymptotic charge of the residual ion.

The effective potential ``\hat{V}(\rr)`` describes the interaction between the ionizing electron and the residual parent ion. We note that here we use the hat notation to indicate that the potential operator is not diagonal in the coordinate space.
Under the framework of the Hartree-Fock method, the effective potential consists of three parts, namely the nuclear Coulomb potential (``V_{\rm{nuc}}``), the direct (``V_{\rm{d}}``) and exchange (``V_{\rm{ex}}``) parts of inter-electron interactions:
```math
\begin{equation}
    \hat{V} = V_{\rm{nuc}} + V_{\rm{d}} + \hat{V}_{\rm{ex}},
\end{equation}
```
with
```math
\begin{equation}
\begin{aligned}
    V_{\rm{nuc}}(\rr) &= - \sum_{A=1}^{N_\rm{atm}} \frac{Z_A}{\abs{\rr-\bm{R}_A}},\\
    V_{\rm{d}}(\rr) &= \quad \sum_{i=1}^N \int \frac{\psi_i^*(\rr') \psi_i(\rr')}{\abs{\rr-\rr'}} \dd \rr', \\
    \hat{V}_{\rm{ex}} \psi_0(\rr) &= -\sum_{i=1}^N \psi_i(\rr) \int \frac{\psi_i^*(\rr') \psi_0(\rr')}{\abs{\rr-\rr'}} \braket{\sigma_i}{\sigma_0} \dd \rr',
\end{aligned}
\end{equation}
```
where ``N`` and ``N_{\rm{atm}}`` denote the number of electrons and atoms, respectively;
``\psi_i(\rr)`` and ``\sigma_i`` denote the molecular orbital and the spin state of the electron with index ``i``, ``\braket{\sigma_i}{\sigma_j}=1`` for electrons ``i`` and ``j`` with the same spin state, and ``\braket{\sigma_i}{\sigma_j}=0`` otherwise;
``Z_A`` and ``\bm{R}_A`` are the nuclear charge and position of atom with index ``A``.

Representing the rotated reference function in Eq.~\eqref{eq:WFAT_G} with a linear combination of spherical harmonics using the Wigner-``D`` matrix allows for efficient numerical evaluation of the structure factor using the coefficients calculated beforehand:
```math
\begin{equation}
    G_\nu(\theta,\chi) = \ee^{-\kappa\mu_F} \sum_{l=\abs{m}}^{\infty} \sum_{m'=-l}^{l} I_{lm'}^\nu d_{mm'}^l(\theta) \ee^{-\ii m' \chi},
\end{equation}
```
where the ``\ee^{-\ii m \phi}`` in the expansion of ``D_{mm'}^l(\phi,\theta,\chi)=\ee^{-\ii m \phi} d_{mm'}^l(\theta) \ee^{-\ii m' \chi}`` is omitted because it doesn't play a part in the final result, and the coefficient ``I_{lm'}^\nu`` has the following expression:
```math
\begin{equation}
    I_{lm'}^\nu = \int \dd\rr\ \Omega_{lm'}^{\nu*}(\rr) \hat{V}_{\rm{c}} \psi_0(\rr).
\end{equation}
```

The original WFAT gives the instantaneous tunneling ionization rate ``w=\dd W/\dd t``, however, without the dependence of ``\kt``.
In order to apply WFAT to give initial condition of the electron, we have to reform the original WFAT to give ``\kt``-dependent rate.
Here we adopt the ``\kt``-dependence in MO-ADK, which gives
```math
\begin{equation}
    \dd W/\dd \kkt \dd t \propto \kt^{2\abs{m}} \ee^{-\kappa \kt^2/F}
\end{equation}
```
under the small-``\kt`` limit.
We modify the field factor ``\mathcal{W}_\nu(F)`` according to the ``\kt``-dependence above, which gives the modified field factor
```math
\begin{equation}
\begin{aligned}
    \mathcal{W}_\nu(F,\kt)
    &= \mathcal{W}_\nu(F) \frac{(\kappa/F)^{\abs{m}+1}}{\abs{m}!} \kt^{2\abs{m}} \ee^{-\kappa \kt^2/F} \\
    &\approx \frac12 \frac{\kappa^{\abs{m}+2}}{F^{\abs{m}+1}\abs{m}!} \left(\frac{4\kappa^2}{F}\right)^{2n^*-2n_\xi-\abs{m}-1} \!\!\!\! \kt^{2\abs{m}} \exp\left[ -\frac23 \frac{(\kt^2+\kappa^2)^{3/2}}{F} \right],
\end{aligned}
\end{equation}
```
where we choose the normalization coefficient so that
```math
\begin{equation}
    \mathcal{W}_\nu(F) = \int_{0}^{\infty} \mathcal{W}_\nu(F,\kt) 2\pi\kt \dd\kt.
\end{equation}
```
In this way we obtain the ``\kt``-dependent rate given by the WFAT:
```math
\begin{equation}
    \frac{\dd W^{\rm{WFAT}}}{\dd \kkt \dd t} = \sum_{\nu} \abs{G_\nu(\theta,\chi)}^2 \mathcal{W}_\nu(F,\kt).
\end{equation}
```
