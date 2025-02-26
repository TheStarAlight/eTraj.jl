# [Theory: Initial Conditions](@id theory_init_cond)

Multiple theoretical frameworks for strong-field ionization provide methodologies to determine the initial conditions of classical electrons in trajectory simulations.
These initial conditions are characterized by three essential properties:
- Initial position ``\rr_0``, i.e., the tunneling exit position;
- Initial momentum ``\kk_0`` [*Note*: Here we use ``\kk`` to denote the momentum/velocity of the electron when the laser is present, to distinguish it from the canonical momentum ``\pp=\kk-\AA``. After the laser vanishes, we have ``\pp=\kk``];
- The corresponding ionization probability ``W`` carried by each electron sample, which depends on the time-dependent laser field and properties of the target atom/molecule.

In this section we briefly revisit the available theories implemented in `eTraj`.
Atomic units (a.u.) are used throughout unless stated otherwise.

```@contents
Pages = ["theory1_init_cond.md"]
Depth = 2
```

---------------------------

## [Strong-Field Approximation with Saddle-Point Approximation (SFA-SPA)](@id SFA)

The *Strong-Field Approximation (SFA)* is originated from the Keldysh theory of strong-field ionization [^Keldysh_1965] [^Faisal_1973] [^Reiss_1980].
In contrast to perturbative methods and adiabatic tunneling theories, the SFA framework encompasses both multi-photon and tunneling processes in laser-atom interactions, as it comprehensively incorporates non-adiabatic effects of the laser-atom coupling.
This comprehensive applicability has facilitated the extensive use of SFA in theoretical studies of strong-field ionization.

[^Keldysh_1965]: L. V. Keldysh, Ionization in the field of a strong electromagnetic wave, *Sov. Phys. JETP* **20**, 1307 (1965).

[^Faisal_1973]: F. H. M. Faisal, Multiple absorption of laser photons by atoms, *J. Phys. B: At. Mol. Phys.* **6**, L89 (1973). DOI: [10.1088/0022-3700/6/4/011](https://doi.org/10.1088/0022-3700/6/4/011)

[^Reiss_1980]:  H. R. Reiss, Effect of an intense electromagnetic field on a weakly bound system, *Phys. Rev. A* **22**, 1786 (1980). DOI: [10.1103/PhysRevA.22.1786](https://doi.org/10.1103/PhysRevA.22.1786)

We consider an electron evolving under a combined field of the Coulomb field ``V(\rr)`` of the parent ion and the laser field ``\FF(t)=-\pd_t \AA(t)``.
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

Here lies the key idea of SFA: when the influence of the Coulomb field to the ionized electrons is weak compared with that of the external laser field, we may neglect the influence of the Coulomb field on the expression of ``M_{\pp}`` by replacing the time-evolution operator with a Coulomb-free one ``U_{\rm{f}}``, and meanwhile replacing the continuum state with the Volkov state ``\ket{\Psi_{\pp}^{\rm{V}}}`` which represents a free electron evolving under the same laser field:
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
Consequently, the expression for ``M_{\pp}`` becomes
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

Given that the variation of the phase factor ``\ee^{\ii\Sp(t)}`` is significantly more sensitive to changes in ``t`` compared to the prefactor, the integrand in Eq. (9) tends to oscillate rapidly in its complex phase. Consequently, contributions tend to cancel out over most of the domain, except near points where the phase ``\Sp(t)`` stabilizes, specifically at the saddle points.
The saddle points ``\ts=\tr+\ii\ti`` are the zeroes of the derivative of the complex function ``\Sp(t)``, which satisfy
```math
\begin{equation}
    -\Sp'(\ts) = \frac12 [\pp+\AA(\ts)]^2 + \Ip = 0.
\end{equation}
```
The integral, representing the second term on the right-hand side of Eq. (9), contributes significantly only in the neighborhood of the endpoints ``t_0``, ``\tf``, and the saddle points ``\ts``. The contributions near the endpoints largely cancel out the first term.
Therefore, the ``M_{\pp}`` is now approximated with the integration around the saddle points:
```math
\begin{equation}
    M_{\pp} \approx \sum_{\ts} \int_{C_{\ts}} \dd\tau \bk{\pp+\AA(\tau)}{\psi_0} \cdot (-\ii\Sp'(\tau)) \ee^{-\ii\Sp(\tau)},
\end{equation}
```
with ``C_{\ts}`` the integration contour following the steepest-descent path related to ``\ts``.

Further evaluation of the prefactor ``\tilde{\psi}_0(\kk) \rvert_{\kk=\pp+\AA(t)} = \bk{\pp+\AA(t)}{\psi_0}`` (i.e., the momentum-space wavefunction) in the vicinity of the saddle points in Eq. (11) is essential before applying the SPA.
We assume the field points towards the ``+ z`` axis, for an atom target at the ``(l,m)`` state with ionization potential ``\Ip``, its wavefunction behaves asymptotically as [^Perelomov_1966]
```math
\begin{equation}
    \psi_0(\rr) \sim 2 C_{\kappa l} \kappa^{3/2} (\kappa r)^{n^*-1} \ee^{-\kappa r} Y_{lm}(\hat{\rr})
\end{equation}
```
for ``\kappa r \gg 1``, with ``\kappa=\sqrt{2\Ip}``, ``n^*=Z/\kappa`` the effective principal quantum number, ``Z`` the charge of the residual ion, ``Y_{lm}`` the spherical harmonics, and ``C_{\kappa l}`` the asymptotic coefficient for atoms, which can be approximated using the Hartree approximation formula [^Hartree_1928]
```math
\begin{equation}
    C_{\kappa l}^2 = \frac{2^{2n^*-2}}{n^* (n^*+l)! (n^*-l-1)!}.
\end{equation}
```
For atomic hydrogen at the ground state we have ``C_{\kappa l} = 1``.
Moreover, for non-integer ``n^*``, the formula can be naturally extended by replacing the factorials ``x!`` with Gamma functions ``\Gamma(x+1)``, i.e.,
```math
\begin{equation}
    C_{\kappa l}^2 = \frac{2^{2n^*-2}}{n^* \Gamma(n^*+l+1) \Gamma(n^*-l)}.
\end{equation}
```
Near the saddle points, corresponding to the limit where ``k^2 \to -\kappa^2``, the form of ``\tilde{\psi}_0(\kk)`` is dictated by the asymptotic behavior of the wavefunction
```math
\begin{equation}
    \tilde{\psi}_0(\kk) = \frac{C_{\kappa l}}{\sqrt{\pi}} \frac{2^{n^*+3/2}\kappa^{2n^*+1/2}\Gamma(n^*+1)}{(k^2+\kappa^2)^{n^*+1}} Y_{lm}(\hat{\kk}),
\end{equation}
```
where ``\Gamma`` is the gamma function [^note2].
By substituting the aforementioned expression into Eq. (11) and using the definition of ``\Sp(t)`` [Eq. (7)], we derive
```math
\begin{equation}
    M_{\pp} = \ii \frac{C_{\kappa l}}{\sqrt{\pi}} 2^{1/2}\kappa^{2n^*+1/2}\Gamma(n^*+1) \sum_{\ts} \int_{C_{\ts}} \frac{Y_{lm}(\hat{\kk}(\tau))}{[\Sp'(\tau)]^{n^*}} \ee^{\ii\Sp(\tau)} \dd\tau,
\end{equation}
```
where ``\hat{\kk}(\tau)`` denotes the complex unit vector along ``\kk(\tau)=\pp+\AA(\tau)``.
The evaluation method for spherical harmonics with complex arguments follows Appendix B of Ref. [^Pisanty_2017]; see also the note [^note3]
To address cases where the integrand exhibits a singularity at ``\ts``, a modified SPA approach can be applied (see Appendix B of Ref. [^Gribakin_1997]):
```math
\begin{equation}
\begin{aligned}
    \int_{C_{\ts}} \frac{Y_{lm}(\hat{\kk}(\tau))}{[\Sp'(\tau)]^{n^*}} \ee^{\ii\Sp(\tau)} \dd\tau
    &\approx \frac{Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{n^*}} \int_{C_{\ts}} \frac{\ee^{\ii\Sp(\tau)}}{(\tau-\ts)^{n^*}} \dd\tau \\
    &\approx \frac{Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{n^*}} \cdot \ii^{n^*} \frac{\Gamma(n^*/2)}{2\Gamma(n^*)} \sqrt{\frac{2\pi}{-\ii\Sp''(\ts)}} [-2\ii\Sp''(\ts)]^{n^*/2} \ee^{\ii\Sp(\ts)}.
\end{aligned}
\end{equation}
```
Through this approach, we arrive at the following expression for the transition amplitude:
```math
\begin{equation}
    M_{\pp} = c_{n^*} C_{\kappa l} \sum_{\ts} \frac{Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{(n^*+1)/2}} \ee^{-\ii\Sp(\ts)},
\end{equation}
```
where ``c_{n^*} = \ii^{(n^*-5)/2} 2^{n^*/2+1} \kappa^{2n^*+1/2} \Gamma(n^*/2+1)`` represents a constant coefficient.

[^note2]: Here ``k`` is actually ``\sqrt{\kk\cdot\kk}`` and is not the conventional "norm" of the complex vector, which is ``\abs{\abs{\kk}}=\sqrt{\kk^*\cdot\kk}``. Here we normalize the complex vector ``\kk`` through ``\hat{\kk}=\kk/\sqrt{\kk\cdot\kk}``, according to Ref. [^Perelomov_1966].

[^note3]: The evaluation of spherical harmonic in the complex domain makes use of the fact that the solid harmonics ``S_{lm}(\rr) := r^l Y_{lm}(\hat{\rr})`` can be expressed in a polynomial of ``x,y,z``, see Ref. [^Caola_1978].

[^Perelomov_1966]: A. Perelomov, V. Popov, and M. Terent’ev, Ionization of atoms in an alternating electric field, *Sov. Phys. JETP* **23**, 924 (1966).

[^Hartree_1928]: D. R. Hartree, The wave mechanics of an atom with a non-coulomb central field. Part I. Theory and methods, *Math. Proc. Cambridge Philos. Soc.* **24**, 89 (1928). DOI: [10.1017/S0305004100011919](https://doi.org/10.1017/S0305004100011919)

[^Pisanty_2017]: E. Pisanty and Á. Jiménez-Galán, Strong-field approximation in a rotating frame: High-order harmonic emission from p states in bicircular fields, *Phys. Rev. A* **96**, 063401 (2017). DOI: [10.1103/PhysRevA.96.063401](https://doi.org/10.1103/PhysRevA.96.063401)

[^Caola_1978]:  M. J. Caola, Solid harmonics and their addition theorems, *J. Phys. A: Math. Gen.* **11**, L23 (1978). DOI: [10.1088/0305-4470/11/2/001](https://doi.org/10.1088/0305-4470/11/2/001)

[^Gribakin_1997]: G. F. Gribakin and M. Yu. Kuchiev, Multiphoton detachment of electrons from negative ions, *Phys. Rev. A* **55**, 3760 (1997). DOI: [10.1103/PhysRevA.55.3760](https://doi.org/10.1103/PhysRevA.55.3760)


The SFA phase ``\Sp(\ts)`` is determined by evaluating the integral
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

where ``S_{\pp,\rm{tun}}`` and ``S_{\pp,\rm{traj}}`` denote the complex phases accumulated during the tunneling process and the subsequent motion in the continuum, respectively.
The phase ``S_{\pp,\rm{tun}}`` corresponds to an imaginary time interval (from ``\ts`` to ``\tr``), during which the electron traverses the potential barrier with an "imaginary" momentum;
its real part signifies the quantum phase, whereas its imaginary part is associated with the ionization probability.

To apply the SFA for preparing initial conditions of photoelectrons, we assume that the electron is emitted at time ``\tr`` from the tunnel exit ``\rr_0`` with momentum ``\kk_0=\kk(\tr)``.
The initial momentum ``\kk_0``, when the Coulomb interaction with the nucleus is neglected, is correlated with the final momentum ``\pp`` via the relationship
```math
\begin{equation}
    \pp = \kk_0 - \int_{\tr}^{\infty} \FF(\tau) \dd\tau = \kk_0 - \AA(\tr).
\end{equation}
```

The initial position ``\rr_0`` corresponding to the tunnel exit, is determined by constructing a quantum tunneling trajectory.
This trajectory starts at a point with vanishing real part, representing the tunnel entrance; the electron tunnels through the barrier over the time interval ``\ts`` to ``\tr`` and emerges as a classical electron at the tunnel exit ``\rr_0`` with real position and momentum.
Thus, the expression for the initial position becomes:
```math
\begin{equation}
    \rr_0^{\rm{SFA-SPA}} = \Re \int_{\ts}^{\tr} [\pp+\AA(\tau)] \dd\tau = \Im \int_0^{\ti} \AA(\tr+\ii\tau) \dd\tau,
\end{equation}
```
with ``\Re`` and ``\Im`` denoting the real and imaginary parts, respectively.

The probability density of the electron sample in the final momentum space ``\pp`` is given by
```math
\begin{equation}
    \dd W^{\rm{SFA-SPA}}/\dd \pp = \sum_{\ts} \abs{\mathcal{P}^{\rm{SFA-SPA}}_{\pp}(\ts)}^2 \exp(-2\Im S_{\pp,\rm{tun}}(\ts)),
\end{equation}
```
where the prefactor encompasses all coefficients:
```math
\begin{equation}
\begin{aligned}
    \mathcal{P}^{\rm{SFA-SPA}}_{\pp}(\ts)
    = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{[\Sp''(\ts)]^{(n^*+1)/2}}
    = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{\left\{[\pp+\AA(\ts)]\cdot\FF(\ts)\right\}^{(n^*+1)/2}}.
\end{aligned}
\end{equation}
```

It should be noted that the ionization probability in Eq. (22) is formulated in terms of the final momentum coordinates ``\pp=(p_x,p_y,p_z)``.
However, in trajectory simulations, initial electrons are sampled in the ``(\tr,\kkt)`` coordinate system, with ``\kkt`` being the initial transverse momentum.
Consequently, if we sample the initial electrons within such a coordinate system, a Jacobian must be introduced as a prefix to the ionization probability.
Assuming laser propagation along the ``z`` axis and polarization in the ``xy`` plane, the transformed expression reads
```math
\begin{equation}
    \dd W^{\rm{SFA-SPA}}/\dd\tr\dd\kkt = \sum_{\ts} J(\kp,\tr) \abs{\mathcal{P}_{\pp}(\ts)}^2 \exp(-2\Im S_{\pp,\rm{tun}}(\ts)),
\end{equation}
```
where ``\kp`` represents the projection of ``\kkt`` onto the polarization plane (i.e., the ``xy`` plane), and the Jacobian is defined as
```math
\begin{equation}
    J(\tr,\kp) = \begin{vmatrix}\frac{\pd(p_x,p_y)}{\pd(\tr,\kp)}\end{vmatrix} =
    \begin{vmatrix}
        \pd p_x/\pd\tr & \pd p_x/\pd\kp \\
        \pd p_y/\pd\tr & \pd p_y/\pd\kp \\
    \end{vmatrix}.
\end{equation}
```

---------------------------

## [SFA-SPA with Non-adiabatic Expansion (SFA-SPANE)](@id SFAAE)

When the Keldysh parameter ``\gamma=\omega\kappa/F_0`` is small (where ``\omega`` denotes the laser angular frequency and ``F_0`` represents the peak field strength), non-adiabatic effects become less significant. Under these conditions, a non-adiabatic expansion scheme can be implemented to develop an approximation based on the SFA-SPA, termed the *SFA-SPA with Non-adiabatic Expansion (SFA-SPANE)* [^Ni_2018] [^Mao_2022] [^Ma_2021] [^Ma_2024].
It captures non-adiabatic effects to a considerable degree and yields results that are comparable to those obtained from the SFA-SPA for relatively small values of the Keldysh parameter.
SFA-SPANE comes with a closed analytical form, avoiding the necessity to solve the saddle-point equation, thereby speeding up the calculation.

[^Ni_2018]: H. Ni, N. Eicke, C. Ruiz, J. Cai, F. Oppermann, N. I. Shvetsov-Shilovski, and L.-W. Pi, Tunneling criteria and a nonadiabatic term for strong-field ionization, *Phys. Rev. A* **98**, 013411 (2018). DOI: [10.1103/PhysRevA.98.013411](https://doi.org/10.1103/PhysRevA.98.013411)

[^Mao_2022]: X. Mao, H. Ni, X. Gong, J. Burgdörfer, and J. Wu, Subcycle-resolved strong-field tunneling ionization: Identification of magnetic dipole and electric quadrupole effects, *Phys. Rev. A* **106**, 063105 (2022). DOI: [10.1103/PhysRevA.106.063105](https://doi.org/10.1103/PhysRevA.106.063105)

[^Ma_2021]: Y. Ma, J. Zhou, P. Lu, H. Ni, and J. Wu, Influence of nonadiabatic, nondipole and quantum effects on the attoclock signal, *J. Phys. B: At. Mol. Opt. Phys.* **54**, 144001 (2021). DOI: [10.1088/1361-6455/ac0d3e](https://doi.org/10.1088/1361-6455/ac0d3e)

[^Ma_2024]: Y. Ma, H. Ni, and J. Wu, Attosecond ionization time delays in strong-field physics, *Chin. Phys. B* **33**, 13201 (2024). DOI: [10.1088/1674-1056/ad0e5d](https://doi.org/10.1088/1674-1056/ad0e5d)


The SFA-SPANE method is applicable when the Keldysh parameter is small, and the non-adiabatic effect is insignificant, which corresponds to the small-``\ti`` case.
We expand the vector potential ``\AA(\ts)=\AA(\tr+\ii\ti)`` in the SFA-SPA around ``\ti=0``, up to the second order of ``\ti``:
```math
\begin{equation}
    \AA(\tr+\ii\ti) = \AA(\tr) - \ii\ti\FF(\tr) + \frac12 \ti^2 \FF'(\tr) + o(\ti^2).
\end{equation}
```
Inserting Eq. (26) into the saddle-point equation in the SFA-SPA [Eq. (10)] leads to
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
which allow for the derivation of analytical expressions of the ionization probability and other quantities.

The initial position ``\rr_0`` in SFA-SPANE is given by
```math
\begin{equation}
    \rr_0^{\rm{SFA-SPANE}} = \Im \int_0^{\ti} \AA(\tr+\ii\tau) \dd\tau = - \frac{\FF}{2} \frac{\kt^2+\kappa^2}{F^2-\kk_0\cdot\FF'}.
\end{equation}
```

The ``\Im S_{\pp,\rm{tun}}`` term, which is related to the ionization probability, within the context of SFA-SPANE, is expressed as
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
where ``\kkt`` corresponds to ``\kk(\tr)`` within the framework of SFA-SPANE due to the vanishing initial longitudinal momentum, as indicated by Eq. (27).
The prefactor reads
```math
\begin{equation}
    \mathcal{P}^{\rm{SFA-SPANE}}_{\pp}(\ts) = c_{n^*} \frac{C_{\kappa l} Y_{lm}(\hat{\kk}(\ts))}{\left[ (\kt^2+\kappa^2)(F^2-\kk_0\cdot\FF') \right]^{(n^*+1)/4}}.
\end{equation}
```

---------------------------

## [Ammosov-Delone-Krainov (ADK)](@id ADK)

The *Ammosov-Delone-Krainov (ADK)* theory [^Ammosov_1986] [^Delone_1998] provides a framework for investigating adiabatic tunneling in strong-field ionization and represents, essentially, the adiabatic limit of the SFA.

[^Ammosov_1986]: M. Ammosov, N. Delone, and V. Krainov, Tunnel ionization of complex atoms and of atomic ions in an alternating electromagnetic field, *Sov. Phys. JETP* **64**, 1191 (1986).

[^Delone_1998]: N. B. Delone and V. P. Krainov, Tunneling and barrier-suppression ionization of atoms and ions in a laser radiation field, *Phys. Usp.* **41**, 469 (1998). DOI: [10.1070/PU1998v041n05ABEH000393](https://doi.org/10.1070/PU1998v041n05ABEH000393)


In the adiabatic limit, the laser field can be considered static; consequently, we have ``\FF'(t)=0``, with higher-order derivatives of ``\FF(t)`` also remaining zero.
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
Expanding Eq. (33) under the small-``\kt`` approximation results in
```math
\begin{equation}
    \dd W^{\rm{ADK}}/\dd\pp \propto \exp \left( -\frac{2\kappa^3}{3F} \right) \exp \left( -\frac{\kappa\kt^2}{F} \right),
\end{equation}
```
which is actually the exponential term of the well-known ADK rate.
It should be noted, however, that the outcome of our approach—specifically, the application of the adiabatic limit within the SFA-SPA framework—differs slightly from the actual ADK rate in terms of the prefactor.
This is because the SFA framework neglects Coulomb potential in the final state, which has been shown to result in a lower ionization rate.

To address this discrepancy, an additional Coulomb-correction (CC) factor, introduced into Eq. (33), helps bridge the gap:
```math
\begin{equation}
    C^{\rm{CC}} = \left(\frac{2\kappa^3}{F}\right)^{n^*} \!\!\!\! \left(1+2\gamma/e\right)^{-2n^*} \left[\Gamma\left(\frac{n^*}{2}+1\right)\right]^{-2}.
\end{equation}
```
We note that this CC factor is implemented in all initial-condition methods that are derived from the SFA.

The tunnel exit is determined using the same methodology:
```math
\begin{equation}
    \rr_0^{\rm{ADK}} = \Im \int_0^{\ti} \AA(\tr+\ii\tau) \dd\tau = - \frac{\FF}{2} \frac{\kt^2+\kappa^2}{F^2},
\end{equation}
```
a result that we associate with the "``\Ip/F``" model. However, there is a subtle distinction in that we replace the ionization potential ``\Ip=\kappa^2/2`` with the effective potential ``\tilde{\Ip}=(\kappa^2+\kt^2)/2`` to account for the initial kinetic energy, ensuring that adiabatic tunneling is accurately described by the condition ``E=\kt^2/2+\rr_0^{\rm{ADK}}\cdot\FF=-\Ip``.


---------------------------

## [Molecular SFA-SPA/SFA-SPANE/ADK](@id MOSFA)

The atomic SFA theory and its adiabatic variants discussed in the previous sections can be systematically extended to molecular systems [^Muth_2000] [^Tong_2002] [^Kjeldsen_2004] [^Kjeldsen_2005].
Under the Born-Oppenheimer approximation and the single-active-electron (SAE) approximation, the strong-field ionization of molecules can be modeled as the interaction between the laser field and the ionizing orbital [often the highest occupied molecular orbital (HOMO)] ``\psi_0(\rr)`` within the effective potential of the parent ion.

[^Muth_2000]: J. Muth-Böhm, A. Becker, and F. H. M. Faisal, Suppressed molecular ionization for a class of diatomics in intense femtosecond laser fields, *Phys. Rev. Lett.* **85**, 2280 (2000). DOI: [10.1103/PhysRevLett.85.2280](https://doi.org/10.1103/PhysRevLett.85.2280)

[^Tong_2002]: X. M. Tong, Z. X. Zhao, and C. D. Lin, Theory of molecular tunneling ionization, *Phys. Rev. A* **66**, 33402 (2002).DOI: [10.1103/PhysRevA.66.033402](https://doi.org/10.1103/PhysRevA.66.033402)

[^Kjeldsen_2004]: T. K. Kjeldsen and L. B. Madsen, Strong-field ionization of N₂ : Length and velocity gauge strong-field approximation and tunnelling theory, *J. Phys. B: At. Mol. Opt. Phys.* **37**, 2033 (2004). DOI: [10.1088/0953-4075/37/10/003](https://doi.org/10.1088/0953-4075/37/10/003)

[^Kjeldsen_2005]: T. K. Kjeldsen, C. Z. Bisgaard, L. B. Madsen, and H. Stapelfeldt, Influence of molecular symmetry on strong-field ionization: Studies on ethylene, benzene, fluorobenzene, and chlorofluorobenzene, *Phys. Rev. A* **71**, 13418 (2005). DOI: [10.1103/PhysRevA.71.042508](https://doi.org/10.1103/PhysRevA.71.042508)

To generalize the atomic SFA to the molecular SFA (MO-SFA), we start from the transition amplitude given by Eq. (11).
In the molecular frame (MF), the asymptotic wavefunction can be expanded into spherical harmonics:
```math
\begin{equation}
    \psi_0^{\rm{MF}}(\rr) \sim \sum_{l,m} 2 C_{lm} \kappa^{3/2} (\kappa r)^{n^*-1} \ee^{-\kappa r} Y_{lm}(\hat{\rr}),
\end{equation}
```
where ``C_{lm}`` are asymptotic coefficients, and we continue to use ``n^*=Z/\kappa`` for simplicity, although it no longer represents the effective principal quantum number.
We assume that in the field frame (FF), the field ``\FF`` points along the ``z`` axis, and the rotation ``\RRh`` from the FF to the MF can be defined via a set of Euler angles ``(\phi,\theta,\chi)`` following the ``z-y'-z''`` convention, which satisfies
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
and the asymptotic behavior of the wavefunction in the FF is obtained by substituting Eq. (40) into Eq. (38), yielding
```math
\begin{equation}
     \psi_0^{\rm{FF}}(\rr) \sim \sum_{l,m,m'} 2 C_{lm} D_{m'm}^l(\phi,\theta,\chi) \kappa^{3/2} (\kappa r)^{n^*-1} \ee^{-\kappa r} Y_{lm'}(\hat{\rr}).
\end{equation}
```

It becomes evident that the primary distinction between the molecular and atomic versions of the theory lies in the formulation of the prefactor ``\mathcal{P}_{\pp}(\ts)``; the expressions for the tunneling exit position and the initial momentum remain unchanged.
Following the same procedure outlined in the previous sections, we derive the prefactor ``\mathcal{P}_{\pp}`` applicable to molecules:
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
Furthermore, it should be noted that upon introducing an additional Coulomb-correction factor [Eq. (36)], the ionization rate conforms to the original MO-ADK theory [^Tong_2002] within the adiabatic and small-``\kt`` limits.

---------------------------

## [Weak-Field Asymptotic Theory (WFAT)](@id WFAT)

The *Weak-Field Asymptotic Theory (WFAT)* generalizes the tunneling ionization from isotropic atomic potentials to arbitrary molecular potentials [^tolstikhin_2011] [^madsen_2013] [^madsen_2017] [^dnestryan_2018].
Compared with the MO-ADK theory, the WFAT naturally accounts for the influence of the permanent dipole moment of the molecule, and can, in its integral representation, calculate the structure factors [a similar concept to the asymptotic coefficients ``C_{lm}`` in Eq. (38)] based on the wavefunction close to the core, rather than using the wavefunction in the asymptotic region, allowing for enhanced accuracy in numerical simulations.

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
The structure factor ``G_\nu(\theta,\chi)`` is found by an integral related to the ionizing orbital and a reference function, which has significant contribution only in the vicinity of the nuclei and is insensitive to the wavefunction's asymptotic behavior:
```math
\begin{equation}
    G_\nu(\theta,\chi) = \ee^{-\kappa\mu_F} \int \dd\rr \ \Omega_\nu^*(\RRh^{-1}\rr) \hat{V}_{\rm{c}}(\rr) \psi_0(\rr),
\end{equation}
```
evaluated within the mean-field framework (MF), where ``\psi_0(\rr)`` represents the wavefunction of the ionizing orbital;
```math
\begin{equation}
    \bm{\mu} = - \int\dd\rr \ \psi_0^*(\rr) \rr \psi_0(\rr)
\end{equation}
```
denotes the orbital dipole moment within the MF, with ``\mu_F`` being its component along the direction of the external field;
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
``\hat{V}_\rm{c}=\hat{V}+Z/r`` denotes the core potential with the Coulomb tail removed, where ``Z`` signifies the asymptotic charge of the residual ion.

The effective potential ``\hat{V}`` characterizes the interaction between the ionizing electron and the residual parent ion.
We note that the hat notation is used to indicate that the potential operator is non-diagonal in the coordinate space.
Under the Hartree-Fock formalism, the effective potential comprises three components, specifically the nuclear Coulomb potential (``V_{\rm{nuc}}``), the direct (``V_{\rm{d}}``) and exchange (``V_{\rm{ex}}``) contributions of inter-electron interactions:
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
    \hat{V}_{\rm{ex}} \psi_0(\rr) &= -\sum_{i=1}^N \psi_i(\rr) \int \frac{\psi_i^*(\rr') \psi_0(\rr')}{\abs{\rr-\rr'}} \bk{\sigma_i}{\sigma_0} \dd \rr',
\end{aligned}
\end{equation}
```
where ``N`` and ``N_{\rm{atm}}`` denote the number of electrons and atoms, respectively;
``\psi_i(\rr)`` and ``\sigma_i`` represent the molecular orbital and the spin state of the electron with index ``i``, ``\bk{\sigma_i}{\sigma_j}=1`` for electrons ``i`` and ``j`` with the same spin state, and ``\bk{\sigma_i}{\sigma_j}=0`` otherwise;
``Z_A`` and ``\bm{R}_A`` corresponds to the nuclear charge and position of atom with index ``A``.

Representing the rotated reference function in Eq. (48) with a linear combination of spherical harmonics using the Wigner-``D`` matrix allows for efficient numerical evaluation of the structure factor using the coefficients calculated beforehand:
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
To adapt WFAT for preparing initial conditions of the electron samples, it is necessary to reformulate the original WFAT to incorporate a ``\kt``-dependent rate.
Here we adopt the ``\kt``-dependence in MO-ADK, which gives
```math
\begin{equation}
    \dd W/\dd t \dd \kkt \propto \kt^{2\abs{m}} \ee^{-\kappa \kt^2/F}
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
where the normalization coefficient is chosen such that
```math
\begin{equation}
    \mathcal{W}_\nu(F) = \int_{0}^{\infty} \mathcal{W}_\nu(F,\kt) 2\pi\kt \dd\kt.
\end{equation}
```
Through this approach, we derive the ``\kt``-dependent rate as provided by the WFAT:
```math
\begin{equation}
    \frac{\dd W^{\rm{WFAT}}}{\dd t \dd \kkt} = \sum_{\nu} \abs{G_\nu(\theta,\chi)}^2 \mathcal{W}_\nu(F,\kt).
\end{equation}
```
