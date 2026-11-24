# Theory
On this page, we aim to provide a basic overview of the important mathematical expressions used within the code. Additionally, we show different formalisms of the Heisenberg Hamiltonian that are used in the codes with an existing interface, as well as the one that OstravaMag uses externally and internally.
  

## Heisenberg Hamiltonian
The energy of a system due to the isotropic exchange interaction is given by the Heisenberg model

$$
E = -\sum_{i\neq j, \boldsymbol{\delta}} J_{ij}(\boldsymbol{\delta})\mathbf{e}_i\cdot \mathbf{e}_j
$$

where $J_{ij}(\boldsymbol{\delta})$ is the exchange interaction between sites $i$ and $j$ with displacement vector $\boldsymbol{\delta}$. The vectors $\boldsymbol{e}_i$ and $\boldsymbol{e}_j$ are unit spin vectors on sites $i$ and $j$. In this formalism, each pair of indices is counted twice.
  

&emsp;There are many different conventions for the Heisenberg Hamiltonian, which affect both the sign and the size of $J_{ij}$. The convention written above is the same as the one used in the TB2J code. OstravaMag also uses it to output the *Jr.dat* file and to read the exchange constants from the input file for the minimum input method.
  

&emsp;Internally, the code uses a different convention, in which the exchange constant is multiplied by 2 ($E$ does not contain the double counting), which is only relevant if one wishes to extend the existing interfaces. This is the same convention that is also used in the GPAW code.
  

&emsp;The remaining two codes, SPR-KKR and jx, use the following transformation of the exchange constants

$$
J'_{ij}(\boldsymbol{\delta}) = J_{ij}(\boldsymbol{\delta})\mathbf{e}_i\cdot \mathbf{e}_j
$$

which means that they already contain the dot product of the unit spin vectors. In the collinear case, this means that sites with opposite spin orientation have an exchange constant of the opposite sign compared with the initial convention.
  

&emsp;Other conventions employed may differ in the global minus sign or use unnormalized spin vectors, in which case the exchange constants must be rescaled accordingly.
## Magnon dispersion
To calculate the dispersion we use the Linear Spin-Wave Theory [5], where we first take the Fourier transformation of the exchange constants in real space

$$
J_{ij}(\boldsymbol{k})=\sum_{\boldsymbol{\delta}}J_{ij}(\boldsymbol{\delta})e^{i\boldsymbol{\delta}\cdot\boldsymbol{k}}.
$$

The magnon spectrum is then calculated as a set of eigenvalues at each k-point of the Brillouin zone. The generalized eigenvalue problem has the following form

$$
\Sigma_z\mathcal{H}_{\mathrm{BdG}}(\mathbf{k})v_n(\mathbf{k})
=\varepsilon_n(\mathbf{k})v_n(\mathbf{k}),
$$

where

$$
\mathcal{H}_{\mathrm{BdG}}(\mathbf{k})=
\begin{pmatrix}A(\mathbf{k})&B(\mathbf{k})\\
B^{\dagger}(\mathbf{k})&A^T(-\mathbf{k})\end{pmatrix},
\qquad
\Sigma_z=\begin{pmatrix}\mathbb{1}&0\\0&-\mathbb{1}\end{pmatrix}
$$

are $2N\times 2N$ matrices. The submatrices $A$ and $B$ have the following form

$$
\begin{array}{l}
A_{\alpha\beta}(\boldsymbol{k}) = \frac{2\delta_{\alpha\beta}}{M_\alpha}\sum_\gamma
\boldsymbol{e}_\alpha\cdot \boldsymbol{e}_\gamma J_{\alpha\gamma}(\mathbf{0})
-\frac{2J_{\alpha\beta}(\mathbf{k})}{\sqrt{M_\alpha M_\beta}}
\frac{1+\boldsymbol{e}_\alpha\cdot \boldsymbol{e}_\beta}{2}, \\
B_{\alpha\beta}(\boldsymbol{k}) = -\frac{2J_{\alpha\beta}(\mathbf{k})}{\sqrt{M_\alpha M_\beta}}
\frac{1-\boldsymbol{e}_\alpha\cdot \boldsymbol{e}_\beta}{2},
\end{array}
$$

where $M_\alpha$ is the size of the magnetic moment on site $\alpha$.
## Thermodynamic quantities
Once the magnon dispersion is calculated, the code determines the magnon density of states $g(\omega)$, which lets us calculate the thermodynamic quantities arising from spin waves. The **occupation** of a certain mode is given by (the branch index is omitted)

$$
n(\omega, T) = \frac{1}{e^{\hbar\omega/k_B T} - 1}
$$

which lets us calculate the **internal energy**

$$
E(T) = \int_0^\infty g(\omega) \hbar\omega \left[\frac{1}{2} + n(\omega, T)\right] d\omega
$$

its derivative with respect to temperature gives us the **constant-volume heat capacity**:

$$
C_V(T) = \int_0^\infty g(\omega) k_B \left(\frac{\hbar\omega}{k_B T}\right)^2 \frac{e^{\hbar\omega/k_B T}}{(e^{\hbar\omega/k_B T} - 1)^2} d\omega
$$

the **free energy** of magnons is given by

$$
F(T) = \frac{1}{2} \int_0^\infty g(\omega) \hbar\omega \, d\omega + k_B T \int_0^\infty g(\omega) \ln\left[1 - e^{-\hbar\omega/k_B T}\right] d\omega
$$

and finally, the **entropy** is calculated as

$$
S(T) = \frac{1}{2T} \int_0^\infty g(\omega) \hbar\omega \coth\left(\frac{\hbar\omega}{2k_B T}\right) d\omega - k_B \int_0^\infty g(\omega) \ln\left[2\sinh\left(\frac{\hbar\omega}{2k_B T}\right)\right] d\omega
$$

## Quasi-harmonic approximation
In the qha, it is necessary to first calculate the free energy as written above for multiple volumes to obtain $F(T, V)$ and to know the total energy of the system for each volume $E_{tot}(V)$. The Gibbs free energy for a given pressure $p$ is then found as

$$
G(T, p)=\min_V{(F(T, V) - F(0, V)+ E_{tot}(V)+ pV)}.
$$

The equilibrium volume for a given pressure and temperature can be obtained in a similar manner

$$
V_0(T,p) = \operatorname*{arg\,min}_{V}\; \left( F(V,T) - F(V, 0)+ E_{tot}(V) + pV \right)
$$

The heat capacity at constant pressure is found as

$$
C_p(T)=\frac{1}{T}\frac{\partial^2 G(T, p)}{\partial T^2}.
$$

## Calculations of $T_c$
To approximate the critical temperature of a given system, two approximations are employed: the Mean-Field Approximation (MFA) and the Random Phase Approximation (RPA). The implementation was done according to Ref. [6].
&emsp;To get $T_c$ in MFA, the eigenvalues $\beta_c^i$ and eigenvectors $\boldsymbol{s}_\alpha^i$ of the following system are found

$$
\sum_\gamma J_{\alpha\gamma}(\boldsymbol{0})\boldsymbol{s}^i_{\alpha} = \beta_c^i\frac{S_\alpha}{(1+S_\alpha)}\boldsymbol{s}^i_{\alpha},
$$

where $S_\alpha$ is the size of the spin moment on site $\alpha$. The critical temperature is then calculated as

$$
T_c^{MFA}=\max_{i}\frac{2\beta_c^i}{3k_B}.
$$

&emsp;For the case of RPA, an iterative scheme is employed. The following matrix is first defined

$$
N_{\alpha\alpha'}(\boldsymbol{q})=\delta_{\alpha\alpha'}\sum_\gamma J(\boldsymbol{0})_{\alpha\gamma}\tilde{\boldsymbol{s}_\gamma^i} - \tilde{\boldsymbol{s}_\alpha^i}J_{\alpha\alpha'}(\boldsymbol{q}).
$$

We then solve a self-consistent problem, in which we start with $\tilde{\boldsymbol{s}_\gamma^i}$ set to the result of the MFA, and it is updated based on an equation given by

$$
(\tilde{\boldsymbol{s}_\alpha^i})_{new}=\frac{2(S_\alpha+1)}{3S_\alpha}\left(\frac{1}{\Omega}\int\text{d}\boldsymbol{q}(N_{old}(\boldsymbol{q}))_{\alpha\alpha}^{-1}\right)^{-1}.
$$

This iterative scheme is then run for as long as two consecutive iterations differ by more than the given tolerance. $T_c^{RPA}$ is given by the size of $\tilde{\boldsymbol{s}_\gamma^i}$ in the last iteration.  
&emsp; To boost the convergence, the code takes advantage of simple mixing with a mixing parameter $\rho\in(0,1]$

$$
(\tilde{\boldsymbol{s}_\alpha^i})'_{new}=\rho(\tilde{\boldsymbol{s}_\alpha^i})_{new}+(1-\rho)(\tilde{\boldsymbol{s}_\alpha^i})_{old},
$$

where this new guess is used in the next iteration.  
**Note:** Sometimes, the prefactor $\frac{S_\alpha + 1}{S_\alpha}$ must be omitted to obtain correct results. Both of those versions are output by the code.
