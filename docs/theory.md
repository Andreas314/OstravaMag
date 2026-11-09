# Theory
On this page, we aim to provide the basic overview of important mathematical expressions used within the code. Additionally, we show different formalisms of the Heisenberg Hamiltonian that are used in the codes with an existing interface as well as the one OstravaMag uses externally and internally. 

## Heisenberg Hamiltonian
The energy of a system due to isotropic exchange interaction is given by the Heisenberg model 
$$
E = -\sum_{i\neq j, \boldsymbol{\delta}} J_{ij}(\boldsymbol{\delta})\mathbf{e}_i\cdot \mathbf{e}_j
$$
where $J_{ij}(\boldsymbol{\delta})$ is the exchange interaction between sites $i$ and $j$ with displacement vector $\boldsymbol{\delta}$. The vectors $\boldsymbol{e}_j$ and $\boldsymbol{e}_j$ are unit spin vectors on sites $i$ and $j$. In this formalism, each pair of the indices is counted twice.  

&emsp;There are many different conventions of the Heisenberg Hamiltonian, which affect both sign and size of $J_{ij}$. The convention written above is the same one as the one being used in the TB2J code. OstravaMag also uses it to output the $\textit{Jr.dat}$ file and to read the exchange constants from the input file for the minimum input method. 

&emsp;Internally, the code uses a different convention, in which the exchange constant is multiplied by 2 ($E$ does not contain the double counting), which is only relevant if one wishes to extend the existing interfaces. This is the same convention, which is also used in the GPAW code.  

&emsp;The remaining 2 codes, SPR-KKR and jx, use the following transformation of the exchange constants
$$
J'_{ij}(\boldsymbol{\delta}) = J_{ij}(\boldsymbol{\delta})\mathbf{e}_i\cdot \mathbf{e}_j
$$
which means that they already contain the dot product of the unit spin vectors. In the collinear case, this means that sites with opposite sign than in the initial convention.

&emsp; Other conventions employed may differ in the global minus sign or use unnormalized spin vectors and the exchange constants must be rescaled accordingly.
## Magnon dispersion
To calculate the dispersion we use the Linear Spin-Wave Theory [5],where we first take the Fourier transformation of the exchange constants in the real space
$$
J_{ij}(\boldsymbol{k})=\sum_{\boldsymbol{\delta}}J_{ij}(\boldsymbol{\delta})e^{i\boldsymbol{\delta}\cdot\boldsymbol{k}}.
$$
The magnon spectrum is then calculated as a set of eigenvalues in each k-point of the Brillouin zone. The generalized eigenvalue problem has the following form
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
are $N\times N$ matrices. The submatrices $A$ and $B$ have the following form
$$
A_{\alpha\beta}(\boldsymbol{k})={}
\frac{2\delta_{\alpha\beta}}{M_\alpha}\sum_\gamma
\boldsymbol{e}_\alpha\cdot \boldsymbol{e}_\gamma
 J_{\alpha\gamma}(\mathbf{0})
-\frac{2J_{\alpha\beta}(\mathbf{k})}{\sqrt{M_\alpha M_\beta}}
\frac{(1+\boldsymbol{e}_\alpha\cdot \boldsymbol{e}_\beta)}{2},\\
B_{\alpha\beta}(\boldsymbol{k})={}-
\frac{2J_{\alpha\beta}(\mathbf{k})}{\sqrt{M_\alpha M_\beta}}
\frac{(1-\boldsymbol{e}_\alpha\cdot \boldsymbol{e}_\beta)}{2},
$$
where $M_\alpha$ is the size of the magnetic moment on site $\alpha$.
## Thermodynamic quantities
Once the magnon dispersion in calculated, the code determines the magnon density of states $g(\omega)$, which lets us calculate thermodynamic quantities coming from spin waves. The **occupation** of a certain mode is given by (the branch index is ommited)
$$
n(\omega, T) = \frac{1}{e^{\hbar\omega/k_B T} - 1}
$$
which lets us calculate the **internal energy**
$$
E(T) = \int_0^\infty g(\omega) \hbar\omega \left[\frac{1}{2} + n(\omega, T)\right] d\omega
$$
its derivative with respect to temperature gives us the **constant-volume heat** capacity:
$$
C_V(T) = \int_0^\infty g(\omega) k_B \left(\frac{\hbar\omega}{k_B T}\right)^2 \frac{e^{\hbar\omega/k_B T}}{(e^{\hbar\omega/k_B T} - 1)^2} d\omega
$$
the **free energy** of magnons is given by 
$$
F(T) = \frac{1}{2} \int_0^\infty g(\omega) \hbar\omega \, d\omega + k_B T \int_0^\infty g(\omega) \ln\left[1 - e^{-\hbar\omega/k_B T}\right] d\omega
$$
and finally, **entropy** is calculated as
$$
S(T) = \frac{1}{2T} \int_0^\infty g(\omega) \hbar\omega \coth\left(\frac{\hbar\omega}{2k_B T}\right) d\omega - k_B \int_0^\infty g(\omega) \ln\left[2\sinh\left(\frac{\hbar\omega}{2k_B T}\right)\right] d\omega
$$
## Quasi-harmonic approximation
In the qha, it is necessary to first calculate the free energy as written above for multiple volumes to obtain $F(T, V)$ and to know the total energy of the system for each volume $E_{tot}(V)$. Gibbs free energy for given pressure $p$ is than found as
$$
G(T, p)=\min_V{(F(T, V) - F(0, V)+ E_{tot}(V)+ pV)}.
$$
The equilibrium volume for given pressure and temperature can be obtained in a similar manner
$$
V_0(T,p) = \operatorname*{arg\,min}_{V}\; \left( F(V,T) - F_0(V, T)+ E_{tot}(V) + pV \right)
$$
The heat capacity at constant pressure is found as
$$
C_p(T)=\frac{1}{T}\frac{\partial^2 G(T, p)}{\partial T^2}.
$$ 
## Calculations of $T_c$ 
To approximate the critical temperature of a given system, two approximation are employed: Mean-Field Approximation (MFA) and Random Phase Approximation (RPA). The implementation was done according to Ref [6]  
&emsp; To get $T_c$ in MFA, the eigenvalues $\beta_c^i$ and eigenvectors $\boldsymbol{s}_\alpha^i$ of the following system are found
$$
\sum_\gamma J_{\alpha\gamma}(\boldsymbol{0})\boldsymbol{s}^i_{\alpha} = \beta_c^i\frac{S_\alpha}{(1+S_\alpha)}\boldsymbol{s}^i_{\alpha},
$$
where $S_\alpha$ is the size of the spin moment on the site $\alpha$. The critical temperature is then calculated as
$$
T_c^{MFA}=\max_{i}\frac{2\beta_c^i}{3k_B}.
$$
&emsp; For the case of RPA, an iterative scheme is employed. The following matrix is first defined
$$
N_{\alpha\alpha'}(\boldsymbol{q})=\delta_{\alpha\alpha'}\sum_\gamma J(\boldsymbol{0})_{\alpha\gamma}\tilde{\boldsymbol{s}_\gamma^i} - \tilde{\boldsymbol{s}_\gamma^i}J_{\alpha\alpha'}(\boldsymbol{q}).
$$
We then solve a self-consistent problem, where we start with $\tilde{\boldsymbol{s}_\gamma^i}$ set to the result of MFA and it is updated based on an equation given by
$$
(\tilde{\boldsymbol{s}_\alpha^i})_{new}=\frac{2(S_\alpha+1)}{3S_\alpha}\left(\frac{1}{\Omega}\int\text{d}\boldsymbol{q}(N_{old}(\boldsymbol{q}))_{\alpha\alpha}^{-1}\right)^{-1} .
$$
This iterative scheme is then run as long as two iteration differ by more than given tolerance. $T_c^{RPA}$ is given by the size of $\tilde{\boldsymbol{s}_\gamma^i}$ in the last iteration.  
**Note:** Sometimes, the prefactor $\frac{S_\alpha + 1}{S_\alpha}$ must be omitted to obtain correct results. Both of those version are output from the code.