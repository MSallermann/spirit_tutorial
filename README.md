## 1.1 Classical Heisenberg Hamiltonian

- A fundamental model in magnetism used to describe the energy of interacting spins. 
- The model treats each magnetic moment as a classical vector.

$$
\begin{aligned}
	\mathcal{H} = & -\sum_{\langle i,j \rangle} \Bigl( J_{ij}\,\mathbf{S}_i\cdot\mathbf{S}_j
	+ \mathbf{D}_{ij}\cdot (\mathbf{S}_i\times\mathbf{S}_j) \Bigr)\\
	&-\sum_{i} K (\mathbf{S}_i\cdot \mathbf{n})^2  -\sum_i \mu_i\, \mathbf{B}\cdot\mathbf{S}_i\\
	&+ \frac{\mu_0}{4\pi} \sum_{i<j} \frac{\mu_i\,\mu_j}{r_{ij}^3}
	\left[\mathbf{S}_i\cdot\mathbf{S}_j - 3\,(\mathbf{S}_i\cdot \hat{r}_{ij})(\mathbf{S}_j\cdot \hat{r}_{ij})\right]\\
	&+ \text{other terms}
\end{aligned}
$$

- Cheap to evaluate. Simulations with $10^6$ spins over nanoseconds are possible.
- Coarse-graining the Heisenberg Hamiltonian is the basis of micromagnetics

- Exchange Interaction:  $-\sum_{\langle i,j \rangle} J_{ij}\, \mathbf{S}_i\cdot\mathbf{S}_j$ 
	- $J_{ij}$: Exchange coupling constant between spins at sites $i$ and $j$.
	- $J_{ij} > 0$ favors ferromagnetic alignment, while $J_{ij} < 0$ favors antiferromagnetic alignment.
    - $\mathbf{S}_i$ and $\mathbf{S}_j$: Spin vectors (unit vectors) at sites $i$ and $j$
    
- Dzyaloshinskii–Moriya Interaction (DMI):  $-\sum_{\langle i,j \rangle} \mathbf{D}_{ij}\cdot (\mathbf{S}_i\times\mathbf{S}_j)$ 
	- $\mathbf{D}_{ij}$: DMI vector for the spin pair at sites $i$ and $j$. 
	-  Favors chiral (twisted or non-collinear) spin configurations. 
	
- Anisotropy: $-\sum_{i} K\, (\mathbf{S}_i\cdot \mathbf{n})^2$ 
	- $K$: Anisotropy constant. 
	- $\mathbf{n}$: Preferred (easy) axis direction. 
	- Encourages spins to align along the direction of $\mathbf{n}$. 
	
- Zeeman Term: $-\sum_i \mu_i\, \mathbf{B}\cdot\mathbf{S}_i$ 
	- $\mu_i$: Magnetic moment at site $i$. 
	- $\mathbf{B}$: External magnetic field.
	- Represents the interaction of the spins with an applied magnetic field.
	
- Dipole–Dipole Interaction (DDI):  $\frac{\mu_0}{4\pi}\sum_{i<j} \frac{\mu_i\,\mu_j}{r_{ij}^3}\left[\mathbf{S}_i\cdot\mathbf{S}_j - 3\, (\mathbf{S}_i\cdot\hat{r}_{ij})(\mathbf{S}_j\cdot\hat{r}_{ij})\right]$ 
	- $\mu_0$: Permeability of free space. 
	- $\mu_i$ and $\mu_j$: Magnetic moments at sites $i$ and $j$. - $r_{ij}$: 
	- Distance between spins at sites $i$ and $j$.
	- $\hat{r}_{ij}$: Unit vector pointing from site $i$ to site $j$. 
	- Accounts for the long-range dipolar interaction between spins.
	
## 1.2 Landau-Lifshitz-Gilbert equation

- The LLG equation governs the time evolution of magnetic moments (spins) in a material under the influence of effective magnetic fields. It can be written as
$$
\frac{d\mathbf{S}}{dt} = -\frac{\gamma}{1+\alpha^2}\,\mathbf{S}\times\Bigl( \mathbf{H}_{\mathrm{eff}} + \boldsymbol{\zeta}(t) \Bigr)
-\frac{\gamma\alpha}{1+\alpha^2}\,\mathbf{S}\times\Bigl[\mathbf{S}\times\Bigl( \mathbf{H}_{\mathrm{eff}} + \boldsymbol{\zeta}(t) \Bigr)\Bigr],
$$

- $\mathbf{S}$: Unit vector representing the spin. 
- $\mathbf{H}_{\mathrm{eff}}$: Effective magnetic field (includes contributions from exchange, anisotropy, Zeeman, etc.). 
- $\gamma$: Gyromagnetic ratio. 
- $\alpha$: Gilbert damping parameter. 
- $\boldsymbol{\zeta}(t)$: Stochastic thermal field representing thermal fluctuations. 
- The Langevin noise-term $\zeta(t)$ can be used to model the time-evolution of a spin system in the canonical (NVT) ensemble. It is a Gaussian noise term with zero mean and the correlation function
$$
\langle \zeta_i(t) \, \zeta_j(t') \rangle = \frac{2\alpha\, k_B T}{\gamma \mu_s} \, \delta_{ij}\, \delta(t-t'),
$$

## 1.3 Geodesic nudged elastic band method
TODO

## 1.4 HTST method
TODO

## 1.5 The Spirit code
TODO

# 2. Exercises
The exercised should be small and self-contained. They should be represented in Jupyter notebooks. TODO.
## 2.1 A small skyrmion