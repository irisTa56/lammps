## Transverse DPD

### Correspondence with existing files

* `pair_dpd_trans.h` - `pair_dpd.h`
* `pair_dpd_trans.cpp` - `pair_dpd.cpp`
* `pair_dpd_trans_tstat.h` - `pair_dpd_tstat.h`
* `pair_dpd_trans_tstat.cpp` - `pair_dpd_tstat.cpp`
* `random_ziggurat.h` - `random_mars.h`
* `random_ziggurat.cpp` - `random_mars.cpp`

### `pair_dpd_trans`

#### Syntax

```
pair_style dpd/trans T cutoff seed
pair_coeff type1 type2 A gamma gamma_t cutoff_pair
```

* `T`: temperature (temperature units)
* `cutoff`: global cutoff for DPD interactions (distance units)
* `seed`: random number seed (positive integer)
* `A`: amplitude of conservative force (force units)
* `gamma`: friction coefficient (force/velocity units)
* `gamma_t`: friction coefficient in transverse direction (force/velocity units)
* `cutoff_pair`: cutoff of each pair (distance units)

#### Expression

$$
\begin{aligned}
\vec{f}_{ij} &= \vec{F}^{\mathrm{C}} + \vec{F}^{\mathrm{D}} + \vec{F}^{\mathrm{R}}, \\
\vec{F}^{\mathrm{C}} &= Aw(r)\vec{r}_{ij}, \\
\vec{F}^{\mathrm{D}} &= -\gamma w^{2}(r)\bm{P}\vec{v}_{ij} -\gamma_{t} w^{2}(r)\bm{Q}\vec{v}_{ij}, \\
\vec{F}^{\mathrm{R}} &= \sigma w(r)\Delta t^{-\frac{1}{2}}\bm{P}\vec{\theta}_{ij} + \sigma_{t} w(r)\Delta t^{-\frac{1}{2}}\bm{Q}\vec{\theta}_{ij},
\end{aligned}
$$

where $\bm{P} = \vec{r}_{ij}\vec{r}_{ij}^{T}$ and $\bm{Q} = \bm{I} - \vec{r}_{ij}\vec{r}_{ij}^{T}$; $\bm{I}$ is an identity matrix. According to the fluctuation-dissipation theorem, $\sigma$ and $\sigma_{t}$ is equal to $\sqrt{2kT\gamma}$ and $\sqrt{2kT\gamma_{t}}$, respectively. $\vec{\theta}_{ij}$ is a three-dimensional vector whose elements are random numbers following the standard normal distribution.


### `pair_dpd_trans_tstat`

#### Syntax

```
pair_style dpd/trans/tstat Tstart Tstop cutoff seed
pair_coeff type1 type2 gamma gamma_t cutoff
```

## Tally computing for pair/hybrid