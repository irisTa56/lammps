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

![expression](./expression.png)

### `pair_dpd_trans_tstat`

#### Syntax

```
pair_style dpd/trans/tstat Tstart Tstop cutoff seed
pair_coeff type1 type2 gamma gamma_t cutoff
```

* `Tstart,Tstop`: desired temperature at start/end of run (temperature units)
* `cutoff`: global cutoff for DPD interactions (distance units)
* `seed`: random number seed (positive integer)
* `gamma`: friction coefficient (force/velocity units)
* `gamma_t`: friction coefficient in transverse direction (force/velocity units)
* `cutoff_pair`: cutoff of each pair (distance units)

#### Expression

Same with `pair_dpd_trans` except for absence of the conservative term.

## Tally computing for pair/hybrid