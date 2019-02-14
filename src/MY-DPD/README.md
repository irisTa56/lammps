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

### `random_ziggurat`

Random number generator:

* Generates a random number following the standard normal distribution.
* Uses [Xorshift method](https://en.wikipedia.org/wiki/Xorshift) to generate a 32-bit integer, and then [Ziggurat method](https://en.wikipedia.org/wiki/Ziggurat_algorithm) to generate (determine) a random number following the standard normal distribution from the integer.
* Faster than [Marsaglia polar method](https://en.wikipedia.org/wiki/Marsaglia_polar_method), which is used for `pair_dpd` and `pair_dpd_tstat`.
  * Even for *standard* DPD, `pair_dpd_trans`/`pair_dpd_trans_tstat` with setting `gamma_t` to zero runs faster than `pair_dpd`/`pair_dpd_tstat`.

## Tally computing for pair/hybrid/overlay

Difference of `pair/hybrid/overlay/tally` from `pair/hybrid/overlay` is how to use `pair_modify`.
For `pair/hybrid/overlay/tally`, the `compute/tally` keyword of `pair_modify` accepts *id* of a `compute` command.
In other words, a Compute instance can be associated with a sub-style of the hybrid pair style by using `compute/tally` keyword with *id* of the instance.
Note that one *Compute* instance should be added to only one sub-style, whereas multiple *Compute* instances can be added to the same sub-style.

**Original Lammps**

```
pair_hybrid tersoff lj/cut/coul/long 12.0
pair_modify pair tersoff compute/tally no
compute 1 lower force/tally upper
compute 2 left pe/tally right
```

The above example computes per-atom forces between the `lower` and `upper` group and per-atom potential energies between the `left` and `right` group for `lj/cut/coul/long` only.
There is no way to computes forces/energies for `tersoff` simultaneously with `lj/cut/coul/long`.
This is because the `compute/tally` keyword of `pair_modify` accepts only `yes` or `no` and affects **all** the *Compute* instances of tally computing.

See also:

* [pair_style hybrid command — LAMMPS documentation](https://lammps.sandia.gov/doc/pair_hybrid.html)
* [compute force/tally command — LAMMPS documentation](https://lammps.sandia.gov/doc/compute_tally.html)
* [pair_modify command — LAMMPS documentation](https://lammps.sandia.gov/doc/pair_modify.html)

**This Package**

```
pair_hybrid tersoff lj/cut/coul/long 12.0
pair_modify pair tersoff compute/tally 1
pair_modify pair lj/cut/coul/long compute/tally 2
compute 1 lower force/tally upper
compute 2 left pe/tally right
```

The above example computes per-atom forces between the `lower` and `upper` group for `tersoff`, and per-atom potential energies between the `left` and `right` group for `lj/cut/coul/long`.

For more examples of usage, please see [examples of irisTa56/wapylmp](https://github.com/irisTa56/wapylmp).