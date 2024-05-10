# Storing twowayfeweights results

In this tutorial, we show how to save the output of `twowayfeweights` in a .tex file and how to integrate its results in some of the most popular post-estimation programs. This vignette is oriented towards Stata users.

## Setup

In this vignette, we use a randomly generated sample from a DGP with differential adoption timing and group specific treatment effects. Under these conditions, the TWFE coefficient will be a weighted average of groups' ATTs, with some negative weights (de Chaisemartin & D'Hautfoeuille, 2020).

``` applescript
clear
set seed 0
local G = 100
local T = 10
set obs `G'
gen g = _n
gen F_g = floor(2 + uniform() * (`T' - 2))
forv j = 1/`T' {
	gen t`j' = `j'
}
reshape long t, i(g F_g) j(j)
drop j
order g t F_g
gen D = t >= F_g
gen D2 = uniform() > 0.5
gen Y = uniform() + D*F_g
```

The output of `twowayfeweights` confirms the presence of negative weights.

```applescript
twowayfeweights Y g t D, type(feTR)
```

## Storing the main table




## Integration of scalars in external regression tables
