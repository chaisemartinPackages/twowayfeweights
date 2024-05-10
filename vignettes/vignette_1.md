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

It is enough to run the following command to store the main table in a .tex file:

```applescript
twowayfeweights_out, saving(filename.tex)
```

The line above will save the main table in *filename.tex* as a TeX tabular. If you want to save the table as an individual TeX file that can be independently run, you could specify the **standalone** argument to the previous line. In this way, the table will be saved in a TeX file with *standalone* document class.

Lastly, the program outputs multiple tables if the *other_treatments* option is specified. In this case `twowayfeweights_out` will concatenate each of these tables in a single TeX tabular (or standalone):

```applescript
twowayfeweights Y G t D, type(feTR) other_treatments(D2)
twowayfeweights_out, saving(filename.tex)
```

## Adding scalars to external regression tables
