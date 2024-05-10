# Storing twowayfeweights results

In this tutorial, we show how to save the output of `twowayfeweights` in a .tex file and how to integrate its results in some of the most popular post-estimation programs. This vignette is oriented towards Stata users.

+ [Setup](#setup)
+ [Storing the main table](#storing-the-main-table)
+ [Adding twowayfeweights output to esttab](#adding-twowayfeweights-output-to-esttab)

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

## Adding twowayfeweights output to esttab

Generally, if we want to store a regression output in .tex table, we would use the `esttab` package. Let's assume that we want to add the sum of the negative weights as a scalar to the regression table. The code below will show an easy procedure to accomplish this. Notice that this technique works for any scalar and it can be extended to multiple scalars.

First, run `twowayfeweights` with the same variables as the regression model and store the results in scalars:

```applescript
est clear
qui twowayfeweights Y g t D, type(feTR)
scalar neg_weights_1 = e(M)[2,2]
qui twowayfeweights Y g t D, type(feTR) controls(X)
scalar neg_weights_2 = e(M)[2,2]
```

Then, run the regressions and store the scalars in `esttab` using `estadd scalar`.

```applescript
reghdfe Y D, abs(g t)
estadd scalar neg_weights neg_weights_1
est sto model_1
reghdfe Y D X, abs(g t)
estadd scalar neg_weights neg_weights_2
est sto model_2
esttab model_* using "filename.tex", replace booktabs se scalar(neg_weights, label("Negative Weights")) standalone
```

The last line saves the regression results in a standalone TeX file. Notice the `s()` (scalar) argument and the corresponding `label` subargument.

---