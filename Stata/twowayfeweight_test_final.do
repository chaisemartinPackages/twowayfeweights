clear
local G = 2
local T = 6
qui do "twowayfeweights.ado"
set obs `=`G' * `T''
gen G = floor((_n-1)/`T')
bys G: gen T = _n
gen D = G == 0 & T >= 4 | G == 1 & T >= 5
gen E1 = G == 0 & T == 4 | G == 1 & T == 5
gen E2 = G == 0 & T == 5 | G == 1 & T == 6
gen E3 = G == 0 & T == 6
gen Y = uniform()
twowayfeweights Y G T E1, other_treatments(E2 E3) type(feTR)

clear
set seed 0
local G = 100
local T = 10
qui do "twowayfeweights.ado"

set obs `G'
gen g = _n
gen F_g = floor(2 + uniform() * (`T' - 2))
forv j = 1/`T' {
	gen t`j' = `j'
}
greshape long t, i(g F_g) j(j)
drop j
order g t F_g

gen D = t >= F_g
gen D2 = uniform() > 0.5
gegen D0 = mean(D), by(g t)
gen X = uniform() > 0.5
gen Y = uniform() + D*F_g
gen weights = uniform()
save test_data, replace

twowayfeweights Y g t D, type(feTR) summary_measures
twowayfeweights Y g t D D0, type(fdTR) summary_measures
twowayfeweights Y g t D, type(feS) summary_measures
twowayfeweights Y g t D, type(fdS) summary_measures
twowayfeweights Y g t D, type(feTR) other_treatments(D2) summary_measures
twowayfeweights Y g t D, type(feTR) controls(X) summary_measures
twowayfeweights Y g t D, type(feTR) test_random_weights(X) summary_measures
twowayfeweights Y g t D, type(feTR) weight(weights) summary_measures
twowayfeweights Y g t D, type(feTR) weight(weight) summary_measures
