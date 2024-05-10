capture program drop twowayfeweights
program twowayfeweights, eclass
	version 12.0
	syntax varlist(min=4 numeric) [if] [in]  [, type(string) test_random_weights(varlist numeric) controls(varlist numeric) path(string) weight(varlist numeric) other_treatments(varlist numeric)  summary_measures]

	foreach p in gtools {
		qui cap which `p'
		if _rc != 0 {
			di as error "twowayfeweights requires the `p' package: " `"{stata ssc install gtools, replace}"'
			exit
		}
	}

	if "`other_treatments'"==""{
	
	if "`type'"==""{
	di as error"Please select the weights you want to estimate using type()"
	}
	
	if "`type'"!=""{
	
/// Preparing the data 	
	
	qui{

	tempvar outcome group time meantreat
	tokenize `varlist'
	gen `outcome'=`1'
	gen `group'=`2'
	gen `time'=`3'
	gen `meantreat'=`4'
	if "`type'"=="fdTR" {
	tempvar treatment
	gen `treatment'=`5'
	}
	
	preserve
	
	*Keeping if sample
	if `"`if'"' != "" {
	keep `if'
	}

	* Account for the presence of weight vars in the data
	if strpos("`weight'", "weight") > 0  { 
		cap rename `weight' var_weight
		local weight "var_weight"
	}
	cap confirm var weight
	if _rc == 0 {
		local p = 1
		foreach v of varlist weight* {
			cap rename `v' weight_OG`p'
			local p = `p' + 1
		}
	}

	if "`type'"!="fdTR" {
	* Keeping only sample used in estimation of regression
	foreach var of varlist `varlist' {
	drop if `var'==.
	}
	if "`controls'"!=""{
	foreach var of varlist `controls' {
	drop if `var'==.
	}
	}
	}
	
	if "`type'"=="fdTR" {
	* Keeping only sample used in estimation of regression & observations with non missing D
	keep if (`time'!=.&`outcome'!=.&`meantreat'!=.)|`treatment'!=.
	if "`controls'"!=""{
	foreach var of varlist `controls' {
	drop if (`time'!=.&`outcome'!=.&`meantreat'!=.)&`var'==.
	}
	}
	}

	*Replacing individual level treatment by (g,t)-level treatment
	capture drop treatment_gt treatment_sd_gt
	gegen treatment_gt=mean(`meantreat'), by(`group' `time')
	gegen treatment_sd_gt=sd(`meantreat'), by(`group' `time')
	sum treatment_sd_gt
	if r(mean)>0&r(mean)!=.{
	noisily di as text "The treatment variable in the regression varies within some group * period cells."
	noisily di as text "The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions" _newline "with a group * period level treatment."
	noisily di as text "The command will replace the treatment by its average value in each group * period."
	noisily di as text "The results below apply to the two-way fixed effects regression with that treatment variable."
	noisily di as text ""
	replace `meantreat'=treatment_gt
	}
	drop treatment_gt treatment_sd_gt	
	
	*Replacing individual level controls by (g,t)-level controls
	if "`controls'"!=""{
	local count=1
	foreach var of varlist `controls' {
	gegen `var'_gt=mean(`var'), by(`group' `time')
	gegen `var'_sd_gt=sd(`var'), by(`group' `time')
	sum `var'_sd_gt
	if r(mean)>0&r(mean)!=.{
	noisily di as text "The control variable " `count' " varies within some group * period cells."
	noisily di as text "The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) on two-way fixed effects regressions" _newline "with controls apply to group * period level controls."
	noisily di as text "The command will replace control variable " `count' " by its average value in each group * period cell."
	noisily di as text "The results below apply to the regression with control variable " `count' " averaged at the group * period level."
	noisily di as text ""
	replace `var'=`var'_gt
	local count=`count'+1
	}
	drop `var'_gt `var'_sd_gt
	}
	}
	
	*Creating the weight variable
	capture drop weight_XX
	if "`weight'"==""{
	gen weight_XX=1
	}
	if "`weight'"!=""{
	gen weight_XX=`weight'
	}
	keep if weight_XX!=.

/// Creating the weights variables 	

/// feTR weights
	
	if "`type'"=="feTR"{
	sum `meantreat' [aweight=weight_XX]
	scalar mean_D=r(mean)
	sum `outcome' [aweight=weight_XX]
	scalar nobs_XX=r(sum_w)
	gegen P_gt=total(weight_XX), by(`group' `time')
	replace P_gt=P_gt/nobs_XX 
	gen nat_weight= P_gt*`meantreat'/mean_D
	
	areg `meantreat' i.`time' `controls' [aweight=weight_XX], absorb(`group')
	predict eps_1, residuals 
	gen eps_1_E_D_gt=eps_1*`meantreat'
	sum eps_1_E_D_gt [aweight=weight_XX]
	scalar denom_W=r(mean)
	gen W=eps_1*mean_D/denom_W
	gen weight=W*nat_weight
	save test_data, replace

	*Computing beta
	
	areg `outcome' i.`time' `meantreat' `controls' [aweight=weight_XX], absorb(`group')
	scalar beta=_b[`meantreat']
	
	* Keeping only one observation in each group * period cell
	bys `group' `time': gen group_period_unit=(_n==1)	
	drop if group_period_unit==0
	drop group_period_unit
	
	}
	
/// fdTR weights

	if "`type'"=="fdTR"{

	sum `treatment' [aweight=weight_XX]
	scalar mean_D=r(mean)
	egen num_obs=total(weight_XX)
	sum num_obs
	scalar nobs_XX=r(mean)	
	gegen P_gt=total(weight_XX), by(`group' `time')
	replace P_gt=P_gt/nobs_XX 
	gen nat_weight= P_gt*`treatment'/mean_D
	reg `meantreat' i.`time' `controls' [aweight=weight_XX]
	predict eps_2, residuals
	// line below sets eps_2=0 for the first period when a group is observed, according to formula in paper
	replace eps_2=0 if eps_2==.
	
	*Computing beta
	
	areg `outcome' `meantreat' `controls' [aweight=weight_XX], absorb(`time')
	scalar beta=_b[`meantreat']

	* Keeping only one observation in each group * period cell
	bys `group' `time': gen group_period_unit=(_n==1)	
	drop if group_period_unit==0
	drop group_period_unit
	
	egen newt=group(`time')
	sort `group' newt
	gen w_tilde_2=eps_2-eps_2[_n+1]*P_gt[_n+1]/P_gt if `group'==`group'[_n+1]&newt+1==newt[_n+1]
	//the condition newt+1==newt[_n+1] above ensures that the computation is right when panel has holes (a group there from period 1 to t0, then disappears, and reappears at t0+k).
	// line below sets w_tilde_2=eps_2 for the last period when a group is observed, according to formula in paper
	replace w_tilde_2=eps_2 if w_tilde_2==.
	gen w_tilde_2_E_D_gt=w_tilde_2*`treatment'
	sum w_tilde_2_E_D_gt [aweight=P_gt]
	scalar denom_W=r(mean)
	gen W=w_tilde_2*mean_D/denom_W
	gen weight=W*nat_weight
	
	}
	
/// feS weights

	if "`type'"=="feS"{
	
	sum `outcome' [aweight=weight_XX]
	scalar nobs_XX=r(sum_w)
	gegen P_gt=total(weight_XX), by(`group' `time')
	replace P_gt=P_gt/nobs_XX

	areg `meantreat' i.`time' `controls' [aweight=weight_XX], absorb(`group')
	predict eps_1, residuals
	gegen newt=group(`time')
	replace newt=newt-1 
	gen eps_1_weight=eps_1*weight_XX
	// Modif. Diego: replace the loop with a reverse sorting cumulative sum //
	gsort `group' -newt
	bys `group': gen E_eps_1_g_geqt_aux = sum(eps_1_weight) 
	bys `group': gen weight_XX_aux = sum(weight_XX) 
	gen E_eps_1_g_geqt = E_eps_1_g_geqt_aux / weight_XX_aux 
	sort `group' newt
	drop eps_1_weight E_eps_1_g_geqt_aux weight_XX_aux

	*Computing beta
	areg `outcome' i.`time' `meantreat' `controls' [aweight=weight_XX], absorb(`group')
	scalar beta=_b[`meantreat']

	* Keeping only one observation in each group * period cell
	
	bys `group' `time': gen group_period_unit=(_n==1)
	drop if group_period_unit==0
	drop group_period_unit
		
	sort `group' `time'
	gen Delta_D=`meantreat'-`meantreat'[_n-1] if `group'==`group'[_n-1]&newt-1==newt[_n-1]
	*NB: the condition newt-1==newt[_n-1] ensures that the computation is right when panel has holes (a group there from period 1 to t0, then disappears, and reappears at t0+k).
	drop if Delta_D==. 
	gen s_gt=(Delta_D>0)-(Delta_D<0)
	gen abs_Delta_D=abs(Delta_D)
	drop Delta_D 
	gen nat_weight= P_gt*abs_Delta_D
	egen P_S=total(nat_weight)
	replace nat_weight=nat_weight/P_S
	
	gen om_tilde_1=s_gt*E_eps_1_g_geqt/P_gt
	sum om_tilde_1 [aweight=nat_weight]
	scalar denom_W=r(mean)
	gen W=om_tilde_1/denom_W
	gen weight=W*nat_weight

	}

//fdS weights

	if "`type'"=="fdS"{

	sum `outcome' [aweight=weight_XX]
	scalar nobs_XX=r(sum_w)
	gegen P_gt=total(weight_XX), by(`group' `time')
	replace P_gt=P_gt/nobs_XX
	
	reg `meantreat' i.`time' `controls' [aweight=weight_XX]
	predict eps_2, residuals
	
	*Computing beta
	areg `outcome' `meantreat' `controls' [aweight=weight_XX], absorb(`time')
	scalar beta=_b[`meantreat']
	
	* Keeping only one observation in each group * period cell	
	bys `group' `time': gen group_period_unit=(_n==1)
	drop if group_period_unit==0
	drop group_period_unit
	
	gen s_gt=(`meantreat'>0)-(`meantreat'<0)
	gen abs_Delta_D=abs(`meantreat')
	gen nat_weight= P_gt*abs_Delta_D
	egen P_S=total(nat_weight)
	replace nat_weight=nat_weight/P_S

	gen W=s_gt*eps_2
	sum W [aweight=nat_weight]
	scalar denom_W=r(mean)
	replace W=W/denom_W
	gen weight=W*nat_weight
	
	}

	// Modif Diego_weights: counting cells with non-zero treatment
	gen count_nat_weights = nat_weight != 0
	sum count_nat_weights
	scalar tot_cells = r(sum)
	drop count_nat_weights
	
/// Results	

	* Adjusting positive and negative weights that are close to 0
	local limit_sensitivity = 10^(-10)
	replace weight = 0 if abs(weight) < `limit_sensitivity' 

	* Computing the sum and the number of positive/negative weights
		
	egen total_weight_plus=total(weight) if weight>0&weight!=.
	egen total_weight_minus=total(weight) if weight<0
	sum total_weight_plus
	scalar nplus=r(N)
	scalar sumplus=r(mean)
    sum total_weight_minus
	scalar nminus=r(N)
	scalar summinus=r(mean)
	
	*Computing the sensitivity measure
	sum W [aweight=nat_weight]
	scalar sensibility=abs(beta)/r(sd)
	
	*Computing the number of weights
	scalar nweights=nplus+nminus
	
	* Regressing the variables in test_random_weights on the weights
	matrix A =0,0,0,0
	if "`test_random_weights'"!=""{
	foreach var of varlist `test_random_weights' {
	reg `var' W [pweight=nat_weight], cluster(`group') robust 
	matrix A =A\_b[W],_se[W],_b[W]/_se[W], ((_b[W]>=0)-(_b[W]<0))*sqrt(e(r2)) 
	}
	matrix B = A[2..., 1...]
	matrix colnames B = Coef SE t-stat Correlation
	matrix rownames B= `test_random_weights'
	}
		
	*Computing the new sensitivity measure
	if summinus<0{
	keep if weight!=0

	// Modif. Diego: change the loop with the cumulative sum
	gsort -W `group' `time' // Replicate the ordering	
	sum W	
	cap drop P_k
	gsort W -`group' -`time'
	gen P_k = sum(nat_weight)	
	cap drop S_k
	gsort W -`group' -`time'
	gen S_k = sum(weight)	
	cap drop T_k
	gsort W -`group' -`time'
	gen sq_weight = nat_weight * W^2
	gen T_k = sum(sq_weight)	
	drop sq_weight
	gsort -W `group' `time' // Replicate the ordering

	gen sens_measure2=abs(beta)/sqrt(T_k+S_k^2/(1-P_k))
	gen ind=(W<-S_k/(1-P_k))
	replace ind=0 if _n==1 
	// Filling holes
	replace ind=max(ind,ind[_n-1])
	// Count
	egen tot_ind=total(ind)
	sum tot_ind
	sum sens_measure2 if _n==r(N)-r(mean)+1
	scalar sensibility2=r(mean)
	}

	*Saving the results in a dataset, if requested

	if "`path'"!=""{
	gen Group_TWFE= `group'
	gen Time_TWFE=`time'
	keep Group_TWFE Time_TWFE weight
	save "`path'", replace
	}
	
	
  restore

*end of quietly condition
}
	
/// Displaying the results and saving them in e()
	if summinus == . {
		scalar summinus = 0
	}
	local ctitle = "ATT"

	local row_1 = ""
	fit_str , str("Treat. var: `4'") len(24) out(row_11) left
	local row_1 = "`row_1'" + r(row_11)
	fit_str , str("# `ctitle's") len(12) out(row_12) left
	local row_1 = "`row_1'" + r(row_12)
	fit_str , str("`=uchar(931)' weights") len(12) out(row_13) left
	local row_1 = "`row_1'" + r(row_13)

	local row_2 = ""
	fit_str , str("Positive weights") len(24) out(row_21) left
	local row_2 = "`row_2'" + r(row_21)
	fit_str , str("`: di %9.0f nplus'") len(12) out(row_22) left
	local row_2 = "`row_2'" + r(row_22)
	fit_str , str("`: di %9.4f sumplus'") len(12) out(row_23) left
	local row_2 = "`row_2'" + r(row_23)
	
	local row_3 = ""
	fit_str , str("Negative weights") len(24) out(row_31) left
	local row_3 = "`row_3'" + r(row_31)
	fit_str , str("`: di %9.0f nminus'") len(12) out(row_32) left
	local row_3 = "`row_3'" + r(row_32)
	fit_str , str("`: di %9.4f summinus'") len(12) out(row_33) left
	local row_3 = "`row_3'" + r(row_33)

	local row_4 = ""
	fit_str , str("Total") len(24) out(row_41) left
	local row_4 = "`row_4'" + r(row_41)
	fit_str , str("`: di %9.0f nweights'") len(12) out(row_42) left
	local row_4 = "`row_4'" + r(row_42)
	fit_str , str("`: di %9.4f `=summinus + sumplus''") len(12) out(row_43) left
	local row_4 = "`row_4'" + r(row_43)

	if inlist("`type'", "feTR", "fdTR") {
		local assumptions "Under the common trends assumption,"
	} 
	else if inlist("`type'", "feS", "fdS") {
		local assumptions "Under the common trends assumption and" _newline "the assumption that groups' treatment effects do not change over time,"
	}


	di ""
	di as text "`assumptions'"
	di as text "the TWFE coefficient beta, equal to `=strtrim("`: di %9.4f beta'")', estimates a weighted sum of " nplus + nminus " `ctitle's. " _newline nplus " `ctitle's receive a positive weight, and " nminus " receive a negative weight."
	if tot_cells > nplus + nminus {
		di as text `=tot_cells' " (g,t) cells receive the treatment, but the `ctitle's of " `=tot_cells - (nplus + nminus)' " cells receive a weight equal to zero."
	}
	di as result "{hline 48}"
	di as result "`row_1'"
	di as result "{hline 48}"
	di as result  "`row_2'"
	di as result  "`row_3'"
	di as text 48 * "-"
	di as result  "`row_4'"
	di as result  "{hline 48}"

	if "`summary_measures'" != "" {
		local subscr = substr("`type'", 1, 2)
		local srow_1 "Summary Measures:"
		local discl "Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)"
		local srow_fe "TWFE coefficient (`=uchar(946)'_`subscr') = `: di %9.4f beta'"
		local srow_2 "min `=uchar(963)'(`=uchar(916)') compatible with `=uchar(946)'_`subscr' and `=uchar(916)'_TR = 0: `: di %9.4f sensibility'"
		di as result ""
		di as result "`srow_1'"
		di as result "`srow_fe'"
		di as result "`srow_2'"
		if summinus<0{
			local srow_3 "min `=uchar(963)'(`=uchar(916)') compatible with treatment effect of opposite sign than `=uchar(946)'_`subscr' in all (g,t) cells: `: di %9.4f sensibility2'"					
			di as result "`srow_3'"
		}
		di as text "`discl'"
	}
	
	ereturn clear 
	matrix mret = J(3, 2, .)
	matrix mret[1,1] = nplus
	matrix mret[2,1] = nminus
	matrix mret[3,1] = nplus + nminus
	matrix mret[1,2] = sumplus
	matrix mret[2,2] = summinus
	matrix mret[3,2] = sumplus + summinus
	matrix coln mret = "N_ATTs" "Sum_weights"
	matrix rown mret = "Pos_Weights" "Neg_Weights" "Tot"
	ereturn matrix M = mret

	ereturn scalar lb_se_te = sensibility
	if summinus<0{
	ereturn scalar lb_se_te2 = sensibility2
	}
	ereturn scalar beta = beta
	ereturn scalar ot = 0

	if "`test_random_weights'"!=""{
	di  as result  _newline "Regression of variables possibly correlated with the treatment effect on the weights"
	matrix list B
	ereturn matrix randomweightstest1 = B
	}

	display _newline
	di as text "The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899)."
	
	}
	
	}

if "`other_treatments'"!=""{

	if "`type'"!="feTR"{
	di as error"When the other_treatments option is specified, you need to specify the type(feTR) option."
	}
	
	if "`type'"=="feTR"{
	
	/// Preparing the data 	
	
	qui{

	tempvar outcome group time meantreat
	tokenize `varlist'
	gen `outcome'=`1'
	gen `group'=`2'
	gen `time'=`3'
	gen `meantreat'=`4'
	
	preserve
	
	*Keeping if sample
	if `"`if'"' != "" {
	keep `if'
	}

	* Account for the presence of weight vars in the data * 
	if strpos("`weight'", "weight") > 0  { 
		cap rename `weight' var_weight
		local weight "var_weight"
	}
	cap confirm var weight
	if _rc == 0 {
		local p = 1
		foreach v of varlist weight* {
			cap rename `v' weight_OG`p'
			local p = `p' + 1
		}
	}
	
	* Keeping only sample used in estimation of regression
	foreach var of varlist `varlist' {
	drop if `var'==.
	}
	if "`controls'"!=""{
	foreach var of varlist `controls' {
	drop if `var'==.
	}
	}
	foreach var of varlist `other_treatments' {
	drop if `var'==.
	}	
	*Replacing individual level treatment by (g,t)-level treatment
	capture drop treatment_gt treatment_sd_gt
	gegen treatment_gt=mean(`meantreat'), by(`group' `time')
	gegen treatment_sd_gt=sd(`meantreat'), by(`group' `time')
	sum treatment_sd_gt
	if r(mean)>0&r(mean)!=.{
	noisily di as text "The treatment variable in the regression varies within some group * period cells."
	noisily di as text "The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions" _newline "with a group * period level treatment."
	noisily di as text "The command will replace the treatment by its average value in each group * period."
	noisily di as text "The results below apply to the two-way fixed effects regression with that treatment variable."
	noisily di as text ""
	replace `meantreat'=treatment_gt
	}
	drop treatment_gt treatment_sd_gt	

	*Replacing individual level other treatments by (g,t)-level treatment
	local count=1
	foreach var of varlist `other_treatments' {
	gegen `var'_gt=mean(`var'), by(`group' `time')
	gegen `var'_sd_gt=sd(`var'), by(`group' `time')
	sum `var'_sd_gt
	if r(mean)>0&r(mean)!=.{
	noisily di as text "The other treatment variable " `count' " varies within some group * period cells."
	noisily di as text "The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) on two-way fixed effects regressions" _newline "with several treatments apply to group * period level treatments."
	noisily di as text "The command will replace other treatment variable " `count' " by its average value in each group * period cell."
	noisily di as text "The results below apply to the regression with other treatment variable " `count' " averaged at the group * period level."
	noisily di as text ""
	replace `var'=`var'_gt
	local count=`count'+1
	}
	drop `var'_gt `var'_sd_gt
	}

	*Replacing individual level controls by (g,t)-level controls
	if "`controls'"!=""{
	local count=1
	foreach var of varlist `controls' {
	gegen `var'_gt=mean(`var'), by(`group' `time')
	gegen `var'_sd_gt=sd(`var'), by(`group' `time')
	sum `var'_sd_gt
	if r(mean)>0&r(mean)!=.{
	noisily di as text "The control variable " `count' " varies within some group * period cells."
	noisily di as text "The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) on two-way fixed effects regressions" _newline "with controls apply to group * period level controls."
	noisily di as text "The command will replace control variable " `count' " by its average value in each group * period cell."
	noisily di as text "The results below apply to the regression with control variable " `count' " averaged at the group * period level."
	noisily di as text ""
	replace `var'=`var'_gt
	local count=`count'+1
	}
	drop `var'_gt `var'_sd_gt
	}
	}
	
	*Creating the weight variable
	capture drop weight_XX
	if "`weight'"==""{
	gen weight_XX=1
	}
	if "`weight'"!=""{
	gen weight_XX=`weight'
	}
	keep if weight_XX!=.

/// Creating the weights variables 	
	
	sum `meantreat' [aweight=weight_XX]
	scalar mean_D=r(mean)
	sum `outcome' [aweight=weight_XX]
	scalar nobs_XX=r(sum_w)
	gegen P_gt=total(weight_XX), by(`group' `time')
	replace P_gt=P_gt/nobs_XX 
	gen nat_weight= P_gt*`meantreat'/mean_D
	areg `meantreat' i.`time' `controls' `other_treatments' [aweight=weight_XX], absorb(`group')
	predict eps_1, residuals 
	gen eps_1_E_D_gt=eps_1*`meantreat'
	sum eps_1_E_D_gt [aweight=weight_XX]
	scalar denom_W=r(mean)
	gen W=eps_1*mean_D/denom_W
	gen weight=W*nat_weight

	// Modif Diego_weights: counting cells with non-zero treatment
	gen count_nat_weights = nat_weight != 0
	sum count_nat_weights
	scalar tot_cells = r(sum)
	drop count_nat_weights

	local j=1
	foreach var of varlist `other_treatments' {
	gen weight_others`j'=W*P_gt*`var'/mean_D

	local limit_sensitivity = 10^(-10)
	replace weight_others`j' = 0 if abs(weight_others`j') < `limit_sensitivity' 

	gen count_`j' = `var' != 0
	sum count_`j'
	scalar tot_cells_`j' = r(sum)
	drop count_`j'

	local j=`j'+1
	}
	
	*Computing beta
	
	areg `outcome' i.`time' `meantreat' `controls' `other_treatments' [aweight=weight_XX], absorb(`group')
	scalar beta=_b[`meantreat']
	
	* Keeping only one observation in each group * period cell
	bys `group' `time': gen group_period_unit=(_n==1)	
	drop if group_period_unit==0
	drop group_period_unit

	* Adjust weights close to zero
	local limit_sensitivity = 10^(-10)
	replace weight = 0 if abs(weight) < `limit_sensitivity' 
	
	* Computing the sum and the number of positive/negative weights
		
	egen total_weight_plus=total(weight) if weight>0&weight!=.
	egen total_weight_minus=total(weight) if weight<0
	sum total_weight_plus
	scalar nplus=r(N)
	scalar sumplus=r(mean)
    sum total_weight_minus
	scalar nminus=r(N)
	scalar summinus=r(mean)
	drop total_weight_plus total_weight_minus
	scalar nweights=nplus+nminus
	
	local j=1
	foreach var of varlist `other_treatments' {
	egen total_weight_plus=total(weight_others`j') if weight_others`j'>0&weight_others`j'!=.
	egen total_weight_minus=total(weight_others`j') if weight_others`j'<0
	sum total_weight_plus
	scalar nplus_others`j'=r(N)
	scalar sumplus_others`j'=r(mean)
    sum total_weight_minus
	scalar nminus_others`j'=r(N)
	scalar summinus_others`j'=r(mean)
	scalar nweights_others`j'=nplus_others`j'+nminus_others`j'
	local j=`j'+1
	drop total_weight_plus total_weight_minus
	}
	
	* Regressing the variables in test_random_weights on the weights
	matrix A =0,0,0,0
	if "`test_random_weights'"!=""{
	foreach var of varlist `test_random_weights' {
	reg `var' W [pweight=nat_weight], cluster(`group') 
	matrix A =A\_b[W],_se[W],_b[W]/_se[W], ((_b[W]>=0)-(_b[W]<0))*sqrt(e(r2)) 
	}
	matrix B = A[2..., 1...]
	matrix colnames B = Coef SE t-stat Correlation
	matrix rownames B= `test_random_weights'
	}
	
	*Saving the results in a dataset, if requested

	if "`path'"!=""{
	capture drop Group 
	capture drop Time
	gen Group= `group'
	gen Time=`time'
	keep Group Time weight weight_others*
	save "`path'", replace
	}

  restore

*end of quietly condition
	}

	/// Displaying the results and saving them in e()
	if summinus == . {
		scalar summinus = 0
	}
	local row_1 = ""
	fit_str , str("Treat. var: `4'") len(24) out(row_11) left
	local row_1 = "`row_1'" + r(row_11)
	fit_str , str("# ATTs") len(12) out(row_12) left
	local row_1 = "`row_1'" + r(row_12)
	fit_str , str("`=uchar(931)' weights") len(12) out(row_13) left
	local row_1 = "`row_1'" + r(row_13)

	local row_2 = ""
	fit_str , str("Positive weights") len(24) out(row_21) left
	local row_2 = "`row_2'" + r(row_21)
	fit_str , str("`: di %9.0f nplus'") len(12) out(row_22) left
	local row_2 = "`row_2'" + r(row_22)
	fit_str , str("`: di %9.4f sumplus'") len(12) out(row_23) left
	local row_2 = "`row_2'" + r(row_23)
	
	local row_3 = ""
	fit_str , str("Negative weights") len(24) out(row_31) left
	local row_3 = "`row_3'" + r(row_31)
	fit_str , str("`: di %9.0f nminus'") len(12) out(row_32) left
	local row_3 = "`row_3'" + r(row_32)
	fit_str , str("`: di %9.4f summinus'") len(12) out(row_33) left
	local row_3 = "`row_3'" + r(row_33)

	local row_4 = ""
	fit_str , str("Total") len(24) out(row_41) left
	local row_4 = "`row_4'" + r(row_41)
	fit_str , str("`: di %9.0f nweights'") len(12) out(row_42) left
	local row_4 = "`row_4'" + r(row_42)
	fit_str , str("`: di %9.4f `=summinus + sumplus''") len(12) out(row_43) left
	local row_4 = "`row_4'" + r(row_43)

	di ""
	di as text "Under the common trends assumption," _newline "the TWFE coefficient beta, equal to `=strtrim("`: di %9.4f beta'")', estimates the sum of several terms."
	di ""
	di as text "The first term is a weighted sum of " nplus + nminus " ATTs of the treatment." _newline nplus " ATTs receive a positive weight, and " nminus " receive a negative weight."
	if tot_cells > nplus + nminus {
		di as text `=tot_cells' " (g,t) cells receive the treatment, but the ATTs of " `=tot_cells - (nplus + nminus)' " cells receive a weight equal to zero."
	}
	di as result "{hline 48}"
	di as result "`row_1'"
	di as result "{hline 48}"
	di as result  "`row_2'"
	di as result  "`row_3'"
	di as text 48 * "-"
	di as result  "`row_4'"
	di as result  "{hline 48}"

	ereturn clear 
	matrix mret = J(3, 2, .)
	matrix mret[1,1] = nplus
	matrix mret[2,1] = nminus
	matrix mret[3,1] = nplus + nminus
	matrix mret[1,2] = sumplus
	matrix mret[2,2] = summinus
	matrix mret[3,2] = sumplus + summinus
	matrix coln mret = "N_ATTs" "Sum_weights"
	matrix rown mret = "Pos_Weights" "Neg_Weights" "Tot"
	ereturn matrix M1 = mret

	local j=1
	foreach var of varlist `other_treatments' {
		foreach sign in plus minus {
			if sum`sign'_others`j' == . {
				scalar sum`sign'_others`j' = 0
			}
		}
		local row_1 = ""
		fit_str , str("Other treat.: `var'") len(24) out(row_11) left
		local row_1 = "`row_1'" + r(row_11)
		fit_str , str("# ATTs") len(12) out(row_12) left
		local row_1 = "`row_1'" + r(row_12)
		fit_str , str("`=uchar(931)' weights") len(12) out(row_13) left
		local row_1 = "`row_1'" + r(row_13)

		local row_2 = ""
		fit_str , str("Positive weights") len(24) out(row_21) left
		local row_2 = "`row_2'" + r(row_21)
		fit_str , str("`: di %9.0f nplus_others`j''") len(12) out(row_22) left
		local row_2 = "`row_2'" + r(row_22)
		fit_str , str("`: di %9.4f sumplus_others`j''") len(12) out(row_23) left
		local row_2 = "`row_2'" + r(row_23)
		
		local row_3 = ""
		fit_str , str("Negative weights") len(24) out(row_31) left
		local row_3 = "`row_3'" + r(row_31)
		fit_str , str("`: di %9.0f nminus_others`j''") len(12) out(row_32) left
		local row_3 = "`row_3'" + r(row_32)
		fit_str , str("`: di %9.4f summinus_others`j''") len(12) out(row_33) left
		local row_3 = "`row_3'" + r(row_33)

		local row_4 = ""
		fit_str , str("Total") len(24) out(row_41) left
		local row_4 = "`row_4'" + r(row_41)
		fit_str , str("`: di %9.0f nweights_others`j''") len(12) out(row_42) left
		local row_4 = "`row_4'" + r(row_42)
		fit_str , str("`: di %9.4f `=summinus_others`j' + sumplus_others`j'''") len(12) out(row_43) left
		local row_4 = "`row_4'" + r(row_43)

		di ""
		di as text "The next term is a weighted sum of " nplus_others`j' + nminus_others`j' " ATTs of treatment " `j' " included in the other_treatments option." _newline nplus_others`j' " ATTs receive a positive weight, and " nminus_others`j' " receive a negative weight."
		if tot_cells_`j' >  nplus_others`j' + nminus_others`j' {
			di as text `=tot_cells_`j'' " (g,t) cells receive the treatment, but the ATTs of " `=tot_cells_`j' - (nplus_others`j' + nminus_others`j')' " cells receive a weight equal to zero."
		}
		di as result "{hline 48}"
		di as result "`row_1'"
		di as result "{hline 48}"
		di as result  "`row_2'"
		di as result  "`row_3'"
		di as text 48 * "-"
		di as result  "`row_4'"
		di as result  "{hline 48}"

		matrix mret = J(3, 2, .)
		matrix mret[1,1] = nplus_others`j'
		matrix mret[2,1] = nminus_others`j'
		matrix mret[3,1] = nplus_others`j' + nminus_others`j'
		matrix mret[1,2] = sumplus_others`j'
		matrix mret[2,2] = summinus_others`j'
		matrix mret[3,2] = sumplus_others`j' + summinus_others`j'
		matrix coln mret = "N_ATTs" "Sum_weights"
		matrix rown mret = "Pos_Weights" "Neg_Weights" "Tot"
		local j=`j'+1
		ereturn matrix M`j' = mret
	}
	
	ereturn scalar beta = beta
	ereturn scalar ot = `j' - 1
	if "`test_random_weights'"!=""{
	di  as result  _newline "Regression of variables possibly correlated with the treatment effect on the weights attached to the treatment"
	matrix list B
	ereturn matrix randomweightstest1 = B
	}

	display _newline
	di as text "The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899)."		
	}	
	}
	
end

cap program drop fit_str
program define fit_str, rclass
syntax , str(string) len(integer) out(string) [left]
{
if "`left'" == "" {
	local n_str = (`len' - length(abbrev("`str'", `len')))* " " + abbrev("`str'", `len')
}
else {
	local n_str = abbrev("`str'", `len') + (`len' - length(abbrev("`str'", `len')))* " " 
}
}
return local `out' = "`n_str'"
end

cap program drop twowayfeweights_out
program define twowayfeweights_out, rclass
version 12.0
syntax , saving(string) [standalone]
cap confirm scalar e(ot)
if _rc != 0 {
	di as err "twowayfeweights output not found in ereturn"
	exit
}
cap file close texout
file open texout using `saving', write replace
local nl = char(10)
if "`standalone'" != "" {
	file write texout "\documentclass{standalone}`nl'"
	file write texout "\begin{document}`nl'"
}
file write texout "\begin{tabular}{lcc}`nl'\hline`nl'"
file write texout "\hline`nl' & N. ATTs & Sum weights \\"	
local ot = e(ot)
if `ot' == 0 {
	file write texout "`nl'\hline`nl'"
	file write texout "Positive weights & `:di %9.0gc e(M)[1,1]' &  `:di %9.4gc e(M)[1,2]' \\ `nl'"
	file write texout "Negative weights & `:di %9.0gc e(M)[2,1]' &  `:di %9.4gc e(M)[2,2]' \\ `nl'"
	file write texout "\hline`nl'Total & `:di %9.0gc e(M)[3,1]' &  `:di %9.4gc e(M)[3,2]' \\ `nl'"
}
if `ot' != 0 {
	forv j = 1/`=`ot'+1' {
		file write texout "\hline `nl'"
		if `j' == 1 file write texout "\multicolumn{3}{l}{\emph{Main Treatment}} \\ `nl'\hline `nl'"
		else file write texout "&& \\ `nl'\multicolumn{3}{l}{\emph{Other Treatment `=`j'-1'}} \\ `nl'\hline `nl'"
		file write texout "Positive weights & `:di %9.0gc e(M`j')[1,1]' &  `:di %9.4gc e(M`j')[1,2]' \\ `nl'"
		file write texout "Negative weights & `:di %9.0gc e(M`j')[2,1]' &  `:di %9.4gc e(M`j')[2,2]' \\ `nl'"
		file write texout "\hline`nl'Total & `:di %9.0gc e(M`j')[3,1]' &  `:di %9.4gc e(M`j')[3,2]' \\ `nl'"
	}
}
file write texout "\hline\hline`nl'\end{tabular}`nl'"
if "`standalone'" != "" {
	file write texout "\end{document}"
}
file close texout
di as text "Output saved in `saving'."
end