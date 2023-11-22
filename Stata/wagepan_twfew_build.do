clear
which bcuse
if _rc != 0 {
	ssc install bcuse, replace
}
bcuse wagepan
sort nr year
forvalue t=1981/1986{
	replace union=0 if year==`t'&union==1&union[_n-1]==0&union[_n+1]==0
	replace union=1 if year==`t'&union==0&union[_n-1]==1&union[_n+1]==1
}
xtset nr year
g diff_lwage = d.lwage
g diff_union = d.union
save "wagepan_twfeweights.dta", replace
