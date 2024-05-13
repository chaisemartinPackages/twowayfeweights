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