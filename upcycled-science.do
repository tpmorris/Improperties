 *! Tim Morris | 14oct2019
*  Simulation study to accompany my brief article for upcycled science
*  requires Ben Jann's hexplot (ssc install heatplot)
*  and Ian White's simsum (ssc install simsum)
version 15
clear *

local run 0  // logical: set to 1 if you want the simulation to run. If you have already run it and have estimates data, set it to anything else
local nrep 10003          // number of repetitions
local numpg 100           // number per group
local numpats = 2*`numpg' // 2-arm trial
local alpha 0             // intercept, arbitrary
local gamma 3             // make effect of x very strong so results are clear
local beta 0              // treatment effect, arbitrary

set seed 89592

tempname elm // Any name. Elms are nice.
postfile `elm' int(repno) str10(rand) str6(xdgm) float(mean0 mean1) using means, replace

set coeftabresults off
if `run' {
quietly {
noi _dots 0, title("Simulation running...")
forval s = 1/`nrep' {
	noi _dots `s' 0
	* Generating data
	clear
	set obs `numpats'
    * Design methods
	egen byte xfixed  = seq(), from(0) to(1) block(`numpg') // two strata, fixed to 86 pts
	gen  byte xrandom = rbinomial(1,.5) // two strata, random with E(x1)=86
	gen  byte z_simple_fixed  = rbinomial(1,.5) // simple randomisation for fixed x
	gen  byte z_simple_random = z_simple_fixed // simple randomisation for random x. I've used exactly the same as above
	egen byte z_stratified_fixed = seq(), from(0) to(1) block(1) by(xfixed) // stratified blocks within strata (length doesn't matter as long as 86 is a multiple)
	egen byte z_stratified_random = seq(), from(0) to(1) block(1) by(xrandom) // stratified blocks within strata (length doesn't matter as long as 86 is a multiple)

	gen float y_fixed  = `alpha' + (`gamma'*xfixed ) + (`beta'*z_stratified_fixed ) + rnormal(0,1)
	gen float y_random = `alpha' + (`gamma'*xrandom) + (`beta'*z_stratified_random) + rnormal(0,1)

	foreach xdgm in fixed random {
		foreach rand in simple stratified {
			mean y_`xdgm' , over(z_`rand'_`xdgm')
			post `elm' (`s') ("`rand'") ("`xdgm'") (_b[y_`xdgm':0]) (_b[y_`xdgm':1])
		}
	}
}
}
postclose `elm'
}

use means, clear
gen float diff = mean1 - mean0
* tidy up labels for plots
lab var rand "Randomisation method"
lab var xdgm "Method for generating x"
	replace xdgm = "Fixed x"  if xdgm == "fixed"
	replace xdgm = "Random x" if xdgm == "random"
lab var mean0 "Mean of y | z=0"
lab var mean1 "Mean of y | z=1"
lab var diff "Difference in mean y between arms"

foreach xdgm in fixed random {
    foreach rand in simple stratified {
        quietly corr mean0 mean1 if xdgm=="`xdgm'" & rand=="`rand'"
        local zlow = 0.5*log( (1+r(rho))/(1-r(rho)) ) + invnormal(0.025)*(sqrt(r(N)-3)^-1)
        local zupp = 0.5*log((1+r(rho))/(1-r(rho))) + invnormal(0.975)*(sqrt(r(N)-3)^-1)
        local rholow = round( (exp(2*`zlow')-1)/(exp(2*`zlow')+1) , 0.01)
        local rhoupp = round( (exp(2*`zupp')-1)/(exp(2*`zupp')+1) , 0.01)
        local corr = round(r(rho), 0.01)
        display as text "`xdgm' x, `rand' randomisation:" _newline "  corr = " as result "`corr'" as text ", 95% Monte Carlo CI:" as result "(`rholow', `rhoupp')"
        local label_`xdgm'_`rand' "Corr = `corr' (95% Monte Carlo CI: `rholow', `rhoupp')"
    }
}


simsum diff, id(repno) by(xdgm rand) empse mcse

* Note: requires Ben Jann's hexplot (ssc install heatplot)
hexplot mean1 mean0, by(xdgm rand, note("")) xlab(,grid) aspect(1) name(corr, replace) ysize(7) xsize(8)
graph export corr.svg, name(corr) replace
graph export corr.pdf, name(corr) replace

twoway hist diff, width(.04) lc(gs8) fc(white) xlab(,grid glw(medium)) ylab(none) ytit("") by(xdgm rand, cols(1) note("")) name(diff, replace) ysize(7) xsize(8)
graph export diff.svg, replace
graph export diff.pdf, replace
