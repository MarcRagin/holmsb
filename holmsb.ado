* Holm SB 1.0.0 by Marc Ragin
* January 2020

cap program drop holmsb
program define holmsb, eclass
version 14
set varabbrev off

* Syntax
syntax varlist(min=1 numeric) [, multiple estpref(name)]
* VARLIST: hypothesis vars only, not ordered
	local hypvars `varlist'
	local ests `estpref'
	di "`hypvars'"
	di "`ests'"

	local nhyps: word count `hypvars'
	if "`estpref'" != "" & !mi("`multiple'") {
		di as error "Stored estimates should not be specified with estpref() when multiple regressions are being tested."
		exit 198
	}

	
	* Syntax 1 (default): Single regression with multiple RHS hypothesis vars
	if mi("`multiple'") {
		* Save original data
		qui tempfile original
		qui save "`original'", replace
	
		* Calculate Bonferroni and Sidak-adjusted p-values
		matrix pbonf_mat = e(b)
		matrix psidak_mat = e(b)
		parmest, norestore
		* Create indicator for hypothesis vars and put at top
		gen hyp = 0
		foreach h in `hypvars' {
			replace hyp = 1 if regexm(parm, "`h'") == 1
		}
		gsort -hyp p
		* Set temporary variable h for Holm's algorithm (WY 1993, p. 64-65)
			* k = num hypotheses
				qui su hyp
				tempname k
				qui gen `k' = `r(sum)'
			* j = current hypothesis in algorithm
				tempname j
				qui gen `j' = _n
			* h = k - j + 1
				tempname h
				qui gen `h' = `k'-`j'+1
				* Bonferroni: pbonf = p*h
				qui gen double 	pbonf = min(p*`h',1) if `j' == 1
				qui replace    	pbonf = min(max(p*`h',pbonf[_n-1]),1) if `j' > 1 & `j' <= `k'
				qui replace		pbonf = p if hyp == 0
				qui replace		pbonf = 99 if p == .
				* Sidak: psidak = 1-[(1-p)^h]
				qui gen double 	psidak = min((1-(1-p)^(`h')),1) if `j' == 1
				qui replace    	psidak = min(max((1-(1-p)^(`h')),psidak[_n-1]),1) if `j' > 1 & `j' <= `k'
				qui replace		psidak = p if hyp == 0
				qui replace		psidak = 99 if p == .
		* Create matrix of adjusted p-values for each adjustment method
		local vars: colnames e(b)
		local pvals "pbonf psidak"
		foreach p in `pvals' {
			foreach v in `vars' {
				qui su `p' if parm == "`v'"
				local m = cond(`r(mean)' == 99, ., `r(mean)')
				matrix `p'_mat[1, colnumb(`p'_mat,"`v'")] = `m'
			}
			* Add adjusted p-values to stored estimates
			estadd matrix `p' = `p'_mat: `ests'
		}
		list parm hyp estimate stderr p pbonf psidak
		use `original', clear
	}
	
	* Syntax 2: Multiple regressions with different outcomes and one RHS hypothesis var
	if !mi("`multiple'") {
		
		* Save original data
		qui tempfile original
		qui save "`original'", replace
		
		* Create an empty dataset to save parameter values
		qui tempfile parms
		parmest, norestore
		qui drop if _n > 0
		qui save `parms', replace
		
		* Save parameter estimates for hypothesis var in each regression
		qui estimates dir
		local regs "`r(names)'"
		qui foreach r in `regs' {
			estimates res `r'
			parmest, norestore
			gen model = "`r'"
			gen depvar = "`e(depvar)'"
			gen hyp = 0
			foreach h in `hypvars' {
				replace hyp = 1 if regexm(parm, "`h'") == 1 & strlen(parm) == strlen("`h'")
			}
			drop if hyp == 0
			append using `parms'
			save `parms', replace
		}
		gsort p
		* Set temporary variable h for Holm's algorithm (WY 1993, p. 64-65)
			* k = num models
				tempname k
				qui gen `k' = _N
			* j = current hypothesis in algorithm
				tempname j
				qui gen `j' = _n
			* h = k - j + 1
				tempname h
				qui gen `h' = `k'-`j'+1
				* Bonferroni: pbonf = p*h
				qui gen double 	pbonf = min(p*`h',1) if `j' == 1
				qui replace    	pbonf = min(max(p*`h',pbonf[_n-1]),1) if `j' > 1 & `j' <= `k'
				qui replace		pbonf = p if hyp == 0
				qui replace		pbonf = 99 if p == .
				* Sidak: psidak = 1-[(1-p)^h]
				qui gen double 	psidak = min((1-(1-p)^(`h')),1) if `j' == 1
				qui replace    	psidak = min(max((1-(1-p)^(`h')),psidak[_n-1]),1) if `j' > 1 & `j' <= `k'
				qui replace		psidak = p if hyp == 0
				qui replace		psidak = 99 if p == .
				qui save `parms', replace
		* Create matrix of adjusted p-values for each adjustment method
		qui foreach r in `regs' {
			use `parms', clear
			drop if model != "`r'"
			local h = parm[1]
			estimates restore `r'
			estimates replay `r'
			mat A = r(table)
			local pvals "pbonf psidak"
			foreach p in `pvals' {
				cap mat drop `p'_mat
				matrix `p'_mat = A[rownumb(A,"pvalue"), 1...]
				su `p'
				local m = `r(min)'
				matrix `p'_mat[1, colnumb(`p'_mat,"`h'")] = `m'
				* Add adjusted p-values to stored estimates
				estadd matrix `p' = `p'_mat
			}
		}
		qui use `parms', clear
		qui sort model
		list model depvar parm hyp estimate stderr p pbonf psidak
		use `original', clear
	}

end	
	/* Syntax 1: one model with multiple outcomes 
	
	
	* Syntax 2: varying models with multiple outcomes/subgroups. 
	if "`outcome_vars'"=="" {
	
		* Perform a full trim to remove leading and trailing spaces from the cmd() option
		mata: st_local("cmd",strtrim(st_local("cmd")))
	
		* If user did NOT use compound double quotes in cmd(), pass through the string asis. This ensures the -tokenize- command below works properly.
		mata: if( substr(st_local("cmd"),1,1)!=char(96) ) stata("syntax, cmd(string asis) *");;
	}

* Determine whether original estimates were stored for later tabulation using eststo command
if length("`e(_estimates_name)'") == 0 {
	local stored = 0
	eststo null
}
else if length("`e(_estimates_name)'") != 0 {
	local stored = 1
	local estname "`e(_estimates_name)'"
	eststo `estname'
}

qui tempfile original
qui save "`original'", replace

* Store original regression command as nullreg
local nullreg "`e(cmdline)'"

* Count number of hypotheses to be tested
local K: word count `varlist'
di "Num hypotheses: `K'"

* Store dependent variable from original regression as depvar
local depvar "`e(depvar)'"

* Set options for bootstrapping
if length(`"`seed'"')!=0 set seed `seed'
local bopts
if length(`"`strata'"')!=0  local bopts `bopts' strata(`strata')
if length(`"`cluster'"')!=0 local bopts `bopts' cluster(`cluster')

* Store whether output window will display RHS variable names (varlab = 0) or labels (varlab = 1)
local varlab = 0
if length(`"`label'"') != 0 local varlab = 1

* Throw errors in certain cases
if length(`"`e(cmd)'"') == 0 {
	dis as error "This is a postestimation command and there are no estimates in memory. Please run a regression."
	exit
}
if "`e(vce)'" == "cluster" & length(`"`cluster'"')==0 {
	dis as error "If you cluster your standard errors, you should also cluster your bootstrap samples."
	dis as error "Please try again, specifying the cluster option."
	exit
}
if strpos("`nullreg'", "i.") != 0 | strpos("`nullreg'", "#") != 0 {
	dis as error "Factor variables, time series operators, and interactions not allowed. Please manually create these variables."
	dis as error "tabulate x, gen(x_f) will create factor variables."
	exit
}
if _rc!=0 {
	dis as error "Your original linear regression does not work."
	dis as error "Please test the regression and try again."
	exit _rc
}

local j=0

* Create temporary dataset to hold adjusted p-values
tempfile bs

* Capture all RHS variables (excluding constant)
local rhs_orig : colfullnames e(b)
local rhs_orig = subinstr("`rhs_orig'", "o.", "", .)
local rhs_orig = subinstr("`rhs_orig'", " _cons", "", .)
local M: word count `rhs_orig'
di "Num RHS vars: `M'"

* Create dummy matrices with variable names to eventually store adjusted p-values
matrix porig_mat = e(b)
matrix pwy_resid_mat = e(b)
matrix pwy_reif_mat = e(b)
matrix pbonf_mat = e(b)
matrix psidak_mat = e(b)

* Define residual regression, replacing Y with e
qui capture drop ehat
qui predict ehat, r

* Determine what's after the RHS variables in the original regression command
* (Could be if/in or just the comma to set the options)
local ifloc = cond(strpos("`nullreg'", " if ") != 0, strpos("`nullreg'", " if ") - 1, 99999)
local inloc = cond(strpos("`nullreg'", " in ") != 0, strpos("`nullreg'", " in ") - 1, 99999)
local comma = cond(strpos("`nullreg'", ",") != 0, strpos("`nullreg'", ",") - 1, 99999)
local min_ifincomma = min(`ifloc', `inloc', `comma') + 1
local nullreg_end = substr("`nullreg'", `min_ifincomma', .)


* Capture matrix of coeffs, p-values, etc.
matrix mat = r(table)
matrix porig_mat = mat[4,1...]
* For variables omitted due to collinearity, assign them a p-value equal to their column number (for sorting later)
foreach var of varlist `rhs_orig' {
	local pval = mat[4, colnumb(mat,"`var'")] 
	if `pval' == . mat mat[4,colnumb(mat,"`var'")] = colnumb(mat,"`var'")
}

********************************************************************************
* Reorder hypothesis variables based on p-values 							   *
********************************************************************************
	* Create matrix with rows as variable names and cols as estimation stats
	matrix mat2 = mat'
	* Create hypothesis variable dummy in column 9 
	* (1 if hypothesis var, 0 if control, -1 if constant term)
		* Set all to 0
		foreach var of varlist `rhs_orig' {
			local `var'1 = 0
			mat mat2[rownumb(mat2,"`var'"), 9] = 0
		}
		* Hypothesis vars
		foreach var of varlist `varlist' {
			* Create dummy for hypothesis vars in results matrix
			mat mat2[rownumb(mat2,"`var'"), 9] = 1
			* Create dummy for hypothesis vars as local
			local `var'1 = 1
		}
		* Constant term
		mat mat2[rownumb(mat2,"_cons"), 9] = -1

	* Sort matrix on p-values by hypothesis vars vs. controls
		* Sort descending on hypothesis dummy (col 9)
		* Then sort ascending on p-value (col 4)
	mata: st_replacematrix("mat2", sort(st_matrix("mat2"), (-9,4)))
	* Determine order of p-values
	foreach var of varlist `rhs_orig' {
		local `var'_ord = 0
		local `var'_b = mat[1,	colnumb(mat,"`var'")]
		local `var'_p = mat[4,	colnumb(mat,"`var'")]
		forvalues i = 1/`M' {
			* Set temp`i' as p-value from column `i'
			local temp`i'_b = mat2[`i', 1]
			local temp`i'_p = mat2[`i', 4]
			* Replace var_ord = `i' if it has the same beta and p-value as in the original matrix
			if `temp`i'_p' != . local `var'_ord = cond(`temp`i'_p' == ``var'_p' & `temp`i'_b' == ``var'_b', `i', ``var'_ord')
		}
		*di "`var'_ord = ``var'_ord'"
		local `var'_p = round(``var'_p', 0.0001)
		*di "`var'_p: ``var'_p'"
	}

	* Create copies of each RHS var, ordered by hyp vs. control and ascending p-value
	local rhs ""
	local cvars ""
	local varlist2 ""
	qui forvalues i = 1/`M' {
		foreach var of varlist `rhs_orig' {
			if ``var'_ord' == `i' & ``var'1' == 1 {
				tempvar x`i'
				gen `x`i'' = `var'
				la var `x`i'' "`var'  (p = ``var'_p')"
				local xname`i' "`var'"
				local xlab`i': var label `var'
				* VARLIST2: hypothesis vars only, ascending order of p-values
				local varlist2 "`varlist2' `var'"
				* RHS: varlist of all RHS vars
				local rhs "`rhs' `x`i''"
			}
			else if ``var'_ord' == `i' & ``var'1' == 0 {
				tempvar c`i'
				gen `c`i'' = `var'
				if ``var'_p' <= 1 la var `c`i'' "`var' (p = ``var'_p')"
				if ``var'_p' > 1 la var `c`i'' "`var' (omitted)"
				* RHS: varlist of all RHS vars
				local rhs "`rhs' `c`i''"
				local cvars "`cvars' `c`i''"
				local cname`i' "`var'"
			}
		}
	}
	*if length(`"`verbose'"')!=0 di "RHS (all rhs vars, renamed as x and c): `rhs'"
	if length(`"`verbose'"')!=0 di "Hypothesis vars, ascending order of p-values: `varlist2'"

	* Setup syntax of residual regression
	local residreg = `"`e(cmd)'"' + " ehat "+ "`rhs'" + "`nullreg_end'"
	if length(`"`verbose'"')!=0 di "Residual regression: `residreg'" _newline(2)

	* Setup syntax of original regression
	local origreg = `"`e(cmd)'"' + " `depvar' "+ "`rhs'" + "`nullreg_end'"
	if length(`"`verbose'"')!=0 di "Original regression: `origreg'"

********************************************************************************


* Conduct original regression using reordered variables, and save regression results to scalars
qui `origreg'
matrix orig = r(table)
* Save estimates for hypothesis vars
	qui forvalues j = `K'(-1)1 {
		* Save beta, se, p, t, and n for each hypothesis variable
		scalar b_`j' 	= orig[rownumb(orig,"b"),		colnumb(orig,"`x`j''")]
		scalar se_`j' 	= orig[rownumb(orig,"se"),		colnumb(orig,"`x`j''")]
		scalar p_`j' 	= orig[rownumb(orig,"pvalue"),	colnumb(orig,"`x`j''")]
		scalar t_`j' 	= abs(orig[rownumb(orig,"t"),	colnumb(orig,"`x`j''")])
		scalar n_`j' 	= e(N)-e(rank)
	}
* Save estimates for control vars
	local k = `K'
	if `M' > `K' + 1 {
		foreach v of varlist `cvars' {
			local ++k
			scalar bc_`k' 	= orig[rownumb(orig,"b"),		colnumb(orig,"`c`k''")]
			scalar sec_`k' 	= orig[rownumb(orig,"se"),		colnumb(orig,"`c`k''")]
			scalar pc_`k' 	= orig[rownumb(orig,"pvalue"),	colnumb(orig,"`c`k''")]
			scalar tc_`k' 	= abs(orig[rownumb(orig,"t"),	colnumb(orig,"`c`k''")])
		}
	}
* Save estimates for constant term
	local k = `k' + 1
	scalar bc_`k' 	= orig[rownumb(orig,"b"),		colnumb(orig,"_cons")]
	scalar sec_`k' 	= orig[rownumb(orig,"se"),		colnumb(orig,"_cons")]
	scalar pc_`k' 	= orig[rownumb(orig,"pvalue"),	colnumb(orig,"_cons")]
	scalar tc_`k' 	= abs(orig[rownumb(orig,"t"),	colnumb(orig,"_cons")])
	local c`k' "_cons"

* Set ctop as total number of vars, including constant
local ctop = `k'
* scalar list
matrix drop orig
matrix drop mat

* Create a dummy for obs included in original regression sample
gen origsample = e(sample)

********************************************************************************
* Run bootstrapped regressions												   *
********************************************************************************
dis "Running `reps' bootstrap replications for each variable.  This may take some time." _newline(1)
forvalues i = 1/`reps' {
	* Display status of bootstrapping replications
		local pct = (`i'/`reps')*100
		if `i' == 1 dis "Bootstrapping:" _continue
		if mod(`pct',10) == 0 & `i' != `reps' dis "...`pct'%" _continue
		if `i' == `reps' dis "...100%. DONE." _newline(1)
	* Preserve original sample
	preserve
	* Take a sample of size N (the default) with replacement
	qui drop if origsample != 1
    bsample, `bopts'
	* Run regression for bootstrap sample
	qui `residreg'
	matrix mat = r(table)
	* Capture each beta, se, p, and n
    forvalues j = `K'(-1)1 {
		local l = `j' + 1
		* WYOUNG1: From book
			* Determine p-value in residual regression
			if `j'==`K' local qstar_`j' = mat[rownumb(mat,"pvalue"), colnumb(mat,"`x`j''")]
			if `j'!=`K' local qstar_`j' = min(`qstar_`l'',mat[rownumb(mat,"pvalue"), colnumb(mat,"`x`j''")])
			* Create local dummy for qstar le p
			local wy_`j' = cond(`qstar_`j'' <= p_`j', 1, 0, .)
	}
	qui `origreg'
	forvalues j = `K'(-1)1 {
		* WYOUNG2: From Julian Reif wyoung Stata command
			* Test bootstrapped beta against null for Westfall-Young
			qui test _b[`x`j''] = b_`j'
			local pstar_`j' = r(p)
			* Create local dummy for pstar le p
			local wy2_`j' = cond(`pstar_`j'' <= p_`j', 1, 0, .)
    }
	matrix drop mat
	* Store results from each bootstrap in dataset bs
	drop _all
	qui set obs `K'
	qui gen i = `i'
	qui gen k = _n
	qui gen qstar = .
	qui gen pstar = .
	qui gen p = .
	qui gen pwy_resid = .
	qui gen pwy_reif = .
	qui gen var = ""
	qui gen varlab = ""
	qui forval k = 1/`K' {
		replace qstar 		= `qstar_`k'' 	if k == `k'
		replace pstar 		= `pstar_`k'' 	if k == `k'
		replace p  			= p_`k' 		if k == `k'
		replace pwy_resid 	= `wy_`k'' 		if k == `k'
		replace pwy_reif 	= `wy2_`k'' 	if k == `k'
		replace var 		= "`xname`k''" 	if k == `k'
		replace varlab 		= "`xlab`k''" 	if k == `k'
	}
	if `i' > 1 append using "`bs'"
	qui save "`bs'", replace
	* Restore original dataset
	restore
}

* After bootstrapping complete, call temp dataset
use "`bs'", clear
* save "C:\Users\mragin\Desktop\bs", replace

*-------------------------------------------------------------------------------
* WESTFALL & YOUNG
*-------------------------------------------------------------------------------

* Steps 3 and 4. Calculate step-down and single-step Westfall-Young adjusted p-value
collapse (mean) pwy_resid pwy_reif p (max) i (firstnm) var varlab, by(k)

* Step 5. Enforce monotonicity using successive maximization.  Include k in the sort to break ties.
sort p k
qui cap replace pwy_resid  	= max(pwy_resid[_n-1], pwy_resid)	if _n > 1 & _n <= `K'
qui cap replace pwy_reif  	= max(pwy_reif[_n-1], pwy_reif)		if _n > 1 & _n <= `K'

* Fill in coeffs and standard errors from saved scalars
qui gen double coef = .
qui gen double stderr = .
qui forval k = 1/`K' {
	replace coef = b_`k'		if k==`k'
	replace stderr = se_`k' 	if k==`k'
}

*-------------------------------------------------------------------------------
* Holm-Bonferroni and Holm-Sidak step-down corrections
* (verbatim from Reif's wyoung program)
*-------------------------------------------------------------------------------
tempname j
qui gen `j' = _N-_n+1

qui gen double pbonf = min(p*`j',1) if _n==1
qui replace    pbonf = min(max(p*`j',pbonf[_n-1]),1) if _n>1 & _n <= `K'

qui gen double psidak = min((1-(1-p)^(`j')),1) if _n==1
qui replace    psidak = min(max((1-(1-p)^(`j')),psidak[_n-1]),1) if _n>1 & _n <= `K'

label var pbonf 		"Bonferroni-Holm p-value"
label var psidak 		"Sidak-Holm p-value"
label var stderr 		"Unadjusted standard error"
label var p  			"Unadjusted p-value"
label var pwy_resid   	"Westfall-Young p-value (residual regression)"
label var pwy_reif   	"Westfall-Young p-value (Reif Wald test)"

// assert psidak<=pbonf+0.00000000001
// foreach v of varlist p* {
// 	assert `v' <= 1 if `v' != 99
// }

* Add in betas, SEs, and p-values for controls
qui set obs `M'
qui replace k = _n if k > `K'
local L = `K' + 1
qui forval k = `L'/`ctop' {
	replace var = "`cname`k''" 		if k==`k'
	replace coef = bc_`k'		if k==`k'
	replace stderr = sec_`k' 	if k==`k'
	replace p = pc_`k' 			if k==`k'
	replace p = 99 				if p == .
	replace pwy_resid = pc_`k'	if k==`k'
	replace pwy_resid = 99		if p == 99
	replace pwy_reif = pc_`k'	if k==`k'
	replace pwy_reif = 99		if p == 99
	replace pbonf = pc_`k'		if k==`k'
	replace pbonf = 99			if p == 99
	replace psidak = pc_`k'		if k==`k'
	replace psidak = 99			if p == 99
}	

* Put back in original order
gen ord = 0
local i = 0
foreach v in `rhs_orig' {
	local ++i
	qui replace ord = `i' if var == "`v'"
}

* Create matrix of adjusted p-values for each adjustment method
* (will eventually output this to estimation matrix)
local pvals "pwy_resid pwy_reif pbonf psidak"
foreach v in `rhs_orig' {
	foreach p in `pvals' {
		qui su `p' if var == "`v'"
		local m = cond(`r(mean)' == 99, ., `r(mean)')
		matrix `p'_mat[1, colnumb(`p'_mat,"`v'")] = `m'
		qui replace `p' = . if var == "`v'" & stderr == .
	}
	qui replace p = . if var == "`v'" & stderr == .
}
qui replace stderr = .a if stderr == .
la de omit .a "(omitted)"
la val stderr omit

* Display matrix of betas, se's, and p-values in output window
order k var varlab coef stderr p pwy_resid pwy_reif pbonf psidak
format coef stderr p pwy_resid pwy_reif pbonf psidak %9.4f
sort ord
drop `j' ord
di "{bf:HYPOTHESIS VARIABLES: ADJUSTED P-VALUES}"
if `varlab' == 0 list var coef stderr p pwy_resid pwy_reif pbonf psidak if k <= `K', sep(10) noobs ab(16) str(20)
else if `varlab' == 1 list varlab coef stderr p pwy_resid pwy_reif pbonf psidak if k <= `K', sep(10) noobs ab(16) str(20)

* Drop variables created
use "`original'", clear
cap drop ehat
cap drop origsample

if `stored' == 0 estimates res null
else if `stored' == 1 estimates res `estname'

local pvals "pwy_resid pwy_reif pbonf psidak porig"
foreach p in `pvals' {
	qui estadd matrix `p' = `p'_mat
}

if `stored' == 0 estimates store null
else if `stored' == 1 estimates store `estname'

end
*/
