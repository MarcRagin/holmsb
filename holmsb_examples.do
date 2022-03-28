***********************************************************************************************
* holmsb: p-value adjustment to control family-wise error rate when testing multiple hypotheses
* Copyright 2020, Marc Ragin, University of Georgia
* 
* 
* EXAMPLE CODE
***********************************************************************************************

*** Example 1: Linear regression, single model with multiple hypotheses on RHS
	sysuse auto, clear
	eststo clear
	
	* Conduct a linear regression
	eststo r1: regress mpg weight displacement foreign

	* Compute adjusted p-values for 2 hypotheses on the right-hand side (weight and displacement).
	holmsb weight displacement, estpref(r1)

	* Adjusted p-values are displayed in the output window and are saved in the active e().
	ereturn list /** Note e(psidak) and e(pbonf) */

*** Example 2: Linear regression, multiple models with the same hypothesis variables on RHS
	sysuse census, clear
	eststo clear
	
	* Conduct linear regressions and store estimates
	eststo r1: regress divorce marriage pop

	eststo r2: regress divorce marriage pop popurban death

	eststo r3: regress divorce marriage pop popurban death i.region

	* Compute adjusted p-values for 1 hypothesis variable on the right-hand side (marriage).
    holmsb pop, multiple
	
	* Adjusted p-values are displayed in the output window and are saved in the active e().
	ereturn list /** Note e(psidak) and e(pbonf) */
	
