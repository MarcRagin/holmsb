{smcl}
{title:Title}

{p 4 4 2}{cmd:holmsb} {hline 2} A postestimation command to control the family-wise error rate using the Sidak and Bonferroni corrections with the Holm step-down algorithm when performing hypothesis tests.{p_end}


{title:Syntax}

{p 8 16 2}{opt holmsb} {varlist} [{cmd:,} {it:options}]

{synoptset}{...}
{synopthdr}
{synoptline}
{synopt :{opt multiple}}specifies that p-values of the coefficients on {cmd:(}{it:{help varlist:varlist}}{cmd:)} should be adjusted across multiple regressions.{p_end}
{synopt :{opt estpref()}}allows the user to specify the stored estimates in which to output adjusted p-values. Default is active estimates.{p_end}
{synoptline}
{p 4 6 2}

{title:Description}

{p 4 4 2}{cmd:holmsb} calculates Sidak-adjusted and Bonferroni-adjusted {it:p}-values using the free step-down methodology of Holm (1979). It outputs the adjusted p-values to estimates stored in memory.

{marker methods}{...}
{title:Methods}

{pstd} All adjustment methods begin with the {it:J} hypotheses sorted by unadjusted p-values so that {it:p(1)<p(2)<...<p(J)}. 

{phang}The Sidak-Holm adjusted {it:p}-values are calculated as {it:{1-(1-p(1))^J, max[p(1),1-(1-p(2))^(J-1)],..., max[p(J-1),p(J)]}}.
If the calculation yields a value larger than 1, then the adjusted {it:p}-value is set equal to 1.

{phang}The Bonferroni-Holm adjusted {it:p}-values are calculated as {it:{p(1)*J, max[p(1),p(2)*(J-1)],..., max[p(J-1),p(J)]}}.

{marker examples}{...}
{title:Example:  linear regression, single model with multiple hypotheses on RHS}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto}{p_end}

{pstd}Conduct a linear regression{p_end}
{phang2}{cmd:. regress mpg weight displacement foreign}{p_end}

{pstd}Compute adjusted p-values for 2 hypotheses on the right-hand side (weight and displacement).{p_end}
{phang2}{cmd:. holmsb weight displacement}{p_end}

{pstd}Adjusted p-values are displayed in the output window and are saved in the active e().{p_end}


{title:Example:  linear regression, multiple models with one hypothesis variable on RHS}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse census}{p_end}

{pstd}Conduct linear regressions and store estimates{p_end}
{phang2}{cmd:. eststo r1: regress divorce marriage pop}{p_end}

{phang2}{cmd:. eststo r2: regress divorce marriage pop popurban death}{p_end}

{phang2}{cmd:. eststo r3: regress divorce marriage pop popurban death i.region}{p_end}

{pstd}Compute adjusted p-values for 1 hypothesis variable on the right-hand side (population).{p_end}
{phang2}{cmd:. holmsb marriage, multiple}{p_end}

{pstd}Adjusted p-values are displayed in the output window and are saved in each of the stored estimates.{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:holmsb} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(psidak)}}the set of Sidak-adjusted p-values.{p_end}
{synopt:{cmd:e(pbonf)}}the set of Bonferroni-adjusted p-values.{p_end}
{p2colreset}{...}

{title:Author}

{p 4 4 2}Marc Ragin, University of Georgia{p_end}
{p 4 4 2}mragin@uga.edu{p_end}


{title:References}

{pstd}Holm, S. (1979). A simple sequentially rejective multiple test procedure. {it:Scandinavian Journal of Statistics}, 65-70.{p_end}
