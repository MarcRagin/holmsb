# holmsb
Stata package for multiple hypothesis testing adjustment using Holm step-down algorithm

**holmsb** calculates Sidak-adjusted and Bonferroni-adjusted p-values using the free step-down methodology of Holm (1979). It outputs the adjusted p-values to estimates stored in memory.

Designed to adjust p-values in two situations:
1. A single linear regression with multiple hypothesis variables on the right-hand side.
2. Multiple linear regressions with the same single hypothesis variable on the right-hand side.

Further details are in the help file.

## Installation
To install the **holmsb** package, type or copy-and-paste the following commands into Stata:

net install holmsb, from("https://raw.github.com/MarcRagin/holmsb/master/") replace
