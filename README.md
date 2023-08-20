
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stoppingrule

<!-- badges: start -->
<!-- badges: end -->

## R package “stoppingrule”: Create and describe early stopping rules for clinical studies

stoppingrule is an R package that provides functionality for
constructing, describing, and evaluating stopping rules. Clinical trials
often include these rules in order to ensure the safety of study
treatments and feasibility of the study.

## Installation

The current version of stoppingrule is available in this repository. If
you have the R package “devtools” installed, stoppingrule can be
installed directly from GitHub with:

``` r
require(devtools)
install_github("mjmartens/stoppingrule")
```

## Example

Suppose a clinical trial of 30 patients is being developed. Study
investigators wish to monitor the occurrence of a specific toxicity
type, which occurs in 20% of patients in the target population based on
historical data. We want to create a stopping rule to test whether the
rate of a given toxicity exceeds 20% using a type I error rate of 10%.
Moreover, we wish to check the stopping rule continuously, after each
patient completes follow-up for the toxicity endpoint.

The function `calc.rule.bin` creates a ‘rule.bin’ object which contains
a matrix with the numbers of evaluable patients and their corresponding
stopping boundary values as well as the rule’s design parameters:

``` r
require(stoppingrule)
#> Loading required package: stoppingrule
poc_rule = calc.rule.bin(ns=1:30,p0=0.20,alpha=0.10,type="Pocock")
print(poc_rule)
#> $Rule
#>       N evaluable Reject bdry
#>  [1,]           1           2
#>  [2,]           2           3
#>  [3,]           3           3
#>  [4,]           4           3
#>  [5,]           5           4
#>  [6,]           6           4
#>  [7,]           7           4
#>  [8,]           8           5
#>  [9,]           9           5
#> [10,]          10           5
#> [11,]          11           6
#> [12,]          12           6
#> [13,]          13           6
#> [14,]          14           7
#> [15,]          15           7
#> [16,]          16           7
#> [17,]          17           8
#> [18,]          18           8
#> [19,]          19           8
#> [20,]          20           8
#> [21,]          21           9
#> [22,]          22           9
#> [23,]          23           9
#> [24,]          24           9
#> [25,]          25          10
#> [26,]          26          10
#> [27,]          27          10
#> [28,]          28          11
#> [29,]          29          11
#> [30,]          30          11
#> 
#> $ns
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> [26] 26 27 28 29 30
#> 
#> $p0
#> [1] 0.2
#> 
#> $type
#> [1] "Pocock"
#> 
#> $alpha
#> [1] 0.1
#> 
#> $param
#> NULL
#> 
#> $cval
#> [1] 0.03766345
#> 
#> attr(,"class")
#> [1] "rule.bin"
```

The function call uses the Pocock-type test proposed by Ivanova et
al. 2005 to construct the stopping boundaries and displays the
boundaries for the first few patients. Note that rejection is impossible
with 1 or 2 evaluable patients because the corresponding boundary values
exceed these numbers.

The function `table.rule.bin` can produce a succinct summary of the
stopping rule for the entire cohort from the rule calculated above:

``` r
table.rule.bin(poc_rule)
#>       N evaluable Reject If N >=
#>  [1,] "3 - 4"     "3"           
#>  [2,] "5 - 7"     "4"           
#>  [3,] "8 - 10"    "5"           
#>  [4,] "11 - 13"   "6"           
#>  [5,] "14 - 16"   "7"           
#>  [6,] "17 - 20"   "8"           
#>  [7,] "21 - 24"   "9"           
#>  [8,] "25 - 27"   "10"          
#>  [9,] "28 - 30"   "11"
```

We can also obtain a graphical summary of the stopping rule using the
`plot` function:

<img src="man/figures/README-plot-1.png" width="100%" />

Lastly, the `OC.rule.bin` function can assess the operating
characteristics of this stopping rule. The rejection probability and
expected numbers of evaluated patients and events are computed at true
toxicity rates of p = 20%, 25%, 30%, 35%, and 40% as follows:

``` r
OC.rule.bin(rule=poc_rule,ps=seq(0.2,0.4,0.05))
#>         p Reject Prob E(evalauted) E(events)
#> [1,] 0.20  0.09915263     28.16703  5.713406
#> [2,] 0.25  0.23312769     26.04001  6.635003
#> [3,] 0.30  0.42603085     22.98251  7.074752
#> [4,] 0.35  0.63398556     19.37012  7.024542
#> [5,] 0.40  0.80576733     15.75017  6.620068
```
