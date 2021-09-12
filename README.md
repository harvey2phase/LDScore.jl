# LDScore

LDSC \[1\] implementation in Julia.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://harvey2phase.github.io/LDScore.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://harvey2phase.github.io/LDScore.jl/dev)
[![Build Status](https://travis-ci.com/harvey2phase/LDScore.jl.svg?branch=master)](https://travis-ci.com/harvey2phase/LDScore.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/harvey2phase/LDScore.jl?svg=true)](https://ci.appveyor.com/project/harvey2phase/LDScore-jl)
[![Build Status](https://api.cirrus-ci.com/github/harvey2phase/LDScore.jl.svg)](https://cirrus-ci.com/github/harvey2phase/LDScore.jl)
[![Coverage](https://codecov.io/gh/harvey2phase/LDScore.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/harvey2phase/LDScore.jl)
[![Coverage](https://coveralls.io/repos/github/harvey2phase/LDScore.jl/badge.svg?branch=master)](https://coveralls.io/github/harvey2phase/LDScore.jl?branch=master)

## To Do
- `ld_score_regression` is not complete
- `estimate_h2`
	- Only has one test case
	- Lots of error checks missing
- Build struct `Gencov` the same way `Hsq` is built to reuse abstract type `LD_Score_Regression` functions

## Estimating SNP-hertiability from summary statistics
To estimate SNP-heritability, we need three sources of data:
1. Summary statistics from a GWAS.
2. List of SNPs.
3. LD scores for each SNP.

LD scores for European and Asian populations can be directly downloaded from the [original LDSC website](https://github.com/bulik/ldsc).

## LD Score Regression
LD Score regression quantifies contribution of polygenic signals and confounding biases in the inflation of the distribution of test statistic in genome-wide association studies (1); it does this by examinig the relationship between test statistics and linkage disequilibrium 

## References
1. Schizophrenia Working Group of the Psychiatric Genomics Consortium, Bulik-Sullivan BK, Loh P-R, Finucane HK, Ripke S, Yang J, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet. 2015 Mar;47(3):291–5. 
