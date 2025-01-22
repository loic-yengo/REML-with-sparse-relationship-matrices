This folders contains a script (`sparseREML_v0.6_toshare.R`) ro run linear mixed model analyses with sparse genomic relationship matrices and an example script (`Example.R`) simulating data for siblings pairs and using the sparseREML() function to re-estimate parameters.

REQUIREMENTS

This code runs under any OS that can run the software package R.
The main function depends on the R package `Matrix`, which handles sparse matrices.
The simulation function also requires the R package `MASS` to simulate multivariate normal distribution. 

RUNNING THE EXAMPLE

To run the example you can simply "source" the file as

R_prompt> `source("Example.R")`

The code should run in less than 30 seconds and print the results of the analysis.
The simulation data consist of N=100,000 sibling pairs with their simulated relatedness (object G: `summary(G)`) ranging between 0.34 and 0.68 and their simulated phenotypes.
Phenotypes are simulated to have an expected mean of 0 and an expected variance of 1. The heritability of the trait is `hsq=0.6` and the variance explained by shared environmental effects is `csq=0.1`.

The result of analysis are in the R object `model_ACE`

R_prompt> `print(model_ACE$h2_g)`

The estimated heritability is 0.676 (s.e. 0.064) and the estimated shared environmental effect is 0.064 (0.032).
