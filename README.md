# Linear Lasso
This repository contains code for the Linear Lasso. The Linear Lasso is a method used in linear regression that finds the important predictors and calculates it coefficients based on least squares estimation. 

The response vector is then seen as the focal point of the space and all other explanatory variables vectors orbit around the response vector. The angles formed between the response vector and the explanatory variables are assumed to be fixed, and will be used as a basis for constructing the method. The information contained in the explanatory variables is projected onto the response vector. The theory of normal linear models allows us to use ordinary least squares (OLS) for the coefficients of the Linear Lasso.

The Linear Lasso (LL) is performed in two steps. First, variables are dropped from the model based on their correlation with the response variable; the number of variables dropped (or ordered) in this step depends on a tuning parameter $\gamma$. Then, an exclusion criterion based on the variance of the distribution of the response variable is introduced to remove (or order) the remaining variables. A repeated cross-validation guides us in the choice of the final model


## Description of the package

This **R** package implements the Linear Lasso which was first used in 
> Fraser, D. A., & Bédard, M. (2022). The linear Lasso: A location model approach. The Canadian Journal of Statistics, 50(2), 437-453.

Some changes have been made to optimize the algorithm.

## Example of Utilization

### Installation

The Linear Lasso algorithm can be run directly from GitHub using the devtools package function ```install_github``` using the following **R** command :

```R
devtools::install_github("yanwatts/linlasso")
```

## Setup 

To use the Lineal Lasso in linear regression, the LL function is available by specifying certain parameters
```R
LL(y, x, gamma = 0.1, K = 10, L = 5, plot = F, french.plot = F)
```

## Input 

The LL function requires minimally the following input :

* ```y``` : the vector of responses;
* ```x``` : the design matrix;
* ```gamma``` : the cutoff parameter (default = 0.1);
* ```K``` : the number of folds for the repeated cross validation (default = 10);
* ```L``` : the number of repetitions for the repeated cross validation (default = 5);
* ```plot``` : plot giving the path of the final model after the one-by-one procedure.
* ```french.plot``` : plot in french.

## Output 

The following is an example of the LL function used for the diabetes dataset introduced by Efron et al. (2004). The output will help us understand how to interpret the Linear Lasso algorithm. We use gamma = 0.2 for illustration purposes. The SEX variable is recoded to 1 and 0.

**First step** : LL eliminates predictors based on their correlation with the response variable

**Second step (One-by-one procedure)** : LL eliminates variables based on the one-by-one procedure

```R
> model.LL = LL(y = diabetes[,11], x = diabetes[,-11], gamma = 0.2, K = 13, L = 50, plot = T)
[1] "Variables left after cutoff = 0.2 : BMI + BP + S1 + S3 + S4 + S5 + S6"
```

The variables left after the first step based on the default cutoff ```gamma = 0.2``` are ```BMI```, ```BP```, ```S1```, ```S3```, ```S4```, ```S5``` and ```S6```. The model gives us the positive correlations between the response variable and the predictors,

```R
> model.LL$c.pos
      BMI        S5        BP        S4        S3        S6        S1       AGE        S2       SEX 
0.5864501 0.5658826 0.4414818 0.4304529 0.3947893 0.3824835 0.2120225 0.1878888 0.1740536 0.0430620 
```

or the table with the cross validated : mean squared errors (MSE.CV), standard deviations (SD.CV) and standard errors (SE.CV) for each rejected variable,

```R
> model.LL$table.MSE
                           
Rejected variables in order Length MSE.CV  SD.CV  SE.CV
                        SEX     10 0.5042 0.1122 0.0311
                        S2       9 0.5199 0.1095 0.0304
                        AGE      8 0.5189 0.1089 0.0302
                        S3       7 0.5176 0.1092 0.0303
                        S6       6 0.5189 0.1091 0.0302
                        S4       5 0.5169 0.1092 0.0303
                        S1       4 0.5199 0.1099 0.0305
                        BP       3 0.5283 0.1088 0.0302
                        S5       2 0.5455 0.1136 0.0315
                        BMI      1 0.6768 0.1446 0.0401
```

or the final model with the smallest cross validated error,

```R
> model.LL$`Variables with minimum MSE`
 [1] "AGE" "SEX" "BMI" "BP"  "S1"  "S2"  "S3"  "S4"  "S5"  "S6" 
```

or the least squares coefficients for the minimum MSE model,

```R
> model.LL$beta.min
            [,1]
AGE  0.009746087
SEX 20.801686649
BMI  5.417986348
BP   0.961537161
S1   1.635507013
S2  -1.604358783
S3  -3.536258970
S4  -7.741258450
S5  -3.598309289
S6   0.080981765
```

or the final model with the 1se standard rule,

```R
> model.LL$`Variables with 1se of minimum MSE`
[1] "BMI" "BP"  "S5" 
```

or the least squares coefficients for the 1se standard rule model,

```R
> model.LL$beta.1se
            [,1]
AGE  0.000000000
SEX  0.000000000
BMI  5.281323909
BP  -0.005810611
S1   0.000000000
S2   0.000000000
S3   0.000000000
S4   0.000000000
S5   3.642028932
S6   0.000000000
```

### Visualization

The following graph plots the One-by-one procedure. The final model is the one with the red line. We can see that basically all variables are kept in the model. A more parcimonious model can be choosen by going at the green point (1se standard rule) and to follow the red line to the right.  


![alt text](diabetes_plot.png)


## Author

[Yan Watts](mailto:yanwatts@hotmail.com?subject=[GitHub]%20Source%20Han%20Sans)

## References
Fraser, D. A., & Bédard, M. (2022). The linear Lasso: A location model approach. The Canadian Journal of Statistics, 50(2), 437-453.

Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). Least Angle Regression. The Annals of Statistics, 32(2), 407–451.
