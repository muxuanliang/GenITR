# GenITR

`GenITR` is an R package for learning and evaluating individualized treatment rules using concordance-based methods, kernel-based estimators, and related value inference procedures. The package is designed for settings where treatment decisions depend on patient-specific covariates and outcomes, with support for estimation of treatment rules as well as inference on the value of a learned rule.

The main exported functions currently include:

- `GenITR()`: estimate an individualized treatment rule
- `GenValueInfer()`: perform inference for the value of an individualized treatment rule

## Package structure

This repository contains the following main components:

- `R/`: core source code for the package
- `man/`: automatically generated `.Rd` documentation files
- `simulation/`: scripts for simulation studies and code used to generate simulation figures
- `real_data/`: code for real data analysis

## Simulation studies

The `simulation/` folder contains scripts used to generate simulated datasets, fit the proposed and competing methods, evaluate treatment rule performance, and produce figures and summary results.

These scripts include code for:

- generating synthetic predictors, treatments, and outcomes under multiple scenarios
- fitting treatment rules using `GenITR()`
- evaluating rule value using `GenValueInfer()`
- comparing concordance-based, kernel-based, and alternative methods
- generating output files used for downstream visualization and figures

## Real data analysis

The `real_data/` folder contains code used for the real data application and analysis pipeline.

The raw real dataset is **not included in this repository** due to privacy and data-use restrictions. Users interested in reproducing the real data analysis should supply their own appropriately formatted dataset and modify the file paths as needed.

## Installation

You can install the development version from the local source:

```r
# install.packages("devtools")
devtools::install()
```

If hosted on GitHub, installation can be done with:

```r
# install.packages("devtools")
devtools::install_github("your_username/GenITR")
```

## Example

Below is a simple example illustrating how to generate simulated data, fit an individualized treatment rule, and perform value inference. This example is adapted from the simulation script in the repository.

```r
library(GenITR)
library(MASS)

set.seed(1)

p <- 4
beta <- c(1, 1, -1, 1) * 0.5
gamma1 <- c(1, -1, 1, 1) * 0.5
gamma2 <- c(1, 0, -1, 0) * 0.5
variance <- diag(p)

generateData <- function(sampleSize, case = 1) {
  x <- MASS::mvrnorm(sampleSize, mu = rep(0, p), Sigma = variance)
  prob <- apply(x, 1, function(t) {
    exp(0.25 * (t[1]^2 + t[2]^2 + t[1] * t[2])) /
      (1 + exp(0.25 * (t[1]^2 + t[2]^2 + t[1] * t[2])))
  })
  trt <- rbinom(sampleSize, 1, prob)

  mu <- 1 + x %*% gamma1
  d <- 2 * x %*% beta
  y <- mu + trt * d + rnorm(sampleSize, 0, 0.5)

  list(
    predictor = x,
    treatment = trt,
    outcome = as.numeric(y)
  )
}

data <- generateData(300)

fit <- GenITR(
  data = data,
  dataRef = NULL,
  compareFun = NULL
)

fit
```

You can also perform inference on the value of the estimated rule:

```r
infer <- GenValueInfer(
  data = data,
  dataRef = NULL,
  compareFun = NULL,
  method = "concordance"
)

infer
```

## Notes

- Some internal helper functions are intentionally not exported.
- Documentation files in `man/` are generated automatically from roxygen comments.
- To regenerate documentation, run:

```r
devtools::document()
```

## License

Please add your preferred license information here.
