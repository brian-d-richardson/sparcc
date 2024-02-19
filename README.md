
# sparcc: semiparametric censored covariate estimation <img id="sparcc_hex" src="man/figures/sparcc_hex.png" align="right" width="125"/>

Brian D. Richardson

## Installation

Installation of `sparcc` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done with the following code:

``` r
# install the package
devtools::install_github(repo = "brian-d-richardson/sparcc", 
                         ref = "main")
```

``` r
# load the package
library(sparcc)
```

The `sparcc` package contains functions to analyze data with a random
right censored covariate using a semiparametric method. The methods
implemented are introduced in the paper, “Doubly robust estimation under
a random right censored covariate,” which is currently in progress.
