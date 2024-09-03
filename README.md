# EcoState

[![Documentation](https://img.shields.io/badge/documentation-EcoState-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/EcoState/)

Package _EcoState_ fits a state-space mass-balance model intended for aquatic ecosystems, using mass-balance equations matching those from Ecopath and dynamical equations matching Ecosim.  Unlike Ecopath with Ecosim (EwE), EcoState fits both biological parameters (e.g., equilibrium biomass and predator-prey vulnerability) and measurement parameters (e.g., catchability coefficients) via fit to time-series data.  EcoState also estimates additional process errors representing nonstationarity in growth efficiency, ecotrophic efficient, migration, or other unmodeled processes.  These process errors allow biomass patterns to closely match available data, so that resulting consumption (and associated productivity and mortality rates) can accurately be conditioned upon any residual patterns.     

## Installation

EcoState can be installed from GitHub using:

``` r
library(remotes)
install_github( "James-Thorson-NOAA/EcoState" )
```

Or to access vignettes from your R session, please instead use:

``` r
remotes::install_github( "James-Thorson-NOAA/EcoState",
                          build_vignettes = TRUE )
browseVignettes("EcoState")
```

# More details 

For more background please read:

Thorson, J.  Kristensen, K., Aydin, K., Gaichas, S., Kimmel, D.G., McHuron, E.A., Nielsen, J.N., Townsend, H., Whitehouse, G.A. EcoState:  Extending Ecopath with Ecosim to estimate biological parameters and process errors using RTMB and time-series data.  Pre-print URL: https://doi.org/10.32942/X2QK81 
