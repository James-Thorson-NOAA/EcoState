# EcoState

[![Documentation](https://img.shields.io/badge/documentation-EcoState-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/EcoState/)

Package _EcoState_ fits a state-space mass-balance model intended for aquatic ecosystems, using mass-balance equations matching those from Ecopath and dynamical equations matching Ecosim.  Unlike Ecopath with Ecosim (EwE), EcoState fits both biological parameters (e.g., equilibrium biomass and predator-prey vulnerability) and measurement parameters (e.g., catchability coefficients) via fit to time-series data.  EcoState also estimates additional process errors representing nonstationarity in growth efficiency, ecotrophic efficient, migration, or other unmodeled processes.  These process errors allow biomass patterns to closely match available data, so that resulting consumption (and associated productivity and mortality rates) can accurately be conditioned upon any residual patterns.     

## Installation

EcoState can be installed from GitHub using:

``` r
library(remotes)
install_github( "James-Thorson-NOAA/EcoState" )
```
