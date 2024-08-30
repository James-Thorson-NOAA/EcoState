
#' @title 
#' eastern Bering Sea ecosystem data
#'
#' @description 
#' Data used to demonstrate a Model of Intermediate Complexity (MICE)
#' for the eastern Bering Sea.  
#' `data(eastern_bering_sea)` loads a list that includes four components:
#' * `Survey` is a long-form data-frame with three columns, providing the Year,
#'    Mass (in relative units for most taxa, and million metric tons for Pollock,
#'    Cod, Arrowtooth, and NFS), and Taxon for each year with available data
#' * `Catch` is a long-form data-frame with three columns, providing the Year,
#'    Mass (in million metric tons), and Taxon for each year with available data
#' * `P_over_B` is a numeric vector with the unitless ratio of biomass production to 
#'   population biomass for each taxon
#' * `Q_over_B` is a numeric vector with the unitless ratio of biomass consumption to 
#'   population biomass for each taxon
#' * `Diet_proportions` is a numeric matrix where each column lists the 
#'   proportion of biomass consumed that is provided by each prey (row)
#'
#' @details
#' The data compiled come from a variety of sources:
#' * Northern fur seal (NFS) survey is an absolute index, corrected for proportion
#'   of time spent in the eastern Bering Sea.  NFS QB is developed from a bioenergetic 
#'   model and also corrected for seasonal residency. Both are provided by 
#'   Elizabeth McHuron. It is post-processed in a variety of ways, and not
#'   to be treated as an index of abundance for NFS for other uses.
#' * Pollock, cod, and arrowtooth surveys are from a bottom trawl survey, and 
#'   cod and arrowtooth are treated as an absolute index.
#' * Copepod and other zooplankton are from an oblique tow bongo net survey, 
#'   with data provided by Dave Kimmel.  It is then post-processed to account
#'   for spatially and seaonally imbalanced data.
#' * Other P_over_B, Q_over_B and Diet_proportions values 
#'   are derived from Rpath models, provided by Andy Whitehouse.
#' * Primary producers is an annual index of relative biomass, developed from monthly
#'   satellite measurements and provided by Jens Nielsen.
#' See Thorson et al. (In review) for more details regarding data standardization 
#'   and sources
#'
#' @name eastern_bering_sea
#' @docType data
#' @usage data(eastern_bering_sea)
#' @keywords data
NULL

#' @title
#' Full rpath inputs for eastern Bering Sea
#'
#' @description
#' All Rpath inputs from Whitehouse et al. 2021
#'
#'
#' @name whitehouse_2021
#' @docType data
#' @usage data(whitehouse_2021)
#' @keywords data
NULL

