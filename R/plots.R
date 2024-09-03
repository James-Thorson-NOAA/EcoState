

#' @title Plot foodweb
#'
#' @description 
#' Plot consumption as a directed graph including all taxa (vertices) 
#' and biomass consumed (arrows).  Taxa are located using tracers,
#' where by default the y-axis is trophic level.  #' 
#'
#' @inheritParams compute_tracer
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
#' @param inverse_method whether to use pseudoinverse or standard inverse
#'
#' @details
#' Trophic level \eqn{l_i} for each predator \eqn{i} is defined as:
#'
#' \deqn{ \mathbf{l - 1 = l Q^*} }
#'
#' where \eqn{\mathbf{Q*}} is the proportion consumption for each predator (column)
#' of different prey (rows).  We identify primary producers as any taxa with no
#' consumption (a column of 0s), and assign them as the first trophic level.
#'
#' @export
plot_foodweb <-
function( Q_ij,
          type_i,
          xtracer_i,
          ytracer_i = rep(1,nrow(Q_ij)),
          B_i = rep( 1, nrow(Q_ij)), 
          taxa_labels = letters[1:nrow(Q_ij)], 
          xloc,
          yloc ){

  #
  if(missing(yloc)){
    yloc = compute_tracer( Q_ij, 
                           inverse_method = "Standard", 
                           tracer_i = ytracer_i,
                           type_i = type_i )[1,]
  }
  if(missing(xloc)){
    xloc = compute_tracer( Q_ij, 
                           inverse_method = "Standard", 
                           tracer_i = xtracer_i,
                           type_i = type_i )[1,]
  }

                
  #
  layout = cbind( x = xloc, 
                  y = yloc )
  dimnames(Q_ij) = list( taxa_labels, taxa_labels )
  #rownames(layout) = taxa_labels                       
  graph = igraph::graph_from_adjacency_matrix( Q_ij, 
                                       weighted = TRUE )
  #plot( graph,
  #           edge.width = log(E(graph)$weight),
  #           vertex.size = log(B_i/min(B_i)) + 10,
  #           type = "width", 
  #           text_size = 4,
  #           #arrow = grid::arrow(type='closed', 18, grid::unit(10,'points')),
  #           layout = layout,
  #           show.legend = FALSE, 
  #           curvature = 0 )
  #axis(1)
  #axis(2)
  #return( list('graph' = graph,
  #             'layout' = layout))

  g = ggnetwork::ggnetwork( x = graph, 
                 layout = layout,
                 scale = FALSE )
  #g = ggnetwork::ggnetwork( Q_ij, 
  #                          layout = layout,
  #                          weighted = TRUE )
  g$log_flow = log(g$weight / min(g$weight,na.rm=TRUE))
  g$log_mass = rep(NA, nrow(g))
  g$log_mass[which(is.na(g$log_flow))] = log(B_i / min(B_i,na.rm=TRUE))
  
  p = ggplot2::ggplot(g) + 
    ggnetwork::geom_edges( ggplot2::aes(x=x, y=y, xend=xend, yend=yend, colour=log_flow) ) +
    ggnetwork::geom_nodes( ggplot2::aes(x=x, y=y, size=log_mass ) ) + 
    ggnetwork::geom_nodetext( ggplot2::aes(x=x, y=y, label=name), fontface="bold", col="red")
  print(p)
  return(invisible(p))
}

