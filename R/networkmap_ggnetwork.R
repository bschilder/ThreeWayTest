#' Network map: make ggraph
#' 
#' Subfunction of \link[ThreeWayTest]{networkmap}.
#' @inheritParams networkmap
#' @returns \link[ggnetwork]{ggnetwork}
#' @keywords internal
networkmap_ggnetwork <- function(tg,
                              layout = "nicely",
                              node_size_range = NULL,
                              edge_alpha=.25, 
                              seed=2023){
    requireNamespace("ggnetwork")
    requireNamespace("ggplot2")
    requireNamespace("pals")
    mean_value <- xend <- yend <- name <-x <- y<-  node_type <- NULL;
    
    set.seed(seed)
    ig <- tidygraph::as.igraph(tg$graph) #|> igraph::simplify()
    # ggnet <- tg$graph |> tidygraph::activate("edges") |>tidygraph::as_tibble()
    ggnet <- suppressWarnings(
        ggnetwork::fortify(ig, list(layout="graphopt"))
    )
    
   ( ggplot2::ggplot(ggnet, 
                    ggplot2::aes(x=x,y=y,
                                        xend=xend,
                                        yend=yend,
                                        label=name)) +
        ggplot2::geom_density2d_filled(adjust=.75) + 
        ggnetwork::geom_edges(
            ggplot2::aes(color=mean_value),
            alpha=.5,
            arrow = ggplot2::arrow(length = ggplot2::unit(5, "pt"),
                                   type = "closed"),
            curvature = 0.2) +
        ggplot2::scale_color_gradientn(colors = pals::plasma(20)) +
        ggnetwork::geom_nodes(ggplot2::aes(color=mean_value,
                                           shape=node_type,
                                           size=as.numeric(node_type))) +
        ggnetwork::geom_nodelabel(ggplot2::aes(color=as.numeric(node_type)), 
                                  fill=ggplot2::alpha("white",.7)) +
        ggplot2::scale_size_continuous(range = node_size_range) + 
        ggnetwork::theme_blank()  )#|>
    
    # plotly::ggplotly()
}



#### GGRAPH APPROACH ####
## IMmplemented this originally but then abandoned when I realized it can't be
## used without importing it.
# networkmap_ggraph <- function(tg,
#                               layout = "nicely",
#                               node_size_range = NULL,
#                               edge_alpha=.25){
#     requireNamespace("ggraph")
#     requireNamespace("ggplot2")
#     requireNamespace("pals")
#     mean_value <- node_type <- name <- value <- NULL;  
#     #### Get to/from for edge bundling #### 
#     con_data <- ggraph::get_con(match(tg$edges$from, tg$nodes$name), 
#                                 match(tg$edges$to, tg$nodes$name),
#                                 value=tg$edges$value)
#     # root <- which((tidygraph::activate(tg$graph,"nodes")|>
#     #                    data.frame())$node_type=="dataset")
#     ggraph::ggraph(graph = tg$graph,
#                    layout = layout
#                    # charge = 0.1,
#                    # spring.length = 1,
#     ) + # "tree" "stress"
#         ggraph::geom_node_voronoi(ggplot2::aes(color=mean_value,
#                                                fill=mean_value),
#                                   show.legend = FALSE) +
#         ggplot2::scale_fill_gradientn(colors=pals::ocean.thermal(1000)[1:200]) + 
#         # ggraph::geom_node_tile(ggplot2::aes(width=as.numeric(node_type)/10, 
#         #                                     height=as.numeric(node_type)/10,
#         #                                     fill=mean_value)) + 
#         # ggraph::geom_edge_density(ggplot2::aes(fill=value)) +
#         # ggraph::geom_edge_diagonal(alpha=.25, 
#         #                            strength = .5,
#         #                            ggplot2::aes(color=value)) +  
#         ggraph::geom_conn_bundle(
#             ggplot2::aes(colour = value),#ggplot2::after_stat(index)),
#             data = con_data,
#             tension = 1,
#             edge_alpha = edge_alpha
#         ) +
#         # ggraph::scale_edge_size(range = c(3,10)) + 
#         # ggraph::scale_edge_colour_distiller(palette = "RdPu") +
#         ggraph::scale_edge_colour_gradientn(
#             colours = ggplot2::alpha(pals::parula(100),alpha = edge_alpha)) +
#         ggraph::geom_node_point(ggplot2::aes(fill=mean_value,
#                                              color=as.numeric(node_type),
#                                              size=as.numeric(node_type)*mean_value,
#                                              # size=mean_value,
#                                              shape=node_type),
#                                 # fill="white",
#                                 stroke=5,
#                                 alpha=0.8) +  
#         ggplot2::scale_size_continuous(range = node_size_range) +
#         ggplot2::scale_shape_manual(
#             values = seq(21,21+length(unique(tg$nodes$node_type)))
#         ) + 
#         # ggraph::geom_node_label(ggplot2::aes(label = name), 
#         #                         fill=ggplot2::alpha("white",alpha = .5)
#         #                         ) + 
#         ggraph::geom_node_text(ggplot2::aes(label = name,
#                                             color=as.numeric(node_type)),
#                                alpha=1, 
#                                family = "mono", fontface = "bold",
#                                # color=ggplot2::alpha("white",alpha = .9)
#         ) +
#         ggraph::geom_node_text(ggplot2::aes(label = name), 
#                                family = "mono",
#                                color="white", 
#                                alpha=.85
#         ) + 
#         ggraph::scale_color_viridis() + 
#         ggplot2::scale_color_gradientn(colors=pals::gnuplot(100)[1:80]) + 
#         # ggraph::scale_fill_viridis() + 
#         ggraph::theme_graph() 
#     # plotly::ggplotly(gg) 
# }

