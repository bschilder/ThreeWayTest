#' Network map: make tidygraph
#' 
#' Subfunction of \link[ThreeWayTest]{networkmap}.
#' @inheritParams networkmap
#' @returns \link[tidygraph]{tbl_graph}
#' @keywords internal
networkmap_tidygraph <- function(dat,
                                 cols,
                                 node_vars){ 
    
    requireNamespace("tidygraph") 
    value <- node_type <- name <- NULL;
    
    dat_melt <- data.table::melt.data.table(dat,
                                            measure.vars = cols, 
                                            variable.name = "dataset", 
                                            na.rm = TRUE)  
    edge_meta_vars <- names(dat_melt)
    edges <- lapply(seq_len(length(node_vars)-1), function(i){
        dat_melt[,c(node_vars[i],node_vars[i+1],
                    unique(c("value", edge_meta_vars))),
                 with=FALSE] |>
            `colnames<-`(c("from","to",unique(c("value",edge_meta_vars))))
    }) |> data.table::rbindlist()
    
    nodes <- data.table::melt.data.table(
        dat_melt,
        measure.vars = node_vars, 
        variable.name = 'node_type',
        value.name = "name")[,list(mean_value=mean(value,na.rm=TRUE),
                                   node_type=unique(node_type)),
                             by=c("name")]|> data.table::setkeyv("name") 
    # tg <- tidygraph::tbl_graph(nodes = nodes, edges = edges)  
    g <- tidygraph::as_tbl_graph(x = edges)
    g <- g |> tidygraph::activate("nodes") |>
        tidygraph::mutate(mean_value=nodes[name,]$mean_value,
                          node_type=nodes[name,]$node_type) 
    return(list(graph=g, 
                nodes=nodes,
                edges=edges))
}