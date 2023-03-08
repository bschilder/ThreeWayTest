#' Network map
#' 
#' Generate an interactive network map of the summary statistics data.
#' @source
#'  \href{https://ggraph.data-imaginist.com/articles/Layouts.html}{
#'  ggraph layouts}
#' @source \href{https://github.com/thomasp85/ggraph/issues/75}{
#' ggraph issue with scale functions}
#' @param node_vars Columns within \code{dat} to set as nodes. 
#' Each node variable will link to the next node in the character vector.
#' You can repeat column names to create more connections between nodes.
#' @param snps A character vector of SNP RSIDs to subset \code{dat} by.
#' @param node_size_range The minimum / maximum size of each node.
#' @inheritParams clustermap
#' @inheritParams ggraph::ggraph
#' @returns A named list containing 
#' a network plot and the data used to create it.
#' 
#' @export
#' @import data.table
#' @examples 
#' dat <- ThreeWayTest::data_matrix_final
#' nm <- networkmap(dat = dat)
networkmap <- function(dat,
                       cols = grep("^w",names(dat), value = TRUE),
                       node_vars = c("dataset","GENE","SNP"),
                       i = seq_len(50),
                       snps = unique(dat$SNP)[i],
                       agg_var = NULL,
                       agg_fun = mean,
                       as_cor = FALSE,
                       k_row = 3,
                       annot_vars = c("TYPE","GENE"),
                       show_plot = TRUE,
                       layout = "nicely",
                       node_size_range = NULL){
    # devoptera::args2vars(networkmap) 
    SNP <- NULL;
    
    if(is.null(agg_var)) agg_var <- "SNP"
    dat <- postprocess_data(dat = dat,
                            cols = cols, 
                            agg_var = agg_var,
                            agg_fun = agg_fun)
    #### Subset data ####
    if(length(snps)>0) dat <- dat[SNP %in% snps,]
    if(length(i)>0) dat <- dat[i,]
    #### Make graph data ####
    tg <- networkmap_tidygraph(dat = dat,
                               cols = cols, 
                               node_vars = node_vars) 
    #### Estimate a reasonable point size range ####
    if(is.null(node_size_range)){
        node_size_range <- c(3,20)*(69/length(tg$graph))
    }
    #### Make plot ####
    # gg <- networkmap_ggraph(tg = tg,
    #                         layout = layout,
    #                         node_size_range = node_size_range)
    gg <- networkmap_ggnetwork(tg = tg,
                               layout = layout,
                               node_size_range = node_size_range) 
    #### Show ####
    if(isTRUE(show_plot)) methods::show(gg)
    #### Return ####
    return(list(plot=gg,
                data=tg)) 
}