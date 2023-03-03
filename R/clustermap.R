#' Cluster map
#' 
#' Generate an interactive clustered heatmap of summary statistics data.
#' @param dat Summary statistics data.
#' @param cols Numeric columns to plot in the heatmap.
#' @param i Indices of rows to include. Set to \code{NULL} to include all rows, 
#' but be warned that this can become very computationally expensive.
#' @param agg_var Variable to aggregate data by. 
#' Set to \code{NULL} to skip this step.
#' @param agg_fun Function to aggregate \code{cols} with.
#' @param as_cor Show the heatmap as a correlation matrix
#'  instead of a feature x sample matrix.
#' @param annot_vars Variables in \code{dat} to include as row-wise annotations.
#' @param show_plot Print the plot.
#' @inheritParams heatmaply::heatmaply
#' @inheritDotParams heatmaply::heatmaply
#' @returns A named list containing 
#' an interactive heatmaply object and the data used to create it.
#' 
#' @export
#' @import data.table
#' @examples 
#' dat <- ThreeWayTest::data_matrix_final
#' cm <- clustermap(dat = dat)
clustermap <- function(dat,
                       cols = grep("^w",names(dat), value = TRUE),
                       i = seq_len(50),
                       agg_var = NULL,
                       agg_fun = mean,
                       as_cor = FALSE,
                       k_row = 3,
                       annot_vars = c("TYPE","GENE"),
                       show_plot = TRUE,
                       ...){
    
    requireNamespace("heatmaply")  
  
    if(is.null(agg_var)) agg_var <- "SNP"
    dat <- postprocess_data(dat = dat,
                            cols = cols, 
                            agg_var = agg_var,
                            agg_fun = agg_fun)
    #### Subset data ####
    if(!is.null(i)) dat <- dat[i,]
    #### Create matrix ####
    X <- as.matrix(dat[,cols,with=FALSE]) |> `rownames<-`(
        dat[[agg_var]]
    )
    #### Convert to correlation matrix ####
    if(isTRUE(as_cor)) X <- stats::cor(X)
    #### Create row annotations ####
    annot_vars <- annot_vars[annot_vars %in% names(dat) &
                                 (!annot_vars %in% agg_var)]
    row_side_colors <- if(length(annot_vars)>0) {
        dat[,annot_vars,with=FALSE]
    } else {
        NULL
    }
    #### Create heatmap ####
    cm <- heatmaply::heatmaply(X, 
                               row_side_colors = row_side_colors, 
                               k_row = min(nrow(X),k_row),
                               ...)
    #### Show plot ####
    if(isTRUE(show_plot)) methods::show(cm)
    return(list(plot=cm,
                data=list(X=X,
                          row_side_colors=row_side_colors)))
}
