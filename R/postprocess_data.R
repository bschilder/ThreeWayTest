#' Post-process results data
#' 
#' Post-process results data.
#' @inheritParams clustermap
#' @returns (Aggregated) \link[data.table]{data.table}
#' 
#' @export
#' @import data.table
#' @examples 
#' dat <- ThreeWayTest::data_matrix_final
#' dat2 <- postprocess_data(dat = dat) 
postprocess_data <- function(dat,
                             cols = grep("^w",names(dat), value = TRUE), 
                             agg_var = NULL,
                             agg_fun = mean){
    
    TYPE <- ANNOT <- GENE <- SNP <- .SD <- n_snps <- NULL;
    
    dat <- data.table::data.table(dat)
    if(!"GENE" %in% names(dat)){
        dat[,GENE:=data.table::tstrsplit(ANNOT,"\\(")[[1]]]    
    }
    if(!"TYPE" %in% names(dat)){
        types <- gsub("\\|","",tolower(
            data.table::fcoalesce(data.table::tstrsplit(dat$ANNOT,"=")[-1])
        )
        )
        dat[,TYPE:=types] 
    }
    if(!is.null(agg_var) &&
       agg_var!='SNP'){
        mean_cols <- c(cols,"MAF")
        dat <- dat[,(mean_cols):=lapply(.SD,agg_fun),
                   .SDcols=mean_cols, 
                   by=agg_var][,n_snps:=(length(unique(SNP))),
                               by=agg_var][,c(agg_var,mean_cols,"n_snps"),
                                           with=FALSE] |>
            unique()
    }
    return(dat)
}