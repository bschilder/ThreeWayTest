#' Get full 1000 Genomes Project genotype data
#' 
#' Download the full version of the \link[ThreeWayTest]{selected_genotype} data, 
#' which contains genotypes for 2,504 individuals from the 
#' \href{https://www.internationalgenome.org/}{1000 Genomes Project}.
#' @source 
#' \code{
#' #### Preprocessing steps ####
#' dat <- get("final_1kg_genotype_correct")
#' dat <- data.table::data.table(dat)
#' tmp <- file.path(tempdir(),"final_1kg_genotype_correct.rds")
#' saveRDS(dat,tmp)
#' piggyback::pb_upload(file = tmp, 
#'                      repo = "bschilder/ThreeWayTest",
#'                      tag = "v0.0.2")
#' }
#' @inheritParams get_data
#' @inheritParams piggyback::pb_download
#' @returns Large \link[data.table]{data.table} of genotype data.
#' 
#' @export
#' @examples 
#' full_genotype <- get_full_genotype()
get_full_genotype <- function(file = "final_1kg_genotype_correct.rds",
                              tag = "latest",
                              overwrite = TRUE,
                              .token = gh::gh_token()){
    
    tmp <- get_data(file = file,
                    overwrite = overwrite,
                    tag = tag,
                    .token = .token)
    #### Load obj ####
    msg <- paste("Loading",file,"into R...")
    message(msg)
    readRDS(tmp)  
}
