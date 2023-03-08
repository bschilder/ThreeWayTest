#' data_matrix_final
#'
#' @description
#' Data matrix of Polyunsaturated fatty acids in real data analysis。
#' CHR is the chromesome number.
#' SNP is the number of SNP.
#' BP is the poisition of the SNP.
#' ANNOT is gene that this SNP is annoted to.
#' Other Allele is the  allele that counts as 0.
#' Effect Allele is the allele that counts as 1.
#' MAF is the allele frequency of other allele.
#' The last 6 columns of wAA, wALA, wDHA, wEDA, wEPA and wLA is the 
#' wald test statistics of 6 phenotypes: arachidonic acid, 
#' alpha-linolenic acid, docosahexanoic acid, eicosadienoic acid
#' and eicosapentanoic acid.
#' For original data and more detailed data description, please visit:
#' https://grasp.nhlbi.nih.gov/FullResults.aspx. The data is in Year 2009 (2).
#' \code{
#' load(file.path("data","data_matrix_final.Rdata"))
#' usethis::use_data(data_matrix_final, overwrite = TRUE)
#' }
#' @usage data("data_matrix_final")
"data_matrix_final"

#' covariance_matrix_data
#'
#' @description 
#' Part of Polyunsaturated fatty acids data after LD pruning to calculate 
#' phenotype covariance matrix. Detailed description is the same as 
#' data_matrix_final. 
#' \code{
#' load(file.path("data","covariance_matrix_estimated.Rdata"))
#' usethis::use_data(covariance_matrix_data, overwrite = TRUE)
#' }
#' @usage data("covariance_matrix_data")
"covariance_matrix_data"

#' selected_genotype
#'
#' @description 
#' 2504 samples in 1000 genomes project for 
#' calculating covariance matrix of genotypes.
#' CHR is the chromosome number.
#' SNP is the number of SNP.
#' BP is the position of the SNP.
#' Counted is the  allele that counts as 1.
#' ALT is the allele that counts as 0.
#' Due to R package limit of data, this is a part-version with 9 SNPs of gene
#' FADS2 to run getting_start.R in the vignettes
#' For full Data, Please go https://github.com/bschilder/ThreeWayTest.
#' \code{
#' load(file.path("data","selected_genotype.Rdata"))
#' usethis::use_data(selected_genotype, overwrite = TRUE)
#' }
#' @usage data("selected_genotype")
"selected_genotype"

#' gene_list
#'
#' @description 
#' List of gene names for performing real data analysis.
#' \code{
#' load(file.path("data","gene_list.Rdata"))
#' usethis::use_data(gene_list, overwrite = TRUE)
#' }
#' @usage data("gene_list")
"gene_list"

#' gene_length_list
#'
#' @description 
#' Number of SNPs contains in each gene.
#' \code{
#' load(file.path("data","gene_length_list.Rdata"))
#' usethis::use_data(gene_length_list, overwrite = TRUE)
#' }
#' @usage data("gene_length_list")
"gene_length_list"

#' PCA_result
#'
#' @description 
#' Principle components for population stratification used as 
#' covariates in regression model. Here we use PCA_result$principal.coordinates
#' as the covariates in our regression model.
#' \code{
#' load(file.path("data","PCA_result.Rdata"))
#' usethis::use_data(PCA_result, overwrite = TRUE)
#' }
#' @usage data("PCA_result")
"PCA_result"
 