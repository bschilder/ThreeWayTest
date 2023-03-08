#' MGAS
#'
#' This is the method proposed by Van_MGAS
#' @param z_vector Column vectorized data matrix with rows represent
#'  phenotype and columns represent genotype.
#' @param est_genetic_cor Estimated phenotype correlation matrix.
#' @param est_pheno_cor Estimated genotype correlation matrix.
#' @returns A numeric value represents the p value of MGAS.
#' 
#' @export
#' @importFrom stats pchisq
#' @examples
#' z_vector<-MASS::mvrnorm(1,mu=rep(0,9),Sigma = diag(nrow = 9, ncol = 9))
#' genotype_covariance<-diag(nrow = 3,ncol = 3)
#' phenotype_covariance<-diag(nrow = 3,ncol = 3)
#' ThreeWayTest::MGAS(z_vector=z_vector,
#'                    est_genetic_cor = genotype_covariance, 
#'                    est_pheno_cor =phenotype_covariance)
MGAS<-function(z_vector,est_genetic_cor,est_pheno_cor){
    X <- kronecker(est_genetic_cor,est_pheno_cor,FUN = "*")
    beta = c(0.3867,0.0021,-0.1347,-0.0104,0.7276,0.0068)
    p_value<-stats::pchisq(z_vector^2,1,lower.tail = FALSE) 
    tmp=sort(p_value,index.return=TRUE) 
    pj=tmp$x 
    iorder=tmp$ix
    length = ncol(est_genetic_cor)*ncol(est_pheno_cor)
    r2=matrix(0,length,length) 
    r2=X[iorder,iorder]
    ro=diag(length)
    for (i1 in 1:length)  {  
        for (i2 in 1:i1) {
            if (i1>i2) {
                er=r2[i1,i2]
                ro[i1,i2]=ro[i2,i1]= beta[1]*er^6+beta[2]*er^5+beta[3]*er^4+beta[4]*
                    er^3+beta[5]*er^2+beta[6]*er
            }}}
    alllam=eigen(ro[1:length,1:length])$values
    qepj=length 
    for (i1 in 1:length) {
        qepj=qepj-(alllam[i1]>1)*(alllam[i1]-1) }
    
    qej=matrix(c(seq(1,length,1)),length,1,byrow=T)
    for (j in 1:length) { 
        sellam=eigen(ro[1:j,1:j])$values
        id=j
        for (i1 in 1:id) {
            qej[j,1]=qej[j,1]-(sellam[i1]>1)*(sellam[i1]-1)
        }
    }
    pg=matrix(0,length,1)
    for (i in 1:length) {
        pg[i,1]=(qepj/qej[i,1])*pj[i]
    }
    pg=pg[iorder]
    p_mgas<-min(pg)
    return(p_mgas)
}