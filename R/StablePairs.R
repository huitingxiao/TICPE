#' @title  get stable pairs from cancer cells
#' @description  the function is used to get highly stable gene pairs from a series of cancer cells.

#' @param cancer_exp a dataset which is used to screen the stable pairs.A matrix with row names as Entrez gene ID and columns as cancer cell samples.
#' @param geneid the rownames of cancer_exp.
#' @param freq a threshold value is used to stipulate that a gene pair maintain the circumstance of > the threshold so that the gene pair should be retained.The threshold should between 0-1, and the default is 0.99.
#'
#' @return a matrix of stable pairs from cancer cells data.
#' @export
#'
#' @examples
#'

StablePairs=function(cancer_exp,geneid,freq){
  cancer_exp=as.matrix(cancer_exp)
  freqs=list()
  for(i in 1:(length(geneid)-1)){
    geneid1=geneid[-c(1:i)]
    pair1=cbind(geneid[i],geneid1)
    coms=cancer_exp[match(pair1[,1],geneid),,drop=F]-cancer_exp[match(pair1[,2],geneid),,drop=F]
    freq1=rowMeans(coms>0)
    freqs[[i]]=freq1
  }
  pairs=t(combn(geneid,2))
  freqs=unlist(freqs)
  stablepairs=rbind(pairs[freqs>freq,],pairs[freqs<(1-freq),c(2,1)])
  return(stablepairs)
}
