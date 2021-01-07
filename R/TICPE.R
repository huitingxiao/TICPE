#' @title  Estimate the proportion of tumor-infiltrating immune cells
#' @description  This package is used to estimate the proportion of tumor-infiltrating immune cells in a specific types of cancer.
#' @param expr the gene expression data set. A matrix with row names as Entrez gene ID and columns as samples.
#' @param select_siggene list of signature genes of each immune cell type.
#' @param stablepairs a matrix of stable pairs from cancer cells with two columns,the first column is the higher expression gene,the second column is the lower expression gene.
#' @param parameter a matrix of transformed parameters for each immune cell type.
#' @param alpha alpha a value to avoid division by 0. Default = 0.5.
#'
#' @return a matrix of the estimated proportion of tumor-infiltrating immune cells.
#' @export
#'
#' @examples
#'
TICPEScores<-function(expr,select_siggene,stablepairs,parameter,alpha = 0.5){
  pair_gid<-unique(c(stable_pairs[,1],stable_pairs[,2]))
  congid<-intersect(rownames(expr),pair_gid)
  loc<-union(which(stable_pairs[,1]%in%setdiff(pair_gid,congid)),which(stable_pairs[,2]%in%setdiff(pair_gid,congid)))
  if (length(loc)==0) {
    stable_pairs1<-stable_pairs
  } else {
    stable_pairs1<-stable_pairs[-loc,]
  }
  
  estimated_prop<-matrix(,length(select_siggene),ncol(expr))
  rownames(estimated_prop)<-names(select_siggene)
  for(i in 1:length(select_siggene)){
    gid<-as.numeric(select_siggene[[i]])
    Nnorm<-expr[which(rownames(expr)%in%gid),]
    Lgene<-length(gid)
    reverse<-matrix(0,ncol(expr),4)
    up_pairs<-stable_pairs[which(stable_pairs[,1]%in%gid),]
    up_exp<-expr[match(up_pairs[,1],rownames(expr)),,drop=F]-expr[match(up_pairs[,2],rownames(expr)),,drop=F]
    up1<-colSums(up_exp>0,na.rm=TRUE)
    down1<-colSums(up_exp<0,na.rm=TRUE)
    down_pairs<-stable_pairs[which(stable_pairs[,2]%in%gid),]
    down_exp<-expr[match(down_pairs[,2],rownames(expr)),,drop=F]-expr[match(down_pairs[,1],rownames(expr)),,drop=F]
    up2<-colSums(down_exp>0,na.rm=TRUE)
    down2<-colSums(down_exp<0,na.rm=TRUE)
    reverse[,1]<-up1+up2
    reverse[,2]<-down1+down2
    reverse[,3]<-nrow(up_pairs)
    reverse[,4]<-nrow(down_pairs)
    up_scores<-apply(reverse,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$estimate)
    estimated_prop[i,]<-((up_scores-min(up_scores))/(as.numeric(parameter[i,2])+alpha))^as.numeric(parameter[i,3])
  }
  return(estimated_prop)
}
