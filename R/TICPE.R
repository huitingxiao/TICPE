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
    C_dat<-matrix(0,Lgene,2)
    N_up<-matrix(0,Lgene,ncol(expr))
    N_down<-matrix(0,Lgene,ncol(expr))
    reverse<-matrix(0,ncol(expr),4)
    for (k in 1:Lgene){
      C_up<-stable_pairs1[which(stable_pairs1[,1]%in%gid[k]),2]
      C_down<-stable_pairs1[which(stable_pairs1[,2]%in%gid[k]),1]
      pairs<-c(C_up,C_down)
      C_dat[k,1]=length(C_up)
      C_dat[k,2]=length(C_down)
      N_tmp<-matrix(rep(Nnorm[k,],length(pairs)),ncol=dim(Nnorm)[2],byrow=T)-expr[match(pairs,rownames(expr)),]##
      N_up[k,]<-colSums(N_tmp>=0)
      N_down[k,]<-colSums(N_tmp<0)
    }
    reverse[,1]<-apply(N_up,2,sum)
    reverse[,2]<-apply(N_down,2,sum)
    reverse[,3]<-rep(sum(C_dat[,1]),ncol(expr))
    reverse[,4]<-rep(sum(C_dat[,2]),ncol(expr))
    up_scores<-apply(reverse,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$estimate)
    estimated_prop[i,]<-((up_scores-min(up_scores))/(as.numeric(parameter[i,2])+alpha))^as.numeric(parameter[i,3])
  }
  return(estimated_prop)
}
