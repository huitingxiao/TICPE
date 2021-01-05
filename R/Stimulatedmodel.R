#' @title  get parameters for upregulated scores transformation
#' @description  the function is used to design a transformation pipeline to transform the cell infiltration scores to the estimating cell proportions.

#' @param combat_edata gene expression data including immune cells,cancer cells,and normal tissues using ComBat to remove batch effects. A matrix with row names as Entrez gene ID and columns as samples.
#' @param select_siggene list of signature genes for each immune cell type.
#' @param stablepairs a matrix of stable pairs from cancer cells with two columns,the first column is the higher expression gene,the second column is the lower expression gene.
#'
#' @return a matrix of transformed parameters for each immune cell type.
#' @import stats

#' @export
#'
#' @examples
#'
Stimulatedmodel<-function(combat_edata,select_siggene,stablepairs){
  prop<-seq(0.008,0.256, by = 0.008)
  parameter<-matrix(,length(select_siggene),3)
  name<-factor(colnames(combat_edata))
  combat_exp<-NULL
  for(i in 1:nrow(combat_edata)){
    b<-t(as.matrix(tapply(as.matrix(as.numeric(combat_edata[i,])),name,simplify = TRUE,median)))
    combat_exp<-rbind(combat_exp,b)
  }
  norm_control<-combat_exp[,which(colnames(combat_exp)%in%"normal")]
  control<-combat_exp[,which(colnames(combat_exp)%in%"cancer cell")]
  mix<-function(x,y,z){
    dat<-NULL
    for(i in 1:32){
      tmp<-x*0.008*i+y*(0.4-0.008*i)+z*0.6
      dat<-cbind(dat,tmp)}
    return(dat)
  }
  for(j in 1:length(select_siggene)){
    mixture<-mix(combat_exp[,which(colnames(combat_exp)%in%names(select_siggene)[j])],norm_control,control)
    rownames(mixture)=rownames(combat_edata)
    gid<-as.numeric(select_siggene[[j]])
    Lgene<-length(which(!gid==''))
    C_dat<-matrix(0,Lgene,2)
    N_up<-matrix(0,Lgene,32)
    N_down<-matrix(0,Lgene,32)
    sig_exp<-mixture[match(gid[which(!gid=='')],rownames(mixture)),]
    reverse<-matrix(0,32,4)
    up_scores<-NULL
    for (k in 1:Lgene){
      C_up<-stablepairs[which(stablepairs[,1]%in%gid[k]),2]
      C_down<-stablepairs[which(stablepairs[,2]%in%gid[k]),1]
      pairs<-c(C_up,C_down)
      C_dat[k,1]<-length(C_up)
      C_dat[k,2]<-length(C_down)
      N_tmp<-matrix(rep(sig_exp[k,],length(pairs)),ncol=dim(sig_exp)[2],byrow=T)-mixture[match(pairs,rownames(mixture)),]
      N_up[k,]<-colSums(N_tmp>0,na.rm=TRUE)
      N_down[k,]<-colSums(N_tmp<0,na.rm=TRUE)
    }
    reverse[,1]<-apply(N_up,2,sum)
    reverse[,2]<-apply(N_down,2,sum)
    reverse[,3]<-rep(sum(C_dat[,1]),ncol(sig_exp))
    reverse[,4]<-rep(sum(C_dat[,2]),ncol(sig_exp))
    temp<-apply(reverse,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$estimate)
    up_scores<-cbind(up_scores,temp)
    x <-prop
    y <- (up_scores-min(up_scores))
    data<-data.frame(x=x,y=y)
    k<-nls(y ~a*x^b,data,start=list(a=1, b=1))
    parameter[j,]<-c(names(select_siggene)[j],coef(k))
  }
  return(parameter)
}
