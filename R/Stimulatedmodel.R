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
Stimulatedmodel<-function(combat_edata,select_siggene,stable_pairs){
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
    reverse<-matrix(0,32,4)
    up_scores<-NULL
    mixture<-mix(combat_exp[,which(colnames(combat_exp)%in%names(select_siggene)[j])],norm_control,control)
    rownames(mixture)=rownames(combat_edata)
    gid<-as.numeric(select_siggene[[j]])
    Lgene<-length(which(!gid==''))
    up_pairs<-stable_pairs[which(stable_pairs[,1]%in%gid),]
    up_exp<-mixture[match(up_pairs[,1],rownames(mixture)),,drop=F]-mixture[match(up_pairs[,2],rownames(mixture)),,drop=F]
    mix_up1<-colSums(up_exp>0,na.rm=TRUE)
    mix_down1<-colSums(up_exp<0,na.rm=TRUE)
    down_pairs<-stable_pairs[which(stable_pairs[,2]%in%gid),]
    down_exp<-mixture[match(down_pairs[,2],rownames(mixture)),,drop=F]-mixture[match(down_pairs[,1],rownames(mixture)),,drop=F]
    mix_up2<-colSums(down_exp>0,na.rm=TRUE)
    mix_down2<-colSums(down_exp<0,na.rm=TRUE)
    reverse[,1]<-mix_up1+mix_up2
    reverse[,2]<-mix_down1+mix_down2
    reverse[,3]<-nrow(up_pairs)
    reverse[,4]<-nrow(down_pairs)
    temp<-apply(reverse,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$estimate)
    up_scores<-cbind(up_scores,temp)
    y <- (up_scores-min(up_scores))
    data<-data.frame(x=seq(0.008,0.256, by = 0.008),y=y)
    k<-nls(y ~a*x^b,data,start=list(a=1, b=1))
    parameter[j,]<-c(names(select_siggene)[j],coef(k))
  }
  return(parameter)
}
