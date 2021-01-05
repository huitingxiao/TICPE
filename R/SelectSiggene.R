#' @title Identify signature genes for each immune cell type compared with cancer cells
#' @description the function is used to curate a specific gene set from overexpressed marker genes as signature genes for each immune cell type.
#'
#' @param cancer_exp the gene expression data set of cancer cell lines. A matrix with row names as Entrez gene ID and columns as cancer cell samples.
#' @param immune_exp the gene expression data set of different types of immune cells. A matrix with row names as Entrez gene ID and columns as immune cell samples.
#' @param marker_geneset list of marker genes of different immune cell subsets collected from literature research and alive tools.
#' @param ratio
#'
#' @return list of signature genes of each immune cell type.
#' @import stats

#' @export
#'
#' @examples
SelectSiggene<-function(cancer_exp,immune_exp,marker_geneset,ratio=0.05){
  sig_list<-list()
  ##identify individual-level differentially upregulated marker genes for each immune cell
  RankComp <-function(expdata,label,gene,freq=0.99){
    Lgene<-length(gene)
    cancer_dat<-expdata[,label==0]
    immune_dat<-expdata[,label==1]
    outlier_dir<-NULL;
    outlier_pvalue<-NULL;
    for (k in 1:Lgene){
      Nnorm<-cancer_dat[k,]
      colN<-dim(cancer_dat)[2]
      colC<-dim(immune_dat)[2]
      N_tmp<-matrix(rep(Nnorm,Lgene),ncol=colN,byrow=T)-cancer_dat
      Nloc_up<-which((rowSums(N_tmp>0)/colN)>freq)
      Nloc_down<-which((rowSums(N_tmp<0)/colN)>freq)
      reverse<-matrix(0,colC,4)
      reverse[,1]<-rep(length(Nloc_up),colC)
      reverse[,2]<-rep(length(Nloc_down),colC)
      Tcanc<-immune_dat[k,]
      if (length(Nloc_up)>0){
        N_tmp<-matrix(rep(Tcanc,length(Nloc_up)),ncol=colC,byrow=T)-immune_dat[Nloc_up,]
        case_p<-colSums(N_tmp<0)
        reverse[,3]<-case_p
      }
      if (length(Nloc_down)>0){
        N_tmpp<-matrix(rep(Tcanc,length(Nloc_down)),ncol=colC,byrow=T)-immune_dat[Nloc_down,]
        case_pp<-colSums(N_tmpp>0)
        reverse[,4]<-case_pp;
      }
      GenePair_sig<-NULL
      GenePair<-rep(0,colC)
      GenePair[which(reverse[,3]>reverse[,4])]<--1
      GenePair[which(reverse[,3]<reverse[,4])]<-1
      tmp=matrix(c(reverse[,1],reverse[,2],reverse[,1]-reverse[,3]+reverse[,4], reverse[,2]-reverse[,4]+reverse[,3]),ncol=4)
      GenePair_sig<-apply(tmp,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$p.value)  #fisher检验
      outlier_dir<-rbind(outlier_dir,GenePair)
      outlier_pvalue<-rbind(outlier_pvalue,GenePair_sig)
    }
    Output<-list(outlier_dir=outlier_dir,outlier_pvalue=outlier_pvalue)
    attr(Output, "class") <- "RCtest"
    return(Output)
  }

  for(j in 1:length(marker_geneset)){
    con_gid<-intersect(intersect(rownames(cancer_exp),rownames(immune_exp)),marker_geneset[[j]])
    loc1<-which(colnames(immune_exp)%in%names(marker_geneset)[j])
    expdata<-as.matrix(cbind(cancer_exp[match(con_gid,rownames(cancer_exp)),],immune_exp[match(con_gid,rownames(immune_exp)),loc1]))
    label<-c(rep(0,dim(cancer_exp)[2]),rep(1,length(loc1)))
    res<-RankComp(expdata,label,con_gid,0.99)
    fdr1<-apply(res$outlier_pvalue,2,function(x){p.adjust(x,method="BH",length(x))})
    sig_p<-matrix(,dim(expdata)[1],2)
    sig_p[,1]<-con_gid
    for(i in 1:dim(expdata)[1]){
      loc2<-which(fdr1[i,]<ratio)
      dir<-res$outlier_dir[i,loc2]
      sig_p[i,2]<-1-pbinom(length(which(dir=="1"))-1,dim(fdr1)[2],0.5)
    }
    fdr2<-p.adjust(sig_p[,2],method = "BH",length(sig_p[,2]))
    siggene<-sig_p[which(fdr2<ratio),1]
    sig_list[[j]]<-siggene
  }

  pair<-t(combn(1:length(marker_geneset),2))
  multgid0<-NULL
  for(i in 1:dim(pair)[1]){
    tmp<-intersect(sig_list[[pair[i,1]]],sig_list[[pair[i,2]]])
    multgid0<-c(multgid0,tmp)
  }
  select_siggene<-lapply(sig_list, function(x){x[-which(x%in%unique(multgid0))]})
  names(select_siggene)<-names(marker_geneset)
  return(select_siggene)
}
