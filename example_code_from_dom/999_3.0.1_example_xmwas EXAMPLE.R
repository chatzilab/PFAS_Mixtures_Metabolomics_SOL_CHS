# # PowerPoint of Tutorial
# https://github.com/kuppal2/xMWAS/blob/master/xMWAS_installation_and_tutorial.pptx
# 
# # install xMWAS
# library(devtools)
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   
#   install.packages("BiocManager"); BiocManager::install("BiocInstaller", version = "3.8")
# BiocManager::install(c("GO.db","graph","RBGL","impute","preprocessCore","mixOmics"),dependencies=TRUE);
# 
# install.packages(c("devtools","WGCNA", "snow","igraph","plyr","plsgenomics", "shiny","shinyBS","visNetwork"),dependencies=TRUE, repos="http://cran.r-project.org")

library(xMWAS)

# load example data
data(exh1n1)
data(classlabels_casecontrol) #example classlabels file for case vs control design
data(classlabels_repeatmeasures) #example classlabels file for repeat measures design

exh1n1 %>% str(2)

# store matrices of omics features 
xMat<-exh1n1$metabolome
yMat<-exh1n1$transcrlptome
zMat<-exh1n1$cytokine
classlabels<-exh1n1$classlabels

# output directory
output<-"/Volumes/GoogleDrive/My Drive/20200601_Postdoc_Envn_Omics_NAFLD/HELIX NAFLD/3 Reports/xMWAS Example Output"

#call the run_xmwas() function:
xmwas_res<-run_xmwas(Xome_data=xMat,Yome_data=yMat,Zome_data=zMat,Wome_data=NA,outloc=output, classlabels=classlabels,class_fname=NA,xmwasmethod="pls",plsmode="regression",max_xvar=10000,max_yvar=10000, max_zvar=10000,max_wvar=10000,rsd.filt.thresh=1,corthresh=0.4,keepX=1000,keepY=1000,keepZ=1000,keepW=1000, pairedanalysis=FALSE,optselect=TRUE,rawPthresh=0.05,numcomps=10,net_edge_colors=c("blue","red"), net_node_colors=c("orange", "green","cyan","pink"),Xname="X",Yname="Y",Zname="Z",Wname="W", net_node_shape=c("square","circle","triangle","star"),all.missing.thresh=0,missing.val=0, seednum=100,label.cex=0.2,vertex.size=6,max_connections=NA, centrality_method="eigenvector",use.X.reference=FALSE,removeRda=TRUE,compare.classes=TRUE,class.comparison.allvar=TRUE) 
suppressWarnings(try(sink(file=NULL),silent=TRUE))
