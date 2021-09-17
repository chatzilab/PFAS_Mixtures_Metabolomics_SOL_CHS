# README.md
# Home / GitHub / xia-lab/MetaboAnalystR / R/general_load_libs.R
# 
# 
# R/general_load_libs.R
# In xia-lab/MetaboAnalystR: An R Package for Comprehensive Analysis of Metabolomics Data
# Defines functions load_biocparallel load_RBGL load_graph load_progress load_msnbase load_Rserve load_RSclient load_httr load_stringr load_grid load_ggplot load_data.table load_rgraphwiz load_kegggraph load_pls load_caret load_rsqlite load_rcolorbrewer load_gplots load_reshape load_igraph load_lattice
# Load lattice, necessary for power analysis
load_lattice <- function(){
  suppressMessages(library(lattice))
}

# Load igraph, necessary for network analysis
load_igraph <- function(){
  suppressMessages(library(igraph))
}

# Load reshape, necessary for graphics
load_reshape <- function(){
  suppressMessages(library(reshape))
}

# Load gplots, necessary for heatmap
load_gplots <- function(){
  suppressMessages(library(gplots))
}

# Load R color brewer, necessary for heatmap
load_rcolorbrewer <- function(){
  suppressMessages(library(RColorBrewer))
}

# Load RSQLite, necessary for network analysis
load_rsqlite <- function(){
  suppressMessages(library(RSQLite))
}

# Load caret, necessary for stats module
load_caret <- function(){
  suppressMessages(library(caret))
}

# Load pls, necessary for stats module
load_pls <- function(){
  suppressMessages(library(pls))
}

# Load KEGGgraph
load_kegggraph <- function(){
  suppressMessages(library(KEGGgraph))
}

# Load RGraphviz
load_rgraphwiz <- function(){
  suppressMessages(library(Rgraphviz))
}

# Load data.table
load_data.table <- function(){
  suppressMessages(library(data.table))
}

# Load ggplot2
load_ggplot <- function(){
  suppressMessages(library(ggplot2))
}

# Load gridExtra
load_grid <- function(){
  suppressMessages(library(gridExtra))
  suppressMessages(library(grid))
}

# Load stringr
load_stringr <- function(){
  suppressMessages(library(stringr))
}

# Load httr
load_httr <- function(){
  suppressMessages(library(httr))
}

# Load RSclient
load_RSclient <- function(){
  installed <- c("RSclient") %in% rownames(installed.packages())
  if(installed){
    suppressMessages(library(RSclient))
  }else{
    print("Please install RSclient R package!")
  }
}

# Load Rserve - need to open only 1 instance for R package users
load_Rserve <- function(){
  
  installed <- c("Rserve") %in% rownames(installed.packages())
  
  if(installed){
    # first need to start up an Rserve instance
    suppressMessages(library(Rserve))
    Rserve::Rserve(args = "--no-save")
  }else{
    print("Please install Rserve R package!")
  }
}

# Load MSnbase
load_msnbase <- function(){
  suppressMessages(library(MSnbase))
}

# Load progress
load_progress <- function(){
  suppressMessages(library(progress))
}

# Load graph
load_graph <- function(){
  suppressMessages(library(graph))
}

# Load RBGL
load_RBGL <- function(){
  suppressMessages(library(RBGL))
}

# Load BiocParallel
load_biocparallel <- function(){
  suppressMessages(library(BiocParallel))
}
# xia-lab/MetaboAnalystR documentation built on Aug. 15, 2021, 8:15 a.m.
# R Package Documentation
# rdrr.io home
# R language documentation
# Run R code online
# Browse R Packages
# CRAN packages
# Bioconductor packages
# R-Forge packages
# GitHub packages
# We want your feedback!
 #  Note that we can't provide technical support on individual packages. You should contact the package authors for that.
 # Tweet to @rdrrHQ
 # GitHub issue tracker
 # ian@mutexlabs.com
 # Personal blog
 # 















# internal variables and functions not to be modified by users

## JAG NOTE: 090221: There was an error on line 1946 which is fixed with one line of code that I added.


# This is only for web version 
.on.public.web <- FALSE; # only TRUE when on metaboanalyst web server

# note, this is usually used at the end of a function
# for local, return itself; for web, push to global environment
.set.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    mSet <<- mSetObj;
    return (1);
  }
  return(mSetObj);
}

.get.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    return(mSet)
  }else{
    return(mSetObj);
  }
}

#'Constructs a dataSet object for storing data 
#'@description This functions handles the construction of a mSetObj object for storing data for further processing and analysis.
#'It is necessary to utilize this function to specify to MetaboAnalystR the type of data and the type of analysis you will perform. 
#'@usage InitDataObjects(data.type, anal.type, paired=FALSE)
#'@param data.type The type of data, either list (Compound lists), conc (Compound concentration data), 
#'specbin (Binned spectra data), pktable (Peak intensity table), nmrpeak (NMR peak lists), mspeak (MS peak lists), 
#'or msspec (MS spectra data)
#'@param anal.type Indicate the analysis module to be performed: stat, pathora, pathqea, msetora, msetssp, msetqea, ts, 
#'cmpdmap, smpmap, or pathinteg
#'@param paired Indicate if the data is paired or not. Logical, default set to FALSE
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import methods

InitDataObjects <- function(data.type, anal.type, paired=FALSE){
  
  if(!.on.public.web){
    if(exists("mSet")){
      mSetObj <- .get.mSet(mSet);
      mSetObj$dataSet$type <- data.type;
      mSetObj$dataSet$paired <- paired;
      mSetObj$analSet$type <- anal.type;
      mSetObj<-CleanDataObjects(mSetObj, anal.type);
      return(.set.mSet(mSetObj));
    }
  }
  
  dataSet <- list();
  dataSet$type <- data.type;
  dataSet$design.type <- "regular"; # one factor to two factor
  dataSet$cls.type <- "disc"; # default until specified otherwise
  dataSet$format <- "rowu";
  dataSet$paired <- paired;
  analSet <- list();
  analSet$type <- anal.type;
  Sys.setenv("OMP_NUM_THREADS" = 2); # to control parallel computing for some packages
  Sys.setenv("OPENBLAS_NUM_THREADS" = 2);
  mSetObj <- list();
  mSetObj$dataSet <- dataSet;
  mSetObj$analSet <- analSet;
  mSetObj$imgSet <- list();
  mSetObj$msgSet <- list(); # store various message during data processing
  mSetObj$msgSet$msg.vec <- vector(mode="character");     # store error messages
  mSetObj$cmdSet <- vector(mode="character"); # store R command
  
  if (anal.type == "mummichog") {
    # Define this parameter set to avoid global variable
    # Author: Zhiqiang
    mSetObj$paramSet$mumRT <- NA;
    mSetObj$paramSet$mumRT.type <- NA;
    mSetObj$paramSet$version <- NA;
    mSetObj$paramSet$mumDataContainsPval <- 1;
    mSetObj$paramSet$mode <- NA;
    mSetObj$paramSet$adducts <- NA;
    mSetObj$paramSet$peakFormat <- "mpt";
  } else if (anal.type == "metapaths") {
    # Define this parameter set to avoid global variable
    # Author: Zhiqiang
    paramSet <- list();
    paramSet$mumRT <- NA;
    paramSet$mumRT.type <- NA;
    paramSet$version <- NA;
    paramSet$mumDataContainsPval <- 1;
    paramSet$mode <- NA;
    paramSet$adducts <- NA;
    paramSet$peakFormat <- "mpt";
    paramSet$metaNum <- 0;
    mSetObj$paramSet <- paramSet;
    # This is an empty paramSet, and will be copied for multiple datasets
    dataNMs <- names(mSetObj)[grepl("MetaData",names(mSetObj))];
    if(length(dataNMs)>0){
      for(n in dataNMs){
        mSetObj[[n]] <- NULL;
      }
    }
  }
  
  .init.global.vars(anal.type);
  print("MetaboAnalyst R objects initialized ...");
  return(.set.mSet(mSetObj));
}

# Clean Data Objects to avoid interference
CleanDataObjects <- function(mSetObj, anal.type){
  if(anal.type == "metapaths") {
    # Define this parameter set to avoid global variable
    # Author: Zhiqiang
    paramSet <- list();
    paramSet$mumRT <- NA;
    paramSet$mumRT.type <- NA;
    paramSet$version <- NA;
    paramSet$mumDataContainsPval <- 1;
    paramSet$mode <- NA;
    paramSet$adducts <- NA;
    paramSet$peakFormat <- "mpt";
    paramSet$metaNum <- 0;
    mSetObj$paramSet <- paramSet;
    # This is an empty paramSet, and will be copied for multiple datasets
    dataNMs <- names(mSetObj)[grepl("MetaData",names(mSetObj))];
    if(length(dataNMs)>0){
      for(n in dataNMs){
        mSetObj[[n]] <- NULL;
      }
    }
  }
  return(mSetObj)
}


# for switching from spec to other modules
PrepareSpec4Switch <- function(){
  InitDataObjects("conc", "stat", FALSE);
  #TableFormatCoerce("metaboanalyst_input.csv", "OptiLCMS", "mummichog");
  Read.TextData(NA, "metaboanalyst_input.csv", "colu", "disc");
}

# this is for switching module
UpdateDataObjects <- function(data.type, anal.type, paired=FALSE){
  
  mSetObj <- .get.mSet(NA);
  mSetObj$dataSet$type <- data.type;
  mSetObj$analSet$type <- anal.type;
  .init.global.vars(anal.type);    
  
  # some specific setup 
  if(anal.type == "mummichog"){
    .set.mSet(mSetObj);
    mSetObj<-.init.MummiMSet(mSetObj);
    load("params.rda");        
    mSetObj<-UpdateInstrumentParameters(mSetObj, peakParams$ppm, peakParams$polarity, "yes", 0.02);
    mSetObj<-.rt.included(mSetObj, "seconds");
    #mSetObj<-Read.TextData(mSetObj, "metaboanalyst_input.csv", "colu", "disc");
    mSetObj<-.get.mSet(NA);
  }
  
  return(.set.mSet(mSetObj));
}

.init.MummiMSet <- function(mSetObj) {  
  mSetObj<-SetPeakFormat(mSetObj, "pvalue");
  #TableFormatCoerce("metaboanalyst_input.csv", "OptiLCMS", "mummichog");
  anal.type <<- "mummichog";
  api.base <<- "http://api.xialab.ca";
  err.vec <<- "";
  return(mSetObj)
}

.init.global.vars <- function(anal.type){
  # other global variables
  msg.vec <<- "";
  err.vec <<- "";
  
  # for network analysis
  module.count <<- 0;
  # counter for naming different json file (pathway viewer)
  smpdbpw.count <<- 0; 
  # for mummichog
  #peakFormat <<- "mpt"  
  # mumRT.type <<- "NA";
  
  # raw data processing
  # rawfilenms.vec <<- vector(); # Disable for now
  
  # for meta-analysis
  mdata.all <<- list(); 
  mdata.siggenes <<- vector("list");
  meta.selected <<- TRUE;
  anal.type <<- anal.type;
  
  if(.on.public.web){
    # disable parallel prcessing for public server
    library(BiocParallel);
    register(SerialParam());
  } else {
    if("stat" %in% anal.type | 
       "msetqea" %in% anal.type | 
       "pathqea" %in% anal.type | 
       "roc" %in% anal.type)
      # start Rserve engine for Rpackage
      load_Rserve();
  }
  
  # plotting required by all
  Cairo::CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",italic="Arial:style=Italic",bolditalic = "Arial:style=Bold Italic",symbol = "Symbol")
  
  # sqlite db path for gene annotation
  if(file.exists("/home/glassfish/sqlite/")){ #.on.public.web
    url.pre <<- "/home/glassfish/sqlite/";
  }else if(file.exists("/home/jasmine/Downloads/sqlite/")){ #jasmine's local
    url.pre <<- "/home/jasmine/Downloads/sqlite/";
    #api.base <<- "localhost:8686"
  }else if(file.exists("/Users/soufanom/Documents/Projects/gene-id-mapping/")){ # soufan laptop
    url.pre <<- "/Users/soufanom/Documents/Projects/gene-id-mapping/";
  }else if(file.exists("~/Documents/Projects/gene-id-mapping/")){
    url.pre <<- "~/Documents/Projects/gene-id-mapping/"
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){ # xia local
    url.pre <<- "/Users/xia/Dropbox/sqlite/";
  }else if(file.exists("/media/zzggyy/disk/sqlite/")){
    url.pre <<-"/media/zzggyy/disk/sqlite/"; #zgy local)
  }else if(file.exists("/home/zgy/sqlite/")){
    url.pre <<-"/home/zgy/sqlite/"; #zgy local)
  } else if(file.exists("/home/le/sqlite/")){# le local
    url.pre <<-"/home/le/sqlite/";
  }else if(file.exists("/home/qiang/Music/")){# qiang local
    url.pre <<-"/home/qiang/sqlite/";
  }else{
    url.pre <<- paste0(dirname(system.file("database", "sqlite/GeneID_25Species_JE/ath_genes.sqlite", package="MetaboAnalystR")), "/")
  }
  
  api.base <<- "http://api.xialab.ca"
  #api.base <<- "132.216.38.6:8987"
  
}
#'For two factor time series only
#'@description For two factor time series only
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param design Input the design type
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SetDesignType <-function(mSetObj=NA, design){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$design.type <- tolower(design);
  return(.set.mSet(mSetObj));
}

#'Record R Commands
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param cmd Commands 
#'@export
RecordRCommand <- function(mSetObj=NA, cmd){
  mSetObj <- .get.mSet(mSetObj); 
  mSetObj$cmdSet <- c(mSetObj$cmdSet, cmd);
  return(.set.mSet(mSetObj));
}

SaveRCommands <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);  
  cmds <- paste(mSetObj$cmdSet, collapse="\n");
  pid.info <- paste0("# PID of current job: ", Sys.getpid());
  cmds <- c(pid.info, cmds);
  write(cmds, file = "Rhistory.R", append = FALSE);
}

#'Export R Command History
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
GetRCommandHistory <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj); 
  if(length(mSetObj$cmdSet) == 0){
    return("No commands found");
  }
  return(mSetObj$cmdSet);
}

#'Constructor to read uploaded CSV or TXT files into the dataSet object
#'@description This function handles reading in CSV or TXT files and filling in the dataSet object 
#'created using "InitDataObjects". 
#'@usage Read.TextData(mSetObj=NA, filePath, format, lbl.type)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param filePath Input the path name for the CSV/TXT files to read.
#'@param format Specify if samples are paired and in rows (rowp), unpaired and in rows (rowu),
#'in columns and paired (colp), or in columns and unpaired (colu).
#'@param lbl.type Specify the data label type, either discrete (disc) or continuous (cont).
#'@param nmdr Boolean. Default set to FALSE (data is uploaded by the user and not fetched
#'through an API call to the Metabolomics Workbench).
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}, Jasmine Chong
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

Read.TextData <- function(mSetObj=NA, filePath, format="rowu", 
                          lbl.type="disc", nmdr = FALSE){
  
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$cls.type <- lbl.type;
  mSetObj$dataSet$format <- format;
  
  if(nmdr){
    dat <- qs::qread("nmdr_study.qs")
  }else{
    dat <- .readDataTable(filePath);
  }
  
  if(class(dat) == "try-error" || ncol(dat) == 1){
    AddErrMsg("Data format error. Failed to read in the data!");
    AddErrMsg("Make sure the data table is saved as comma separated values (.csv) format!");
    AddErrMsg("Please also check the followings: ");
    AddErrMsg("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.");
    AddErrMsg("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.");
    AddErrMsg("Make sure sample names and feature (peak, compound) names are unique.");
    AddErrMsg("Missing values should be blank or NA without quote.");
    AddErrMsg("Make sure the file delimeters are commas.");
    return(0);
  }
  
  msg <- NULL;
  
  if(substring(format,4,5)=="ts"){
    # two factor time series data
    if(substring(format,1,3)=="row"){ # sample in row
      msg<-c(msg, "Samples are in rows and features in columns");
      smpl.nms <-dat[,1];
      all.nms <- colnames(dat);
      facA.lbl <- all.nms[2];
      cls.lbl <- facA <- dat[,2]; # default assign facA to cls.lbl in order for one-factor analysis
      facB.lbl <- all.nms[3];
      facB <- dat[,3];
      conc <- dat[,-c(1:3)];
      var.nms <- colnames(conc);
    }else{ # sample in col
      msg<-c(msg, "Samples are in columns and features in rows.");
      all.nms <- dat[,1];
      facA.lbl <- all.nms[1];
      cls.lbl <- facA <- dat[1,-1];
      facB.lbl <- all.nms[2];
      facB <- dat[2,-1];
      var.nms <- dat[-c(1:2),1];
      conc<-t(dat[-c(1:2),-1]);
      smpl.nms <- rownames(conc);
    }
    
    if(mSetObj$dataSet$design.type =="time" | mSetObj$dataSet$design.type =="time0"){
      # determine time factor
      if(!(tolower(facA.lbl) == "time" | tolower(facB.lbl) == "time")){
        AddErrMsg("No time points found in your data");
        AddErrMsg("The time points group must be labeled as <b>Time</b>");
        return(0);
      }
    }
  }else{
    
    if(substring(format,1,3)=="row"){ # sample in row
      msg <- c(msg, "Samples are in rows and features in columns");
      smpl.nms <-dat[,1];
      dat[,1] <- NULL; #remove sample names
      if(lbl.type == "qc"){
        rownames(dat) <- smpl.nms;
        #mSetObj$dataSet$orig <- dat;
        qs::qsave(dat, file="data_orig.qs");
        mSetObj$dataSet$cmpd <- colnames(dat);
        return(1);
      }
      cls.lbl <- dat[,1];
      conc <- dat[,-1, drop=FALSE];
      var.nms <- colnames(conc);
      if(lbl.type == "no"){ #no class label
        cls.lbl <- rep(1, nrow(dat));
        conc <- dat[,, drop=FALSE]; 
        var.nms <- colnames(conc);
      }
    }else{ # sample in col
      msg<-c(msg, "Samples are in columns and features in rows.");
      if(lbl.type == "no"){
        cls.lbl <- rep(1, ncol(dat));
        conc <- t(dat[,-1]);
        var.nms <- dat[,1];
        smpl.nms <- colnames(dat[,-1]);
      }else{
        var.nms <- dat[-1,1];
        dat[,1] <- NULL;
        smpl.nms <- colnames(dat);
        cls.lbl <- dat[1,];
        conc <- t(dat[-1,]);
      }
    }
  }
  
  mSetObj$dataSet$type.cls.lbl <- class(cls.lbl);
  
  msg <- c(msg, "The uploaded file is in comma separated values (.csv) format.");
  
  # try to remove empty line if present
  # identified if no sample names provided
  
  empty.inx <- is.na(smpl.nms) | smpl.nms == ""
  if(sum(empty.inx) > 0){
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), "empty rows</font> were detected and excluded from your data."));
    smpl.nms <- smpl.nms[!empty.inx];
    cls.lbl <-  cls.lbl[!empty.inx];
    conc <- conc[!empty.inx, ];
  }
  
  # try to check & remove empty lines if class label is empty
  # Added by B. Han
  empty.inx <- is.na(cls.lbl) | cls.lbl == ""
  if(sum(empty.inx) > 0){
    if(mSetObj$analSet$type != "roc"){
      msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), "empty labels</font> were detected and excluded from your data."));
      smpl.nms <- smpl.nms[!empty.inx];
      cls.lbl <-  cls.lbl[!empty.inx];
      conc <- conc[!empty.inx, ];
    }else{
      # force all NA to empty string, otherwise NA will become "NA" class label
      cls.lbl[is.na(cls.lbl)] <- "";
      msg <- c(msg, paste("<font color=\"orange\">", sum(empty.inx), "new samples</font> were detected from your data."));
    }
  }
  
  if(anal.type == "roc"){
    if(length(unique(cls.lbl[!empty.inx])) > 2){
      AddErrMsg("ROC analysis is only defined for two-group comparisions!");
      return(0);
    }
  }
  
  # check for uniqueness of dimension name
  if(length(unique(smpl.nms))!=length(smpl.nms)){
    dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse=" ");
    AddErrMsg("Duplicate sample names are not allowed!");
    AddErrMsg(dup.nm);
    return(0);
  }
  
  # try to remove check & remove empty line if feature name is empty
  empty.inx <- is.na(var.nms) | var.nms == "";
  if(sum(empty.inx) > 0){
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), "empty features</font> were detected and excluded from your data."));
    var.nms <- var.nms[!empty.inx];
    conc <- conc[,!empty.inx];
  }
  
  if(length(unique(var.nms))!=length(var.nms)){
    dup.nm <- paste(var.nms[duplicated(var.nms)], collapse=" ");
    AddErrMsg("Duplicate feature names are not allowed!");
    AddErrMsg(dup.nm);
    return(0);
  }
  
  if(anal.type == "mummichog"){
    is.rt <- mSetObj$paramSet$mumRT;
    if(!is.rt){
      mzs <- as.numeric(var.nms);
      if(sum(is.na(mzs) > 0)){
        AddErrMsg("Make sure that feature names are numeric values (mass or m/z)!");
        return(0);
      }
    }
  }
  
  # now check for special characters in the data labels
  if(sum(is.na(iconv(smpl.nms)))>0){
    na.inx <- is.na(iconv(smpl.nms));
    nms <- paste(smpl.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", nms, collapse=" "));
    return(0);
  }
  
  if(sum(is.na(iconv(var.nms)))>0){
    na.inx <- is.na(iconv(var.nms));
    nms <- paste(var.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in feature names!", nms, collapse=" "));
    return(0);
  }
  
  # only keep alphabets, numbers, ",", "." "_", "-" "/"
  url.smp.nms <- CleanNames(smpl.nms);
  names(url.smp.nms) <- smpl.nms;
  
  url.var.nms <- CleanNames(var.nms); # allow space, comma and period
  names(url.var.nms) <- var.nms;
  
  cls.lbl <- ClearStrings(as.vector(cls.lbl));
  
  # now assgin the dimension names
  rownames(conc) <- smpl.nms;
  colnames(conc) <- var.nms;
  
  # check if paired or not
  if(mSetObj$dataSet$paired){
    # save as it is and process in sanity check step
    mSetObj$dataSet$orig.cls <- mSetObj$dataSet$pairs <- cls.lbl;
  } else {
    if(lbl.type == "disc"){
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- as.factor(as.character(cls.lbl));
      
      if(substring(format,4,5)=="ts"){
        
        mSetObj$dataSet$facA.type <- is.numeric(facA);
        mSetObj$dataSet$orig.facA <- mSetObj$dataSet$facA <- as.factor(as.character(facA));
        mSetObj$dataSet$facA.lbl <- facA.lbl;
        
        mSetObj$dataSet$facB.type <- is.numeric(facB);
        mSetObj$dataSet$orig.facB <- mSetObj$dataSet$facB <- as.factor(as.character(facB));
        mSetObj$dataSet$facB.lbl <- facB.lbl;
      }
      
    } else { # continuous
      
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- tryCatch({
        as.numeric(cls.lbl);
      },warning=function(na) {
        print("Class labels must be numeric and continuous!");
        return(0);
      })
      
      if(mSetObj$dataSet$cls == 0){
        AddErrMsg("Class labels must be numeric and continuous!");
        return(0)
      }
    }
  }
  
  # for the current being to support MSEA and MetPA
  if(mSetObj$dataSet$type == "conc"){
    mSetObj$dataSet$cmpd <- var.nms;
  }
  
  mSetObj$dataSet$mumType <- "table";
  mSetObj$dataSet$url.var.nms <- url.var.nms;
  mSetObj$dataSet$url.smp.nms <- url.smp.nms;
  #mSetObj$dataSet$orig <- conc; # copy to be processed in the downstream
  qs::qsave(conc, file="data_orig.qs");
  mSetObj$msgSet$read.msg <- c(msg, paste("The uploaded data file contains ", nrow(conc),
                                          " (samples) by ", ncol(conc), " (", tolower(GetVariableLabel(mSetObj$dataSet$type)), ") data matrix.", sep=""));
  
  return(.set.mSet(mSetObj));
}


#'Read paired peak or spectra files
#'@description This function reads paired peak lists or spectra files. The pair information
#'is stored in a file where each line is a pair and names are separated by ":".
#'@param filePath Set file path
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
ReadPairFile <- function(filePath="pairs.txt"){
  all.pairs <- scan(filePath, what='character', strip.white = T);
  labels <- as.vector(rbind(1:length(all.pairs), -(1:length(all.pairs))));
  all.names <- NULL;
  for(i in 1:length(all.pairs)){
    all.names=c(all.names, unlist(strsplit(all.pairs[i],":", fixed=TRUE), use.names=FALSE));
  }
  names(labels) <- all.names;
  return(labels);
}

#'Save the processed data with class names
#'@description This function saves the processed data with class names as CSV files. 
#'Several files may be generated, the original data, processed data, peak normalized, and/or normalized data. 
#'@param mSetObj Input name of the created mSet Object
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SaveTransformedData <- function(mSetObj=NA){
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.save.data")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_savedata.Rc");    
    }
    return(my.save.data(mSetObj));
  }else{
    return(my.save.data(mSetObj));
  }
}

#' Read an mzTab tab separated file from the passed in file.
#' Adapted from: https://github.com/lifs-tools/rmzTab-m/blob/master/R/MzTabReader.r
#' @param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#' @param filename The name of the mzTab file to parse.
#' @param identifier The identifier to be used when the table is parsed. Use "name"
#' to use the chemical_name, "mass" to use the theoretical_neutral_mass and "sml_id"
#' to use the SML_ID. If the number of missing name and mass entries is greater than 90%,
#' then the SML_ID will be used.
#' @export
Read.mzTab <- function(mSetObj=NA, filename, identifier = "name") {
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.parse.mztab")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_mztab.Rc");    
    }
    return(my.parse.mztab(mSetObj, filename, identifier));
  }else{
    return(my.parse.mztab(mSetObj, filename, identifier));
  }
}

#'Read peak list files
#'@description This function reads peak list files and fills the data into a dataSet object.  
#'For NMR peak lists, the input should be formatted as two-columns containing numeric values (ppm, int).
#'Further, this function will change ppm to mz, and add a dummy 'rt'.
#'For MS peak data, the lists can be formatted as two-columns (mz, int), in which case the function will add a dummy 'rt', or
#'the lists can be formatted as three-columns (mz, rt, int).
#'@usage Read.PeakList(mSetObj=NA, foldername)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param foldername Name of the folder containing the NMR or MS peak list files to read.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#'@export

Read.PeakList<-function(mSetObj=NA, foldername="upload"){
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.parse.peaklist")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_peaks.Rc");    
    }
    return(my.parse.peaklist(mSetObj, foldername));
  }else{
    return(my.parse.peaklist(mSetObj, foldername));
  }
}

#' Adds an error message
#'@description The error message will be printed in all cases.
#'Used in higher functions. 
#'@param msg Error message to print 
#'@export
AddErrMsg <- function(msg){
  err.vec <<- c(err.vec, msg);
  print(msg);
}

# general message only print when running local
AddMsg <- function(msg){
  msg.vec <<- c(msg.vec, msg);
  if(!.on.public.web){
    print(msg);
  }
}

# return the latest message
GetCurrentMsg <- function(){
  return(msg.vec[length(msg.vec)]);
}

GetMetaCheckMsg <- function(mSetObj=NA){
  return(current.msg);
}

#'Plot compound summary
#'change to use dataSet$proc instead of dataSet$orig in
#'case of too many NAs
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param cmpdNm Input the name of the compound to plot
#'@param format Input the format of the image to create
#'@param dpi Input the dpi of the image to create
#'@param width Input the width of the image to create
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotCmpdSummary <- function(mSetObj=NA, cmpdNm, version, format="png", dpi=72, width=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(.on.public.web){
    load_ggplot()
    load_grid()
  }
  
  imgName <- mSetObj$dataSet$url.var.nms[cmpdNm];
  imgName <- paste(imgName, "_", version, "_summary_dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 7.5;
  }else{
    w <- width;
  }
  
  if(substring(mSetObj$dataSet$format,4,5)!="ts"){
    
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height= w*0.65, type=format, bg="white");
    
    # need to consider norm data were edited, different from proc
    smpl.nms <- rownames(mSetObj$dataSet$norm);
    
    mns <- by(as.numeric(mSetObj$dataSet$proc[smpl.nms, cmpdNm]), mSetObj$dataSet$cls, mean, na.rm=T);
    sds <- by(as.numeric(mSetObj$dataSet$proc[smpl.nms, cmpdNm]), mSetObj$dataSet$cls, sd, na.rm=T);
    
    ups <- mns + sds;
    dns <- mns - sds;
    
    # all concentration need start from 0
    y <- c(0, dns, mns, ups);
    
    rg <- range(y) + 0.05 * diff(range(y)) * c(-1, 1)
    pt <- pretty(y)
    
    axp=c(min(pt), max(pt[pt <= max(rg)]),length(pt[pt <= max(rg)]) - 1);
    
    # ggplot alternative
    col <- unique(GetColorSchema(mSetObj$dataSet$cls));
    
    df.orig <- data.frame(value = as.vector(mns), name = levels(mSetObj$dataSet$cls), up = as.vector(ups), down = as.vector(dns))
    p.orig <- ggplot(df.orig, aes(x = name, y = value, fill = name)) + geom_bar(stat = "identity", colour = "black") + theme_bw()
    p.orig <- p.orig + scale_y_continuous(breaks=pt, limits = range(pt)) + ggtitle("Original Conc.")
    p.orig <- p.orig + theme(plot.title = element_text(size = 11, hjust=0.5)) + theme(axis.text.x = element_text(angle=90, hjust=1))
    p.orig <- p.orig + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
    p.orig <- p.orig + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
    p.orig <- p.orig + geom_segment(aes(xend=name, y=up, yend=dns)) + scale_fill_manual(values=col)
    p.orig <- p.orig + theme(plot.margin = margin(t=0.35, r=0.5, b=0.15, l=0.15, "cm"), axis.text = element_text(size=10))
    
    df.norm <- data.frame(value=mSetObj$dataSet$norm[, cmpdNm], name = mSetObj$dataSet$cls)
    p.norm <- ggplot2::ggplot(df.norm, aes(x=name, y=value, fill=name)) + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA) + theme_bw() + geom_jitter(size=1)
    p.norm <- p.norm + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
    p.norm <- p.norm + stat_summary(fun.y=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
    p.norm <- p.norm + scale_fill_manual(values=col) + ggtitle(cmpdNm) + theme(axis.text.x = element_text(angle=90, hjust=1))
    p.norm <- p.norm + ggtitle("Normalized Conc.") + theme(plot.title = element_text(size = 11, hjust=0.5)) 
    p.norm <- p.norm + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
    p.norm <- p.norm + theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.5, "cm"), axis.text = element_text(size=10))
    
    gridExtra::grid.arrange(p.orig, p.norm, ncol=2, top = grid::textGrob(paste(cmpdNm), gp=grid::gpar(fontsize=14, fontface="bold")))
    
    dev.off();
    
  }else if(mSetObj$dataSet$design.type =="time0"){
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=8, height= 6, type=format, bg="white");
    plotProfile(mSetObj, cmpdNm);
    dev.off();
    
  }else{
    if(mSetObj$dataSet$design.type =="time"){ # time trend within phenotype
      out.fac <- mSetObj$dataSet$exp.fac;
      in.fac <- mSetObj$dataSet$time.fac;
      xlab = "Time";
    }else{ # factor a split within factor b
      out.fac <- mSetObj$dataSet$facB;
      in.fac <- mSetObj$dataSet$facA;
      xlab = mSetObj$dataSet$facA.lbl;
    }
    
    # two images per row
    img.num <- length(levels(out.fac));
    row.num <- ceiling(img.num/2)
    
    if(row.num == 1){
      h <- w*5/9;
    }else{
      h <- w*0.5*row.num;
    }
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    
    col <- unique(GetColorSchema(in.fac));
    
    p_all <- list()
    
    for(lv in levels(out.fac)){
      inx <- out.fac == lv;
      df.orig <- data.frame(facA = lv, value = mSetObj$dataSet$norm[inx, cmpdNm], name = in.fac[inx])
      p_all[[lv]] <- df.orig
    }
    
    alldata <- do.call(rbind, p_all)
    
    p.time <- ggplot2::ggplot(alldata, aes(x=name, y=value, fill=name)) + geom_boxplot(outlier.shape = NA, outlier.colour=NA) + theme_bw() + geom_jitter(size=1) 
    p.time <- p.time + facet_wrap(~facA, nrow = row.num) + theme(axis.title.x = element_blank(), legend.position = "none")
    p.time <- p.time + scale_fill_manual(values=col) + theme(axis.text.x = element_text(angle=90, hjust=1))
    p.time <- p.time + ggtitle(cmpdNm) + theme(plot.title = element_text(size = 11, hjust=0.5, face = "bold")) + ylab("Abundance")
    p.time <- p.time + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
    p.time <- p.time + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10)) 
    
    print(p.time)
    dev.off()
  }
  
  if(.on.public.web){
    return(imgName);
  }else{
    return(.set.mSet(mSetObj));
  }
}


##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

GetErrMsg<-function(){
  return(err.vec);
}

#'Save compound name for mapping
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param qvec Input the vector to query
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Setup.MapData <- function(mSetObj=NA, qvec){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$cmpd <- qvec;
  return(.set.mSet(mSetObj));
}

GetMetaInfo <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  if(mSetObj$dataSet$design.type == "regular"){
    return("Group");
  }else{
    return(c(mSetObj$dataSet$facA.lbl, mSetObj$dataSet$facB.lbl));
  }
}

#'Get all group names
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
# 
GetGroupNames <- function(mSetObj=NA, exp.fac=NA){
  mSetObj <- .get.mSet(mSetObj);  
  if(mSetObj$dataSet$design.type == "regular"){
    cls.lbl <- mSetObj$dataSet$prenorm.cls;
    if(mSetObj$analSet$type=="roc"){
      empty.inx <- is.na(cls.lbl) | cls.lbl == "";
      # make sure re-factor to drop level
      my.cls <- factor(cls.lbl[!empty.inx]);
    }else{
      my.cls <- cls.lbl;
    }
  }else{
    if(exp.fac == mSetObj$dataSet$facA.lbl){
      my.cls <- mSetObj$dataSet$facA;  
    }else{
      my.cls <- mSetObj$dataSet$facB;
    }
    my.fac <- exp.fac;
  }
  current.cls <<- my.cls;
  
  if(.on.public.web){
    return(levels(my.cls));
  }else{
    return(.set.mSet(mSetObj));
  }
}

# groups entering analysis
GetNormGroupNames <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  levels(mSetObj$dataSet$cls);
}

GetOrigSmplGroupNames <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  as.character(mSetObj$dataSet$orig.cls);
}

GetPrenormSmplGroupNames <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  as.character(mSetObj$dataSet$prenorm.cls);
}

GetFilesToBeSaved <-function(naviString){
  partialToBeSaved <- c("Rload.RData");
  return(unique(partialToBeSaved));
}

GetExampleDataPath<-function(naviString){
  return(url.pre);
}

Read.TextDataMumMixed <- function(mSetObj=NA, filePath,filePath2, format="rowu", lbl.type="disc", order=NULL){
  
  temp_mSet1 <- Read.TextData(NA, filePath, format, lbl.type)
  temp_mSet1 <- .get.mSet(mSetObj);
  qs::qsave(temp_mSet1, file="data_orig_temp1.qs");
  
  temp_mSet1 <- qs::qread("data_orig_temp1.qs");
  orig1 <- qs::qread("data_orig.qs");
  
  temp_mSet2 <- Read.TextData(NA, filePath2, format, lbl.type)
  temp_mSet2 <-.get.mSet(mSetObj);
  orig2 <- qs::qread("data_orig.qs");
  mSetObj <- temp_mSet1
  mSetObj$dataSet$mumType <- "table";
  mSetObj$dataSet$mode <- "mixed"
  orig <- cbind(orig1, orig2);
  mSetObj$dataSet$pos_inx <- c(rep(TRUE, ncol(orig1)), rep(FALSE, ncol(orig2)))
  qs::qsave(orig, file="data_orig.qs");
  .set.mSet(mSetObj)
  
  return(1);
}


TableFormatCoerce <- function(oriFile = NULL, oriFormat = "Unknown", targetModule = "mummichog", sampleIn = "col"){
  
  df <- read.csv(oriFile);
  df[is.na(df)] <- df[df == ""] <- 0;
  
  if (targetModule == "mummichog") {
    ## This is used to convert the data format from Spectral Processing to Mummichog
    if (oriFormat == "OptiLCMS") {
      if (sampleIn == "col") {
        newFea <- strsplit(df[, 1], "@")
        
        newFea <- lapply(newFea, function(x) {
          return(paste0(x[1], "__", x[2]))
        })
      } else {
        # cannot be in row for OptiLCMS source
      }
      
      df[, 1] <- unlist(newFea)
      df[1, 1] <- "Label"  
    }
  }
  
  write.csv(df, file = oriFile, quote = FALSE, row.names = FALSE)
  
}

################################################################################################
################# First create table of all studies with named metabolites #####################
################################################################################################

#'Function to retrieve all available datasets from the Metabolomics Workbench.
#'@description This function uses the httr R package to make an API
#'call to the Metabolomics Workbench to retrieve a table of
#'all compatible datasets.
#'@usage ListNMDRStudies(mSetObj=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}, Jasmine Chong
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
ListNMDRStudies <- function(mSetObj=NA){
  
  # The Metabolomics Workbench API url
  call <- "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/named_metabolites"
  
  # Use httr::GET to send the request to the Metabolomics Workbench API
  # The response will be saved in query_results
  query_results <- httr::GET(call, encode = "json")
  
  # Check if response is ok
  # 200 is ok! 401 means an error has occured on the user's end.
  if(query_results$status_code!=200){
    AddErrMsg("REST url to Metabolomics Workbench failed!")
    return(0)
  }
  
  # Parse the response into a table
  query_results_text <- content(query_results, "text", encoding = "UTF-8")
  query_results_json <- RJSONIO::fromJSON(query_results_text, flatten = TRUE)
  query_results_table <- t(rbind.data.frame(query_results_json))
  rownames(query_results_table) <- query_results_table[,1]
  
  # Keep studies with > 10 metabolites
  num_metabolites <- as.numeric(query_results_table[,"num_metabolites"])
  keep.inx <- num_metabolites > 10
  query_results_table_keep <- query_results_table[keep.inx, ]
  query_results_table_keep <- subset(query_results_table_keep, select = - c(analysis_id, analysis_type, 
                                                                            num_metabolites,
                                                                            ms_type, units))
  
  query_results_table_keep <- data.frame(query_results_table_keep)
  
  load_stringr()
  query_results_table_keep$lipid <- stringr::str_detect(query_results_table_keep$study_title, fixed("lipid", ignore_case=TRUE))
  
  mSetObj$dataSet$NMDR_studies <- query_results_table_keep
  
  if(!.on.public.web){
    fast.write.csv(query_results_table_keep, "all_metabolomics_workbench_studies.csv")
  }else{
    qs::qsave(query_results_table_keep, "metabolomics_workbench_studies.qs")
  }
  
  return(.set.mSet(mSetObj));
}

##################################################################################################
################# Now using the StudyID download the study dataset as matrix #####################
##################################################################################################

#'Function to retrieve dataset from the Metabolomics Workbench.
#'@description This function uses the httr R package to make an API
#'call to the Metabolomics Workbench to download and save a dataset
#'based on the Study ID into the current working directory.
#'@usage GetNMDRStudy(mSetObj=NA, StudyID)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param StudyID Input the StudyID of the study from the Metabolomics Workbench. 
#'Use the ListNMDRStudies function to obtain a list of all available studies
#'from the Metabolomics Workbench.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}, Jasmine Chong
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetNMDRStudy <- function(mSetObj=NA, StudyID){
  
  mSetObj <- .get.mSet(mSetObj);
  
  check.inx <-  nchar(StudyID) == 8 & grep("ST", StudyID)
  
  if(!check.inx){
    AddErrMsg("Error! Invalid Study ID selected!")
    return(0)
  }
  
  mSetObj$dataSet$NMDR_id <- StudyID
  
  load_httr()
  call <- paste0("https://www.metabolomicsworkbench.org/rest/study/study_id/", StudyID, "/metaboanalyst")
  
  # set timeout to 45 seconds
  query_results <- httr::GET(call, encode = "json", timeout = 45)
  
  if(query_results$status_code!=200){
    AddErrMsg("Error! REST url to Metabolomics Workbench failed!")
    return(0)
  }
  
  # Parse the response into a table
  query_results_text <- content(query_results, "text", encoding = "UTF-8")
  query_study_dataset <- read.delim(text = query_results_text, header = T, check.names = F)
  
  if(nrow(query_study_dataset) < 3){
    AddErrMsg("Error! Selected study does not have enough samples for downstream statistical analyses!")
    return(0)
  }
  
  if(ncol(query_study_dataset) < 3){
    AddErrMsg("Error! Selected study does not have enough features for downstream statistical analyses!")
    return(0)
  }
  
  qs::qsave(query_study_dataset, "nmdr_study.qs")
  
  return(.set.mSet(mSetObj));
}

























### An R package for pathway enrichment analysis for untargeted metabolomics
### based on high-resolution LC-MS platform
### This is based on the mummichog algorithm implemented in python (http://mummichog.org/)
### The goals of developing MummichogR are
### 1) to make this available to the R user community
### 2) high-performance (similar or faster compared to python)
### 3) broader pathways support - by adding support for 21 common organisms based on KEGG pathways
### 4) companion web interface on MetaboAnalyst - the "MS Peaks to Pathways" module
### @authors J. Chong \email{jasmine.chong@mail.mcgill.ca}, J. Xia \email{jeff.xia@mcgill.ca}
### McGill University, Canada
### License: GNU GPL (>= 2)

#'Save adduct names for mapping
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param qvec Input the vector to query
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Setup.AdductData <- function(mSetObj=NA, qvec){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$adduct.list <- qvec;
  mSetObj$dataSet$adduct.custom <- TRUE
  return(.set.mSet(mSetObj));
}

#'Set the peak enrichment method for the MS Peaks to Paths module
#'@description This function sets the peak enrichment method.
#'@param mSetObj Input the name of the created mSetObj.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SetPeakEnrichMethod <- function(mSetObj=NA, algOpt, version="v2"){
  
  mSetObj <- .get.mSet(mSetObj);
  #mSetObj$peaks.alg <- algOpt
  
  mSetObj$paramSet$version <- version
  
  if(algOpt == "gsea"){
    #anal.type <<- "gsea_peaks"
    mSetObj$paramSet$anal.type <- "gsea_peaks";
  }else if(algOpt == "mum"){
    #anal.type <<- "mummichog"
    mSetObj$paramSet$anal.type <- "mummichog";
  }else{
    #anal.type <<- "integ_peaks"
    mSetObj$paramSet$anal.type <- "integ_peaks";
  }
  return(.set.mSet(mSetObj));
}

#'Constructor to read uploaded user files into the mummichog object
#'@description This function handles reading in CSV or TXT files and filling in the mSet object
#' for mummichog analysis. It makes sure that all necessary columns are present.
#'@usage Read.PeakListData(mSetObj=NA, filename = NA)
#'@param mSetObj Input the name of the created mSetObj.
#'@param filename Input the path name for the CSV/TXT files to read.
#'@param meta.anal Logical, TRUE if data will be used for meta-analysis.
#'@param method Input the type of statistical scores included in the
#'mummichog input. "pvalue" for p-values, "es" for effect-sizes, and
#'"both" for both p-values and effect-sizes.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#'@export
Read.PeakListData <- function(mSetObj=NA, filename = NA, meta.anal = FALSE,
                              method = "pvalue") {
  
  mSetObj <- .get.mSet(mSetObj);
  
  file_name <- tools::file_path_sans_ext(basename(filename)) 
  mumDataContainsPval = 1; #whether initial data contains pval or not
  input <- as.data.frame(.readDataTable(filename));
  user_cols <- gsub("[^[:alnum:]]", "", colnames(input));
  mummi.cols <- c("m.z", "p.value", "t.score", "r.t")
  
  if(meta.anal & method %in% c("es", "both")){
    mummi.cols <- c(mummi.cols, "effect.size", "lower.ci", "upper.ci")
  }
  
  # first check if mode included
  mumDataContainsMode <- "mode" %in% user_cols;
  
  if(mumDataContainsMode){
    mode.info <- input$mode  
    input <- subset(input, select=-mode)
    user_cols <- gsub("[^[:alnum:]]", "", colnames(input))
  }
  
  # next check what column names are there
  hit <- "mz" %in% user_cols;
  
  if(sum(hit) < 1){
    AddErrMsg("Missing information, data must contain a 'm.z' column!");
    return(0);
  }
  
  if(length(colnames(input) %in% mummi.cols) == 1){
    peakFormat <- mSetObj$paramSet$peakFormat;
  }else{
    # subset to what's needed for ms peaks
    # then rename columns
    hits2 <- match(gsub("[^[:alnum:]]", "", mummi.cols), user_cols)
    input <- input[, na.omit(hits2)]  
    user_cols <- user_cols[na.omit(hits2)]
    hits.colnames <- match(user_cols, gsub("[^[:alnum:]]", "", mummi.cols))
    user.cols <- mummi.cols[na.omit(hits.colnames)]
    peakFormat <- paste0(substr(sort(user.cols), 1, 1), collapse = "")
    colnames(input) <- user.cols
  }
  
  rt.hit <- "r.t" %in% colnames(input)
  
  if(rt.hit > 0){
    rt = TRUE
  }else{
    rt = FALSE
  }
  
  qs::qsave(input, "mum_raw.qs");
  
  if(!"p.value" %in% colnames(input)){
    mumDataContainsPval <- 0;
    input[,'p.value'] <- rep(0, length=nrow(input))
  }
  
  if(!"t.score" %in% colnames(input)){
    input[,'t.score'] <- rep(0, length=nrow(input))
  }
  
  if(rt){
    mSetObj$dataSet$mummi.orig <- cbind(input$p.value, input$m.z, input$t.score, input$r.t);
    colnames(mSetObj$dataSet$mummi.orig) = c("p.value", "m.z", "t.score", "r.t")
  }else{
    mSetObj$dataSet$mummi.orig <- cbind(input$p.value, input$m.z, input$t.score);
    colnames(mSetObj$dataSet$mummi.orig) = c("p.value", "m.z", "t.score")
  }
  
  if(meta.anal & method %in% c("es", "both")){
    mSetObj$dataSet$mummi.orig <- cbind(mSetObj$dataSet$mummi.orig, effect.size=input$effect.size, 
                                        lower.ci=input$lower.ci, upper.ci=input$upper.ci);
  }
  
  if(mSetObj$dataSet$mode == "positive"){
    mSetObj$dataSet$pos_inx <- rep(TRUE, nrow(mSetObj$dataSet$mummi.orig))
  }else if(mSetObj$dataSet$mode == "negative"){
    mSetObj$dataSet$pos_inx <- rep(FALSE, nrow(mSetObj$dataSet$mummi.orig) )
  }else{ # mixed
    mSetObj$dataSet$pos_inx <- mode.info == "positive"
  }
  
  mSetObj$paramSet$mumRT = rt
  mSetObj$dataSet$mumType = "list";
  mSetObj$msgSet$read.msg <- paste("A total of", length(input$p.value), 
                                   "m/z features were found in your uploaded data.");
  mSetObj$dataSet$fileName <- file_name;
  mSetObj$paramSet$mumDataContainsPval <- mumDataContainsPval;
  mSetObj$paramSet$peakFormat <- peakFormat;
  
  return(.set.mSet(mSetObj));
}

# function to format peak table from mzmine to right format for metaboanalyst
# @format "row" for features in rows, "col" for features in columns.
.format_mzmine_pktable <- function(mSetObj=NA, fileName, format="rowu", group1=NA, group2=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  mzmine_table <- .readDataTable(fileName)
  
  if(class(mzmine_table) == "try-error" || ncol(mzmine_table) == 1){
    AddErrMsg("Data format error. Failed to read in the data!");
    AddErrMsg("Make sure the data table is saved as comma separated values (.csv) format!");
    AddErrMsg("Please also check the followings: ");
    AddErrMsg("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.");
    AddErrMsg("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.");
    AddErrMsg("Make sure sample names and feature (peak, compound) names are unique.");
    AddErrMsg("Missing values should be blank or NA without quote.");
    AddErrMsg("Make sure the file delimeters are commas.");
    return(0);
  }
  
  rownames(mzmine_table) <- mzmine_table[,1]
  mzmine_table <- mzmine_table[,-1]
  
  if(!is.na(group1) & !is.na(group2)){
    if(format == "row"){
      keep.inx <- mzmine_table[1, ] %in% c(group1, group2)
      mzmine_table <- mzmine_table[, keep.inx]
    }else{
      keep.inx <- mzmine_table[, 1] %in% c(group1, group2)
      mzmine_table <- mzmine_table[keep.inx, ]
    }
    newnames <- paste0(tools::file_path_sans_ext(basename(fileName)), "_", group1, "_", group2, ".csv")
  }else{
    
    mzmine_table[1, ] <- tolower(mzmine_table[1, ])
    groups <- unique(unlist(mzmine_table[1, ]))
    
    if(length(groups) > 2){
      AddErrMsg("Only two groups permitted!");
      if("qc" %in% groups){
        AddErrMsg("Remove QCs prior to uploading the mzmine table!");
      }
      return(0);
    }
    
    if(.on.public.web){
      newnames <- "mzmine_peaktable_metaboanalyst.csv"
    }else{
      newnames <- paste0(tools::file_path_sans_ext(basename(fileName)), "_formatted.csv")
    }
  }
  
  if(format == "colu"){ # samples in columns so features are in row
    feats <- gsub("/", "__", rownames(mzmine_table) )
    feats <- sub(".*?__", "", feats )
    feats <- sub("min", "", feats )
    feats <- sub("mz", "", feats )
    feats <- make.unique(feats, sep="")
    rownames(mzmine_table) <- feats
  }else{ # features in columns
    feats <- gsub("/", "__", colnames(mzmine_table) )
    feats <- sub(".*?__", "", feats )
    feats <- sub("min", "", feats )
    feats <- sub("mz", "", feats )
    feats <- make.unique(feats, sep="")
    colnames(mzmine_table) <- feats
  }
  write.csv(mzmine_table, newnames, row.names = TRUE)
  return(.set.mSet(mSetObj))
}

#'Set the peak format for the mummichog analysis
#'@description Set the peak format for mummichog analysis.
#'@param type Input the name of mummichog analysis type, usually 'mpt'.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SetPeakFormat <-function(mSetObj=NA, type = "mpt"){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$paramSet$peakFormat <- type;
  return(.set.mSet(mSetObj))
}

GetPeakFormat <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$paramSet$peakFormat)
}

#'Convert mSetObj to proper format for MS Peaks to Pathways module
#'@description Following t-test analysis or effect size calculation, 
#'this functions converts the results from the mSetObj 
#'to the proper format for mummichog analysis. 
#'@param mSetObj Input the name of the created mSetObj.
#'@param rt Logical, whether or not to include retention time information.
#'@param rds.file Logical, if true, the "annotated_peaklist.rds"
#'must be in the current working directory to get corresponding retention time
#'information for the features. If not, the retention time information
#'will be taken from the feature names. Feature names must be formatted
#'so that the mz and retention time for a single peak is separated by two
#'underscores. For instance, m/z of 410.2148 and retention time of 42.46914 seconds
#'must be formatted as 410.2148__42.46914.
#'@param rt.type Character, input whether retention time is in seconds (default as RT using
#'MetaboAnalystR is seconds) or minutes (as from MZmine).
#'@param test Character, input what statistical values to include in the mummichog input. 
#'For p-values and t-scores only from t-test, use "tt".
#'For log2FC from the fold-change analsis, use "fc".
#'For effect-sizes, use "es".
#'For, p-values, fold-changes and effect sizes, use "all". 
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

Convert2Mummichog <- function(mSetObj=NA, rt=FALSE, rds.file=FALSE, rt.type="seconds", 
                              test="tt", mode=NA){
  
  if(test=="tt"|test=="all"){    
    if(is.null(mSetObj$analSet$tt)){
      AddErrMsg("T-test was not performed!")
      return(0)
    }
    
    tt.pval <- sort(mSetObj$analSet$tt$p.value);
    fdr <- p.adjust(tt.pval, "fdr")
    mz.pval <- names(tt.pval)
    pvals <- cbind(mz.pval, as.numeric(fdr))
    colnames(pvals) <- c("m.z", "p.value")
    
    tt.tsc <- sort(mSetObj$analSet$tt$t.score);
    mz.tsc <- names(tt.tsc)
    tscores <- cbind(mz.tsc, as.numeric(tt.tsc))
    colnames(tscores) <- c("m.z", "t.score")
  }
  
  if(test=="es"|test=="all"){
    
    effect.size <- mSetObj$analSet$effect.size
    
    if(is.null(effect.size)){
      AddErrMsg("Effect size was not calculated!")
      return(0)
    }
    
    mz <- rownames(effect.size)
    esize <- cbind(mz, effect.size)
    colnames(esize) <- c("m.z", "effect.size", "st.dev", "lower.ci", "upper.ci")
  }
  
  if(test=="fc"|test=="all"){
    
    log2fc <- mSetObj$analSet$fc$fc.log
    
    if(is.null(log2fc)){
      AddErrMsg("Fold-change was not calculated!")
      return(0)
    }
    
    mz.fc <- names(log2fc)
    fcs <- cbind(mz.fc, as.numeric(log2fc))
    colnames(fcs) <- c("m.z", "log2.fc")
  }
  
  if(rt & rds.file){
    
    if(!file.exists("annotated_peaklist.rds")){
      AddErrMsg("annotated_peaklist.rds not found in current working directory!")
      return(0)
    }
    
    camera_output <- readRDS("annotated_peaklist.rds")
    mz.cam <- round(camera_output$mz, 5) 
    rt.cam <- round(camera_output$rt, 5) 
    camera <- cbind(mz.cam, rt.cam)
    colnames(camera) <- c("m.z", "r.t")
    
    if(test == "tt"){
      mummi_new <- Reduce(function(x,y) merge(x,y,by="m.z", all = TRUE), list(pvals, tscores, camera))
      complete.inx <- complete.cases(mummi_new[,c("p.value", "t.score", "r.t")]) # filter out m/zs without pval and tscore
    }else if(test == "es"){
      mummi_new <- Reduce(function(x,y) merge(x,y,by="m.z", all = TRUE), list(esize, camera))
      complete.inx <- complete.cases(mummi_new[,c("effect.size", "r.t")]) 
    }else if(test == "all"){
      mummi_new <- Reduce(function(x,y) merge(x,y,by="m.z", all = TRUE), list(pvals, tscores, fcs, esize, camera))
      complete.inx <- complete.cases(mummi_new[,c("p.value", "t.score", "log2.fc", "effect.size", "r.t")]) # filter out m/zs without pval and tscore
    }else if(test=="fc"){
      mummi_new <- Reduce(function(x,y) merge(x,y,by="m.z", all = TRUE), list(fcs, camera))
      complete.inx <- complete.cases(mummi_new[,c("log2.fc", "r.t")]) 
    }
    
    mummi_new <- mummi_new[complete.inx,]
    
  }else{
    
    if(test=="tt"){
      mummi_new <- merge(pvals, tscores)
    }else if(test=="es"){
      mummi_new <- esize
    }else if(test=="all"){
      mummi_new <- Reduce(merge, list(pvals, tscores, fcs, esize))
      mummi_new[] <- lapply(mummi_new, as.character)
    }else if(test=="fc"){
      mummi_new <- fcs
    }
    
    if(rt){ # taking retention time information from feature name itself
      feat_info <- mummi_new[,1]
      feat_info_split <- matrix(unlist(strsplit(feat_info, "__", fixed=TRUE)), ncol=2, byrow=T)
      colnames(feat_info_split) <- c("m.z", "r.t")
      
      if(rt.type == "minutes"){
        rtime <- as.numeric(feat_info_split[,2])
        rtime <- rtime * 60
        feat_info_split[,2] <- rtime
      }
      
      mummi_new <- cbind(feat_info_split, mummi_new[,-1])
    }
  }
  
  if(!is.na(mode)){
    if(mode=="positive"){
      mode <- rep("positive", nrow(mummi_new))
    }else{
      mode <- rep("negative", nrow(mummi_new))
    }
    mummi_new <- cbind(mummi_new, mode)
  }
  
  mummi_new[,1] <- as.numeric(make.unique(as.character(mummi_new[,1]), sep=""))
  
  filename <- paste0("mummichog_input_", Sys.Date(), ".txt")
  input_filename <<- filename;
  write.table(mummi_new, filename, row.names = FALSE)
  
  return(.set.mSet(mSetObj))
}

#'Update the mSetObj with user-selected parameters for MS Peaks to Pathways.
#'@description This functions handles updating the mSet object for mummichog analysis. It is necessary to utilize this function
#'to specify to the organism's pathways to use (libOpt), the mass-spec mode (msModeOpt) and mass-spec instrument (instrumentOpt).
#'@usage UpdateInstrumentParameters(mSetObj=NA, instrumentOpt, msModeOpt, custom=FALSE)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param instrumentOpt Numeric. Define the mass-spec instrument used to perform untargeted metabolomics.
#'@param msModeOpt  Character. Define the mass-spec mode of the instrument used to perform untargeted metabolomics.
#'@param custom Logical, select adducts for mummichog to consider.
#'@param force_primary_ion Character, if "yes", only mz features that match compounds with a primary ion are kept.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

UpdateInstrumentParameters <- function(mSetObj=NA, instrumentOpt, msModeOpt, 
                                       force_primary_ion = "yes", rt_frac = 0.02, 
                                       rt_tol = NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(!is.numeric(instrumentOpt)){
    AddErrMsg("Mass accuracy must be numeric!")
  }else{
    mSetObj$dataSet$instrument <- instrumentOpt;
  }
  
  mSetObj$dataSet$mode <- msModeOpt;
  mSetObj$dataSet$primary_ion <- force_primary_ion;
  mSetObj$dataSet$rt_tol <- rt_tol;
  mSetObj$dataSet$rt_frac <- rt_frac;
  
  return(.set.mSet(mSetObj));
}

#'Sanity Check Data
#'@description SanityCheckData is used for data processing, and performs a basic sanity
#'check of the uploaded data, ensuring that the data is suitable for further analysis.
#'The function ensure that all parameters are properly set based on updated parameters.
#'@usage SanityCheckMummichogData(mSetObj=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#'@export

SanityCheckMummichogData <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  if(mSetObj$dataSet$mumType == "table"){
    mSetObj$paramSet$mumDataContainsPval <- 1;
    orig.data<- qs::qread("data_orig.qs");
    l = sapply(colnames(orig.data),function(x) return(unname(strsplit(x,"/", fixed=TRUE)[[1]][1])))
    colnames(orig.data) <- l;
    qs::qsave(orig.data, file="data_orig.qs");
    
    if(.on.public.web){
      return(SanityCheckData(NA));
    }else{
      mSetObj <- SanityCheckData(mSetObj)
      return(.set.mSet(mSetObj));
    }
  }
  
  msg.vec <- NULL;
  mSetObj$mum_nm <- "mummichog_query.json"
  mSetObj$mum_nm_csv <- "mummichog_pathway_enrichment.csv"
  ndat <- mSetObj$dataSet$mummi.orig;
  pos_inx = mSetObj$dataSet$pos_inx
  ndat <- data.frame(cbind(ndat, pos_inx), stringsAsFactors = FALSE)
  
  rawdat <- qs::qread("mum_raw.qs");
  
  if(mSetObj$paramSet$mumRT){
    na.num <- sum(is.na(ndat$r.t))
    # filter out any reads w. NA RT
    if(na.num>0){
      na.inx <- which(is.na(ndat$r.t))
      ndat <- ndat[-na.inx,]
      msg.vec <- c(msg.vec, paste("A total of <b>", na.num, "</b> mz features with missing retention times were removed."));
    }
  }
  
  # another check to ensure no missing or NA values
  missing.inx <- apply(ndat, 2, function(x) any(is.na(x)))
  
  if(any(missing.inx)){
    AddErrMsg("NA values found in the uploaded data!")
    return(0)
  }
  
  read.msg <- mSetObj$msgSet$read.msg
  
  # sort mzs by p-value
  ord.inx <- order(ndat[,1]);
  ndat <- ndat[ord.inx,]; # order by p-vals
  
  # filter based on mz
  mznew <- ndat[,2];
  
  # trim to fit within 50 - 2000
  my.inx <- mznew > 50 & mznew < 2001;
  trim.num <- sum(!my.inx);
  range = paste0(min(ndat[,1], "-", max(ndat[,1])))
  if(length(unique(mSetObj[["dataSet"]][["pos_inx"]])) > 1){
    colNMadd <- "and mode";
    colnumadd <- 1;
  } else {
    colnumadd <- 0;
    colNMadd <- NULL;
  }
  
  msg.vec <- c(msg.vec, paste("The instrument's mass accuracy is <b>", mSetObj$dataSet$instrument , "</b> ppm."));
  msg.vec <- c(msg.vec, paste("The instrument's analytical mode is <b>", mSetObj$dataSet$mode , "</b>."));
  msg.vec <- c(msg.vec, paste("The uploaded data contains <b>", length(colnames(rawdat)) + colnumadd, "</b> columns."));
  
  peakFormat <- mSetObj$paramSet$peakFormat;
  
  if (peakFormat == "rmp") {
    msg.vec <- c(msg.vec, paste("The peaks are ranked by <b>p-values</b>."));
  } else if (peakFormat == "rmt") {
    msg.vec <- c(msg.vec, paste("The peaks are ranked by <b>t-scores</b>."));
  }
  
  msg.vec <- c(msg.vec, paste("The column headers of uploaded data are <b>", paste(colnames(rawdat),collapse=", "), colNMadd ,"</b>."));
  msg.vec <- c(msg.vec, paste("The range of m/z peaks is trimmed to 50-2000. <b>", trim.num, "</b> features have been trimmed."));
  
  if(trim.num > 0){
    ndat <- ndat[my.inx,]
    #msg.vec <- c(msg.vec, paste("A total of", trim.num, "were excluded to fit within mz range of 50-2000"));
  }
  
  # remove duplicated mzs (make unique)
  dup.inx <- duplicated(ndat);
  dup.num <- sum(dup.inx);
  
  if(dup.num > 0){
    ndat <- ndat[!dup.inx,];
    msg.vec <- c(msg.vec, paste("A total of <b>", dup.num, "</b> duplicated mz features were removed."));
  }
  
  # make mzs unique
  mzs <- ndat[,2]
  # ensure features are unique
  mzs_unq <- mzs[duplicated(mzs)]
  set.seed(123);
  while(length(mzs_unq)>0){
    mzs[duplicated(mzs)] <- sapply(mzs_unq, function(x) paste0(x, sample(1:9, 1, replace = TRUE)));
    mzs_unq <- mzs[duplicated(mzs)]
  }
  
  ndat[,2] <- mzs
  ref_mzlist <- ndat[,2];
  
  # set up expression (up/dn)
  tscores <- as.numeric(ndat[,3]);
  names(tscores) <- ref_mzlist;
  
  # set up rt
  if(mSetObj$paramSet$mumRT){
    
    retention_time <- as.numeric(ndat[,4]);
    names(retention_time) <- ref_mzlist;
    mSetObj$dataSet$pos_inx <- as.numeric(ndat$pos_inx) == 1;
    mSetObj$dataSet$ret_time <- retention_time;
    
    if(is.na(mSetObj$dataSet$rt_tol)){
      rt_tol <- max(mSetObj$dataSet$ret_time) * mSetObj$dataSet$rt_frac 
      print(paste0("Retention time tolerance is ", rt_tol))
      mSetObj$dataSet$rt_tol <- rt_tol
    }
    
  }else{
    mSetObj$dataSet$pos_inx <- as.numeric(ndat$pos_inx) == 1;
  }
  
  ref.size <- length(ref_mzlist);
  
  msg.vec <- c(msg.vec, paste("A total of ", ref.size, "input mz features were retained for further analysis."));
  
  if(ref.size > 20000){
    msg.vec <- c(msg.vec, "There are too many input features, the performance may be too slow.");
  }
  
  if(ref.size < 100){
    AddErrMsg("There are too few m/z features. Ensure that all of your m/z features have been uploaded!");
    return(0)
  }
  
  if(min(ndat[,"p.value"])<0 || max(ndat[,"p.value"])>1){
    msg.vec <- c(msg.vec, "Please make sure the p-values are between 0 and 1.");
  }
  
  mSetObj$msgSet$check.msg <- c(mSetObj$msgSet$check.msg, read.msg, msg.vec);
  mSetObj$dataSet$mummi.proc <- ndat;
  mSetObj$dataSet$expr_dic <- tscores;
  mSetObj$dataSet$ref_mzlist <- ref_mzlist;
  
  return(.set.mSet(mSetObj));
}

#'Set the cutoff for mummichog analysis
#'@description Set the p-value cutoff for mummichog analysis.
#'@param mSetObj Input the name of the created mSetObj.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

SetMummichogPvalFromPercent <- function(mSetObj=NA, fraction){
  
  mSetObj <- .get.mSet(mSetObj);
  peakFormat <- mSetObj$paramSet$peakFormat;
  
  if(peakFormat %in% c("rmp", "rmt")){
    maxp <- 0;
  }else{
    pvals <- c(0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
    ndat <- mSetObj$dataSet$mummi.proc;
    n <- floor(fraction*length(ndat[,"p.value"]))
    cutoff <- ndat[n+1,1]
    if(!any(pvals <= cutoff)){
      maxp <- 0.00001
    }else{
      maxp <- max(pvals[pvals <= cutoff])
    }
  }
  mSetObj$dataSet$cutoff <- maxp
  .set.mSet(mSetObj);
  return(SetMummichogPval("NA", maxp));
}

#'Set the cutoff for mummichog analysis
#'@description Set the p-value cutoff for mummichog analysis.
#'@param mSetObj Input the name of the created mSetObj.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

SetMummichogPval <- function(mSetObj=NA, cutoff){
  
  mSetObj <- .get.mSet(mSetObj);
  
  msg.vec <- NULL;
  mSetObj$dataSet$cutoff <- cutoff;
  ndat <- mSetObj$dataSet$mummi.proc;
  msg.vec <- c(msg.vec, "Only a small percentage (below 10%) peaks in your input peaks should be significant.");
  msg.vec <- c(msg.vec, "The algorithm works best for <u>200~500 significant peaks in 3000~7000 total peaks</u>.");
  
  ref_mzlist <- ndat[,2];
  if(mSetObj$paramSet$mumDataContainsPval == 1){
    my.inx <- ndat[,1] <= cutoff;
    input_mzlist <- ref_mzlist[my.inx];
    # note, most of peaks are assumed to be not changed significantly, more than 25% should be warned
    
  }else{
    inxToKeep = length(ref_mzlist)/10
    if(inxToKeep > 500){
      inxToKeep = 500;
    }
    input_mzlist <- ref_mzlist[1:inxToKeep];
  }
  
  sig.size <- length(input_mzlist);
  sig.part <- round(100*sig.size/length(ref_mzlist),2);
  
  if(sig.part > 25){
    msg.vec <- c(msg.vec, paste("<font color=\"orange\">Warning: over", sig.part, "percent were significant based on your cutoff</font>."));
    msg.vec <- c(msg.vec, "You should adjust p-value cutoff to control the percentage");
  }else{
    msg.vec <- c(msg.vec, paste("A total of", sig.size, "or", sig.part, "percent signficant mz features were found based on the selected p-value cutoff:", cutoff));
  }
  
  if(sig.size > 2000){
    msg.vec <- c(msg.vec, "There are too many significant features based on the current cutoff, possibly too slow.");
  }else if(sig.size == 0){
    AddErrMsg("No significant features were found based on the current cutoff! Failed to perform analysis.");
    return(0);
  }else if(sig.size < 30){
    msg.vec <- c(msg.vec, "The number of significant features is small based on the current cutoff, possibly not accurate.");
  }
  
  mSetObj$dataSet$input_mzlist <- input_mzlist;
  mSetObj$dataSet$N <- sig.size;
  
  if(.on.public.web){
    mSet <<- mSetObj;
    return(round(sig.part, 0))
  } else {
    return(.set.mSet(mSetObj));
  }
}

#' For meta-analysis, set the p-value
#' cutoff used to define the new
#' input_mzlist
SetMetaPeaksPvals <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$cutoff <- cutoff;
  return(.set.mSet(mSetObj));
}

#'Function to perform peak set enrichment analysis
#'@description This is the main function that performs either the mummichog
#'algorithm, GSEA, or both for peak set enrichment analysis. 
#'@usage PerformPSEA(mSetObj=NA, lib, libVersion, minLib, permNum = 100)
#'@param mSetObj Input the name of the created mSetObj object. 
#'@param lib Input the name of the organism library, default is hsa_mfn. 
#'@param libVersion Input the version of the KEGG pathway libraries ("current" or "old").
#'@param minLib Numeric, input the minimum number of metabolites needed to consider the pathway 
#'or metabolite set. 
#'@param permNum Numeric, input the number of permutations to perform. Default is 100.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import qs

PerformPSEA <- function(mSetObj=NA, lib, libVersion, minLib = 3, permNum = 100){
  require("plyr");
  mSetObj <- .get.mSet(mSetObj);
  mSetObj <- .setup.psea.library(mSetObj, lib, libVersion, minLib);
  version <- mSetObj$paramSet$version;
  mSetObj$dataSet$paramSet <- mSetObj$paramSet;
  
  if(mSetObj$paramSet$mumRT & version=="v2"){
    mSetObj <- .init.RT.Permutations(mSetObj, permNum)
  } else {
    mSetObj <- .init.Permutations(mSetObj, permNum)
  }
  
  if(class(mSetObj) != "list"){
    if(mSetObj == 0){
      AddErrMsg("MS Peaks to Paths analysis failed! Likely not enough m/z to compound hits for pathway analysis!")
      return(0)
    }
  }
  return(.set.mSet(mSetObj));
}

#' Internal function to perform PSEA, no retention time
#' @importFrom data.table setDF
.init.Permutations <- function(mSetObj, permNum){
  anal.type0 <- mSetObj$paramSet$anal.type;
  
  if(anal.type0 == "mummichog"){
    mSetObj <- .perform.mummichogPermutations(mSetObj, permNum);
    mSetObj <- .compute.mummichogSigPvals(mSetObj);
  } else if(anal.type0 == "gsea_peaks") {
    mSetObj <- .compute.mummichog.fgsea(mSetObj, permNum);
  } else {
    # need to perform mummichog + gsea then combine p-values
    mSetObj <- .compute.mummichog.fgsea(mSetObj, permNum);
    mSetObj <- .perform.mummichogPermutations(mSetObj, permNum);
    mSetObj <- .compute.mummichogSigPvals(mSetObj);
    pathResults <- vector("list")
    
    pathResults$ora.all <- mSetObj$mummi.resmat
    pathResults$gsea.all <- mSetObj$mummi.gsea.resmat
    
    path.names.all <- lapply(pathResults, rownames)
    path.intersect <- Reduce(intersect, path.names.all)
    
    # remove paths not found by all three - complete.cases
    path.intersected <- lapply(pathResults, function(x) x[row.names(x) %in% path.intersect,])
    #li_2 <- lapply(seq_along(path.intersected), function(i) {colnames(path.intersected[[i]]) <- paste0(colnames(path.intersected[[i]]), names(path.intersected)[[i]]) ; path.intersected[[i]] } )
    path2 <- data.table::setDF(Reduce(merge, 
                                      lapply(path.intersected, 
                                             data.table::data.table, 
                                             keep.rownames = TRUE)))
    
    mum.df <- path2[, c("rn", "Pathway total", "Hits.total", "Hits.sig", "FET", "P_val")]
    rownames(mum.df) <- mum.df$rn;
    # combine p-values
    # replace 1 with 0.999 and 0 with low number
    mum.df[,5:6][mum.df[,5:6]==1] <- 0.999
    mum.df[,5:6][mum.df[,5:6]==0] <- .Machine$double.xmin
    
    combo.all <- apply(mum.df[,5:6], 1, function(x) sumlog(x))
    
    #extract p-values
    all.ps <- unlist(lapply(combo.all, function(z) z["p"]))
    df.combo <- as.matrix(mum.df[,2:6])
    dfcombo <- round(cbind(df.combo, all.ps), 5)
    colnames(dfcombo) <- c("Total_Size", "Hits", "Sig_Hits", "Mummichog_Pvals", "GSEA_Pvals", "Combined_Pvals")
    ord.inx <- order(dfcombo[,6]);
    dfcombo <- signif(as.matrix(dfcombo[ord.inx, ]), 4);
    
    write.csv(dfcombo, "mummichog_integ_pathway_enrichment.csv", row.names = TRUE)
    mSetObj$integ.resmat <- dfcombo;
    matched_cpds <- names(mSetObj$cpd_exp)
    colnames(mum.df)[1] = "pathways";
    inx2<-na.omit(match(mum.df$pathways[ord.inx], mSetObj$pathways$name))
    filt_cpds <- lapply(inx2, function(f) { mSetObj$pathways$cpds[f] })
    
    cpds <- lapply(filt_cpds, function(x) intersect(unlist(x), matched_cpds))
    cpds_sig <- lapply(filt_cpds, function(x) intersect(unlist(x), mSetObj$input_cpdlist))
    
    mSetObj$path.nms <- rownames(dfcombo)
    mSetObj$path.hits <- convert2JsonList(cpds)
    mSetObj$path.pval <- as.numeric(dfcombo[,6])
    mSetObj <<- mSetObj;
    matched_res <- qs::qread("mum_res.qs");
    json.res <- list(
      cmpd.exp = mSetObj$cpd_exp,
      path.nms = rownames(dfcombo),
      hits.all = convert2JsonList(cpds),
      hits.all.size = as.numeric(dfcombo[,2]),
      hits.sig = convert2JsonList(cpds_sig),
      hits.sig.size = as.numeric(dfcombo[,3]),
      mum.p = as.numeric(dfcombo[,4]),
      gsea.p = as.numeric(dfcombo[,5]),
      comb.p = as.numeric(dfcombo[,6]),
      peakToMet = mSetObj$cpd_form_dict,
      peakTable = matched_res
    );
    
    json.mat <- RJSONIO::toJSON(json.res, .na='null');
    sink(mSetObj$mum_nm);
    cat(json.mat);
    sink();
  }
  return(mSetObj);
}

#' Internal function to perform PSEA, with RT
.init.RT.Permutations <- function(mSetObj, permNum){
  
  anal.type0 <- mSetObj$paramSet$anal.type;
  
  if(anal.type0 == "mummichog"){
    mSetObj <- .perform.mummichogRTPermutations(mSetObj, permNum);
    mSetObj <- .compute.mummichogRTSigPvals(mSetObj);
  }else if(anal.type0 == "gsea_peaks"){
    mSetObj <- .compute.mummichog.RT.fgsea(mSetObj, permNum);
  }else{
    # need to perform mummichog + gsea then combine p-values
    mSetObj <- .compute.mummichog.RT.fgsea(mSetObj, permNum); # OK (jg)
    mSetObj <- .perform.mummichogRTPermutations(mSetObj, permNum); # OK (jg)
    mSetObj <- .compute.mummichogRTSigPvals(mSetObj); # OK (jg)
    
    pathResults <- vector("list")
    
    pathResults$ora.all <- mSetObj$mummi.resmat
    pathResults$gsea.all <- mSetObj$mummi.gsea.resmat
    
    path.names.all <- lapply(pathResults, rownames)
    path.intersect <- Reduce(intersect, path.names.all)
    
    # remove paths not found by all three - complete.cases
    path.intersected <- lapply(pathResults, function(x) x[row.names(x) %in% path.intersect,])
        #li_2 <- lapply(seq_along(path.intersected), function(i) {colnames(path.intersected[[i]]) <- paste0(colnames(path.intersected[[i]]), names(path.intersected)[[i]]) ; path.intersected[[i]] } )
    path2 <- data.table::setDF(Reduce(merge, lapply(path.intersected, data.table::data.table, keep.rownames = TRUE)))
    path2 <- unique(path2) #(Added JAG)
    
    mum.df <- path2[, c(1:4, 6, 12)]
    rownames(mum.df) <- mum.df$rn 
    mum.df <- mum.df[,-1] 
    colnames(mum.df) <- c("Total_Size", "Hits", "Sig_Hits", "Mummichog_Pvals", "GSEA_Pvals") 
    
    # combine p-values
    combo.all <- apply(mum.df[,c("Mummichog_Pvals", "GSEA_Pvals")], 1, function(x) sumlog(x))
    
    #extract p-values
    all.ps <- unlist(lapply(combo.all, function(z) z["p"]))
    df.combo <- as.matrix(mum.df)
    df.combo <- round(cbind(df.combo, all.ps), 5)
    colnames(df.combo) <- c("Total_Size", "Hits", "Sig_Hits", "Mummichog_Pvals", "GSEA_Pvals", "Combined_Pvals")
    ord.inx <- order(df.combo[,6]);
    df.combo <- signif(as.matrix(df.combo[ord.inx, ]), 4);
    write.csv(df.combo, "mummichog_integ_pathway_enrichment.csv", row.names = TRUE)
    mSetObj$integ.resmat <- df.combo
    
    ## transform ecpd to cpd for json files
    # for all cpds
    total_ecpds <- unique(mSetObj$total_matched_ecpds) #all matched compounds
    current.mset <- mSetObj$pathways$emp_cpds[match(rownames(df.combo), mSetObj$pathways$name)]
    ecpds <- lapply(current.mset, function(x) intersect(x, total_ecpds)); #pathways & all ref ecpds
    cpds <- lapply(ecpds, function(x) unique(unlist(mSetObj$ecpd_cpd_dict[match(x, names(mSetObj$ecpd_cpd_dict))])) )
    
    # for sig ecpds
    qset <- unique(unlist(mSetObj$input_ecpdlist));
    current.mset <- mSetObj$pathways$emp_cpds[match(rownames(df.combo), mSetObj$pathways$name)]
    feats <- lapply(current.mset, function(x) intersect(x, qset));
    cpds_feats <- lapply(feats, function(x) unique(unlist(mSetObj$ecpd_cpd_dict[match(x, names(mSetObj$ecpd_cpd_dict))])) )  
    
    # now make exp vec for all compounds
    cpds2ec <- mSetObj$cpd_ecpd_dict
    cpds.all <- unique(unlist(mSetObj$ecpd_cpd_dict[match(total_ecpds, names(mSetObj$ecpd_cpd_dict))]))
    cpds.exp <- sapply(cpds.all, function(x) sapply(seq_along(x), function(i) mean(mSetObj$ec_exp[match(unique(unlist(cpds2ec[match(x[[i]], names(cpds2ec))])), names(mSetObj$ec_exp))]) ) )
    
    mSetObj$path.nms <- rownames(df.combo)
    mSetObj$path.hits <- convert2JsonList(cpds)
    mSetObj$path.pval <- as.numeric(df.combo[,6])
    
    mSetObj <<- mSetObj;
    matched_res <- qs::qread("mum_res.qs");
    
    json.res <- list(
      cmpd.exp = cpds.exp,
      path.nms = rownames(df.combo),
      hits.all = convert2JsonList(cpds),
      hits.all.size = as.numeric(df.combo[,2]),
      hits.sig = convert2JsonList(cpds_feats),
      hits.sig.size = as.numeric(df.combo[,3]),
      mum.p = as.numeric(df.combo[,4]),
      gsea.p = as.numeric(df.combo[,5]),
      comb.p = as.numeric(df.combo[,6]),
      peakToMet = mSetObj$cpd_form_dict,
      peakTable = matched_res
    );
    
    json.mat <- RJSONIO::toJSON(json.res, .na='null');
    sink(mSetObj$mum_nm);
    cat(json.mat);
    sink();
  }
  return(mSetObj);
}

# Internal function to set up library for PSEA
.setup.psea.library <- function(mSetObj = NA, 
                                lib, 
                                libVersion, 
                                minLib=3, 
                                metaAnalysis = FALSE,
                                metaLevel = "pathway", 
                                combine.level,
                                pval.method, 
                                es.method, 
                                rank.metric, 
                                mutual.feats = TRUE){
  
  version <- mSetObj$paramSet$version;
  filenm <- paste(lib, ".qs", sep="")
  biocyc <- grepl("biocyc", lib);
  
  if(!is.null(mSetObj$curr.cust)){
    
    if(biocyc){
      user.curr <- mSetObj$curr.map$BioCyc;
    }else{
      user.curr <- mSetObj$curr.map$KEGG;
    }
    
    currency <<- user.curr;
    
    if(length(currency)>0){
      mSetObj$mummi$anal.msg <- c("Currency metabolites were successfully uploaded!")
    }else{
      mSetObj$mummi$anal.msg <- c("Errors in currency metabolites uploading!")
    }
  }
  
  if(.on.public.web){
    mum.url <- paste("../../libs/mummichog/", filenm, sep="");
    print(paste("Adding mummichog library:", mum.url));
    mummichog.lib <- try(qs::qread(mum.url), silent = TRUE)
    if(class(mummichog.lib) == "try-error"){
      AddErrMsg("The database you have selected is not matched/found!")
      return(0)
    }
    
  }else{
    if(!file.exists(filenm)){
      mum.url <- paste("https://www.metaboanalyst.ca/resources/libs/mummichog/", filenm, sep="");
      download.file(mum.url, destfile = filenm, method="libcurl", mode = "wb")
      mummichog.lib <- qs::qread(filenm);
    }else{
      mummichog.lib <- qs::qread(filenm);
    }
  }
  
  if(!is.null(mSetObj$dataSet$adduct.custom)){
    mw <- mummichog.lib$cpd.lib$mw;
    new_adducts <- new_adduct_mzlist(mSetObj, mw);
    ms0 <- new_adducts$pos;
    
    cpd.lib <- list(
      mz.matp = new_adducts$pos,
      mz.matn = new_adducts$neg,
      mw = mummichog.lib$cpd.lib$mw,
      id = mummichog.lib$cpd.lib$id,
      name = mummichog.lib$cpd.lib$name
    );
    
  }else{
    ms <- mummichog.lib$cpd.lib$adducts[["positive"]];
    cpd.lib <- list(
      mz.matp = mummichog.lib$cpd.lib$adducts[["positive"]],
      mz.matn = mummichog.lib$cpd.lib$adducts[["negative"]],
      mw = mummichog.lib$cpd.lib$mw,
      id = mummichog.lib$cpd.lib$id,
      name = mummichog.lib$cpd.lib$name
    );
  }
  
  cpd.treep <- mummichog.lib$cpd.tree[["positive"]];
  cpd.treen <- mummichog.lib$cpd.tree[["negative"]];
  
  # build empirical compound library after
  path.length <- sapply(mummichog.lib$pathways$cpds, length)
  
  min.inx <- which(path.length >= minLib)
  
  cleaned.pathways <- vector("list")
  cleaned.pathways$cpds <- mummichog.lib$pathways$cpds[min.inx]
  cleaned.pathways$id <- mummichog.lib$pathways$id[min.inx]
  cleaned.pathways$name <- mummichog.lib$pathways$name[min.inx]
  
  mSetObj$pathways <- cleaned.pathways;
  
  if(metaAnalysis & metaLevel %in% c("cpd", "ec")){
    mSetObj <- .search.compoundLibMeta(mSetObj, cpd.lib, cpd.treep, cpd.treen, metaLevel, combine.level,
                                       pval.method, es.method, rank.metric, mutual.feats);
  }else{
    mSetObj <- .search.compoundLib(mSetObj, cpd.lib, cpd.treep, cpd.treen);
  }
  
  if(mSetObj$paramSet$mumRT & version=="v2"){
    # only for empirical compounds
    if(metaLevel %in% c("ec", "pathway", "pooled")){
      # map cpds to empirical cpds
      cleaned.pathways$emp_cpds <- lapply(cleaned.pathways$cpds, 
                                          function(x) {
                                            unique(unlist(mSetObj$cpd_ecpd_dict[na.omit(match(x, names(mSetObj$cpd_ecpd_dict)))]))
                                          })
      
      # delete emp_cpds, cpds and names with no emp_cpds
      null.inx <- sapply(cleaned.pathways$emp_cpds, is.null)
      
      new_pathways <- vector("list");
      
      new_pathways$cpds <- cleaned.pathways$cpds[!null.inx]
      new_pathways$name <- cleaned.pathways$name[!null.inx]
      new_pathways$emp_cpds <- cleaned.pathways$emp_cpds[!null.inx]
      
      mSetObj$pathways <- new_pathways;
    }
  }
  
  mSetObj$lib.organism <- lib; #keep track of lib organism for sweave report
  return(mSetObj);
}

# version 2
# internal function for searching compound library
.search.compoundLib <- function(mSetObj, 
                                cpd.lib, 
                                cpd.treep, 
                                cpd.treen){
  
  ref_mzlist <- as.numeric(mSetObj$dataSet$ref_mzlist);
  print(paste0("Got ", length(ref_mzlist), " mass features."))
  pos_inx <- mSetObj$dataSet$pos_inx;
  ref_mzlistp <- ref_mzlist[pos_inx];
  ref_mzlistn <- ref_mzlist[!pos_inx];
  version <- mSetObj$paramSet$version;
  
  # for empirical compounds
  if(mSetObj$paramSet$mumRT & version=="v2"){
    ord_rt <- rank(mSetObj$dataSet$ret_time, ties.method = "random")
    ret_time_pos <- mSetObj$dataSet$ret_time[pos_inx];
    ret_time_rank_pos <- ord_rt[pos_inx]
    ret_time_neg <- mSetObj$dataSet$ret_time[!pos_inx]
    ret_time_rank_neg <- ord_rt[!pos_inx]
    rt_tol <- mSetObj$dataSet$rt_tol
    rt_tol_rank <- length(ref_mzlist) * mSetObj$dataSet$rt_frac
  } else {
    # add fake RT
    ret_time_pos <- rep(1, length(ref_mzlistp))
    ret_time_rank_pos <- rep(1, length(ref_mzlistp))
    ret_time_neg <- rep(1, length(ref_mzlistn))
    ret_time_rank_neg <- rep(1, length(ref_mzlistn))
  }
  
  modified.statesp <- colnames(cpd.lib$mz.matp);
  modified.statesn <- colnames(cpd.lib$mz.matn);
  my.tolsp <- mz_tolerance(ref_mzlistp, mSetObj$dataSet$instrument);
  my.tolsn <- mz_tolerance(ref_mzlistn, mSetObj$dataSet$instrument);
  
  # get mz ladder (pos index)
  self.mzsp <- floor(ref_mzlistp);
  all.mzsp <- cbind(self.mzsp-1, self.mzsp, self.mzsp+1);
  
  self.mzsn <- floor(ref_mzlistn);
  all.mzsn <- cbind(self.mzsn-1, self.mzsn, self.mzsn+1);
  
  # matched_res will contain detailed result (cmpd.id. query.mass, mass.diff) for all mz;
  # use a high-performance variant of list
  matched_resp <- myFastList();
  matched_resn <- myFastList();
  
  if(mSetObj$dataSet$mode != "negative"){
    for(i in 1:length(ref_mzlistp)){
      mz <- ref_mzlistp[i];
      rt <- ret_time_pos[i];
      rt_rank <- ret_time_rank_pos[i];
      my.tol <- my.tolsp[i];
      all.mz <- all.mzsp[i,];
      pos.all <- as.numeric(unique(unlist(cpd.treep[all.mz])));
      
      for(pos in pos.all){
        id <- cpd.lib$id[pos];
        mw.all <- cpd.lib$mz.matp[pos,]; #get modified mzs
        diffs <- abs(mw.all - mz); #modified mzs - mz original
        hit.inx <- which(diffs < my.tol);
        if(length(hit.inx)>0){
          for(spot in 1:length(hit.inx)){
            hit.pos <- hit.inx[spot];# need to match all
            index <- paste(mz, id, rt, hit.pos, sep = "___");
            matched_resp$add(index, c(i, id, mz, rt, rt_rank, mw.all[hit.pos], modified.statesp[hit.pos], diffs[hit.pos])); #replaces previous when hit.inx>1
          }
        }
      }
    }
  }
  
  all.mzsn <<- all.mzsn
  
  if (mSetObj$dataSet$mode != "positive") {
    for(i in 1:length(ref_mzlistn)){
      mz <- ref_mzlistn[i];
      rt <- ret_time_neg[i];
      rt_rank <- ret_time_rank_neg[i];
      my.tol <- my.tolsn[i];
      all.mz <- all.mzsn[i,];
      pos.all <- as.numeric(unique(unlist(cpd.treen[all.mz])));
      
      for(pos in pos.all){
        id <- cpd.lib$id[pos]; # position of compound in cpd.tree
        mw.all <- cpd.lib$mz.matn[pos,]; #get modified mzs
        diffs <- abs(mw.all - mz); #modified mzs - mz original
        hit.inx <- which(diffs < my.tol);
        if(length(hit.inx)>0){
          for(spot in 1:length(hit.inx)){
            hit.pos <- hit.inx[spot];# need to match all
            index <- paste(mz, id, rt, hit.pos, sep = "___"); #name in fast_list
            matched_resn$add(index, c(i, id, mz, rt, rt_rank, mw.all[hit.pos], modified.statesn[hit.pos], diffs[hit.pos])); #replaces previous when hit.inx>1
          }
        }
      }
    }
  }
  
  # convert to regular list
  if (mSetObj$dataSet$mode == "mixed") {
    
    matched_resn <- matched_resn$as.list();
    matched_resp <- matched_resp$as.list();
    
    neg_matches <- length(matched_resn) > 0
    pos_matches <- length(matched_resp) > 0
    
    if(!neg_matches & !pos_matches){
      msg.vec <<- "No compound matches from upload peak list!"
      return(0)
    }
    
    if(neg_matches){
      matched_resn <- data.frame(matrix(unlist(matched_resn), nrow=length(matched_resn), byrow=T), stringsAsFactors = FALSE);
    }
    
    if(pos_matches){
      matched_resp <- data.frame(matrix(unlist(matched_resp), nrow=length(matched_resp), byrow=T), stringsAsFactors = FALSE);
    }
    
    if(neg_matches & pos_matches){ # both w. matches
      matched_res <- rbind(matched_resp, matched_resn)
    }else if(neg_matches & !pos_matches){ # only neg w. matches
      matched_res <- matched_resn
    }else{ # only pos w. matches
      matched_res <- matched_resp
    }
    
  } else if(mSetObj$dataSet$mode == "positive") {
    matched_resp <- matched_resp$as.list();
    
    if(is.null(unlist(matched_resp))){
      msg.vec <<- "No compound matches from upload peak list!"
      return(0)
    }
    
    matched_resp <- data.frame(matrix(unlist(matched_resp), nrow=length(matched_resp), byrow=T), stringsAsFactors = FALSE);
    matched_res <- matched_resp;
    
  } else {
    matched_resn <- matched_resn$as.list();
    if(is.null(unlist(matched_resn))){
      msg.vec <<- "No compound matches from upload peak list!"
      return(0)
    }
    
    matched_resn <- data.frame(matrix(unlist(matched_resn), nrow=length(matched_resn), byrow=T), stringsAsFactors = FALSE);
    matched_res <- matched_resn
  }
  
  # re-order columns for output
  matched_res <- matched_res[, c(3,2,7,8,4,5)];
  colnames(matched_res) <- c("Query.Mass", "Matched.Compound", "Matched.Form", "Mass.Diff", "Retention.Time", "RT.Rank");
  
  if(!mSetObj$paramSet$mumRT & version=="v2"){
    matched_res <- matched_res[,-(5:6)]
  }
  
  #print(paste0(length(unique(matched_res[,2])), " matched compounds! cpd2mz"))
  
  # now create empirical compounds if necessary!
  # 1 compound matches to multiple m/z, filter by RT 
  if(mSetObj$paramSet$mumRT & version=="v2"){
    start <- Sys.time()
    # mz, ion
    empirical.cpd.list <- split(matched_res[,c(1,3,5,6)], matched_res[,2]); # split mz, ion and rt by compound
    empirical.cpds2cpds <- vector(length=(length(empirical.cpd.list)), "list")
    names(empirical.cpds2cpds) <- names(empirical.cpd.list)
    
    # for each compound, if multiple matches, split into ECpds if > RT tolerance - rt_tol
    for(i in 1:length(empirical.cpd.list)){
      
      mzs <- empirical.cpd.list[[i]]$Query.Mass
      ions <- empirical.cpd.list[[i]]$Matched.Form
      rts <- empirical.cpd.list[[i]]$Retention.Time
      rt.rank <- empirical.cpd.list[[i]]$RT.Rank
      cpds <- names(empirical.cpd.list)[i]
      
      # first, for each compound, determine ECs among matched ions
      if(length(mzs)>1){ # if multiple ECs per compound
        
        # first group together to create empirical cpds by rt
        rts <- as.numeric(rts)
        names(rts) <- paste0(mzs, ";", ions, ";", rts, ";", cpds)
        rts <- sort(rts)
        
        # second, group together to create empirical cpds by rt rank
        rt.ranks <- as.numeric(rt.rank)
        names(rt.ranks) <- paste0(mzs, ";", ions, ";", rts, ";", cpds)
        rt.ranks <- sort(rt.ranks)
        
        split.inx <- c(0, cumsum(Reduce("&", list(abs(diff(rts)) > rt_tol, abs(diff(rt.ranks)) > rt_tol_rank) )))
        
        # need to deal w. multiple rts but only 1 EC
        if(length(unique(split.inx)) > 1){
          e.cpds <- split(rts, split.inx)
          empirical.cpds2cpds[[i]] <- lapply(e.cpds, names)
        }else{
          empirical.cpds2cpds[[i]] <- paste0(names(rts), collapse="__")
        }
        
      }else{ # if only 1 EC per compound
        empirical.cpds2cpds[[i]] <- paste0(mzs, ";", ions, ";", rts, ";", cpds)
      }
    }
    
    initial_ecs <- unlist(empirical.cpds2cpds, recursive=FALSE)
    names(initial_ecs) <- paste0("EC", 1:length(initial_ecs))
    print(paste0(length(initial_ecs), " inital ECs created!"))
    
    # second, merge ECs if same m/z and form - append compounds
    try <- melt(initial_ecs)
    try2 <- strsplit(as.character(try[,1]), split="__", fixed=TRUE) # deals with multiple rts belonging to 1 EC
    try2.df <- data.frame(value=unlist(try2), L1 = rep(try$L1, sapply(try2, length)))
    
    info <- strsplit(as.character(try2.df[,1]), split=";")
    df_ecs <- data.frame(ec=as.character(try2.df[,2]), mz=sapply(info, `[[`, 1), form=sapply(info, `[[`, 2), rt = sapply(info, `[[`, 3), cpd = sapply(info, `[[`, 4), stringsAsFactors = F)
    df_ecs$str_row_inx <- paste(df_ecs$mz, df_ecs$form, df_ecs$rt, sep = "___")
    qs::qsave(df_ecs, "initial_ecs.qs")
    merged_ecs <- aggregate(. ~ str_row_inx, df_ecs, paste, collapse=";")
    
    # cleaning the df
    # merged_ecs$ec <- sapply(strsplit(merged_ecs$ec, ";", fixed=TRUE), function(x) unlist(x)[1]) - keep as long name
    merged_ecs$mz <- sapply(strsplit(merged_ecs$mz, ";", fixed=TRUE), function(x) unique(unlist(x)))
    merged_ecs$form <- sapply(strsplit(merged_ecs$form, ";", fixed=TRUE), function(x) unique(unlist(x)))
    merged_ecs$rt <- sapply(strsplit(merged_ecs$rt, ";", fixed=TRUE), function(x) unique(unlist(x)))
    print(paste0(length(unique(merged_ecs$ec)), " merged ECs identified!"))
    
    # third, check if primary ion is present
    # needs to be per EC!
    if(mSetObj$dataSet$primary_ion=="yes"){
      
      ecs <- unique(merged_ecs$ec);
      
      # function to group ECs and verify if contains primary ion
      new_info <- lapply(ecs, function(x) { 
        new_info <- merged_ecs[which(merged_ecs$ec == x),] # subset merged_ecs to rows containing ECx
        primary.inx <- length(intersect(new_info$form, primary_ions))
        
        if(primary.inx>0){
          new_info <- new_info
        }else{
          new_info <- NULL
        }
        new_info
      })  
      
      final_ecs <- do.call(args=new_info, what=rbind)[,-1]
      
    }else{
      final_ecs <- merged_ecs[,-1]
    }
    
    colnames(final_ecs) <- c("Empirical.Compound", "Query.Mass", "Matched.Form", "Retention.Time", "Matched.Compound")
    
    # transform to long format
    cpd_split <- strsplit(as.character(final_ecs$Matched.Compound), ";", fixed=TRUE)
    reps <- pmax(lengths(cpd_split))
    df2 <- final_ecs[rep(1:nrow(final_ecs), reps), 1:4]
    df2$Matched.Compound <- unlist(mapply(function(x,y) c(x, rep(NA, y)), cpd_split, reps-lengths(cpd_split)))
    
    matched_res <- merge(matched_res, df2)
    matched_res <- matched_res[,-6] #rm rt rank
    matched_res[,6] <- as.character(matched_res[,6])
    
    # now deal with the fact that if at least one EC overlap, need to count as same EC per compound...
    my_final_cpds <- aggregate(. ~ Matched.Compound, matched_res, paste, collapse="___")
    my_final_cpds_list <- lapply(split(my_final_cpds$Empirical.Compound, my_final_cpds$Matched.Compound), unlist)
    
    cpd2ec1 <- lapply(seq_along(my_final_cpds_list), function(x) { # function used to make grouping of ecs per cpd
      
      ecs <- unlist(strsplit(my_final_cpds_list[[x]], "___", fixed=TRUE))
      
      if(length(ecs) > 1){
        ecs.list <- as.list(strsplit(ecs, ";", fixed=TRUE))
        library(igraph)
        m = sapply(ecs.list, function(x) sapply(ecs.list, function(y) length(intersect(x,y))>0))
        g = igraph::groups(components(graph_from_adjacency_matrix(m)))
        ecs <- paste0(sapply(g, function(z) paste0(ecs[z], collapse = "|") ), collapse = "___")
      }
      ecs
    })
    
    names(cpd2ec1) <- names(my_final_cpds_list)
    
    update_ecs <- lapply(seq_along(cpd2ec1), function(z) {
      
      ecs.old <- unlist(strsplit(my_final_cpds_list[[z]], "___", fixed=TRUE))
      ecs.new <- unlist(strsplit(cpd2ec1[[z]], "___", fixed=TRUE))
      
      for(i in seq_along(ecs.new)){
        pattern <- ecs.new[i]
        pattern_vec <- unlist(strsplit(pattern, "\\|"))
        up.pattern <- paste0(unique(pattern_vec), collapse = "|")
        ecs.old[ ecs.old %in% pattern_vec  ] <- up.pattern
      }
      
      ecs.old <- paste0(ecs.old, collapse = "___")
      ecs.old
    })
    
    updated_ecs <- do.call(rbind, update_ecs)
    my_final_cpds$Empirical.Compound <- updated_ecs
    
    new_dt <- data.table::data.table(my_final_cpds)
    new_dt <- new_dt[, list(Query.Mass = unlist(strsplit(as.character(Query.Mass), "___", fixed=TRUE)), 
                            Matched.Form = unlist(strsplit(as.character(Matched.Form), "___", fixed=TRUE)),
                            Retention.Time = unlist(strsplit(as.character(Retention.Time), "___", fixed=TRUE)),
                            Mass.Diff = unlist(strsplit(as.character(Mass.Diff), "___", fixed=TRUE)),
                            Empirical.Compound = unlist(strsplit(as.character(Empirical.Compound), "___", fixed=TRUE))),
                     by = Matched.Compound]
    
    matched_res <- data.frame(Query.Mass = new_dt$Query.Mass, Matched.Compound = new_dt$Matched.Compound, Matched.Form = new_dt$Matched.Form,
                              Retention.Time = new_dt$Retention.Time, Mass.Diff = new_dt$Mass.Diff, Empirical.Compound = new_dt$Empirical.Compound, stringsAsFactors = FALSE)
    
    # make EC names
    ec <- matched_res$Empirical.Compound
    ec.unique <- unique(matched_res$Empirical.Compound)
    
    for(i in seq_along(ec.unique)){
      ec <- replace(ec, grep(paste0("\\b", ec.unique[i], "\\b"), ec, perl=TRUE), paste0("EC000", i))
    }
    
    #for(i in seq_along(ec.unique)){
    #  ec <- gsub(paste0("\\b", ec.unique[i], "\\b"), paste0("EC000", i), ec)
    #}
    
    matched_res$Empirical.Compound <- gsub("\\|.*", "", ec)
    end <- Sys.time()
    totaltime <- end-start
    print(paste0(length(unique(matched_res$Empirical.Compound)), " empirical compounds identified in ", totaltime, " seconds."))
  }
  
  fast.write.csv(matched_res, file="mummichog_matched_compound_all.csv", row.names=FALSE);
  qs::qsave(matched_res, "mum_res.qs");
  
  # now update expr. profile
  matched_mz <- matched_res[,1];
  matched_ts <- mSetObj$dataSet$expr_dic[matched_mz];
  
  if(mSetObj$paramSet$mumRT & version=="v2"){ # RT need to be in EC space
    # first create ecpd to expression dict
    ec.exp.mat <- data.frame(key=matched_res[,6], value=as.numeric(matched_ts), stringsAsFactors = F)
    ec_exp_dict <- Convert2Dictionary(ec.exp.mat);
    ec.exp.vec <- unlist(lapply(ec_exp_dict, max));
    
    # also need to make cpd_exp_dict for KEGG network view
    exp.mat <- data.frame(key=matched_res[,2], value=as.numeric(matched_ts));
    cpd_exp_dict <- Convert2Dictionary(exp.mat);
    
    # ecpd to cpd dict
    cpd_ecpd_dict <- Convert2Dictionary(matched_res[,c(2,6)])
    ecpd_cpd_dict <- Convert2Dictionary(matched_res[,c(6,2)])
    
    # now mz 2 ecpd dict
    mz2cpd_dict <- Convert2Dictionary(matched_res[,c(1,2)]); #indexed/named by mz
    mz2ec_dict <- Convert2Dictionary(matched_res[,c(1,6)])
    ec2mz_dict <- Convert2Dictionary(matched_res[,c(6,1)])
    
    # save to mSetObj
    mSetObj$ec_exp_dict <- ec_exp_dict
    mSetObj$cpd_exp_dict <- cpd_exp_dict;
    mSetObj$ec_exp <- ec.exp.vec
    mSetObj$mz2cpd_dict <- mz2cpd_dict;
    mSetObj$mz2ec_dict <- mz2ec_dict
    mSetObj$ec2mz_dict <- ec2mz_dict
    mSetObj$ecpd_cpd_dict <- ecpd_cpd_dict
    mSetObj$cpd_ecpd_dict <- cpd_ecpd_dict
    mSetObj$cpd_ecpd_counts <- cpd2ec1
    
    # now do matching to identify significant input_ecpdlist
    refmz <- names(mz2ec_dict)
    hits.index <- which(refmz %in% as.character(mSetObj$dataSet$input_mzlist));
    ec1 <- unique(unlist(mz2ec_dict[hits.index]));
    mSetObj$input_ecpdlist <- ec1;
    mSetObj$total_matched_ecpds <- unique(as.vector(matched_res$Empirical.Compound));
  } else {
    # get the expression profile for each 
    exp.mat <- data.frame(key=matched_res[,2], value=as.numeric(matched_ts));
    cpd_exp_dict <- Convert2Dictionary(exp.mat);
    # create average exp
    exp.vec <- unlist(lapply(cpd_exp_dict, mean));
    
    # now need to get the mapping from mz to compound id (one mz can have 0, 1, or more id hits)
    mz2cpd_dict <- Convert2Dictionary(matched_res[,c(1,2)]); #indexed/named by mz
    cpd2mz_dict <- Convert2Dictionary(matched_res[,c(2,1)]); # indexed/named by id
    
    # now do matching to identify significant input_cpdlist
    refmz <- names(mz2cpd_dict)
    hits.index <- which(refmz %in% as.character(mSetObj$dataSet$input_mzlist));
    cpd1 <- unique(unlist(mz2cpd_dict[hits.index]));
    cpd1 <- cpd1[!(cpd1 %in% currency)];
    
    mSetObj$mz2cpd_dict <- mz2cpd_dict;
    mSetObj$cpd_exp_dict <- cpd_exp_dict;
    mSetObj$cpd_exp <- exp.vec;
    mSetObj$cpd2mz_dict <- cpd2mz_dict;
    mSetObj$input_cpdlist <- cpd1;
    mSetObj$total_matched_cpds <- unique(as.vector(matched_res$Matched.Compound));
  }
  
  form.mat <- cbind(matched_res[,2], matched_res[,3]);
  cpd_form_dict <- Convert2Dictionary(form.mat);
  mSetObj$cpd_form_dict <- cpd_form_dict;
  
  return(mSetObj);
}

# version 2
# internal function for searching compound library
.search.compoundLibMeta <- function(mSetObjMeta, cpd.lib, cpd.treep, cpd.treen, metaLevel = "cpd",
                                    combine.level = "pvalue", pval.method = "fisher", es.method = "fixed",
                                    rank.metric = "mean", mutual.feats = TRUE){
  
  metaFiles <- unique(metaFiles)
  
  metaMsetObj <- vector("list")
  version <- mum.version ##TODO: There must be an error bc mum.version is not global variable any more
  
  # first do compound mapping
  for(meta_file in seq_along(metaFiles)){
    
    mSetObj <- qs::qread(metaFiles[meta_file])
    ref_mzlist <- as.numeric(mSetObj$dataSet$ref_mzlist);
    print(paste0("Got ", length(ref_mzlist), " mass features."))
    pos_inx <- mSetObj$dataSet$pos_inx;
    ref_mzlistp <- ref_mzlist[pos_inx];
    ref_mzlistn <- ref_mzlist[!pos_inx];
    
    # for empirical compounds
    if(mSetObj$paramSet$mumRT){ 
      ret_time_pos <- mSetObj$dataSet$ret_time[pos_inx];
      ret_time_neg <- mSetObj$dataSet$ret_time[!pos_inx]
      rt_tol <- mSetObj$dataSet$rt_tol
    }else{
      # add fake RT
      ret_time_pos <- rep(1, length(ref_mzlistp))
      ret_time_neg <- rep(1, length(ref_mzlistn))
    }
    
    modified.statesp <- colnames(cpd.lib$mz.matp);
    modified.statesn <- colnames(cpd.lib$mz.matn);
    my.tolsp <- mz_tolerance(ref_mzlistp, mSetObj$dataSet$instrument);
    my.tolsn <- mz_tolerance(ref_mzlistn, mSetObj$dataSet$instrument);
    
    # get mz ladder (pos index)
    self.mzsp <- floor(ref_mzlistp);
    all.mzsp <- cbind(self.mzsp-1, self.mzsp, self.mzsp+1);
    
    self.mzsn <- floor(ref_mzlistn);
    all.mzsn <- cbind(self.mzsn-1, self.mzsn, self.mzsn+1);
    
    # matched_res will contain detailed result (cmpd.id. query.mass, mass.diff) for all mz;
    # use a high-performance variant of list
    matched_resp <- myFastList();
    matched_resn <- myFastList();
    
    if(mSetObj$dataSet$mode != "negative"){
      for(i in seq_along(ref_mzlistp)){
        mz <- ref_mzlistp[i];
        rt <- ret_time_pos[i];
        my.tol <- my.tolsp[i];
        all.mz <- all.mzsp[i,];
        pos.all <- as.numeric(unique(unlist(cpd.treep[all.mz])));
        
        for(pos in pos.all){
          id <- cpd.lib$id[pos];
          mw.all <- cpd.lib$mz.matp[pos,]; #get modified mzs
          diffs <- abs(mw.all - mz); #modified mzs - mz original
          hit.inx <- which(diffs < my.tol);
          if(length(hit.inx)>0){
            for(spot in seq_along(hit.inx)){
              hit.pos <- hit.inx[spot];# need to match all
              index <- paste(mz, id, rt, hit.pos, sep = "_");
              matched_resp$add(index, c(i, id, mz, rt, mw.all[hit.pos], modified.statesp[hit.pos], diffs[hit.pos])); #replaces previous when hit.inx>1
            }
          }
        }
      }
    }
    
    all.mzsn <<- all.mzsn
    
    if(mSetObj$dataSet$mode != "positive"){
      for(i in seq_along(ref_mzlistn)){
        mz <- ref_mzlistn[i];
        rt <- ret_time_neg[i];
        my.tol <- my.tolsn[i];
        all.mz <- all.mzsn[i,];
        pos.all <- as.numeric(unique(unlist(cpd.treen[all.mz])));
        
        for(pos in pos.all){
          id <- cpd.lib$id[pos]; # position of compound in cpd.tree
          mw.all <- cpd.lib$mz.matn[pos,]; #get modified mzs
          diffs <- abs(mw.all - mz); #modified mzs - mz original
          hit.inx <- which(diffs < my.tol);
          if(length(hit.inx)>0){
            for(spot in seq_along(hit.inx)){
              hit.pos <- hit.inx[spot];# need to match all
              index <- paste(mz, id, rt, hit.pos, sep = "_"); #name in fast_list
              matched_resn$add(index, c(i, id, mz, rt, mw.all[hit.pos], modified.statesn[hit.pos], diffs[hit.pos])); #replaces previous when hit.inx>1
            }
          }
        }
      }
    }
    
    # convert to regular list
    if(mSetObj$dataSet$mode == "mixed"){
      
      matched_resn <- matched_resn$as.list();
      matched_resp <- matched_resp$as.list();
      
      neg_matches <- length(matched_resn) > 0
      pos_matches <- length(matched_resp) > 0
      
      if(!neg_matches & !pos_matches){
        msg.vec <<- "No compound matches from upload peak list!"
        return(0)
      }
      
      if(neg_matches){
        matched_resn <- data.frame(matrix(unlist(matched_resn), nrow=length(matched_resn), byrow=T), stringsAsFactors = FALSE);
      }
      
      if(pos_matches){
        matched_resp <- data.frame(matrix(unlist(matched_resp), nrow=length(matched_resp), byrow=T), stringsAsFactors = FALSE);
      }
      
      if(neg_matches & pos_matches){ # both w. matches
        matched_res <- rbind(matched_resp, matched_resn)
      }else if(neg_matches & !pos_matches){ # only neg w. matches
        matched_res <- matched_resn
      }else{ # only pos w. matches
        matched_res <- matched_resp
      }
      
    }else if(mSetObj$dataSet$mode == "positive"){
      matched_resp <- matched_resp$as.list();
      
      if(is.null(unlist(matched_resp))){
        msg.vec <<- "No compound matches from upload peak list!"
        return(0)
      }
      
      matched_resp <- data.frame(matrix(unlist(matched_resp), nrow=length(matched_resp), byrow=T), stringsAsFactors = FALSE);
      matched_res <- matched_resp
    }else{
      matched_resn <- matched_resn$as.list();
      
      if(is.null(unlist(matched_resn))){
        msg.vec <<- "No compound matches from upload peak list!"
        return(0)
      }
      
      matched_resn <- data.frame(matrix(unlist(matched_resn), nrow=length(matched_resn), byrow=T), stringsAsFactors = FALSE);
      matched_res <- matched_resn
    }
    
    # re-order columns for output
    matched_res <- matched_res[, c(3,2,6,7,4)];
    colnames(matched_res) <- c("Query.Mass", "Matched.Compound", "Matched.Form", "Mass.Diff", "Retention.Time");
    
    if(metaLevel == "cpd"){
      fileName = paste0("mummichog_matched_compound_", mSetObj$dataSet$fileName, "_all.csv")
      fast.write.csv(matched_res, file=fileName, row.names=FALSE);
    }
    
    # objects save to meta msetObj
    metafile <- tools::file_path_sans_ext(metaFiles[meta_file])
    metaMsetObj[[metafile]] <- vector("list", 5)
    
    # now update expr. profile
    matched_mz <- matched_res[,1];
    matched_ts <- mSetObj$dataSet$expr_dic[matched_mz];
    
    if(combine.level %in% c("pvalue", "both", "pool")){
      pvals <- mSetObj$dataSet$mummi.proc$p.value
      names(pvals) <- mSetObj$dataSet$mummi.proc$m.z
      matched_pvals <- pvals[matched_mz]
      metaMsetObj[[metafile]]$matched_pvals <- matched_pvals;
    } 
    
    if(combine.level %in% c("es", "both")){
      es <- mSetObj$dataSet$mummi.proc[,c("effect.size", "lower.ci", "upper.ci")]
      rownames(es) <- mSetObj$dataSet$mummi.proc$m.z
      match.inx <- match(matched_mz, rownames(es))
      matched_es <- es[match.inx,]
      metaMsetObj[[metafile]]$matched_es <- matched_es;
    }
    
    metaMsetObj[[metafile]]$mumResTable <- matched_res;
    metaMsetObj[[metafile]]$ref_mzlist <- mSetObj$dataSet$ref_mzlist 
    metaMsetObj[[metafile]]$input_mzlist <- mSetObj$dataSet$input_mzlist
    metaMsetObj[[metafile]]$matched_ts <- matched_ts;
    metaMsetObj[[metafile]]$mumRT <-mSetObj$paramSet$mumRT
  }
  
  # second fill in p-value and effect size information
  mSetObj <- mSetObjMeta;
  
  matched_res <- lapply(metaMsetObj, "[[", "mumResTable")
  matched_res <- data.table::rbindlist(matched_res, idcol = TRUE)
  matched_ts <- unlist(lapply(metaMsetObj, "[[", "matched_ts"))
  matched_res <- cbind(matched_res, matched_ts)
  
  if(combine.level %in% c("pvalue", "both", "pool")){
    matched_pval <- unlist(lapply(metaMsetObj, "[[", "matched_pvals"))
  }else{ # initialize empty
    matched_pval <- rep(NA, length(matched_ts))
  }
  
  if(combine.level %in% c("es", "both")){
    matched_es <- lapply(metaMsetObj, "[[", "matched_es")
    matched_es <- Reduce(rbind, matched_es)
  }else{ # initialize empty
    matched_es <- matrix("1", nrow = length(matched_ts), ncol=3)
    colnames(matched_es) <- c("effect.size", "lower.ci", "upper.ci") 
  }
  
  matched_res <- cbind(matched_res, Matched.Pval = matched_pval, matched_es)
  
  # combine at compound level
  if(metaLevel == "cpd"){
    
    if(mutual.feats){
      matched_res_ag <- aggregate(. ~ Matched.Compound, matched_res, paste, collapse=";")
      
      # keep compounds that only match across all files
      matched <- strsplit(matched_res_ag$.id, ";", fixed=TRUE)
      matched.inx <- vapply(matched, function(x) length(unique(x))==length(metaFiles), logical(1))
      
      if(sum(matched.inx) == 0){
        AddErrMsg("No compounds found across all files!")
        return(0)
      }
      
      matched_res_ag <- matched_res_ag[matched.inx,]
      # undo aggregate
      
      matched_res <- splitstackshape::cSplit(matched_res_ag, c(".id", "Query.Mass", "Matched.Form", "Mass.Diff", "Retention.Time", "matched_ts", 
                                                               "Matched.Pval", "effect.size", "lower.ci", "upper.ci"), 
                                             ";", "long", makeEqual = FALSE, type.convert = "as.character")
      matched_res <- data.table::setDF(matched_res)
      colnames(matched_res)[2] <- "File.Name"
      colnames(matched_res)[7:11] <- c("Matched.Scores", "Matched.Pvalue", "Matched.ES", "Matched.L.CI", "Matched.U.CI")
      
    }else{
      matched_res <- matched_res[, c(3,1,2,4:11)]
      matched_res <- data.table::setDF(matched_res)
      colnames(matched_res)[2] <- "File.Name"
      colnames(matched_res)[7:11] <- c("Matched.Scores", "Matched.Pvalue", "Matched.ES", "Matched.L.CI", "Matched.U.CI")
    }
    
    fast.write.csv(matched_res, file="mummichog_matched_compound_all.csv", row.names=FALSE);
    # or empirical compound combining here!  
  }else{
    
    all_ref_mz <- unlist(lapply(metaMsetObj, "[[", "ref_mzlist"))
    
    matched_res$RT.Rank <- rank(as.numeric(matched_res$Retention.Time), ties.method = "random")
    rt_tol_rank <- length(all_ref_mz) * mSetObj$dataSet$rt_frac
    
    # now create empirical compounds if necessary!
    # 1 compound matches to multiple m/z, filter by RT 
    if(mSetObj$paramSet$mumRT & version=="v2"){
      
      start <- Sys.time()
      empirical.cpd.list <- data.table:::split.data.table(matched_res, by="Matched.Compound", sorted = TRUE); # split all info by compound
      empirical.cpds2cpds <- vector(length=(length(empirical.cpd.list)), "list")
      names(empirical.cpds2cpds) <- names(empirical.cpd.list)
      
      # for each compound, if multiple matches, split into ECpds if > RT tolerance - rt_tol
      for(i in seq_along(empirical.cpd.list)){
        
        id <- empirical.cpd.list[[i]]$.id
        mzs <- empirical.cpd.list[[i]]$Query.Mass
        ions <- empirical.cpd.list[[i]]$Matched.Form
        rts <- as.numeric(empirical.cpd.list[[i]]$Retention.Time)
        rt.rank <- as.numeric(empirical.cpd.list[[i]]$RT.Rank)
        score <- empirical.cpd.list[[i]]$matched_ts
        mass.diff <- empirical.cpd.list[[i]]$Mass.Diff
        p.val <- empirical.cpd.list[[i]]$Matched.Pval
        es <- empirical.cpd.list[[i]]$effect.size
        l.ci <- empirical.cpd.list[[i]]$lower.ci
        u.ci <- empirical.cpd.list[[i]]$upper.ci
        cpds <- names(empirical.cpd.list)[i]
        
        # first, for each compound, determine ECs among matched ions
        if(length(mzs)>1){ # if multiple ECs per compound
          
          # first group together to create empirical cpds by rt
          names(rts) <- paste0(mzs, ";", ions, ";", rts, ";", cpds, ";", id, ";", score,
                               ";", mass.diff, ";", p.val,";", es, ";", l.ci, ";", u.ci)
          rts <- sort(rts)
          
          # second, group together to create empirical cpds by rt rank
          names(rt.rank) <- paste0(mzs, ";", ions, ";", rts, ";", cpds, ";", id, ";", score,
                                   ";", mass.diff, ";", p.val,";", es, ";", l.ci, ";", u.ci)
          rt.rank <- sort(rt.rank)
          
          split.inx <- c(0, cumsum(Reduce("&", list(abs(diff(rts)) > rt_tol, abs(diff(rt.rank)) > rt_tol_rank) )))
          
          # need to deal w. multiple rts but only 1 EC
          if(length(unique(split.inx)) > 1){
            e.cpds <- split(rts, split.inx)
            empirical.cpds2cpds[[i]] <- lapply(e.cpds, names)
          }else{
            empirical.cpds2cpds[[i]] <- paste0(names(rts), collapse="__")
          }
          
        }else{ # if only 1 EC per compound
          empirical.cpds2cpds[[i]] <- paste0(mzs, ";", ions, ";", rts, ";", cpds, ";", id, ";", score,
                                             ";", mass.diff, ";", p.val,";", es, ";", l.ci, ";", u.ci)
        }
      }
      
      initial_ecs <- unlist(empirical.cpds2cpds, recursive=FALSE)
      names(initial_ecs) <- paste0("EC", seq_along(initial_ecs))
      print(paste0(length(initial_ecs), " inital ECs created!"))
      
      # second, merge ECs if same m/z and form - append compounds
      try <- melt(initial_ecs)
      try2 <- strsplit(as.character(try[,1]), split="__", fixed=TRUE) # deals with multiple rts belonging to 1 EC
      try2 <- data.frame(value=unlist(try2), L1 = rep(try$L1, sapply(try2, length)))
      
      info <- strsplit(as.character(try2[,1]), split=";", fixed=TRUE)
      
      df_ecs <- data.frame(ec = as.character(try2[,2]), mz = sapply(info, `[[`, 1), form = sapply(info, `[[`, 2), rt = sapply(info, `[[`, 3), 
                           cpd = sapply(info, `[[`, 4), id = sapply(info, `[[`, 5), score = sapply(info, `[[`, 6), 
                           mass_diff = sapply(info, `[[`, 7), pval = sapply(info, `[[`, 8), es = sapply(info, `[[`, 9),
                           lci = sapply(info, `[[`, 10), uci = sapply(info, `[[`, 11), stringsAsFactors = F)
      
      df_ecs$str_row_inx <- paste(df_ecs$mz, df_ecs$form, df_ecs$rt, sep = "___")
      merged_ecs <- aggregate(. ~ str_row_inx, df_ecs, paste, collapse=";")
      
      # cleaning the df
      # merged_ecs$ec <- sapply(strsplit(merged_ecs$ec, ";"), function(x) unlist(x)[1]) - keep as long name
      merged_ecs$mz <- sapply(strsplit(merged_ecs$mz, ";", fixed=TRUE), function(x) unique(unlist(x)))
      merged_ecs$form <- sapply(strsplit(merged_ecs$form, ";", fixed=TRUE), function(x) unique(unlist(x)))
      merged_ecs$rt <- sapply(strsplit(merged_ecs$rt, ";", fixed=TRUE), function(x) unique(unlist(x)))
      merged_ecs$id <- sapply(strsplit(merged_ecs$id, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";"))
      merged_ecs$score <- sapply(strsplit(merged_ecs$score, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";") )
      merged_ecs$mass_diff <- sapply(strsplit(merged_ecs$mass_diff, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";") )
      merged_ecs$pval <- sapply(strsplit(merged_ecs$pval, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";") )
      merged_ecs$es <- sapply(strsplit(merged_ecs$es, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";") )
      merged_ecs$lci <- sapply(strsplit(merged_ecs$lci, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";") )
      merged_ecs$uci <- sapply(strsplit(merged_ecs$uci, ";", fixed=TRUE), function(x) paste0(unique(unlist(x)), collapse=";") )
      print(paste0(length(unique(merged_ecs$ec)), " merged ECs identified!"))
      
      # third, check if primary ion is present
      # needs to be per EC!
      if(mSetObj$dataSet$primary_ion=="yes"){
        
        ecs <- unique(merged_ecs$ec)
        
        # function to group ECs and verify if contains primary ion
        new_info <- lapply(ecs, function(x) { 
          new_ec_info <- merged_ecs[which(merged_ecs$ec == x),] # subset merged_ecs to rows containing ECx
          primary.inx <- length(intersect(new_ec_info$form, primary_ions))
          
          if(primary.inx>0){
            new_ec_info <- new_ec_info
          }else{
            new_ec_info <- NULL
          }
          new_ec_info
        })  
        
        final_ecs <- do.call(args=new_info, what=rbind)[,-1]
        
      }else{
        final_ecs <- merged_ecs[,-1]
      }
      
      colnames(final_ecs) <- c("Empirical.Compound", "Query.Mass", "Matched.Form", "Retention.Time", "Matched.Compound", "FileName", "Matched.Scores",
                               "Mass.Diff", "Matched.Pvalue", "Matched.ES", "Matched.L.CI", "Matched.U.CI")
      
      # transform to long format
      cpd_split <- strsplit(as.character(final_ecs$Matched.Compound), ";", fixed=TRUE)
      reps <- pmax(lengths(cpd_split))
      df2 <- final_ecs[rep(1:nrow(final_ecs), reps), c(1:4, 6:12)]
      df2$Matched.Compound <- unlist(mapply(function(x,y) c(x, rep(NA, y)), cpd_split, reps-lengths(cpd_split)))
      df2 <- unique(df2)
      
      # now deal with the fact that if at least one EC overlap, need to count as same EC per compound...
      my_final_cpds <- aggregate(. ~ Matched.Compound, df2, paste, collapse="___")
      my_final_cpds_list <- lapply(split(my_final_cpds$Empirical.Compound, my_final_cpds$Matched.Compound), unlist)
      
      cpd2ec1 <- lapply(seq_along(my_final_cpds_list), function(x) { # function used to make grouping of ecs per cpd
        
        ecs <- unlist(strsplit(my_final_cpds_list[[x]], "___", fixed=TRUE))
        
        if(length(ecs) > 1){
          ecs.list <- as.list(strsplit(ecs, ";", fixed=TRUE))
          library(igraph)
          m = sapply(ecs.list, function(x) sapply(ecs.list, function(y) length(intersect(x,y))>0))
          g = igraph::groups(components(graph_from_adjacency_matrix(m)))
          ecs <- paste0(sapply(g, function(z) paste0(ecs[z], collapse = "|") ), collapse = "___")
        }
        ecs
      })
      
      names(cpd2ec1) <- names(my_final_cpds_list)
      
      update_ecs <- lapply(seq_along(cpd2ec1), function(z) {
        
        ecs.old <- unlist(strsplit(my_final_cpds_list[[z]], "___", fixed=TRUE))
        ecs.new <- unlist(strsplit(cpd2ec1[[z]], "___", fixed=TRUE))
        
        for(i in seq_along(ecs.new)){
          pattern <- ecs.new[i]
          pattern_vec <- unlist(strsplit(pattern, "\\|", fixed=TRUE))
          up.pattern <- paste0(unique(pattern_vec), collapse = "|")
          ecs.old[ ecs.old %in% pattern_vec  ] <- up.pattern
        }
        
        ecs.old <- paste0(ecs.old, collapse = "___")
        ecs.old
      })
      
      updated_ecs <- do.call(rbind, update_ecs)
      my_final_cpds$Empirical.Compound <- updated_ecs
      
      new_dt <- data.table::data.table(my_final_cpds)
      new_dt <- new_dt[, list(Query.Mass = unlist(strsplit(as.character(Query.Mass), "___", fixed=TRUE)), 
                              Matched.Form = unlist(strsplit(as.character(Matched.Form), "___", fixed=TRUE)),
                              Retention.Time = unlist(strsplit(as.character(Retention.Time), "___", fixed=TRUE)),
                              Empirical.Compound = unlist(strsplit(as.character(Empirical.Compound), "___", fixed=TRUE)),
                              File.Name = unlist(strsplit(as.character(FileName), "___", fixed=TRUE)),
                              Matched.Scores = unlist(strsplit(as.character(Matched.Scores), "___", fixed=TRUE)),
                              Mass.Diff = unlist(strsplit(as.character(Mass.Diff), "___", fixed=TRUE)),
                              Matched.Pvalue = unlist(strsplit(as.character(Matched.Pvalue), "___", fixed=TRUE)),
                              Matched.ES = unlist(strsplit(as.character(Matched.ES), "___", fixed=TRUE)),
                              Matched.L.CI = unlist(strsplit(as.character(Matched.L.CI), "___", fixed=TRUE)),
                              Matched.U.CI = unlist(strsplit(as.character(Matched.U.CI), "___", fixed=TRUE))
      ),
      by = Matched.Compound]
      
      matched_res <- data.frame(Query.Mass = new_dt$Query.Mass, Matched.Compound = new_dt$Matched.Compound, Matched.Form = new_dt$Matched.Form, Mass.Diff = new_dt$Mass.Diff,
                                Retention.Time = new_dt$Retention.Time, Matched.Scores = new_dt$Matched.Scores, Matched.Pvalue = new_dt$Matched.Pvalue, 
                                Matched.ES = new_dt$Matched.ES, Matched.L.CI = new_dt$Matched.L.CI, Matched.U.CI = new_dt$Matched.U.CI,
                                Empirical.Compound = new_dt$Empirical.Compound, File.Name = new_dt$File.Name, stringsAsFactors = FALSE)
      
      # make new EC names
      ec <- as.list(matched_res$Empirical.Compound)
      ec.unique <- unique(matched_res$Empirical.Compound)
      ec.new <- paste0("EC000", seq_along(ec.unique))
      
      ec.new <- vapply(seq_along(ec), function(i) { 
        
        inx <- match(ec[[i]], ec.unique)
        ec <- ec.new[inx]
        ec
        
      }, character(1))
      
      matched_res$Empirical.Compound <- gsub("\\|.*", "", ec.new)
      end <- Sys.time()
      totaltime <- end-start
      print(paste0(length(unique(matched_res$Empirical.Compound)), " empirical compounds identified in ", totaltime, " seconds."))
      
      if(mutual.feats){
        # keep empirical compounds that only match across all files
        matched_res <- aggregate(. ~ Empirical.Compound, matched_res, paste, collapse="___")
        matched <- strsplit(matched_res$File.Name, ";|___")
        matched.inx <- vapply(matched, function(x) length(unique(x))==length(metaFiles), logical(1))
        
        if(sum(matched.inx)==0){
          AddErrMsg("No empirical compounds found across all studies!")
          return(0)
        }else if(sum(matched.inx) < 50){
          AddErrMsg("Not enough empirical compounds matched across all studies! Try meta-analysis at a higher level (compound or pathway).")
          return(0)
        }
        
        print(paste0(sum(matched.inx), "matched empirical compounds identified across all studies!"))
        
        matched_res <- matched_res[matched.inx,]
        matched_res <- splitstackshape::cSplit(matched_res, c("Query.Mass", "Matched.Compound", "Matched.Form", "Mass.Diff", "Retention.Time", "Matched.Scores", 
                                                              "Matched.Pvalue", "Matched.ES", "Matched.L.CI", "Matched.U.CI", "File.Name"), 
                                               "___", "long", makeEqual = FALSE, type.convert = "as.character")
        matched_res <- data.table::setDF(matched_res)
      }else{
        matched_res <- matched_res[,c(11,1:10,12)]
        matched_res <- data.table::setDF(matched_res)
      }
      fast.write.csv(matched_res, file="mummichog_matched_compound_postmerge.csv", row.names=FALSE);
    }else{
      AddErrMsg("Meta-analysis at empirical compound level is invalid!")
      return(0)
    }
  }
  
  ref_mzlist <- lapply(metaMsetObj, "[[", "ref_mzlist") 
  ref_mzlist <- unlist(unique(ref_mzlist))
  mSetObj$dataSet$ref_mzlist <- ref_mzlist
  
  mumRT <- unlist(lapply(metaMsetObj, "[[", "mumRT")) 
  qs::qsave(matched_res, "mum_res.qs");
  
  ##############################################
  # COMBINE EITHER P-VALUES
  # EFFECT-SIZES
  # OR BOTH (only for mummichog or integ_peaks)
  # then re-make input_cpdlist and input_ecpdlist
  # gsea uses rank metric only
  
  #  if(anal.type %in% c("mummichog", "integ_peaks")){
  
  if(combine.level %in% c("pvalue", "both")){
    
    if(metaLevel %in% "ec"){
      # merge to empirical compounds
      all_p <- aggregate(. ~ Empirical.Compound, matched_res, paste, collapse="___")
      old_p <- strsplit(all_p$Matched.Pvalue, "___", fixed=TRUE)
      scores <- strsplit(all_p$Matched.Pvalue, ";|___")
      scores <- lapply(scores, function(x) {
        if(length(x) == 1){
          x <- rep(x, 2)
        }
        x;}) 
      
    }else{
      # merge to compounds
      all_p <- aggregate(. ~ Matched.Compound, matched_res, paste, collapse="___")
      old_p <- strsplit(all_p$Matched.Pvalue, "___", fixed=TRUE)
      scores <- strsplit(all_p$Matched.Pvalue, ";|___")
      scores <- lapply(scores, function(x) {
        if(length(x) == 1){
          x <- rep(x, 2)
        }
        x;}) 
    }
    
    # combine p-values
    if(pval.method=="fisher"){
      meta.pvals <- lapply(scores, function(x) sumlog(as.numeric(x)))
    }else if(pval.method=="edgington"){ 
      meta.pvals <- lapply(scores, function(x) sump(as.numeric(x)))
    }else if(pval.method=="stouffer"){
      meta.pvals <- lapply(scores, function(x) sumz(x))
    }else if(pval.method=="vote"){
      meta.pvals <- lapply(scores, function(x) votep(x))
    }else if(pval.method=="min"){
      Meta.P <- lapply(scores, function(x) min(x) )
    }else if(pval.method=="max") {
      Meta.P <- lapply(scores, function(x) max(x) )
    }else{
      AddErrMsg("Invalid meta-analysis method!")
      return(0)
    }
    
    #extract p-values
    if(exists("meta.pvals")){
      Meta.P <- unlist(lapply(meta.pvals, function(z) z["p"]))
    }
    
    Meta.P2 <- rep(Meta.P, vapply(old_p, length, numeric(1)))
    matched_res$Meta.P <- Meta.P2
    #    }
    
    # now create input mzlist - used to create
    # input cpd/ec list
    cutoff <- mSetObj$dataSet$cutoff
    
    if(combine.level == "both"){
      my.inx <- matched_res[,"Meta.P.Both"] < cutoff
    }else if(combine.level == "pvalue"){
      my.inx <- matched_res[,"Meta.P"] < cutoff
    }else{
      my.inx <- matched_res[,"Meta.ES.Pval"] < cutoff
    }
    
    input_mzlist <- unlist(unique(matched_res[as.vector(my.inx), "Query.Mass"]))
    
  }else{ # gsea
    input_mzlist <- lapply(metaMsetObj, "[[", "input_mzlist") 
    input_mzlist <- unlist(unique(input_mzlist)) # will be updated later
  }
  
  sig.size <- length(input_mzlist);
  mSetObj$dataSet$N <- sig.size;
  mSetObj$dataSet$input_mzlist <- input_mzlist
  
  if(all(mumRT) & version=="v2"){ # RT need to be in EC space
    
    # for GSEA
    # need to merge t-scores if same m/z in the data
    if(rank.metric == "mean"){     # default using the mean
      matched_res$Matched.Scores <- vapply(matched_res$Matched.Scores, function(x) mean(as.numeric(unlist(strsplit(as.character(x), ";", fixed=TRUE)))), numeric(1))
    }else if(rank.metric == "min"){
      matched_res$Matched.Scores <- vapply(matched_res$Matched.Scores, function(x) min(as.numeric(unlist(strsplit(as.character(x), ";", fixed=TRUE)))), numeric(1))
    }else if(rank.metric == "max"){
      matched_res$Matched.Scores <- vapply(matched_res$Matched.Scores, function(x) max(as.numeric(unlist(strsplit(as.character(x), ";", fixed=TRUE)))), numeric(1))
    }else if(rank.metric == "median"){
      matched_res$Matched.Scores <- vapply(matched_res$Matched.Scores, function(x) median(as.numeric(unlist(strsplit(as.character(x), ";", fixed=TRUE)))), numeric(1))
    }else{
      AddErrMsg("Invalid method selected for merging scores!")
      return(0)
    }
    
    ec.exp.mat <- data.frame(key=matched_res$Empirical.Compound, value=as.numeric(matched_res$Matched.Scores), stringsAsFactors = F)
    ec_exp_dict <- Convert2Dictionary(ec.exp.mat);
    ec.exp.vec <- unlist(lapply(ec_exp_dict, max));
    
    # also need to make cpd_exp_dict for KEGG network view
    exp.mat <- data.frame(key=matched_res$Matched.Compound, value=as.numeric(matched_res$Matched.Scores), stringsAsFactors = F);
    cpd_exp_dict <- Convert2Dictionary(exp.mat);
    
    # ecpd to cpd dict
    cpd_ecpd_dict <- Convert2Dictionary(matched_res[,c(3,1)])
    ecpd_cpd_dict <- Convert2Dictionary(matched_res[,c(1,3)])
    
    # now mz 2 ecpd dict
    mz2cpd_dict <- Convert2Dictionary(matched_res[,c(2,3)]); #indexed/named by mz
    mz2ec_dict <- Convert2Dictionary(matched_res[,c(2,1)])
    ec2mz_dict <- Convert2Dictionary(matched_res[,c(1,2)])
    
    # save to mSetObj
    mSetObj$ec_exp_dict <- ec_exp_dict
    mSetObj$cpd_exp_dict <- cpd_exp_dict;
    mSetObj$ec_exp <- ec.exp.vec
    mSetObj$mz2cpd_dict <- mz2cpd_dict;
    mSetObj$mz2ec_dict <- mz2ec_dict
    mSetObj$ec2mz_dict <- ec2mz_dict
    mSetObj$ecpd_cpd_dict <- ecpd_cpd_dict
    mSetObj$cpd_ecpd_dict <- cpd_ecpd_dict
    mSetObj$cpd_ecpd_counts <- cpd2ec1
    
    # now do matching to identify significant input_ecpdlist
    # trio.list <- data.frame(mz = names(mz2ec_dict), ec = sapply(mz2ec_dict, paste, collapse="; "), cpd = sapply(mz2cpd_dict, paste, collapse="; "))
    refmz <- names(mz2ec_dict)
    hits.index <- which(refmz %in% as.character(input_mzlist));
    ec1 <- unique(unlist(mz2ec_dict[hits.index]));
    mSetObj$input_ecpdlist <- ec1;
    mSetObj$total_matched_ecpds <- unique(as.vector(matched_res$Empirical.Compound));
    form.mat <- cbind(matched_res[,2], matched_res[,4]);
    
  }else{ # compound level
    
    # get the expression profile for each 
    exp.mat <- data.frame(key=matched_res$Matched.Compound, value=as.numeric(matched_res$Matched.Scores), stringsAsFactors = F);
    cpd_exp_dict <- Convert2Dictionary(exp.mat);
    # create average exp
    exp.vec <- unlist(lapply(cpd_exp_dict, function(x) mean(unlist(x))));
    
    # now need to get the mapping from mz to compound id (one mz can have 0, 1, or more id hits)
    mz2cpd_dict <- Convert2Dictionary(matched_res[ , c("Query.Mass", "Matched.Compound")]) #indexed/named by mz
    cpd2mz_dict <- Convert2Dictionary(matched_res[ , c("Matched.Compound", "Query.Mass")]) # indexed/named by id
    
    # now do matching to identify significant input_cpdlist
    refmz <- names(mz2cpd_dict)
    hits.index <- which(refmz %in% as.character(input_mzlist));
    cpd1 <- unique(unlist(mz2cpd_dict[hits.index]));
    cpd1 <- cpd1[!(cpd1 %in% currency)];
    form.mat <- cbind(matched_res$Query.Mass, matched_res$Matched.Form);
    
    mSetObj$mz2cpd_dict <- mz2cpd_dict;
    mSetObj$cpd_exp_dict <- cpd_exp_dict;
    mSetObj$cpd_exp <- exp.vec;
    mSetObj$cpd2mz_dict <- cpd2mz_dict;
    mSetObj$input_cpdlist <- cpd1;
    mSetObj$dataSet$N <- length(input_mzlist);
    mSetObj$total_matched_cpds <- unique(as.vector(matched_res$Matched.Compound));
  }
  
  cpd_form_dict <- Convert2Dictionary(form.mat);
  mSetObj$cpd_form_dict <- cpd_form_dict;
  return(mSetObj);
}

#' @param method If "cohen", computes Cohen's d, if "hedges",
#' computes Hegdes' g effect size statistics.
#' @noRd
CalculateEffectSize <- function(mSetObj=NA, paired=FALSE, method="cohen"){
  
  mSetObj <- .get.mSet(mSetObj);
  
  inx1 <- which(mSetObj$dataSet$cls==levels(mSetObj$dataSet$cls)[1])
  inx2 <- which(mSetObj$dataSet$cls==levels(mSetObj$dataSet$cls)[2])
  
  # samples in row, features in columns
  x <- mSetObj$dataSet$norm[inx1,]
  y <- mSetObj$dataSet$norm[inx2,]
  
  library(effsize)
  
  my.fun <- function(x) {
    
    if(method == "cohen"){
      tmp <- try(cohen.d(x[inx1], x[inx2], paired=paired, hedges.correction=FALSE));
    }else{
      tmp <- try(cohen.d(x[inx1], x[inx2], paired=paired, hedges.correction=TRUE));
    }
    
    if(class(tmp) == "try-error") {
      return(c(NA, NA, NA, NA));
    }else{
      return(c(tmp$estimate, tmp$sd, tmp$conf.int));
    }
  }
  
  results <- apply(as.matrix(mSetObj$dataSet$norm), 2, my.fun)
  rownames(results) <- c("effect_size", "win_group_stdev", "lower_ci", "upper_ci")
  
  mSetObj$analSet$effect.size <- t(results);
  return(.set.mSet(mSetObj));
}

# Internal function for permutation, no RT
.perform.mummichogPermutations <- function(mSetObj, permNum){
  
  num_perm <- permNum;
  print(paste('Resampling, ', num_perm, 'permutations to estimate background ...'));
  permutation_hits <- permutation_record <- vector("list", num_perm);
  matched_res <- qs::qread("mum_res.qs");
  set.seed(123)
  for(i in 1:num_perm){ # for each permutation, create list of input compounds and calculate pvalues for each pathway
    input_mzlist <- sample(mSetObj$dataSet$ref_mzlist, mSetObj$dataSet$N)
    t <- make_cpdlist(mSetObj, input_mzlist);
    perm <- ComputeMummichogPermPvals(t, mSetObj$total_matched_cpds, mSetObj$pathways, matched_res, input_mzlist, mSetObj$cpd2mz_dict);
    permutation_record[[i]] <- perm[1]
    permutation_hits[[i]] <- perm[2]
  }
  
  mSetObj$perm_record <- permutation_record;
  mSetObj$perm_hits <- permutation_hits;
  
  return(mSetObj);
}

# Calculate p-values for each Lperm
# Used in higher mummichogR functions w.out RT
ComputeMummichogPermPvals <- function(input_cpdlist, total_matched_cpds, pathways, matches.res, input_mzlist, cpd2mz_dict){
  
  ora.vec <- input_cpdlist; #Lperm
  query_set_size <- length(ora.vec)
  current.mset <- pathways$cpds; #all
  total_cpds <- unique(total_matched_cpds) #matched compounds
  total_feature_num <- length(total_cpds)
  
  size <- negneg <- vector(mode="list", length=length(current.mset));
  
  cpds <- lapply(current.mset, function(x) intersect(x, total_cpds)); # pathways & all ref cpds
  feats <- lapply(current.mset, function(x) intersect(x, ora.vec)); #pathways & lsig
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  set.num <- unlist(lapply(cpds, length)); #cpdnum
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset));
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- min(feat_len[i], count_cpd2mz(cpd2mz_dict, unlist(feats[i]), input_mzlist))
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size;
  }
  
  unsize <- as.integer(unlist(sizes))
  res.mat <- matrix(0, nrow=length(current.mset), ncol=1)
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg)), query_set_size)
  res.mat[,1] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  perm_records <- list(res.mat, as.matrix(unsize));
  return(perm_records);
}

# Internal function for permutation
.perform.mummichogRTPermutations <- function(mSetObj, permNum){
  num_perm <- permNum;
  print(paste('Resampling, ', num_perm, 'permutations to estimate background ...'));
  permutation_hits <- permutation_record <- vector("list", num_perm);
  matched_res <- qs::qread("mum_res.qs");
  set.seed(123)
  for(i in 1:num_perm){ # for each permutation, create list of input emp compounds and calculate pvalues for each pathway
    input_mzlist <- sample(mSetObj$dataSet$ref_mzlist, mSetObj$dataSet$N)
    t <- make_ecpdlist(mSetObj, input_mzlist);
    perm <- ComputeMummichogRTPermPvals(t, mSetObj$total_matched_ecpds, mSetObj$pathways, matched_res, input_mzlist);
    permutation_record[[i]] <- perm[1]
    permutation_hits[[i]] <- perm[2]
  }
  
  mSetObj$perm_record <- permutation_record;
  mSetObj$perm_hits <- permutation_hits
  
  return(mSetObj);
}

# Calculate p-values for each Lperm
# Used in higher mummichogR functions w. RT
ComputeMummichogRTPermPvals <- function(input_ecpdlist, total_matched_ecpds, pathways, matches.res, input_mzlist){
  
  ora.vec <- input_ecpdlist; #Lperm
  query_set_size <- length(ora.vec) # query set size
  current.mset <- pathways$emp_cpds; #all
  total_ecpds <- unique(total_matched_ecpds) # matched number of empirical compounds
  total_feature_num <- length(total_ecpds)
  
  size <- negneg <- vector(mode="list", length=length(current.mset));
  
  ecpds <- lapply(current.mset, function(x) intersect(x, total_ecpds)); # pathways & all ref ecpds
  feats <- lapply(current.mset, function(x) intersect(x, ora.vec)); #pathways & query ecpds (perm lsig)
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  set.num <- unlist(lapply(ecpds, length)); #cpdnum
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset));
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- feat_len[i] # for ecs, just use length of overlap feats - overlap_size
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size;
  }
  
  unsize <- as.integer(unlist(sizes))
  res.mat <- matrix(0, nrow=length(current.mset), ncol=1)
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg)), query_set_size)
  res.mat[,1] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  perm_records <- list(res.mat, as.matrix(unsize));
  return(perm_records);
}

# Internal function for significant p value 
.compute.mummichogSigPvals <- function(mSetObj){
  
  qset <- unique(unlist(mSetObj$input_cpdlist)); #Lsig ora.vec
  query_set_size <- length(qset); #q.size
  
  total_cpds <- unique(mSetObj$total_matched_cpds) #all matched compounds
  total_feature_num <- length(total_cpds)
  
  current.mset <- mSetObj$pathways$cpds; #all compounds per pathway
  path.num <- unlist(lapply(current.mset, length));
  
  cpds <- lapply(current.mset, function(x) intersect(x, total_cpds)); #pathways & all ref cpds
  set.num <- unlist(lapply(cpds, length)); #cpdnum
  
  feats <- lapply(current.mset, function(x) intersect(x, qset)); #pathways & lsig
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  feat_vec <- sapply(feats, function(x) paste(x, collapse=";"))
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset)); #empty lists
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- min(feat_len[i], count_cpd2mz(mSetObj$cpd2mz_dict, unlist(feats[i]), mSetObj$dataSet$input_mzlist)) #min overlap or mz hits
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size; # failure in left part
  }
  
  #error fixing for negatives, problem occurs when total_feat_num and query_set_size too close (lsig too close to lall)
  negneg <- rapply(negneg, function(x) ifelse(x<0,0,x), how = "replace") 
  
  unsize <- as.integer(unlist(sizes));
  
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));
  
  # prepare for the result table
  res.mat <- matrix(0, nrow=length(current.mset), ncol=6);
  
  #fishermatrix for phyper
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg) - unsize), query_set_size); 
  first <- unlist(lapply(sizes, function(x) max(0, x-1)));
  easematrix <- cbind(first, (set.num - unsize + 1), (query_set_size - unsize), unlist(negneg)); 
  
  res.mat[,1] <- path.num;  
  res.mat[,2] <- set.num;
  res.mat[,3] <- unsize;
  res.mat[,4] <- query_set_size*(path.num/uniq.count); #expected
  res.mat[,6] <- apply(easematrix, 1, function(x) fisher.test(matrix(x, nrow=2), alternative = "greater")$p.value);
  res.mat[,5] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  colnames(res.mat) <- c("Pathway total", "Hits.total", "Hits.sig", "Expected", "FET", "EASE");
  rownames(res.mat) <- mSetObj$pathways$name
  
  mSetObj$pvals <- res.mat;
  permutations_hits <- matrix(unlist(mSetObj$perm_hits), nrow=length(mSetObj$perm_hits), byrow=TRUE);
  sig_hits <- res.mat[,3]; # sighits
  sigpvalue <- res.mat[,6]; # EASE scores
  
  perm_record <- unlist(mSetObj$perm_record);
  perm_minus <- abs(0.9999999999 - perm_record);
  
  if(length(sig_hits[sig_hits!=0]) < round(length(sig_hits)*0.05)){ # too few hits that can't calculate gamma dist!
    if(!exists("adjustedp")){
      adjustedp <- rep(NA, length = length(res.mat[,1]))
    }
    res.mat <- cbind(res.mat, Gamma=adjustedp);
  }else{
    tryCatch({
      fit.gamma <- fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1));
      rawpval <- as.numeric(sigpvalue);
      adjustedp <- 1 - (pgamma(1-rawpval, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["scale"]));
    }, error = function(e){
      if(mSetObj$dataSet$mumType == "table"){
        if(!exists("adjustedp")){
          adjustedp <- rep(NA, length = length(res.mat[,1]))
        }
        res.mat <- cbind(res.mat, Gamma=adjustedp);
      }
      print(e)   
    }, finally = {
      if(!exists("adjustedp")){
        adjustedp <- rep(NA, length = length(res.mat[,1]))
      }
      res.mat <- cbind(res.mat, Gamma=adjustedp);
    })
  }
  
  #calculate empirical p-values
  record <- mSetObj$perm_record
  fisher.p <- as.numeric(res.mat[,5])
  
  #pathway in rows, perms in columns
  record_matrix <- do.call(cbind, do.call(cbind, record))
  num_perm <- ncol(record_matrix)
  
  #number of better hits for web
  better.hits <- sapply(seq_along(record_matrix[,1]), function(i) sum(record_matrix[i,] <= fisher.p[i])  )
  
  #account for a bias due to finite sampling - Davison and Hinkley (1997)
  emp.p <- sapply(seq_along(record_matrix[,1]), function(i) (sum(record_matrix[i,] <= fisher.p[i])/num_perm) )
  
  res.mat <- cbind(res.mat, Emp.Hits=better.hits, Empirical=emp.p, Cpd.Hits = feat_vec)
  
  # remove those no hits
  hit.inx <- as.numeric(as.character(res.mat[,3])) > 0;
  res.mat <- res.mat[hit.inx, , drop=FALSE];
  
  if(nrow(res.mat) <= 1){
    AddErrMsg("Not enough m/z to compound hits for pathway analysis!")
    return(0)
  }
  
  # prepare json element for network
  hits.all <- cpds[hit.inx];
  hits.sig <- feats[hit.inx];  
  path.nms <- mSetObj$pathways$name[hit.inx];
  
  # order by p-values
  ord.inx <- order(res.mat[,7]);
  
  Cpd.Hits <- res.mat[ord.inx, 10]
  res.mat <- signif(apply(as.matrix(res.mat[ord.inx, 1:9, drop=FALSE]), 2, as.numeric), 5);
  rownames(res.mat) <- path.nms[ord.inx]
  mSetObj$mummi.resmat <- res.mat[,-9];
  
  mSetObj$path.nms <- path.nms[ord.inx]
  mSetObj$path.hits <- convert2JsonList(hits.all[ord.inx])
  mSetObj$path.pval <- as.numeric(res.mat[,5])
  matched_res <- qs::qread("mum_res.qs");
  
  json.res <- list(
    cmpd.exp = mSetObj$cpd_exp,
    path.nms = path.nms[ord.inx],
    hits.all = convert2JsonList(hits.all[ord.inx]),
    hits.all.size = as.numeric(res.mat[,2]),
    hits.sig = convert2JsonList(hits.sig[ord.inx]),
    hits.sig.size = as.numeric(res.mat[,3]),
    fisher.p = as.numeric(res.mat[,7]),
    peakToMet = mSetObj$cpd_form_dict,
    peakTable = matched_res
  );
  
  matri = res.mat[,-8]
  matri = cbind(res.mat, paste0("P", seq.int(1, nrow(res.mat))))
  colnames(matri)[ncol(matri)] = "Pathway Number"
  matri <- cbind(matri, Cpd.Hits)
  fast.write.csv(matri, file=mSetObj$mum_nm_csv, row.names=TRUE);
  json.mat <- RJSONIO::toJSON(json.res, .na='null');
  sink(mSetObj$mum_nm);
  cat(json.mat);
  sink();
  return(mSetObj);
}

# Internal function for significant p value with RT
.compute.mummichogRTSigPvals <- function(mSetObj){
  
  qset <- unique(unlist(mSetObj$input_ecpdlist)); #Lsig ora.vec
  query_set_size <- length(qset); #q.size
  
  total_ecpds <- unique(mSetObj$total_matched_ecpds) #all matched compounds
  total_feature_num <- length(total_ecpds)
  
  current.mset <- mSetObj$pathways$emp_cpds; #all compounds per pathway
  path.num <- unlist(lapply(current.mset, length));
  
  ecpds <- lapply(current.mset, function(x) intersect(x, total_ecpds)); #pathways & all ref ecpds
  set.num <- unlist(lapply(ecpds, length)); # total ecpd num in pathway
  
  feats <- lapply(current.mset, function(x) intersect(x, qset)); #pathways & lsig
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  feat_vec <- sapply(feats, function(x) paste(x, collapse=";"))
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset)); #empty lists
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- feat_len[i] # overlap size
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size; # failure in left part
  }
  
  #error fixing for negatives, problem occurs when total_feat_num and query_set_size too close (lsig too close to lall)
  negneg <- rapply(negneg, function(x) ifelse(x<0,0,x), how = "replace") 
  
  unsize <- as.integer(unlist(sizes));
  
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));
  
  # prepare for the result table
  res.mat <- matrix(0, nrow=length(current.mset), ncol=6);
  
  #fishermatrix for phyper
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg) - unsize), query_set_size); 
  first <- unlist(lapply(sizes, function(x) max(0, x-1)));
  easematrix <- cbind(first, (set.num - unsize + 1), (query_set_size - unsize), unlist(negneg)); 
  
  res.mat[,1] <- path.num;  
  res.mat[,2] <- set.num;
  res.mat[,3] <- unsize;
  res.mat[,4] <- query_set_size*(path.num/uniq.count); #expected
  res.mat[,5] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  res.mat[,6] <- apply(easematrix, 1, function(x) fisher.test(matrix(x, nrow=2), alternative = "greater")$p.value);
  
  colnames(res.mat) <- c("Pathway total", "Hits.total", "Hits.sig", "Expected", "FET", "EASE");
  rownames(res.mat) <- mSetObj$pathways$name
  
  mSetObj$pvals <- res.mat;
  permutations_hits <- matrix(unlist(mSetObj$perm_hits), nrow=length(mSetObj$perm_hits), byrow=TRUE);
  sig_hits <- res.mat[,3]; # sighits
  sigpvalue <- res.mat[,5]; # EASE scores
  
  perm_record <- unlist(mSetObj$perm_record);
  perm_minus <- abs(0.9999999999 - perm_record);
  
  if(length(sig_hits[sig_hits!=0]) < round(length(sig_hits)*0.05)){ # too few hits that can't calculate gamma dist!
    if(!exists("adjustedp")){
      adjustedp <- rep(NA, length = length(res.mat[,1]))
    }
    res.mat <- cbind(res.mat, Gamma=adjustedp);
  }else{
    tryCatch({
      fit.gamma <- fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1));
      rawpval <- as.numeric(sigpvalue);
      adjustedp <- 1 - (pgamma(1-rawpval, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["scale"]));
    }, error = function(e){
      if(mSetObj$dataSet$mumType == "table"){
        if(!exists("adjustedp")){
          adjustedp <- rep(NA, length = length(res.mat[,1]))
        }
        res.mat <- cbind(res.mat, Gamma=adjustedp);
      }
      print(e)   
    }, finally = {
      if(!exists("adjustedp")){
        adjustedp <- rep(NA, length = length(res.mat[,1]))
      }
      res.mat <- cbind(res.mat, Gamma=adjustedp);
    })
  }
  
  #calculate empirical p-values
  record <- mSetObj$perm_record
  fisher.p <- as.numeric(res.mat[,5])
  
  #pathway in rows, perms in columns
  record_matrix <- do.call(cbind, do.call(cbind, record))
  num_perm <- ncol(record_matrix)
  
  #number of better hits for web
  better.hits <- sapply(seq_along(record_matrix[,1]), function(i) sum(record_matrix[i,] <= fisher.p[i])  )
  
  #account for a bias due to finite sampling - Davison and Hinkley (1997)
  emp.p <- sapply(seq_along(record_matrix[,1]), function(i) (sum(record_matrix[i,] <= fisher.p[i])/num_perm) )
  
  res.mat <- cbind(res.mat, Emp.Hits = better.hits, Empirical = emp.p, EC.Hits = feat_vec)
  
  # remove pathways with no hits
  hit.inx <- as.numeric(as.character(res.mat[,3])) > 0;
  res.mat <- res.mat[hit.inx, , drop=FALSE];
  
  if(nrow(res.mat) <= 1){
    AddErrMsg("Not enough m/z to compound hits for pathway analysis! Try Version 1 (no RT considerations)!")
    return(0)
  }
  
  # prepare json element for network
  # need to convert ecpds to cpds
  # and get average expression based on ec
  cpds <- lapply(ecpds, function(x) unique(unlist(mSetObj$ecpd_cpd_dict[match(x, names(mSetObj$ecpd_cpd_dict))])) )  
  cpd.exp.vec <- sapply(ecpds, function(x) mean(mSetObj$ec_exp[match(x, names(mSetObj$ec_exp))]) )
  cpds_feats <- lapply(feats, function(x) unique(unlist(mSetObj$ecpd_cpd_dict[match(x, names(mSetObj$ecpd_cpd_dict))])) )  
  
  # now make exp vec for all compounds
  cpds2ec <- mSetObj$cpd_ecpd_dict
  cpds.all <- unique(unlist(mSetObj$ecpd_cpd_dict[match(total_ecpds, names(mSetObj$ecpd_cpd_dict))]))
  cpd.exp.vec <- sapply(cpds.all, function(x) sapply(seq_along(x), function(i) mean(mSetObj$ec_exp[match(unique(unlist(cpds2ec[match(x[[i]], names(cpds2ec))])), names(mSetObj$ec_exp))]) ) )
  
  hits.all <- cpds[hit.inx];
  hits.sig <- cpds_feats[hit.inx];  
  path.nms <- mSetObj$pathways$name[hit.inx];
  
  # order by p-values
  if(length(na.omit(res.mat[,7])) == 0){
    ord.inx <- order(res.mat[,5]); # order by FET if gamma not able to be calc
  }else{
    ord.inx <- order(res.mat[,7]); # order by gamma
  }
  
  EC.Hits = res.mat[ord.inx, 10]
  res.mat <- signif(apply(as.matrix(res.mat[ord.inx, 1:9]), 2, as.numeric), 5); # loop through columns and keep rownames
  rownames(res.mat) <- path.nms[ord.inx]
  
  mSetObj$mummi.resmat <- res.mat[,-9];
  mSetObj$path.nms <- path.nms[ord.inx]
  mSetObj$path.hits <- convert2JsonList(hits.all[ord.inx])
  mSetObj$path.pval <- as.numeric(res.mat[,5])
  matched_res <- qs::qread("mum_res.qs");
  
  json.res <- list(
    cmpd.exp = cpd.exp.vec,
    path.nms = path.nms[ord.inx],
    hits.all = convert2JsonList(hits.all[ord.inx]),
    hits.all.size = as.numeric(res.mat[,2]),
    hits.sig = convert2JsonList(hits.sig[ord.inx]),
    hits.sig.size = as.numeric(res.mat[,3]),
    fisher.p = as.numeric(res.mat[,7]),
    peakToMet = mSetObj$cpd_form_dict,
    peakTable = matched_res
  );
  
  matri = res.mat[,-8]
  matri = cbind(res.mat, paste0("P", seq.int(1, nrow(res.mat))))
  colnames(matri)[ncol(matri)] = "Pathway Number"
  matri <- cbind(matri, EC.Hits)
  fast.write.csv(matri, file=mSetObj$mum_nm_csv, row.names=TRUE);
  json.mat <- RJSONIO::toJSON(json.res, .na='null');
  sink(mSetObj$mum_nm);
  cat(json.mat);
  sink();
  return(mSetObj);
}

#' Internal function for calculating GSEA, no RT
.compute.mummichog.fgsea <- function(mSetObj, permNum){
  
  num_perm <- permNum;
  total_cpds <- mSetObj$cpd_exp #scores from all matched compounds
  
  current.mset <- mSetObj$pathways$cpds; #all compounds per pathway
  names(current.mset) <- mSetObj$pathways$name
  path.size <- unlist(lapply(mSetObj$pathways$cpds, length)) #total size of pathways
  
  df.scores <- data.frame(id=names(total_cpds), scores=total_cpds)
  ag.scores <- aggregate(id ~ scores, data = df.scores, paste, collapse = "; ")
  
  ag.sorted <- ag.scores[order(-ag.scores$scores),]
  row.names(ag.sorted) <- NULL
  
  dt.scores <- data.table::data.table(ag.sorted)
  dt.scores.out <- dt.scores[, list(scores=scores, id = unlist(strsplit(id, "; ", fixed = TRUE))), by=1:nrow(dt.scores)]
  
  rank.vec <- as.numeric(dt.scores.out$nrow)
  names(rank.vec) <- as.character(dt.scores.out$id)
  
  scores.vec <- as.numeric(ag.sorted$scores)
  names(scores.vec) <- as.character(ag.sorted$id)
  
  # run fgsea
  if(mSetObj$paramSet$mumDataContainsPval == 0){
    rank.vec = seq.int(1, length(mSetObj$cpd_exp))
    names(rank.vec) <- names(mSetObj$cpd_exp)
    scores.vec = seq.int(1, length(mSetObj$cpd_exp))
    names(scores.vec) <- names(mSetObj$cpd_exp)
  }
  
  fgseaRes <- fgsea2(mSetObj, current.mset, scores.vec, rank.vec, num_perm)
  
  res.mat <- matrix(0, nrow=length(fgseaRes$pathway), ncol=5)
  
  path.size <- unlist(lapply(current.mset, length))
  matched.size <- path.size[match(fgseaRes$pathway, names(path.size))]
  
  # create result table
  res.mat[,1] <- matched.size
  res.mat[,2] <- fgseaRes$size
  res.mat[,3] <- fgseaRes$pval
  res.mat[,4] <- fgseaRes$padj
  res.mat[,5] <- fgseaRes$NES
  
  rownames(res.mat) <- fgseaRes$pathway
  colnames(res.mat) <- c("Pathway_Total", "Hits", "P_val", "P_adj", "NES")
  
  # order by p-values
  ord.inx <- order(res.mat[,3]);
  res.mat <- signif(as.matrix(res.mat[ord.inx, ]), 4);
  
  mSetObj$mummi.gsea.resmat <- res.mat;
  
  Cpd.Hits <- qs::qread("pathwaysFiltered.qs")
  Cpd.Hits <- unlist(lapply(seq_along(Cpd.Hits), function(i) paste(names(Cpd.Hits[[i]]), collapse = ";")))
  Cpd.Hits <- Cpd.Hits[Cpd.Hits != ""]
  
  res.mat <- cbind(res.mat, Cpd.Hits[ord.inx])
  fast.write.csv(res.mat, file="mummichog_fgsea_pathway_enrichment.csv", row.names=TRUE);
  
  matched_cpds <- names(mSetObj$cpd_exp)
  inx2<-na.omit(match(rownames(res.mat), mSetObj$pathways$name))
  filt_cpds <- lapply(inx2, function(f) { mSetObj$pathways$cpds[f] })
  
  cpds <- lapply(filt_cpds, function(x) intersect(unlist(x), matched_cpds))
  mSetObj$path.nms <- rownames(res.mat)
  mSetObj$path.hits<- convert2JsonList(cpds)
  mSetObj$path.pval <- as.numeric(res.mat[,3])
  json.res <- list(cmpd.exp = total_cpds,
                   path.nms = rownames(res.mat),
                   hits.all = convert2JsonList(cpds),
                   hits.all.size = as.numeric(res.mat[,2]),
                   nes = fgseaRes$NES,
                   fisher.p = as.numeric(res.mat[,3]))
  
  json.mat <- RJSONIO::toJSON(json.res, .na='null');
  sink(mSetObj$mum_nm);
  cat(json.mat);
  sink();
  
  return(mSetObj);
}

#' Internal function for calculating GSEA, with RT
#' @importFrom stats aggregate
.compute.mummichog.RT.fgsea <- function(mSetObj, permNum){
  
  #Declare variable
  scores <- NULL;
  
  # Need to perform in EC space
  num_perm <- permNum;
  total_ecpds <- mSetObj$ec_exp #scores from all matched compounds
  
  current.mset <- mSetObj$pathways$emp_cpds; #all compounds per pathway
  names(current.mset) <- mSetObj$pathways$name
  path.size <- unlist(lapply(mSetObj$pathways$ecpds, length)) #total size of pathways
  
  df.scores <- data.frame(id=names(total_ecpds), scores=total_ecpds)
  ag.scores <- aggregate(id ~ scores, data = df.scores, paste, collapse = "; ")
  
  ag.sorted <- ag.scores[order(-ag.scores$scores),]
  row.names(ag.sorted) <- NULL
  
  dt.scores <- data.table::data.table(ag.sorted)
  dt.scores.out <- dt.scores[, list(scores=scores, id = unlist(strsplit(id, "; ", fixed = TRUE))), by=1:nrow(dt.scores)]
  
  rank.vec <- as.numeric(dt.scores.out$nrow)
  names(rank.vec) <- as.character(dt.scores.out$id)
  
  scores.vec <- as.numeric(ag.sorted$scores)
  names(scores.vec) <- as.character(ag.sorted$id)
  
  # run fgsea
  if(mSetObj$paramSet$mumDataContainsPval == 0){
    rank.vec = seq.int(1, length(mSetObj$ec_exp))
    names(rank.vec) <- names(mSetObj$ec_exp)
    scores.vec = seq.int(1, length(mSetObj$ec_exp))
    names(scores.vec) <- names(mSetObj$ec_exp)
  }
  
  fgseaRes <- fgsea2(mSetObj, current.mset, scores.vec, rank.vec, num_perm)
  
  res.mat <- matrix(0, nrow=length(fgseaRes$pathway), ncol=5)
  
  path.size <- unlist(lapply(current.mset, length))
  matched.size <- path.size[match(fgseaRes$pathway, names(path.size))]
  
  # create result table
  res.mat[,1] <- matched.size
  res.mat[,2] <- fgseaRes$size
  res.mat[,3] <- fgseaRes$pval
  res.mat[,4] <- fgseaRes$padj
  res.mat[,5] <- fgseaRes$NES
  
  rownames(res.mat) <- fgseaRes$pathway
  colnames(res.mat) <- c("Pathway_Total", "Hits", "P_val", "P_adj", "NES")
  
  # order by p-values
  ord.inx <- order(res.mat[,3]);
  res.mat <- signif(as.matrix(res.mat[ord.inx, ]), 4);
  
  mSetObj$mummi.gsea.resmat <- res.mat;
  
  EC.Hits <- qs::qread("pathwaysFiltered.qs")
  EC.Hits <- lapply(seq_along(EC.Hits), function(i) paste(names(EC.Hits[[i]]), collapse = ";"))
  res.mat <- cbind(res.mat, EC.Hits)
  fast.write.csv(res.mat, file="mummichog_fgsea_pathway_enrichment.csv", row.names=TRUE);
  
  # need to convert ECs to compounds for json
  total_ecpds <- unique(mSetObj$total_matched_ecpds) #all matched compounds
  current.mset <- current.mset[match(rownames(res.mat), mSetObj$pathways$name)]
  ecpds <- lapply(current.mset, function(x) intersect(x, total_ecpds)); #pathways & all ref ecpds
  cpds <- lapply(ecpds, function(x) unique(unlist(mSetObj$ecpd_cpd_dict[match(x, names(mSetObj$ecpd_cpd_dict))])) )
  
  # now make exp vec for all compounds
  cpds2ec <- mSetObj$cpd_ecpd_dict
  cpds.all <- unique(unlist(mSetObj$ecpd_cpd_dict[match(total_ecpds, names(mSetObj$ecpd_cpd_dict))]))
  cpds.exp <- sapply(cpds.all, function(x) sapply(seq_along(x), function(i) mean(mSetObj$ec_exp[match(unique(unlist(cpds2ec[match(x[[i]], names(cpds2ec))])), names(mSetObj$ec_exp))]) ) )
  
  mSetObj$path.nms <- rownames(res.mat)
  mSetObj$path.hits <- convert2JsonList(cpds)
  mSetObj$path.pval <- as.numeric(res.mat[,3])
  
  json.res <- list(cmpd.exp = cpds.exp,
                   path.nms = rownames(res.mat),
                   hits.all = convert2JsonList(cpds),
                   hits.all.size = as.numeric(res.mat[,2]),
                   nes = fgseaRes$NES,
                   fisher.p = as.numeric(res.mat[,3]))
  
  json.mat <- RJSONIO::toJSON(json.res, .na='null');
  sink(mSetObj$mum_nm);
  cat(json.mat);
  sink();
  
  return(mSetObj);
}

#####################################################################

#'Map currency metabolites to KEGG & BioCyc
#'@description This function maps the user selected list
#'of compounds to its corresponding KEGG IDs and BioCyc IDs
#'@param mSetObj Input the name of the created mSetObj object 
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PerformCurrencyMapping <- function(mSetObj = NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  qvec <- mSetObj$dataSet$cmpd;
  curr_db <- .get.my.lib("currency_cmpd.qs");
  hit.inx <- match(tolower(qvec), tolower(curr_db$DisplayName));
  
  num_hits <- length(na.omit(hit.inx))
  
  if(num_hits == 0){
    mSetObj$mummi$curr.msg <- c("No currency metabolites were selected or mapped!")
    print(mSetObj$mummi$curr.msg)
    return(0)
  }
  
  match.values <- curr_db[hit.inx,];
  curr.met <- nrow(match.values)
  
  mSetObj$curr.map <- match.values
  
  if(curr.met > 0){
    mSetObj$mummi$curr.msg <- paste("A total of ", curr.met ," currency metabolites were successfully uploaded!", sep = "")
  }
  
  mSetObj$curr.cust <- TRUE;
  return(.set.mSet(mSetObj));
}

#'Read Adduct List
#'@description This function reads in the user's adduct list and 
#'saves it as a matrix. 
#'@usage Read.AdductData(mSetObj=NA, adductList)
#'@param mSetObj Input the name of the created mSetObj object 
#'@param adductList Input the name of the adduct list
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PerformAdductMapping <- function(mSetObj=NA, add.mode){
  
  mSetObj <- .get.mSet(mSetObj);
  
  adducts <- mSetObj$dataSet$adduct.list
  
  if(add.mode == "positive"){
    add_db <- .get.my.lib("pos_adduct.qs");
  }else if(add.mode == "negative"){
    add_db <- .get.my.lib("neg_adduct.qs");
  }else if(add.mode == "mixed"){
    add_db <- .get.my.lib("mixed_adduct.qs");
  }else{
    msg <- c("Adduct mode is not valid")
  }
  
  hit.inx <- match(tolower(adducts), tolower(add_db$Ion_Name));  
  hits <- length(na.omit(hit.inx))
  
  if(hits == 0){
    mSetObj$mummi$add.msg <- c("No adducts were selected!");
    return(0)
  }
  
  match.values <- add_db[na.omit(hit.inx),];
  sel.add <- nrow(match.values);
  
  if(sel.add > 0){
    mSetObj$mummi$add.msg <- paste("A total of ", sel.add ," adducts were successfully selected!", sep = "")
  }
  
  mSetObj$dataSet$adduct.custom <- TRUE
  mSetObj$dataSet$add.map <- match.values
  
  return(.set.mSet(mSetObj));
}

# internal function to create new mz matrix from user-curated list of adducts
new_adduct_mzlist <- function(mSetObj=NA, mw){
  
  # if(!exists("metaFiles")){
  #   mSetObj <- .get.mSet(mSetObj);
  # }
  
  mode <- mSetObj$dataSet$mode;
  
  ion.name <- mSetObj$dataSet$add.map$Ion_Name
  ion.mass <- mSetObj$dataSet$add.map$Ion_Mass
  
  mw_modified <- NULL;
  
  if(mode!="mixed"){ #pos or neg
    
    mass.list <- as.list(ion.mass)
    mass.user <- lapply(mass.list, function(x) eval(parse(text=paste(gsub("PROTON", 1.00727646677, x)))) )
    mw_modified <- cbind(mw, do.call(cbind, mass.user));
    
    if(mode == "positive"){
      mw_modified.pos <- mw_modified[,-1, drop = FALSE]
      mw_modified.neg <- as.matrix(mw_modified[,1, drop = FALSE])
      colnames(mw_modified.pos) <- ion.name;
      colnames(mw_modified.neg) <- "M"
    }else{ #negative
      mw_modified.neg <- mw_modified[,-1, drop = FALSE]
      mw_modified.pos <- as.matrix(mw_modified[,1, drop = FALSE])
      colnames(mw_modified.neg) <- ion.name;
      colnames(mw_modified.pos) <- "M"
    }
    
    mw_modified <- list(mw_modified.neg, mw_modified.pos)
    names(mw_modified) <- c("neg", "pos")
    
  } else {
    #deal w. mixed ion mode, need to return pos and neg 
    
    neg.ions <- c("M-H [1-]", "M-2H [2-]", "M-3H [3-]", "M-H2O-H [1-]", "M-H+O [1-]", "M+K-2H [1-]", "M+Na-2H [1- ]", "M+Cl [1-]", "M+Cl37 [1-]",   
                  "M+K-2H [1-]", "M+FA-H [1-]", "M+Hac-H [1-]", "M+Br [1-]", "M+Br81 [1-]", "M+TFA-H [1-]", "M+ACN-H [1-]", "M+HCOO [1-]", "M+CH3COO [1-]", 
                  "2M-H [1-]", "2M+FA-H [1-]", "2M+Hac-H [1-]", "3M-H [1-]", "M(C13)-H [1-]", "M(S34)-H [1-]", "M(Cl37)-H [1-]")
    
    ion.name.neg <- intersect(ion.name, neg.ions)
    ion.mass.neg <- ion.mass[which(ion.name %in% neg.ions)] 
    
    ion.name.pos <- setdiff(ion.name, neg.ions)
    ion.mass.pos <- ion.mass[which(ion.name %in% ion.name.pos)] 
    
    mass.list.neg <- as.list(ion.mass.neg)
    mass.user.neg <- lapply(mass.list.neg, function(x) eval(parse(text=paste(gsub("PROTON", 1.00727646677, x)))) )
    mw_modified.neg <- do.call(cbind, mass.user.neg);
    colnames(mw_modified.neg) <- ion.name.neg;
    
    mass.list.pos <- as.list(ion.mass.pos)
    mass.user.pos <- lapply(mass.list.pos, function(x) eval(parse(text=paste(gsub("PROTON", 1.00727646677, x)))) )
    mw_modified.pos <- do.call(cbind, mass.user.pos);
    colnames(mw_modified.pos) <- ion.name.pos;
    
    mw_modified <- list(mw_modified.neg, mw_modified.pos)
    names(mw_modified) <- c("neg", "pos")
  }
  return(mw_modified);
}

#'Update the mSetObj with user-selected parameters for MS Peaks to Pathways.
#'@description This functions handles updating the mSet object for mummichog analysis. 
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param force_primary_ion Character, if "yes", only mz features that match compounds with a primary ion are kept.
#'@param rt_tol Numeric. Input the retention time tolerance used for determining ECs (in seconds).
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
UpdateEC_Rules <- function(mSetObj = NA, force_primary_ion, rt_tol){
  
  mSetObj <- .get.mSet(mSetObj);
  
  mSetObj$dataSet$primary_ion <- force_primary_ion;
  
  ok <- is.numeric(rt_tol)
  
  if(ok){
    mSetObj$dataSet$rt_tol <- rt_tol;
  }else{
    AddErrMsg("Retention time tolerance must be numeric!")
    return(0)
  }
  
  msg.vec <- "EC Rules successfully updated."
  mSetObj$mummi$ec.msg <- msg.vec
  
  return(.set.mSet(mSetObj));
}

##############################
##### Plotting Functions #####
##############################

#'Plot MS Peaks to Paths mummichog permutation p-values
#'@description Plots the mummichog permutation p-values
#'@param mSetObj Input name of the created mSet Object
#'@param pathway Input the name of the pathway
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width. 
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotMSPeaksPerm <- function(mSetObj=NA, pathway, imgName, format="png", dpi=72, width=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  bw.vec <- unlist(mSetObj$perm_record);
  
  len <- length(bw.vec);
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 8;
  }else if(width == 0){
    w <- 7;
  }else{
    w <- width;
  }
  h <- w*6/8;
  
  mSetObj$imgSet$mspeaks.permut <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(mar=c(5,5,2,4));
  hst <- hist(bw.vec, breaks = "FD", freq=T,
              ylab="Frequency", xlab= 'Permutation test statistics', col="lightblue", main="");
  
  # add the indicator using original label
  h <- max(hst$counts)
  inx <- which(rownames(mSetObj$mummi.resmat) == pathway)
  raw.p <- mSetObj$mummi.resmat[inx,5]
  arrows(raw.p, h/5, raw.p, 0, col="red", lwd=2);
  text(raw.p, h/3.5, paste('Raw \n statistic \n', raw.p), xpd=T);
  dev.off();
  return(.set.mSet(mSetObj));
}

#' PlotPeaks2Paths
#' @description Plots either the original mummichog or GSEA results.
#' @param mSetObj Input the name of the created mSetObj object
#' @param imgName Input a name for the plot
#' @param format Character, input the format of the image to create.
#' @param dpi Numeric, input the dpi of the image to create.
#' @param width Numeric, input the width of the image to create.
#' @param Labels Character, indicate if the plot should be labeled. By default
#' it is set to "default", and the 5 top-ranked pathways per each algorithm will be plotted.
#' Users can adjust the number of pathways to be annotated per pathway using the "num_annot" 
#' parameter.
#' @author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export

PlotPeaks2Paths <- function(mSetObj=NA, imgName, format = "png", dpi = 72, width = 9, labels = "default",
                            num_annot = 5){
  mSetObj <- .get.mSet(mSetObj);
  #cat("mSetObj$paramSet$anal.type: --->", mSetObj$paramSet$anal.type, "\n")
  #cat("anal.type: --->", anal.type, "\n")
  anal.type0 <- mSetObj$paramSet$anal.type
  if(anal.type0 == "mummichog"){
    mummi.mat <- mSetObj$mummi.resmat
    y <- -log10(mummi.mat[,5]);
    x <- mummi.mat[,3]/mummi.mat[,4]
    pathnames <- rownames(mummi.mat)
  }else{
    gsea.mat <- mSetObj$mummi.gsea.resmat
    y <- -log10(gsea.mat[,3])
    x <- gsea.mat[,2]/gsea.mat[,1]
    pathnames <- rownames(gsea.mat)
  }
  
  inx <- order(y, decreasing= T);
  
  y <- y[inx];
  x <- x[inx]; 
  path.nms <- pathnames[inx];
  
  # set circle size based on enrichment factor
  sqx <- sqrt(x);
  min.x <- min(sqx, na.rm = TRUE);
  max.x <- max(sqx, na.rm = TRUE);
  
  if(min.x == max.x){ # only 1 value
    max.x = 1.5*max.x;
    min.x = 0.5*min.x;
  }
  
  maxR <- (max.x - min.x)/40;
  minR <- (max.x - min.x)/160;
  radi.vec <- minR+(maxR-minR)*(sqx-min.x)/(max.x-min.x);
  
  # set background color according to combo.p
  bg.vec <- heat.colors(length(y));
  
  if(format == "png"){
    bg = "transparent";
  }else{
    bg="white";
  }
  
  if(is.na(width)){
    w <- 7;
  }else if(width == 0){
    w <- 7;
  }else{
    w <- width;
  }
  h <- w;
  
  df <- data.frame(path.nms, x, y)
  
  if(labels == "default"){
    pk.inx <- GetTopInx(df$y, num_annot, T)
  }
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  mSetObj$imgSet$mummi.plot <- imgName
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg=bg);
  op <- par(mar=c(6,5,2,3));
  plot(x, y, type="n", axes=F, xlab="Enrichment Factor", ylab="-log10(p)", bty = "l");
  axis(1);
  axis(2);
  symbols(x, y, add = TRUE, inches = F, circles = radi.vec, bg = bg.vec, xpd=T);
  
  if(labels=="default"){
    text(x[pk.inx], y[pk.inx], labels = path.nms[pk.inx], pos=3, xpd=T, cex=0.8)
  }else if(labels == "all"){
    text(x, y, labels = path.nms, pos=3, xpd=T, cex=0.8)
  }
  
  par(op);
  dev.off();
  if(anal.type0 == "mummichog"){
    df <- list(pval=unname(y), enr=unname(x), pathnames=path.nms);
    sink("scattermum.json");
    cat(RJSONIO::toJSON(df));
    sink();
  }else{
    df <- list(pval=unname(y), enr=unname(gsea.mat[,5]), pathnames=path.nms);
    sink("scattergsea.json");
    cat(RJSONIO::toJSON(df));
    sink();
  }
  
  return(.set.mSet(mSetObj));
  
}

#' PlotIntegPaths
#' @description Plots both the original mummichog and the GSEA results by combining p-values
#' using the Fisher's method (sumlog). 
#' @param mSetObj Input the name of the created mSetObj object
#' @param imgName Input a name for the plot
#' @param format Character, input the format of the image to create.
#' @param dpi Numeric, input the dpi of the image to create.
#' @param width Numeric, input the width of the image to create.
#' @param Labels Character, indicate if the plot should be labeled. By default
#' it is set to "default", and the 5 top-ranked pathways per each algorithm will be plotted.
#' Users can adjust the number of pathways to be annotated per pathway using the "labels.x" 
#' and "labels.y" parameters.
#' Users can set this to "none" for no annotations, or "all" to annotate all pathways. 
#' @param labels.x Numeric, indicate the number of top-ranked pathways using the fGSEA algorithm 
#'  to annotate on the plot. 
#' @param labels.y Numeric, indicate the number of top-ranked pathways using the original 
#' mummichog algorithm to annotate on the plot. 
#' @author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import scales

PlotIntegPaths <- function(mSetObj=NA, imgName, format = "png", dpi = 72, width = 9, labels = "default", 
                           labels.x = 5, labels.y = 5, scale.axis = TRUE){
  
  mSetObj <- .get.mSet(mSetObj);
  
  # check if mummichog + gsea was performed
  if(is.null(mSetObj$mummi.resmat) | is.null(mSetObj$mummi.gsea.resmat)){
    print("Both mummichog and fGSEA must be performed!")
    return(0)
  }
  
  combo.resmat <- mSetObj$integ.resmat
  pathnames <- rownames(combo.resmat)
  # Sort values based on combined pvalues
  y <- -log10(combo.resmat[,4]);
  x <- -log10(combo.resmat[,5]);
  combo.p <- -log10(combo.resmat[,6])
  
  if(scale.axis){
    y <- scales::rescale(y, c(0,4))
    x <- scales::rescale(x, c(0,4))
    combo.p <- scales::rescale(combo.p, c(0,4))
  }
  
  inx <- order(combo.p, decreasing= T);
  
  combo.p <- combo.p[inx]
  x <- x[inx]; 
  y <- y[inx];
  path.nms <- pathnames[inx];
  
  # set circle size based on combined pvalues
  min.x <- min(combo.p, na.rm = TRUE);
  max.x <- max(combo.p, na.rm = TRUE);
  
  if(min.x == max.x){ # only 1 value
    max.x = 1.5*max.x;
    min.x = 0.5*min.x;
  }
  
  maxR <- (max.x - min.x)/40;
  minR <- (max.x - min.x)/160;
  radi.vec <- minR+(maxR-minR)*(combo.p-min.x)/(max.x-min.x);
  
  # set background color according to combo.p
  bg.vec <- heat.colors(length(combo.p));
  
  if(format == "png"){
    bg = "transparent";
  }else{
    bg="white";
  }
  
  if(is.na(width)){
    w <- 7;
  }else if(width == 0){
    w <- 7;
  }else{
    w <- width;
  }
  h <- w;
  
  df <- data.frame(path.nms, x, y)
  
  if(labels == "default"){
    mummi.inx <- GetTopInx(df$y, labels.y, T)
    gsea.inx <- GetTopInx(df$x, labels.x, T)
    all.inx <- mummi.inx | gsea.inx;
  }
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  mSetObj$imgSet$integpks.plot <- imgName
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg=bg);
  op <- par(mar=c(6,5,2,3));
  
  # color blocks only make sense if scaled...
  if(scale.axis){
    plot(x, y, type="n", axes=F, xlab="GSEA -log10(p)", ylab="Mummichog -log10(p)", bty = "l");
    axis(1);
    axis(2);
    symbols(x, y, add = TRUE, inches = F, circles = radi.vec, bg = bg.vec, xpd=T);
    
    axis.lims <- par("usr")
    
    # mummichog sig
    mum.x <- c(axis.lims[1], axis.lims[1], axis.lims[2], axis.lims[2])
    mum.y <- c(2, axis.lims[4], axis.lims[4], 2)
    polygon(mum.x, mum.y, col=rgb(82/255,193/255,188/255,0.3), border = NA)
    
    # gsea sig
    gsea.x <- c(2,2,axis.lims[4],axis.lims[4])
    gsea.y <- c(axis.lims[1],axis.lims[4],axis.lims[4],axis.lims[1])
    polygon(gsea.x, gsea.y, col=rgb(216/255,126/255,178/255,0.3), border = NA)
  }else{
    plot(x, y, type="n", xlim=c( 0, round(max(x)) ), ylim=c(0, round(max(y)) ), xlab="GSEA -log10(p)", ylab="Mummichog -log10(p)", bty = "l");
    symbols(x, y, add = TRUE, inches = F, circles = radi.vec, bg = bg.vec, xpd=T);
  }
  
  if(labels=="default"){
    text(x[all.inx], y[all.inx], labels = path.nms[all.inx], pos=3, xpd=T, cex=0.8)
  }else if(labels == "all"){
    text(x, y, labels = path.nms, pos=3, xpd=T, cex=0.8)
  }
  
  par(op);
  dev.off();
  
  df <- list(pval=unname(y), enr=unname(x), metap= unname(combo.p), pathnames=pathnames);
  sink("scatterinteg.json");
  cat(RJSONIO::toJSON(df));
  sink();
  
  return(.set.mSet(mSetObj));
}

#' Plot m/z hits in a pathway
#' @description Function to create a boxplot of m/z features
#' within a specific pathway. m/z features used by the original
#' mummichog algorithm are highlighted with an asterisk. 
#' @param mSetObj Input the name of the created mSetObj object.
#' @param msetNM Character, input the name of the pathway. 
#' @param format Character, input the format of the image to create. 
#' @param dpi Numeric, input the dpi of the image to create. Default 
#' is set to 300.
#' @param width Numeric, input the width of the image to create.
#' Default is set to 10.
#' @author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import ggplot2
#' @importFrom scales rescale

PlotPathwayMZHits <- function(mSetObj=NA, msetNM, format="png", dpi=300, width=10){
  
  mSetObj <- .get.mSet(mSetObj);
  mum.version <- mSetObj$paramSet$version;
  inx <- which(mSetObj$pathways$name == msetNM)
  
  # Brief sanity check
  if(length(inx) == 0){
    print(paste0(msetNM, " is not found in the pathway library! Please verify the spelling."))
    AddErrMsg("Invalid pathway selected!")
    return(0)
  }
  
  if(mum.version=="v2" & mSetObj$paramSet$mumRT){
    mset <- mSetObj$pathways$emp_cpds[[inx]];
    mzs <- as.numeric(unique(unlist(mSetObj$ec2mz_dict[mset])))
  }else{
    mset <- mSetObj$pathways$cpds[[inx]];
    mzs <- as.numeric(unique(unlist(mSetObj$cpd2mz_dict[mset])))
  }
  
  if(anal.type == "integ_peaks"){
    # check if mummichog + gsea was performed
    if(is.null(mSetObj$mummi.resmat) | is.null(mSetObj$mummi.gsea.resmat)){
      AddErrMsg("Both mummichog and fGSEA must be performed!")
      return(0)
    }
  }
  
  pvals <- mSetObj$dataSet$mummi.proc[mSetObj$dataSet$mummi.proc[, 2] %in% mzs, ]
  
  pval.cutoff <- mSetObj$dataSet$cutoff
  
  if(anal.type %in% c("integ_peaks", "mummichog")){
    used.inx <- pvals[,1] < pval.cutoff
    mummi <- which(used.inx)
    mummi_mzs <- pvals[,2]
    
    # add astericks if used by mummichog
    mummi_mzs_star <- mummi_mzs
    mummi_mzs_star[mummi] <- paste(mummi_mzs_star[mummi], "*",sep="");
  }else{
    mummi_mzs_star <- pvals[,2]
  }
  
  # create boxdata
  data <- data.frame(mSetObj$dataSet$proc)
  
  boxdata <- data[,as.character(mummi_mzs)]
  colnames(boxdata) <- mummi_mzs_star
  boxdata$class <- mSetObj$dataSet$cls
  
  num.feats <- length(result)
  
  boxdata.m <- melt(boxdata, id.vars="class")
  boxdata.m$value <- scales::rescale(boxdata.m$value, c(0,1))
  
  boxplotName <- paste(msetNM, ".", format, sep="");
  
  if(num.feats == 1){
    
    w = width
    h = width
    p <- ggplot(data = boxdata.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=class), outlier.shape = NA, outlier.colour=NA)
    p <- p + ggtitle(msetNM) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Group"))
    p <- p + xlab("m/z feature") + ylab("Intensity")
    
    ggsave(p, filename = boxplotName, dpi=300, width=w, height=h, limitsize = FALSE)
    
    return(mSetObj)
    
  }else if(num.feats<10){
    w = width
    h <- num.feats
    cols = 3
  }else if(num.feats<5){
    w = width
    h = width
  }else{
    w = width
    h <- num.feats * 0.35
    cols = 6
  }
  
  p <- ggplot(data = boxdata.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=class), outlier.shape = NA, outlier.colour=NA)
  p <- p + facet_wrap( ~ variable, scales="free", ncol=cols) + xlab("m/z features") + ylab("Intensity")
  p <- p + ggtitle(msetNM) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Group"))
  
  ggsave(p, filename = boxplotName, dpi=300, width=w, height=h, limitsize = FALSE)
  
  return(.set.mSet(mSetObj));
}

###############################
######## For R Package ########
###############################

#' Function to get compound details from a specified pathway
#' @description Function to get compound details from a specified pathway.
#' The results will be both printed in the console as well as saved
#' as a csv file. Note that performing this function multiple times will
#' overwrite previous queries. Significant compounds will be indicated with an asterisk.
#' @param mSetObj Input the name of the created mSetObj object.
#' @param msetNm Input the name of the pathway
#' @export
GetMummichogPathSetDetails <- function(mSetObj=NA, msetNm){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(!is.null(mSetObj$api$lib)){
    
    # get file from api.metaboanalyst.ca
    toSend = list(pathName = msetNm)
    
    load_httr()
    base <- api.base
    endpoint <- paste0("/pathsetdetails/", mSetObj$api$guestName)
    call <- paste(base, endpoint, sep="")
    query_results <- httr::POST(call, body = toSend, encode= "json")
    
    if(query_results$status_code == 200){
      filename <- httr::content(query_results, "text")
    }
    
    endpointfile <- paste0("/getFile", "/", mSetObj$api$guestName, "/", filename)
    callfile <- paste(base, endpointfile, sep="")
    download.file(callfile, destfile = basename(filename))
    print(paste0(filename, " saved to current working directory!"))
    return(.set.mSet(mSetObj));
  }
  
  version <- mum.version <- mSetObj$paramSet$version;
  inx <- which(mSetObj$pathways$name == msetNm)
  
  if(is.na(inx)){
    AddErrMsg("Invalid pathway name!")
    return(0)
  }
  
  if(version=="v2" & mSetObj$paramSet$mumRT){
    mset <- mSetObj$pathways$emp_cpds[[inx]];
    mset_cpds <- mSetObj$pathways$cpds[[inx]];
    
    hits.all <- unique(mSetObj$total_matched_ecpds)
    hits.sig <- mSetObj$input_ecpdlist;
    
    refs <- mset %in% hits.all;
    sigs <- mset %in% hits.sig;
    
    ref.ecpds <- mset[which(refs & !sigs)]
    sig.ecpds <- mset[sigs]
    
    ref.mzs <- lapply(ref.ecpds, function(x) paste(as.numeric(unique(unlist(mSetObj$ec2mz_dict[x]))), collapse = "; ")) 
    sig.mzs <- lapply(sig.ecpds, function(x) paste(as.numeric(unique(unlist(mSetObj$ec2mz_dict[x]))), collapse = "; "))  
    
    ref.cpds <- lapply(ref.ecpds, function(x) paste(unique(unlist(mSetObj$ecpd_cpd_dict[x])), collapse = "; "))
    sig.cpds <- lapply(sig.ecpds, function(x) paste(unique(unlist(mSetObj$ecpd_cpd_dict[x])), collapse = "; "))
    
    path.results <- matrix(c(unlist(sig.mzs), unlist(ref.mzs), unlist(sig.cpds), unlist(ref.cpds)), ncol=2) 
    colnames(path.results) <- c("mzs", "cpds")
    rownames(path.results) <- c(paste0(sig.ecpds, "*"), ref.ecpds)
    
    name <- paste0(gsub(" ", "_", msetNm), "_ecpd_mz_info.csv")
    fast.write.csv(path.results, name)
  }else{
    mset <- mSetObj$pathways$cpds[[inx]];
    
    hits.all <- unique(mSetObj$total_matched_cpds)
    hits.sig <- mSetObj$input_cpdlist;
    
    refs <- mset %in% hits.all;
    sigs <- mset %in% hits.sig;
    
    ref.cpds <- mset[which(refs & !sigs)]
    sig.cpds <- mset[sigs]
    
    ref.mzs <- lapply(ref.cpds, function(x) paste(as.numeric(unique(unlist(mSetObj$cpd2mz_dict[x]))), collapse = "; ")) 
    sig.mzs <- lapply(sig.cpds, function(x) paste(as.numeric(unique(unlist(mSetObj$cpd2mz_dict[x]))), collapse = "; "))  
    
    path.results <- matrix(c(unlist(sig.mzs), unlist(ref.mzs)), ncol=1) 
    colnames(path.results) <- "mzs"
    rownames(path.results) <- c(paste0(sig.cpds, "*"), ref.cpds)
    
    name <- paste0(gsub(" ", "_", msetNm), "_cpd_mz_info.csv")
    fast.write.csv(path.results, name)
  }
  
  return(.set.mSet(mSetObj));
}

#' Function to get adduct details from a specified compound
#' @description Function to get adduct details from a specified compound.
#' The results will be both printed in the console as well as saved
#' as a csv file. Note that performing this function multiple times will
#' overwrite previous queries.
#' @param mSetObj Input the name of the created mSetObj object.
#' @param cmpd.id Input the name of the selected compound.
#'@import qs
#' @export
GetCompoundDetails <- function(mSetObj=NA, cmpd.id){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(!is.null(mSetObj$api$lib)){
    
    # get file from api.metaboanalyst.ca
    toSend = list(cmpdId = cmpd.id)
    
    load_httr()
    base <- api.base
    endpoint <- paste0("/compounddetails/", mSetObj$api$guestName)
    call <- paste(base, endpoint, sep="")
    query_results <- httr::POST(call, body = toSend, encode= "json")
    
    if(query_results$status_code == 200){
      filename <- httr::content(query_results, "text")
    }
    
    endpointfile <- paste0("/getFile", "/", mSetObj$api$guestName, "/", filename)
    callfile <- paste(base, endpointfile, sep="")
    download.file(callfile, destfile = basename(filename))
    print(paste0(filename, " saved to current working directory!"))
    return(.set.mSet(mSetObj));
  }
  
  forms <- mSetObj$cpd_form_dict[[cmpd.id]];
  
  if(is.null(forms)){
    print("This compound is not valid!")
    return(0)
  }
  
  matched_res <- qs::qread("mum_res.qs");
  mz <- matched_res[which(matched_res$Matched.Compound == cmpd.id), 1] 
  mass.diff <- matched_res[which(matched_res$Matched.Compound == cmpd.id), 4]
  tscores <- mSetObj$cpd_exp_dict[[cmpd.id]];
  
  res <- cbind(rep(cmpd.id, length(mz)), mz, forms, mass.diff, tscores) 
  colnames(res) <- c("Matched.Compound", "m.z", "Matched.Form", "Mass.Diff", "T.Scores")
  fast.write.csv(res, "mummichog_compound_details.csv")
  return(.set.mSet(mSetObj));
}

# Function to return the unique m/zs from the selected pathways 
# based on the compounds

GetMummichogMZHits <- function(mSetObj=NA, msetNm){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(!is.null(mSetObj$api$lib)){
    
    # get file from api.metaboanalyst.ca
    toSend = list(pathName = msetNm)
    
    load_httr()
    base <- api.base
    endpoint <- paste0("/mzhits/", mSetObj$api$guestName)
    call <- paste(base, endpoint, sep="")
    query_results <- httr::POST(call, body = toSend, encode= "json")
    
    if(query_results$status_code == 200){
      result <- httr::content(query_results, "text")
    }
    mSetObj$mz.hits <- result
    print(paste0("Unique m/z features in ", msetNm, ": ", result))
    return(.set.mSet(mSetObj));
  }
  
  inx <- which(mSetObj$pathways$name == msetNm)
  mset <- mSetObj$pathways$cpds[[inx]];
  mzs <- as.numeric(unique(unlist(mSetObj$cpd2mz_dict[mset]))) 
  result <- intersect(mzs, mSetObj$dataSet$input_mzlist)
  mSetObj$mz.hits <- result
  
  return(.set.mSet(mSetObj));
}

###############################
####### Getters For Web #######
###############################

GetMatchingDetails <- function(mSetObj=NA, cmpd.id){
  
  mSetObj <- .get.mSet(mSetObj);
  forms <- mSetObj$cpd_form_dict[[cmpd.id]];
  tscores <- mSetObj$cpd_exp_dict[[cmpd.id]];
  # create html table
  res <- paste("<li>", "<b>", forms, "</b>: ", tscores, "</li>",sep="", collapse="");
  return(res);
}

GetMummichogHTMLPathSet <- function(mSetObj=NA, msetNm){
  
  mSetObj <- .get.mSet(mSetObj);
  inx <- which(mSetObj$pathways$name == msetNm)
  
  mset <- mSetObj$pathways$cpds[[inx]];
  mum.version <- mSetObj$paramSet$version;
  
  # all matched compounds
  if(mum.version == "v2" & mSetObj$paramSet$mumRT){
    hits.all <- unique(unlist(mSetObj$ecpd_cpd_dict))
  }else{
    hits.all <- unique(mSetObj$total_matched_cpds) 
  }
  
  if(anal.type == "mummichog"|anal.type == "integ_peaks"){
    
    # get the sig compounds
    if(mum.version == "v2" & mSetObj$paramSet$mumRT){
      hits.sig <- mSetObj$input_ecpdlist;
      hits.sig.inx <- match(hits.sig, names(mSetObj$ecpd_cpd_dict))
      hits.sig <- unlist(mSetObj$ecpd_cpd_dict[hits.sig.inx])
    }else{
      hits.sig <- mSetObj$input_cpdlist;
    }
    
    # highlighting with different colors
    refs <- mset %in% hits.all;
    sigs <- mset %in% hits.sig;
    
    red.inx <- which(sigs);
    blue.inx <- which(refs & !sigs);
    
    # use actual cmpd names
    nms <- mset;
    nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");
    nms[blue.inx] <- paste("<font color=\"blue\">", "<b>", nms[blue.inx], "</b>", "</font>",sep="");
  }else{
    refs <- mset %in% hits.all;
    red.inx <- which(refs);
    nms <- mset;
    nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");
  }
  
  # for large number, return only hits (context already enough)
  if(length(nms) > 200){
    nms <- nms[refs];
  }
  return(cbind(msetNm, paste(unique(nms), collapse="; ")));
}

GetMummiResMatrix <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  anal.type0 <- mSetObj$paramSet$anal.type;
  if(anal.type0 == "mummichog"){    
    return(mSetObj$mummi.resmat);
  }else if(anal.type0 == "gsea_peaks"){
    return(mSetObj$mummi.gsea.resmat);
  }else{
    return(mSetObj$integ.resmat);
  }
}

GetMummiResRowNames <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  anal.type0 <- mSetObj$paramSet$anal.type;
  if(anal.type0 == "mummichog"){
    return(rownames(mSetObj$mummi.resmat));
  }else if(anal.type0 == "gsea_peaks"){
    return(rownames(mSetObj$mummi.gsea.resmat));
  }else{
    return(rownames(mSetObj$integ.resmat));
  }
}

GetMummiResColNames <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(colnames(mSetObj$mummi.resmat));
}

GetCurrencyMsg <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$mummi$curr.msg)
}

GetAdductMsg <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$mummi$add.msg)
}

GetECMsg <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$mummi$ec.msg)
}

GetDefaultPvalCutoff <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet();
  peakFormat <- mSetObj$paramSet$peakFormat;
  
  if(peakFormat %in% c("rmp", "rmt")){
    maxp <- 0;
  }else{
    pvals <- c(0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
    ndat <- mSetObj$dataSet$mummi.proc;
    
    ### Handle something very wrong for mass table
    if(is.null(ndat)){
      res <- PerformFastUnivTests(mSetObj$dataSet$norm, mSetObj$dataSet$cls, var.equal=TRUE);
      ndat <- res[order(res[,2]),];
      n <- floor(0.1*length(ndat[,2]))
      cutoff <- ndat[n+1,2]
    } else {
      n <- floor(0.1*length(ndat[,"p.value"]))
      cutoff <- ndat[n+1,1]
    }
    
    if(!any(pvals <= cutoff)){
      maxp <- 0.00001
    }else{
      maxp <- max(pvals[pvals <= cutoff])
    }
  }
  return(maxp)
}

GetDefaultRTTol <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  rt_tol <- mSetObj$dataSet$rt_tol;
  if(is.na(rt_tol)){
    rt_tol <- 0
  }
  return(rt_tol)
}

GetMummiMode <- function(mSetObj){
  mSetObj <- .get.mSet(mSetObj);
  mode <- mSetObj$dataSet$mode
  return(mode);
}

GetMummiDataType <- function(mSetObj){
  mSetObj <- .get.mSet(mSetObj);
  type <- mSetObj$dataSet$type
  return(type)
}

# Replicate because do not want to have to read in stats_univariate to perform MS Peaks
GetTopInx <- function(vec, n, dec=T){
  inx <- order(vec, decreasing = dec)[1:n];
  # convert to T/F vec
  vec<-rep(F, length=length(vec));
  vec[inx] <- T;
  return (vec);
}

#########################################
########### Utility Functions ###########
#########################################

# Global variables define currency compounds
currency <- c('C00001', 'C00080', 'C00007', 'C00006', 'C00005', 'C00003',
              'C00004', 'C00002', 'C00013', 'C00008', 'C00009', 'C00011',
              'G11113', '', 'H2O', 'H+', 'Oxygen', 'NADP+', 
              'NADPH', 'NAD+', 'NADH', 'ATP', 
              'Pyrophosphate', 'ADP', 'CO2');

all_currency <- c('C00001', 'C00080', 'C00007', 'C00006', 'C00005', 'C00003',
                  'C00004', 'C00002', 'C00013', 'C00008', 'C00009', 'C00011',
                  'G11113', '', 'H2O', 'Water', 'H+', 'Hydron', 'O2', 'Oxygen', 'NADP+', 
                  'NADP', 'NADPH', 'NAD+', 'NAD', 'NADH', 'ATP', 'Diphosphate',
                  'Pyrophosphate', 'ADP','Orthophosphate', 'CO2', 'Carbon dioxide');

primary_ions <- c('M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]',
                  'M+H [1+]', 'M+Na [1+]', 'M-H2O+H [1+]', 'M-H [1-]', 'M-2H [2-]', 'M-H2O-H [1-]')

# mz tolerance based on instrument type
# input: a vector of mz,
# output: a vector of distance tolerance
# Review on mass accuracy by Fiehn: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1464138/

mz_tolerance <- function(mz, ms.type){
  return(ms.type*1e-06*mz)
}

#'Utility function to create compound lists for permutation analysis
#'@description From a vector of m/z features, this function outputs a vector of compounds.
#'@usage make_cpdlist(mSetObj=NA, input_mzs)
#'@param mSetObj Input the name of the created mSetObj
#'@param input_mzs The vector of randomly drawn m/z features.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
make_cpdlist <- function(mSetObj=NA, input_mzs){
  cpd <- unique(unlist(mSetObj$mz2cpd_dict[input_mzs]));
  cpd <- cpd[!is.null(cpd)];
  return(cpd);
}

#'Utility function to create compound lists for permutation analysis
#'@description From a vector of m/z features, this function outputs a vector of compounds.
#'@usage make_cpdlist(mSetObj=NA, input_mzs)
#'@param mSetObj Input the name of the created mSetObj
#'@param input_mzs The vector of randomly drawn m/z features.
#'@author Jasmine Chong, Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
make_ecpdlist <- function(mSetObj=NA, input_mzs){
  ecpd <- unique(unlist(mSetObj$mz2ec_dict[input_mzs]));
  ecpd <- ecpd[!is.null(ecpd)];
  return(ecpd);
}

# Utility function to adjust for the fact that a single m/z feature can match to several compound identifiers
# input: a vector of compound ids
# output: a length of unique mzs corresponding to those compounds

count_cpd2mz <- function(cpd2mz_dict, cpd.ids,  inputmzlist){ # inputmz is either t or input cpd_list and cpd.ids are overlap features
  
  if(length(cpd.ids)==0){
    return(0);
  }
  mzs <- as.numeric(unique(unlist(cpd2mz_dict[cpd.ids])));
  if(length(mzs)==0){
    return(0);
  }else{
    result <- intersect(mzs, inputmzlist); #intersect to only mzs from the input mz list
    return(length(result));
  }
}

# convert single element vector in list to matrix
# b/c single element vector will convert to scalar in javascript, force to matrix
convert2JsonList <- function(my.list){
  lapply(my.list, function(x){
    if(length(x) == 1){
      matrix(x);
    }else{
      x;
    }
  });
}

# input: a two-col (id, val) data with potential duplicates (same id may be associated with 1 or more values
# output: a list named by unique id, with multiple values will be merged to vector
Convert2Dictionary <- function(data, quiet=T){
  
  all.ids <- data[,1];
  dup.inx <- duplicated(all.ids);
  if(sum(dup.inx) > 0){
    uniq.ids <- all.ids[!dup.inx];
    uniq.vals <- data[!dup.inx,2];
    
    # convert two-col data it to list (vals as list values, ids as list names)
    uniq.list <- split(uniq.vals, uniq.ids)
    
    # the list element orde will be sorted by the names alphabetically, need to get updated ones
    uniq.id.list <- names(uniq.list)
    
    dup.ids <- all.ids[dup.inx];
    uniq.dupids <- unique(dup.ids);
    uniq.duplen <- length(uniq.dupids);
    
    for(id in uniq.dupids){ # only update those with more than one hits
      hit.inx.all <- which(all.ids == id);
      hit.inx.uniq <- which(uniq.id.list == id);
      uniq.list[[hit.inx.uniq]]<- data[hit.inx.all,2];
    }
    
    AddMsg(paste("A total of ", sum(dup.inx), " of duplicates were merged.", sep=""));
    return(uniq.list);
  }else{
    AddMsg("All IDs are unique.");
    uniq.list <- split(data[,2], data[,1]);
    return(uniq.list);
  }
}

# utility function for fast list expanding (dynamic length)
# We need to repeatedly add an element to a list. With normal list concatenation
# or element setting this would lead to a large number of memory copies and a
# quadratic runtime. To prevent that, this function implements a bare bones
# expanding array, in which list appends are (amortized) constant time.
# https://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time-o1

myFastList <- function(capacity = 100) {
  buffer <- vector('list', capacity)
  names <- character(capacity)
  length <- 0
  methods <- list()
  
  methods$double.size <- function() {
    buffer <<- c(buffer, vector('list', capacity))
    names <<- c(names, character(capacity))
    capacity <<- capacity * 2
  }
  
  methods$add <- function(name, val) {
    if(length == capacity) {
      methods$double.size()
    }
    
    length <<- length + 1
    buffer[[length]] <<- val
    names[length] <<- name
  }
  
  methods$as.list <- function() {
    b <- buffer[0:length]
    names(b) <- names[0:length]
    return(b)
  }
  
  methods
}

####
####
#### from the fgsea R package, minor edits to adapt to untargeted metabolomics

#'Pre-ranked gsea adapted for untargeted metabolomics
#'@export
#'@import fgsea

fgsea2 <- function(mSetObj, pathways, stats, ranks,
                   nperm,
                   minSize=1, maxSize=Inf,
                   nproc=0,
                   gseaParam=1,
                   BPPARAM=NULL) {
  
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.fgsea")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_fgsea.Rc");    
    }
    return(my.fgsea(mSetObj, pathways, stats, ranks, nperm,
                    minSize, maxSize, nproc, gseaParam, BPPARAM));
  }else{
    return(my.fgsea(mSetObj, pathways, stats, ranks, nperm, 
                    minSize, maxSize, nproc, gseaParam, BPPARAM));
  }
  
  
}

calcGseaStat2 <- function(stats, selectedStats, gseaParam=1,
                          returnAllExtremes=FALSE,
                          returnLeadingEdge=FALSE) {
  
  S <- selectedStats
  r <- stats
  p <- gseaParam
  
  S <- sort(S)
  
  # account for 1 mz can be multiple cpds
  S.scores <- r[S]
  u.S <- S[!duplicated(S.scores)]
  scores <- unique(S.scores)
  
  m <- length(scores)
  N <- length(r)
  
  if (m == N) {
    stop("GSEA statistic is not defined when all genes are selected")
  }
  
  NR <- (sum(abs(scores)^p))
  rAdj <- abs(scores)^p
  if (NR == 0) {
    # this is equivalent to rAdj being rep(eps, m)
    rCumSum <- seq_along(rAdj) / length(rAdj)
  } else {
    rCumSum <- cumsum(rAdj) / NR
  }
  
  tops <- rCumSum - (u.S - seq_along(u.S)) / (N - m)
  if (NR == 0) {
    # this is equivalent to rAdj being rep(eps, m)
    bottoms <- tops - 1 / m
  } else {
    bottoms <- tops - rAdj / NR
  }
  
  maxP <- max(tops)
  minP <- min(bottoms)
  
  if(maxP > -minP) {
    geneSetStatistic <- maxP
  } else if (maxP < -minP) {
    geneSetStatistic <- minP
  } else {
    geneSetStatistic <- 0
  }
  
  if (!returnAllExtremes && !returnLeadingEdge) {
    return(geneSetStatistic)
  }
  
  res <- list(res=geneSetStatistic)
  if (returnAllExtremes) {
    res <- c(res, list(tops=tops, bottoms=bottoms))
  }
  if (returnLeadingEdge) {
    leadingEdge <- if (maxP > -minP) {
      u.S[seq_along(u.S) <= which.max(bottoms)]
    } else if (maxP < -minP) {
      rev(u.S[seq_along(u.S) >= which.min(bottoms)])
    } else {
      NULL
    }
    
    res <- c(res, list(leadingEdge=leadingEdge))
  }
  res
}

####
####
#### for heatmap view (online only)

#' SetRTincluded
#'
#' @param mSetObj mSetObj
#' @param rt retention time types, "minutes", "seconds" or "no"
#'
#' @return mSetObj
#' @export
#' 
SetRTincluded <- function(mSetObj = NA, rt = "no") {
  return(.rt.included(mSetObj = mSetObj, rt))
}
.rt.included <- function(mSetObj = NA, rt){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(rt %in% c("minutes", "seconds")){
    #is.rt <<- TRUE
    mSetObj$paramSet$mumRT <- TRUE;
    #mumRT.type <<- rt
    mSetObj$paramSet$mumRT.type <- rt;
  } else {
    #is.rt <<- FALSE
    mSetObj$paramSet$mumRT <- FALSE
    #mumRT.type <<- rt
    mSetObj$paramSet$mumRT.type <- rt;
  }
  
  return(.set.mSet(mSetObj));
}

CreateHeatmapJson <- function(mSetObj=NA, libOpt, libVersion, minLib, 
                              fileNm, filtOpt, version="v1"){
  
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.heatmap.json")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_heatmap.Rc");    
    }
    return(my.heatmap.json(mSetObj, libOpt, libVersion, minLib, fileNm, filtOpt, version));
  }else{
    return(my.heatmap.json(mSetObj, libOpt, libVersion, minLib, fileNm, filtOpt, version));
  }
}

DoPeakConversion <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  res <- Ttests.Anal(mSetObj, F, 1, FALSE, TRUE)
  
  ####### NOTE: need to use Ttests.Anal because Convert2Mummichog function takes as input result list from Ttests.anal: mSetObj$analSet$tt
  
  if(.on.public.web){
    
    is.rt <- mSetObj$paramSet$mumRT;
    mSetObj <- .get.mSet();
    .on.public.web <<- F;
    
    mSetObj<-Convert2Mummichog(mSetObj, is.rt, F, mSetObj$paramSet$mumRT.type, "tt", mSetObj$dataSet$mode);
    
    SetPeakFormat(mSetObj, "mpt")
    filename <- paste0("mummichog_input_", Sys.Date(), ".txt");
    
    mSetObj<-Read.PeakListData(mSetObj, filename);
    mSetObj<-SanityCheckMummichogData(mSetObj);
    .on.public.web <<- T;
    .set.mSet(mSetObj);
    return(1)
    
  } else {
    
    mSetObj <- res;
    is.rt <- mSetObj$paramSet$mumRT;
    
    if(mSetObj$analSet$tt$sig.num == 0){
      AddErrMsg("T-test failed, please consider trying different data normalization option!");
      return(0);
    }
    
    mSetObj <- Convert2Mummichog(mSetObj, is.rt, F, mSetObj$paramSet$mumRT.type, "tt", mSetObj$dataSet$mode);
    mSetObj$paramSet$peakFormat <- "mpt"; # SetPeakFormat("mpt")
    filename <- paste0("mummichog_input_", Sys.Date(), ".txt")
    mSetObj <- Read.PeakListData(mSetObj, filename);
    mSetObj <- SanityCheckMummichogData(mSetObj)
  }
  return(.set.mSet(mSetObj));
}

CreateListHeatmapJson <- function(mSetObj=NA, libOpt, libVersion, 
                                  minLib, fileNm, filtOpt, version="v1"){
  
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.list.heatmap")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_listheatmap.Rc");    
    }
    return(my.list.heatmap(mSetObj, libOpt, libVersion, minLib, fileNm, filtOpt, version));
  }else{
    return(my.list.heatmap(mSetObj, libOpt, libVersion, minLib, fileNm, filtOpt, version));
  }
}

#'Set organism for further analysis
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param org Set organism ID
#'@export
SetOrganism <- function(mSetObj=NA, org){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$org <- org;
  return(.set.mSet(mSetObj))
}


###### Functions got from metap package ######
sumlog <-
  function(p) {
    keep <- (p > 0) & (p <= 1)
    invalid <- sum(1L * keep) < 2
    if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(chisq = NA_real_, df = NA_integer_,
                  p = NA_real_, validp = p[keep])
    } else {
      lnp <- log(p[keep])
      chisq <- (-2) * sum(lnp)
      df <- 2 * length(lnp)
      if(length(lnp) != length(p)) {
        warning("Some studies omitted")
      }
      res <- list(chisq = chisq, df = df,
                  p = pchisq(chisq, df, lower.tail = FALSE), validp = p[keep])
    }
    class(res) <- c("sumlog", "metap")
    res
  }

sump <-
  function(p)  {
    keep <- (p >= 0) & (p <= 1)
    invalid <- sum(1L * keep) < 2
    if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(p = NA_real_, conservativep = NA_real_, validp = p[keep])
    } else {
      sigmap <- sum(p[keep])
      k <- length(p[keep])
      conservativep <- exp( k * log(sigmap) - lgamma(k + 1))
      nterm <- floor(sigmap) + 1 # how many values of sump
      denom <- lfactorial(k)
      psum <- 0
      terms <- vector("numeric", nterm)
      for (i in 1:nterm) {
        terms[i] <- lchoose(k, i - 1) + k * log(sigmap - i + 1) - denom
        pm <- 2 * (i %% 2) - 1
        psum <- psum + pm * exp(terms[i])
      }
      if(k != length(p)) {
        warning("Some studies omitted")
      }
      if(sigmap > 20) {
        warning("Likely to be unreliable, check with another method")
      }
      res <- list(p = psum, conservativep = conservativep, validp = p[keep])
    }
    class(res) <- c("sump", "metap")
    res
  }

sumz <-
  function(p, weights = NULL, data = NULL, subset = NULL, na.action = na.fail)  {
    if(is.null(data)) data <- sys.frame(sys.parent())
    mf <- match.call()
    mf$data <- NULL
    mf$subset <- NULL
    mf$na.action <- NULL
    mf[[1]] <- as.name("data.frame")
    mf <- eval(mf, data)
    if(!is.null(subset)) mf <- mf[subset,]
    mf <- na.action(mf)
    p <- as.numeric(mf$p)
    weights <- mf$weights
    noweights <- is.null(weights)
    if(noweights) weights <- rep(1, length(p))
    if(length(p) != length(weights)) warning("Length of p and weights differ")
    keep <- (p > 0) & (p < 1)
    invalid <- sum(1L * keep) < 2
    if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(z = NA_real_, p = NA_real_,
                  validp = p[keep], weights = weights)
    } else {
      if(sum(1L * keep) != length(p)) {
        warning("Some studies omitted")
        omitw <- weights[!keep]
        if((sum(1L * omitw) > 0) & !noweights)
          warning("Weights omitted too")
      }
      zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep]) /
        sqrt(sum(weights[keep]^2))
      res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE),
                  validp = p[keep], weights = weights)
    }
    class(res) <- c("sumz", "metap")
    res
  }

votep <-
  function(p, alpha = 0.5) {
    alpha <- ifelse(alpha > 1, alpha / 100, alpha) # if percent
    stopifnot(alpha > 0, alpha < 1)
    keep <- (p >= 0) & (p <= 1)
    alp <- vector("numeric", 2)
    if(alpha <= 0.5) {
      alp[1] <- alpha
      alp[2] <- 1 - alpha
    } else {
      alp[2] <- alpha
      alp[1] <- 1 - alpha
    }
    invalid <- sum(1L * keep) < 2
    if(invalid) {
      warning("Must have at least two valid p values")
      res = list(p = NA_real_, pos = NA_integer_, neg = NA_integer_,
                 alpha = alpha, validp = p[keep])
    } else {
      pi <- p[keep]
      k <- length(pi)
      pos <- sum(1L * (pi < alp[1]))
      neg <- sum(1L * (pi > alp[2]))
      if(k != length(p)) {
        warning("Some studies omitted")
      }
      if((pos + neg) <= 0) {
        warning("All p values are within specified limits of alpha")
        p <- 1
      } else {
        p = binom.test(pos, pos + neg, 0.5, alternative = "greater")$p.value
      }
      res = list(p = p, pos = pos, neg = neg, alpha = alpha, validp = pi)
    }
    class(res) <- c("votep", "metap")
    res
  }




### Perform miscellaneous tasks
### Perform misc tasks
### Jeff Xia\email{jeff.xia@mcgill.ca}
### McGill University, Canada
###License: GNU GPL (>= 2)

# Limit of detection (1/5 of min for each var)
.replace.by.lod <- function(x){
  lod <- min(x[x>0], na.rm=T)/5;
  x[x==0|is.na(x)] <- lod;
  return(x);
}

ReplaceMissingByLoD <- function(int.mat){
  int.mat <- as.matrix(int.mat);
  
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  int.mat <- apply(int.mat, 2, .replace.by.lod);
  rownames(int.mat) <- rowNms;
  colnames(int.mat) <- colNms;
  return (int.mat);
}

#'Given a data with duplicates, remove duplicates
#'@description Dups is the one with duplicates
#'@param data Input data to remove duplicates
#'@param lvlOpt Set options, default is mean
#'@param quiet Set to quiet, logical, default is T
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
RemoveDuplicates <- function(data, lvlOpt="mean", quiet=T){
  
  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];
    
    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);
    
    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);
      
      # average the whole sub matrix 
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
    AddMsg(paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""));
    return(uniq.data);
  }else{
    AddMsg("All IDs are unique.");
    return(data);
  }
} 

# in public web, this is done by microservice
.perform.computing <- function(){
  dat.in <- qs::qread("dat.in.qs"); 
  dat.in$my.res <- dat.in$my.fun();
  qs::qsave(dat.in, file="dat.in.qs");    
}

#' Read data table
#' @description Function to read in a data table. First, it will try to use fread, however, it has issues with 
#' some windows 10 files. In such case, use the slower read.table method.
#' @param fileName Input filename
#' @author Jeff Xia\email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @importFrom data.table fread
#' @export

.readDataTable <- function(fileName){
  
  dat <- tryCatch(
    {
      data.table::fread(fileName, header=TRUE, check.names=FALSE, blank.lines.skip=TRUE, data.table=FALSE);
    }, error=function(e){
      print(e);
      return(.my.slowreaders(fileName));    
    }, warning=function(w){
      print(w);
      return(.my.slowreaders(fileName));
    });
  
  if(any(dim(dat) == 0)){
    dat <- .my.slowreaders(fileName);
  }
  return(dat);
}

.my.slowreaders <- function(fileName){
  print("Using slower file reader ...");
  formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
  if(formatStr == "txt"){
    dat <- try(read.table(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
  }else{ # note, read.csv is more than read.table with sep=","
    dat <- try(read.csv(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
  }  
  return(dat);
}

.get.sqlite.con <- function(sqlite.path){
  load_rsqlite();
  return(dbConnect(SQLite(), sqlite.path)); 
}

#'Transform two column text to data matrix
#'@description Transform two column input text to data matrix (single column data frame)
#'@param txtInput Input text
#'@param sep.type Indicate the seperator type for input text. Default set to "space"
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
getDataFromTextArea <- function(txtInput, sep.type="space"){
  
  lines <- unlist(strsplit(txtInput, "\r|\n|\r\n")[1]);
  if(substring(lines[1],1,1)=="#"){
    lines <- lines[-1];
  }
  
  # separated by tab 
  if(sep.type=="tab"){
    my.lists <- strsplit(lines, "\\t");
  }else{ # from any space
    my.lists <- strsplit(lines, "\\s+");
  }
  my.mat <- do.call(rbind, my.lists);
  
  if(dim(my.mat)[2] == 1){ # add 0
    my.mat <- cbind(my.mat, rep(0, nrow(my.mat)));
  }else if(dim(my.mat)[2] > 2){
    my.mat <- my.mat[,1:2];
    msg <- "More than two columns found in the list. Only first two columns will be used."
    AddErrMsg(msg);
  }
  rownames(my.mat) <- data.matrix(my.mat[,1]);
  my.mat <- my.mat[,-1, drop=F];
  return(my.mat);
}

#'Permutation
#'@description Perform permutation, options to change number of cores used
#'@param perm.num Numeric, input the number of permutations to perform
#'@param fun Dummy function
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@usage Perform.permutation(perm.num, fun)
#'@export
#'
Perform.permutation <- function(perm.num, fun){
  
  # for public server, perm.num is not always followed to make sure loop will not continue for very long time
  # before the adventure, see how long it takes for 10 permutations
  # if it is extremely slow (>60 sec) => max 20 (<0.05)
  # if it is very slow (30-60 sec) => max 100 (<0.01)
  
  start.num <- 1; 
  perm.res <- NULL;
  if(.on.public.web & perm.num > 20){
    start.time <- Sys.time();
    perm.res <- lapply(1:10, fun);
    end.time <- Sys.time();
    
    time.taken <- end.time - start.time;
    print(paste("time taken for 10 permutations: ", time.taken));
    
    if(time.taken > 60){
      perm.num <- 20;
    }else if(time.taken > 30){
      perm.num <- 100;
    }
    start.num <- 11;
  }
  print(paste("performing", perm.num, "permutations ..."));
  perm.res <- c(perm.res, lapply(start.num:perm.num, fun));
  return(list(perm.res=perm.res, perm.num = perm.num));
}

#'Unzip .zip files
#'@description Unzips uploaded .zip files, removes the uploaded file, checks for success
#'@param inPath Input the path of the zipped files
#'@param outPath Input the path to directory where the unzipped files will be deposited
#'@param rmFile Logical, input whether or not to remove files. Default set to T
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
UnzipUploadedFile<-function(inPath, outPath, rmFile=T){
  sys.cmd <- paste("unzip",  "-o", inPath, "-d", outPath);
  #print(sys.cmd);
  sys.res <- try(system(sys.cmd));
  if(sys.res == 0){ # success code for system call
    return (1);
  }else{  # success code for system call
    print(sys.res);
    r.res <- unzip(inPath, exdir=outPath);
    return(length(r.res)>0);
  }
}

#'Perform data cleaning
#'@description Cleans data and removes -Inf, Inf, NA, negative and 0s.
#'@param bdata Input data to clean
#'@param removeNA Logical, T to remove NAs, F to not. 
#'@param removeNeg Logical, T to remove negative numbers, F to not. 
#'@param removeConst Logical, T to remove samples/features with 0s, F to not. 
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'

CleanData <-function(bdata, removeNA=T, removeNeg=T, removeConst=T){
  
  if(sum(bdata==Inf, na.rm=TRUE)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- max(bdata, na.rm=T)*2
  }
  if(sum(bdata==-Inf, na.rm=TRUE)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- min(bdata, na.rm=T)/2
  }
  if(removeNA){
    if(sum(is.na(bdata))>0){
      bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
    }
  }
  if(removeNeg){
    if(sum(as.numeric(bdata<=0)) > 0){
      inx <- bdata <= 0;
      bdata[inx] <- NA;
      bdata[inx] <- min(bdata, na.rm=T)/2
    }
  }
  if(removeConst){
    varCol <- apply(data.frame(bdata), 2, var, na.rm=T); # getting an error of dim(X) must have a positive length, fixed by data.frame
    constCol <- (varCol == 0 | is.na(varCol));
    constNum <- sum(constCol, na.rm=T);
    if(constNum > 0){
      bdata <- data.frame(bdata[,!constCol, drop=FALSE], check.names = F); # got an error of incorrect number of dimensions, added drop=FALSE to avoid vector conversion
    }
  }
  bdata;
}

#'Replace infinite numbers
#'@description Replace -Inf, Inf to 99999 and -99999
#'@param bdata Input matrix to clean numbers
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'
CleanNumber <-function(bdata){
  if(sum(bdata==Inf)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- 999999;
  }
  if(sum(bdata==-Inf)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- -999999;
  }
  bdata;
}

# only keep alphabets, numbers, "." "_", "-" and @
CleanNames <- function(query){
  query <- gsub("[^[:alnum:].@_-]", "", query);
  return(make.unique(query));
}

#'Remove spaces
#'@description Remove from, within, leading and trailing spaces
#'@param query Input the query to clear
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)

ClearStrings<-function(query){
  # kill multiple white space
  query <- gsub(" +"," ",query);
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  return(query);
}

# Remove HTML tag
PrepareLatex <- function(stringVec){
  stringVec <- gsub("<(.|\n)*?>","",stringVec);
  stringVec <- gsub("%", "\\\\%", stringVec);
  stringVec;
}

#'Determine value label for plotting
#'@description Concentration or intensity data type
#'@param data.type Input concentration or intensity data
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
GetAbundanceLabel<-function(data.type){
  if(data.type=="conc"){
    return("Concentration");
  }else {
    return("Intensity");
  }
}

#'Determine variable label for plotting
#'@description Determine data type, binned spectra, nmr peak, or ms peak
#'@param data.type Input the data type
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
GetVariableLabel<-function(data.type){
  if(data.type=="conc"){
    return("Compounds");
  }else if(data.type=="specbin"){
    return("Spectra Bins");
  }else if(data.type=="nmrpeak"){
    return("Peaks (ppm)");
  }else if(data.type=="mspeak"){
    return("Peaks (mass)");
  }else{
    return("Peaks(mz/rt)");
  }
}

#'Create Latex table
#'@description generate Latex table
#'@param mat Input matrix
#'@param method Input method to create table
#'@param data.type Input the data type
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

GetSigTable<-function(mat, method, data.type){
  if(!isEmptyMatrix(mat)){ # test if empty
    cap<-"Important features identified by";
    if(nrow(mat)>50){
      smat<-as.matrix(mat[1:50,]); # only print top 50 if too many
      colnames(smat)<-colnames(mat); # make sure column names are also copied
      mat<-smat;
      cap<-"Top 50 features identified by";
    }
    # change the rowname to first column
    col1<-rownames(mat);
    cname<-colnames(mat);
    cname<-c(GetVariableLabel(data.type), cname);
    mat<-cbind(col1, mat);
    rownames(mat)<-NULL;
    colnames(mat)<-cname;
    print(xtable::xtable(mat, caption=paste(cap, method)), caption.placement="top", size="\\scriptsize");
  }else{
    print(paste("No significant features were found using the given threshold for", method));
  }
}

#'Sig table matrix is empty
#'@description Test if a sig table matrix is empty
#'@param mat Matrix to test if empty
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)

isEmptyMatrix <- function(mat){
  if(is.null(mat) | length(mat)==0){
    return(TRUE);
  }
  if(nrow(mat)==0 | ncol(mat)==0){
    return(TRUE);
  }
  if(is.na(mat[1,1])){
    return(TRUE);
  }
  return(FALSE);
}

# from color names or code 1, 2, 3 the rbg strings
my.col2rgb <- function(cols){
  rgbcols <- col2rgb(cols);
  return(apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")}));
}

my.col2rgba <- function(cols, alpha){
  rgbcols <- col2rgb(cols);
  rgbcols <- rbind(rgbcols, alpha);
  return(as.vector(apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ")", sep="")})));
}

# #FFFFFF to rgb(1, 0, 0)
hex2rgba <- function(cols, alpha=0.8){
  my.cols <- apply(sapply(cols, col2rgb), 2, function(x){paste("rgba(", x[1], ",", x[2], ",", x[3], ",",  alpha, ")", sep="")});
  return(as.vector(my.cols));
}

hex2rgb <- function(cols){
  return(apply(sapply(cols, col2rgb), 2, function(x){paste("rgb(", x[1], ",", x[2], ",", x[3], ")", sep="")}));
}

# List of objects
# Improved list of objects
# Jeff Xia\email{jeff.xia@mcgill.ca}
# McGill University, Canada
# License: GNU GPL (>= 2)

.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  mSetObj <- .get.mSet(mSetObj);
  print(lapply(mSetObj$dataSet, object.size));
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

#'Extend axis
#'@description Extends the axis range to both ends
#'vec is the values for that axis
#'unit is the width to extend, 10 will increase by 1/10 of the range
#'@param vec Input the vector
#'@param unit Numeric
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
GetExtendRange<-function(vec, unit=10){
  var.max <- max(vec, na.rm=T);
  var.min <- min(vec, na.rm=T);
  exts <- (var.max - var.min)/unit;
  c(var.min-exts, var.max+exts);
}

getVennCounts <- function(x, include="both") {
  x <- as.matrix(x)
  include <- match.arg(include,c("both","up","down"))
  x <- sign(switch(include,
                   both = abs(x),
                   up = x > 0,
                   down = x < 0
  ))
  nprobes <- nrow(x)
  ncontrasts <- ncol(x)
  names <- colnames(x)
  if(is.null(names)) names <- paste("Group",1:ncontrasts)
  noutcomes <- 2^ncontrasts
  outcomes <- matrix(0,noutcomes,ncontrasts)
  colnames(outcomes) <- names
  for (j in 1:ncontrasts)
    outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
  xlist <- list()
  for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
  counts <- as.vector(table(xlist))
  structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

# Perform utilities for MetPa
# borrowed from Hmisc
# Jeff Xia\email{jeff.xia@mcgill.ca}
# McGill University, Canada
# License: GNU GPL (>= 2)
all.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")){
  what <- match.arg(what)
  old <- options(warn = -1)
  on.exit(options(old));
  x <- sub("[[:space:]]+$", "", x);
  x <- sub("^[[:space:]]+", "", x);
  inx <- x %in% c("", extras);
  xs <- x[!inx];
  isnum <- !any(is.na(as.numeric(xs)))
  if (what == "test") 
    isnum
  else if (isnum) 
    as.numeric(x)
  else x
}

ClearNumerics <-function(dat.mat){
  dat.mat[is.na(dat.mat)] <- -777;
  dat.mat[dat.mat == Inf] <- -999;
  dat.mat[dat.mat == -Inf] <- -111;
  dat.mat;
}

#'Calculate Pairwise Differences
#'@description Mat are log normalized, diff will be ratio. Used in higher functions. 
#'@param mat Input matrix of data to calculate pair-wise differences.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'
CalculatePairwiseDiff <- function(mat){
  f <- function(i, mat) {
    z <- mat[, i-1] - mat[, i:ncol(mat), drop = FALSE]
    colnames(z) <- paste(colnames(mat)[i-1], colnames(z), sep = "/")
    z
  }
  res <- do.call("cbind", sapply(2:ncol(mat), f, mat));
  round(res,5);
}

##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

#'Update graph settings
#'@description Function to update the graph settings.
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
# col.vec should already been created
UpdateGraphSettings <- function(mSetObj=NA, colVec, shapeVec){
  mSetObj <- .get.mSet(mSetObj);
  grpnms <- levels(current.cls);
  
  # default styles
  grp.num <- length(grpnms);
  if(grp.num <= 18){ 
    cols <- pal_18[1:grp.num];
  }else{
    cols <- colorRampPalette(pal_18)(grp.num);
  }
  
  # make sure the NA
  na.inx <- colVec == "#NA";
  colVec[na.inx] <- cols[na.inx];
  shapeVec[shapeVec == 0] <- 21;
  
  names(colVec) <- grpnms;
  names(shapeVec) <- grpnms;
  
  colVec <<- colVec;
  shapeVec <<- shapeVec;
  return(.set.mSet(mSetObj));
}

GetShapeSchema <- function(mSetObj=NA, show.name, grey.scale){
  mSetObj <- .get.mSet(mSetObj);
  if(exists("shapeVec") && all(shapeVec >= 0)){
    sps <- rep(0, length=length(mSetObj$dataSet$cls));
    clsVec <- as.character(mSetObj$dataSet$cls)
    grpnms <- names(shapeVec);
    for(i in 1:length(grpnms)){
      sps[clsVec == grpnms[i]] <- shapeVec[i];
    }
    shapes <- sps;
  }else{
    if(show.name | grey.scale){
      shapes <- as.numeric(mSetObj$dataSet$cls)+1;
    }else{
      shapes <- rep(21, length(mSetObj$dataSet$cls));
    }
  }
  return(shapes);
}

pal_18 <- c("#e6194B", "#3cb44b", "#4363d8", "#42d4f4", "#f032e6", "#ffe119", "#911eb4", "#f58231", "#bfef45",
            "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075");
cb_pal_18 <- c("#E69F00", "#b12c6f", "#56B4E9", "#009E73", "#F0E442", "#004488", 
               "#D55E00", "#EE6677", "#CCBB44", "#A95AA1", "#DCB69F", "#661100", 
               "#63ACBE", "#332288", "#EE7733", "#EE3377", "#0072B2", "#999933");

# return a gradient color vec based on value 
GetRGBColorGradient <- function(vals){
  library(RColorBrewer);
  #seed.cols <- brewer.pal(3, "YlOrRd");
  #seed.cols <- brewer.pal(9, "Oranges")[c(2,5,7)]
  seed.cols <- c("#FCF5DF", "#FFEDA0", "#F03B20")
  cols <- colorRampPalette(seed.cols)(length(vals));
  
  # set alpha for 
  my.alpha <- signif(seq(from=0.3, to=0.8, length.out=length(vals)),2);
  rgb.cols <- my.col2rgba(cols, alpha=my.alpha);
  
  # now need make sure values and colors are matched using names
  nms.orig <- names(vals);
  names(rgb.cols) <- names(sort(vals));
  ord.cols <- rgb.cols[nms.orig];
  return(as.vector(ord.cols)); # note remove names
}


GetSizeGradient <- function(vals){
  
  my.sizes <- round(seq(from=4, to=6, length.out=length(vals)));
  
  # now need make sure values and colors are matched using names
  nms.orig <- names(vals);
  names(my.sizes) <- names(sort(vals));
  ord.sizes <- my.sizes[nms.orig];
  return(as.vector(ord.sizes)); # note remove names
}

GetColorSchema <- function(my.cls, grayscale=F){
  
  lvs <- levels(my.cls); 
  grp.num <- length(lvs);
  if(grayscale){
    dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num);
  }else if(exists("colVec") && !any(colVec =="#NA") && length(colVec) == length(levels(my.cls))){
    dist.cols <- colVec;
  }else{             
    if(grp.num <= 18){ # update color and respect default
      dist.cols <- pal_18[1:grp.num];
    }else{
      dist.cols <- colorRampPalette(pal_18)(grp.num);
    }
  }
  
  colors <- vector(mode="character", length=length(my.cls));
  for(i in 1:length(lvs)){
    colors[my.cls == lvs[i]] <- dist.cols[i];
  }
  return (colors);
}

#'Remove folder
#'@description Remove folder
#'@param folderName Input name of folder to remove
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'
RemoveFolder<-function(folderName){
  a<-system(paste("rm",  "-r", folderName), intern=T);
  if(!length(a)>0){
    AddErrMsg(paste("Could not remove file -", folderName));
    return (0);
  }
  return(1);
}

#'Remove file
#'@description Remove file
#'@param fileName Input name of file to remove
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
RemoveFile<-function(fileName){
  if(file.exists(fileName)){
    file.remove(fileName);
  }
}

# do memory cleaning after removing many objects
cleanMem <- function(n=10) { for (i in 1:n) gc() }

#'Retrieve last command from the Rhistory.R file
#'@description Fetches the last command from the Rhistory.R file
#'@param regexp Retrieve last command from Rhistory file
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
GetCMD<-function(regexp){
  # store all lines into a list object
  all.lines<-readLines("Rhistory.R");
  
  all.matches<-grep(regexp, all.lines, value=T);
  if(length(all.matches)==0){
    return(NULL);
  }else{
    # only return the last command
    return(all.matches[length(all.matches)]);
  }
}

# Memory functions
ShowMemoryUse <- function(..., n=40) {
  library(pryr);
  sink(); # make sure print to screen
  print(mem_used());
  print(sessionInfo());
  print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
  print(warnings());
}

#'Perform utilities for cropping images
#'@description Obtain the full path to convert (from imagemagik)
#'for cropping images
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
GetConvertFullPath<-function(){
  path <- system("which convert", intern=TRUE);
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find convert in the PATH!");
    return("NA");
  }
  return(path);
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}

#'Converts xset object from XCMS to mSet object for MetaboAnalyst
#'@description This function converts processed raw LC/MS data from XCMS 
#'to a usable data object (mSet) for MetaboAnalyst. The immediate next step following using this 
#'function is to perform a SanityCheck, and then further data processing and analysis can continue.
#'@param xset The name of the xcmsSet object created.
#'@param dataType The type of data, either list (Compound lists), conc (Compound concentration data), 
#'specbin (Binned spectra data), pktable (Peak intensity table), nmrpeak (NMR peak lists), mspeak (MS peak lists), 
#'or msspec (MS spectra data).
#'@param analType Indicate the analysis module to be performed: stat, pathora, pathqea, msetora, msetssp, msetqea, ts, 
#'cmpdmap, smpmap, or pathinteg.
#'@param paired Logical, is data paired (T) or not (F).
#'@param format Specify if samples are paired and in rows (rowp), unpaired and in rows (rowu),
#'in columns and paired (colp), or in columns and unpaired (colu).
#'@param lbl.type Specify the data label type, either discrete (disc) or continuous (cont).
#'@export

XSet2MSet <- function(xset, dataType, analType, paired=F, format, lbl.type){
  
  # data <- xcms::groupval(xset, "medret", "into");
  # data2 <- rbind(class= as.character(phenoData(xset)$class), data);
  # rownames(data2) <- c("group", paste(round(groups(xset)[,"mzmed"], 3), round(groups(xset)[,"rtmed"]/60, 1), sep="/"));
  # fast.write.csv(data2, file="PeakTable.csv");
  # mSet <- InitDataObjects("dataType", "analType", paired)
  # mSet <- Read.TextData(mSet, "PeakTable.csv", "format", "lbl.type")
  # print("mSet successfully created...")
  # return(.set.mSet(mSetObj));
}

#'Get fisher p-values
#'@param numSigMembers Number of significant members
#'@param numSigAll Number of all significant features
#'@param numMembers Number of members
#'@param numAllMembers Number of all members
#'@export
GetFisherPvalue <- function(numSigMembers, numSigAll, numMembers, numAllMembers){
  z <- cbind(numSigMembers, numSigAll-numSigMembers, numMembers-numSigMembers, numAllMembers-numMembers-numSigAll+numSigMembers);
  z <- lapply(split(z, 1:nrow(z)), matrix, ncol=2);
  z <- lapply(z, fisher.test, alternative = 'greater');
  p.values <- as.numeric(unlist(lapply(z, "[[", "p.value"), use.names=FALSE));
  return(p.values);
}

saveNetworkInSIF <- function(network, name){
  edges <- .graph.sif(network=network, file=name);
  sif.nm <- paste(name, ".sif", sep="");
  if(length(list.edge.attributes(network))!=0){
    edge.nms <- .graph.eda(network=network, file=name, edgelist.names=edges);
    sif.nm <- c(sif.nm, edge.nms);
  }
  if(length(list.vertex.attributes(network))!=0){
    node.nms <- .graph.noa(network=network, file=name);
    sif.nm <- c(sif.nm, node.nms);
  }
  # need to save all sif and associated attribute files into a zip file for download
  zip(paste(name,"_sif",".zip", sep=""), sif.nm);
}

.graph.sif <- function(network, file){
  edgelist.names <- igraph::get.edgelist(network, names=TRUE)
  edgelist.names <- cbind(edgelist.names[,1], rep("pp", length(E(network))), edgelist.names[,2]);
  write.table(edgelist.names, row.names=FALSE, col.names=FALSE, file=paste(file, ".sif", sep=""), sep="\t", quote=FALSE)
  return(edgelist.names) 
}

# internal method to write cytoscape node attribute files
.graph.noa <- function(network, file){
  all.nms <- c();
  attrib <- list.vertex.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    noa <- cbind(V(network)$name, rep("=", length(V(network))), get.vertex.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    file.nm <- paste(file, "_", attrib[i], ".NA", sep="");
    write(first.line, file=file.nm, ncolumns = 1, append=FALSE, sep=" ")
    write.table(noa, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

# internal method to write cytoscape edge attribute files
.graph.eda <- function(network, file, edgelist.names){
  all.nms <- c();
  attrib <- list.edge.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.edge.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    eda <- cbind(cbind(edgelist.names[,1], rep("(pp)", length(E(network))), edgelist.names[,3]), rep("=", length(E(network))), get.edge.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="");
    file.nm <- paste(file, "_", attrib[i], ".EA", sep="");
    write(first.line, file=file.nm, ncolumns=1, append=FALSE, sep =" ")
    write.table(eda, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

PlotLoadBoxplot <- function(mSetObj=NA, cmpd){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(.on.public.web){
    load_ggplot()
  }
  
  cls.lbls <- mSetObj$dataSet$cls;
  y.label <- GetAbundanceLabel(mSetObj$dataSet$type);
  cmpd.name = paste0("Met_", cmpd, ".png")
  
  Cairo::Cairo(file=cmpd.name, width=240, height=400, bg = "transparent", type="png");
  
  col <- unique(GetColorSchema(cls.lbls))
  df <- data.frame(conc = mSetObj$dataSet$norm[, cmpd], class = cls.lbls)
  p <- ggplot2::ggplot(df, aes(x=class, y=conc, fill=class)) + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA) + theme_bw() + geom_jitter(size=1)
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
  p <- p + stat_summary(fun.y=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
  p <- p + theme(text = element_text(size=15), plot.margin = margin(t=0.45, r=0.25, b=1.5, l=0.25, "cm"))
  p <- p + scale_fill_manual(values=col) + ggtitle(cmpd) + theme(axis.text.x = element_text(angle=90, hjust=1), axis.text = element_text(size=10))
  p <- p + theme(plot.title = element_text(size = 14, hjust=0.5, face="bold", vjust=2))
  print(p)
  
  dev.off()
}

#'Compute within group and between group sum of squares
#'(BSS/WSS) for each row of a matrix which may have NA
#'@description Columns have labels, x is a numeric vector,
#'cl is consecutive integers
#'@param x Numeric vector
#'@param cl Columns
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)

Get.bwss<-function(x, cl){
  K <- max(cl) - min(cl) + 1
  tvar <- var.na(x);
  tn <- sum(!is.na(x));
  wvar <- wn <- numeric(K);
  
  for(i in (1:K)) {
    if(sum(cl == (i + min(cl) - 1)) == 1){
      wvar[i] <- 0;
      wn[i] <- 1;
    }
    
    if(sum(cl == (i + min(cl) - 1)) > 1) {
      wvar[i] <- var.na(x[cl == (i + min(cl) - 1)]);
      wn[i] <- sum(!is.na(x[cl == (i + min(cl) - 1)]));
    }
  }
  
  WSS <- sum.na(wvar * (wn - 1));
  TSS <- tvar * (tn - 1)
  (TSS - WSS)/WSS;
}

# compute the distance to the centroid of the given data
# note each col is a x,y,z
# since this is centered, the centroid is origin
GetDist3D <-function(mat, target=c(0,0,0)){
  dist.vec <- apply(mat, 2, function(x) dist(rbind(x, target)));
  return(dist.vec);
}

sum.na <- function(x,...){
  res <- NA
  tmp <- !(is.na(x) | is.infinite(x))
  if(sum(tmp) > 0)
    res <- sum(x[tmp])
  res
}

var.na <- function(x){
  res <- NA
  tmp <- !(is.na(x) | is.infinite(x))
  if(sum(tmp) > 1){
    res <- var(as.numeric(x[tmp]))
  }
  res
}

end.with <- function(bigTxt, endTxt){
  return(substr(bigTxt, nchar(bigTxt)-nchar(endTxt)+1, nchar(bigTxt)) == endTxt);
}

## fast T-tests/F-tests, and use cache to avoid redundant computing
PerformFastUnivTests <- function(data, cls, var.equal=TRUE){
  if(!exists("mem.univ")){
    require("memoise");
    mem.univ <<- memoise(.perform.fast.univ.tests);
  }
  return(mem.univ(data, cls, var.equal));
}

.perform.fast.univ.tests <- function(data, cls, var.equal=TRUE){
  
  print("Performing fast univariate tests ....");
  # note, feature in rows for gene expression
  data <- t(as.matrix(data));
  if(length(levels(cls)) > 2){
    res <- try(rowcolFt(data, cls, var.equal = var.equal));
  }else{
    res <- try(rowcoltt(data, cls, FALSE, 1L, FALSE));
  }  
  
  if(class(res) == "try-error") {
    res <- cbind(NA, NA);
  }else{
    # res <- cbind(res$statistic, res$p.value);
    # make sure row names are kept
    res <- res[, c("statistic", "p.value")];
  }
  
  return(res);
}

GetCurrentPathForScheduler <- function(){
  path <-getwd();
  path <- gsub(basename(path), "", path)
  return(path)
}

fast.write.csv <- function(dat, file, row.names=TRUE){
  tryCatch(
    {
      if(is.data.frame(dat)){
        # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
        data.table::fwrite(dat, file, row.names=row.names);
      }else{
        write.csv(dat, file, row.names=row.names);  
      }
    }, error=function(e){
      print(e);
      write.csv(dat, file, row.names=row.names);   
    }, warning=function(w){
      print(w);
      write.csv(dat, file, row.names=row.names); 
    });
}

rowcolFt =  function(x, fac, var.equal, which = 1L) {
  
  if(!(which %in% c(1L, 2L)))
    stop(sQuote("which"), " must be 1L or 2L.")
  
  if(which==2L)
    x = t(x)
  
  if (typeof(x) == "integer")
    x[] <- as.numeric(x)
  
  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]
  
  ## Number of levels (groups)
  k <- nlevels(fac)
  
  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- matrix(
    sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
    nrow = nrow(x),
    ncol = nlevels(fac))
  
  ## x1: a matrix of group means, with as many rows as x, columns correspond to groups 
  x1 <- xm[,fac, drop=FALSE]
  
  ## degree of freedom 1
  dff    <- k - 1
  
  if(var.equal){
    ## x0: a matrix of same size as x with overall means
    x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
    
    ## degree of freedom 2
    dfr    <- ncol(x) - dff - 1
    
    ## mean sum of squares
    mssf   <- rowSums(sqr(x1 - x0)) / dff
    mssr   <- rowSums(sqr( x - x1)) / dfr
    
    ## F statistic
    fstat  <- mssf/mssr
    
  } else{
    
    ## a nrow(x) x nlevels(fac) matrix with the group size  of each factor
    ## level
    ni <- t(matrix(tapply(fac,fac,length),ncol=nrow(x),nrow=k))
    
    ## wi: a nrow(x) x nlevels(fac) matrix with the variance * group size of each factor
    ## level
    sss <- sqr(x-x1)
    x5 <- matrix(
      sapply(levels(fac), function(fl) rowSums(sss[,which(fac==fl), drop=FALSE])),
      nrow = nrow(sss),
      ncol = nlevels(fac))          
    wi <- ni*(ni-1) /x5
    
    ## u : Sum of wi
    u  <- rowSums(wi)
    
    ## F statistic
    MR <- rowSums(sqr((1 - wi/u)) * 1/(ni-1))*1/(sqr(k)-1)
    fsno <- 1/dff * rowSums(sqr(xm - rowSums(wi*xm)/u) * wi)
    fsdeno <- 1+ 2* (k-2)*MR
    fstat <- fsno/fsdeno
    
    ## degree of freedom 2: Vector with length nrow(x)
    dfr <- 1/(3 * MR)
    
  }
  
  res = data.frame(statistic = fstat,
                   p.value   = pf(fstat, dff, dfr, lower.tail=FALSE),
                   row.names = rownames(x))
  
  attr(res, "df") = c(dff=dff, dfr=dfr)
  return(res)
}

rowcoltt =  function(x, fac, tstatOnly, which, na.rm) {
  
  #if(.on.public.web){
  #  dyn.load(.getDynLoadPath());
  #}
  
  if (!missing(tstatOnly) && (!is.logical(tstatOnly) || is.na(tstatOnly)))
    stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
  
  f = checkfac(fac)
  if ((f$nrgrp > 2) || (f$nrgrp <= 0))
    stop("Number of groups is ", f$nrgrp, ", but must be >0 and <=2 for 'rowttests'.")
  
  if (typeof(x) == "integer")
    x[] <- as.numeric(x)
  
  #cc = .Call("rowcolttests", x, f$fac, f$nrgrp, which-1L, na.rm)
  cc = XiaLabCppLib::rowcolttestsR(x, f$fac, f$nrgrp, which-1L, na.rm)
  
  res = data.frame(statistic = cc$statistic,
                   dm        = cc$dm,
                   row.names = dimnames(x)[[which]])
  
  if (!tstatOnly)
    res = cbind(res, p.value = 2*pt(abs(res$statistic), cc$df, lower.tail=FALSE))
  
  attr(res, "df") = cc$df    
  return(res)
}

checkfac = function(fac) {
  
  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac, na.rm=TRUE)+1)
    fac   = as.integer(fac)
  }
  ## this must precede the factor test
  if(is.character(fac))
    fac = factor(fac)
  
  if (is.factor(fac)) {
    nrgrp = nlevels(fac)
    fac   = as.integer(as.integer(fac)-1)
  } 
  if(!is.integer(fac))
    stop("'fac' must be factor, character, numeric, or integer.")
  
  if(any(fac<0, na.rm=TRUE))
    stop("'fac' must not be negative.")
  
  return(list(fac=fac, nrgrp=nrgrp))
}

# to convert all rds files to qs file for faster access
convert.rds2qs <- function(){
  rds.files <- list.files(".", pattern=".rds$");
  for(rds in rds.files){
    lib.rds <- readRDS(rds);
    nm <- substr(rds, 1, nchar(rds)-4);
    qs::qsave(lib.rds,paste0(nm, ".qs"));
  }
}

# to convert all rda files to qs file for faster access
convert.rda2qs <- function(){
  rda.files <- list.files(".", pattern=".rda$");
  for(rda in rda.files){
    nm <- substr(rda, 1, nchar(rda)-4);
    lib.rda <- load(rda);
    # here the name is inmexpa (jointpa)
    # qs::qsave(inmexpa,paste0(nm, ".qs"));
    # here the name is metpa (metpa)
    qs::qsave(metpa,paste0(nm, ".qs"));
    # here the name is current.msetlib (msets)
    #qs::qsave(current.msetlib,paste0(nm, ".qs"));
  }
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

# a utility function to get pheatmap image size (before saving to PNG)
# https://stackoverflow.com/questions/61874876/get-size-of-plot-in-pixels-in-r
get_pheatmap_dims <- function(dat, cellheight = 15, cellwidth = 15){
  png("NUL"); # trick to avoid open device in server 
  heat_map <- pheatmap::pheatmap(dat, cellheight = cellheight, cellwidth = cellwidth);
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"));
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"));
  dev.off();
  return(list(height = plot_height, width = plot_width));
}




#' my.fgsea
#'
#' @param mSetObj mSetObj
#' @param pathways pathways
#' @param stats stats
#' @param ranks ranks
#' @param nperm nperm
#' @param minSize minSize
#' @param maxSize maxSize
#' @param nproc nproc
#' @param gseaParam gseaParam
#' @param BPPARAM BPPARAM
#' @noRd
#' @importFrom data.table data.table rbindlist
#' @return

my.fgsea <- function(mSetObj, pathways, stats, ranks,
                     nperm,
                     minSize=1, maxSize=Inf,
                     nproc=0,
                     gseaParam=1,
                     BPPARAM=NULL) {
  
  # Warning message for ties in stats
  ties <- sum(duplicated(stats[stats != 0]))
  if (ties != 0) {
    warning("There are ties in the preranked stats (",
            paste(round(ties * 100 / length(stats), digits = 2)),
            "% of the list).\n",
            "The order of those tied m/z features will be arbitrary, which may produce unexpected results.")
  }
  
  # Warning message for duplicate gene names
  if (any(duplicated(names(stats)))) {
    warning("There are duplicate m/z feature names, fgsea may produce unexpected results")
  }
  
  granularity <- 1000
  permPerProc <- rep(granularity, floor(nperm / granularity))
  if (nperm - sum(permPerProc) > 0) {
    permPerProc <- c(permPerProc, nperm - sum(permPerProc))
  }
  set.seed(123)
  seeds <- sample.int(10^9, length(permPerProc))
  
  if (is.null(BPPARAM)) {
    if (nproc != 0) {
      if (.Platform$OS.type == "windows") {
        # windows doesn't support multicore, using snow instead
        BPPARAM <- BiocParallel::SnowParam(workers = nproc)
      } else {
        BPPARAM <- BiocParallel::MulticoreParam(workers = nproc)
      }
    } else {
      BPPARAM <- BiocParallel::bpparam()
    }
  }
  
  minSize <- max(minSize, 1)
  stats <- abs(stats) ^ gseaParam
  
  # returns list of indexs of matches between pathways and rank names
  pathwaysPos <- lapply(pathways, function(p) { as.vector(na.omit(fmatch_r(p, names(ranks)))) })
  pathwaysFiltered <- lapply(pathwaysPos, function(s) { ranks[s] })
  qs::qsave(pathwaysFiltered, "pathwaysFiltered.qs")
  
  # adjust for the fact that a single m/z feature can match to several compound identifiers (not when in EC space)
  # subsets m/z features responsible for a compound and matches it to total set of matched m/z features
  # returns the length
  
  matched_res <- qs::qread("mum_res.qs");
  
  if(mSetObj$paramSet$mumRT){
    pathwaysSizes <- sapply(pathwaysFiltered, length)
  }else{
    pathway2mzSizes <- sapply(pathways, function(z) { length(intersect(as.numeric(unique(unlist(mSetObj$cpd2mz_dict[z]))), unique(matched_res[,1])))} )
    oldpathwaysSizes <- sapply(pathwaysFiltered, length)
    pathwaysSizes <- pmin(pathway2mzSizes, oldpathwaysSizes)
  }
  
  toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
  m <- length(toKeep)
  
  if (m == 0) {
    return(data.table::data.table(pathway=character(),
                                  pval=numeric(),
                                  padj=numeric(),
                                  ES=numeric(),
                                  NES=numeric(),
                                  nMoreExtreme=numeric(),
                                  size=integer(),
                                  leadingEdge=list()))
  }
  
  pathwaysFiltered <- pathwaysFiltered[toKeep]
  pathwaysSizes <- pathwaysSizes[toKeep]
  
  K <- max(pathwaysSizes)
  
  #perform gsea
  gseaStatRes <- do.call(rbind,
                         lapply(pathwaysFiltered, fgsea::calcGseaStat,
                                stats=stats,
                                returnLeadingEdge=TRUE))
  
  leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
  pathwayScores <- unlist(gseaStatRes[, "res"])
  
  #perform permutations
  universe <- seq_along(stats)
  set.seed(123)
  counts <- BiocParallel::bplapply(seq_along(permPerProc), function(i) {
    nperm1 <- permPerProc[i]
    leEs <- rep(0, m)
    geEs <- rep(0, m)
    leZero <- rep(0, m)
    geZero <- rep(0, m)
    leZeroSum <- rep(0, m)
    geZeroSum <- rep(0, m)
    if (m == 1) {
      for (i in seq_len(nperm1)) {
        randSample <- sample.int(length(universe), K)
        randEsP <- fgsea::calcGseaStat(
          stats = stats,
          selectedStats = randSample,
          gseaParam = 1)
        leEs <- leEs + (randEsP <= pathwayScores)
        geEs <- geEs + (randEsP >= pathwayScores)
        leZero <- leZero + (randEsP <= 0)
        geZero <- geZero + (randEsP >= 0)
        leZeroSum <- leZeroSum + pmin(randEsP, 0)
        geZeroSum <- geZeroSum + pmax(randEsP, 0)
      }
    } else {
      if (packageVersion("fgsea") > "1.12.0"){
        aux <- fgsea:::calcGseaStatCumulativeBatch(
          stats = stats,
          gseaParam = 1,
          pathwayScores = pathwayScores,
          pathwaysSizes = pathwaysSizes,
          iterations = nperm1,
          seed = seeds[i],
          scoreType = "std")} else {
            aux <- fgsea:::calcGseaStatCumulativeBatch(
              stats = stats,
              gseaParam = 1,
              pathwayScores = pathwayScores,
              pathwaysSizes = pathwaysSizes,
              iterations = nperm1,
              seed = seeds[i])
          }
      leEs = get("leEs", aux)
      geEs = get("geEs", aux)
      leZero = get("leZero", aux)
      geZero = get("geZero", aux)
      leZeroSum = get("leZeroSum", aux)
      geZeroSum = get("geZeroSum", aux)
    }
    data.table::data.table(pathway=seq_len(m),
                           leEs=leEs, geEs=geEs,
                           leZero=leZero, geZero=geZero,
                           leZeroSum=leZeroSum, geZeroSum=geZeroSum
    )
  }, BPPARAM=BPPARAM)
  
  counts <- data.table::rbindlist(counts)
  
  # Getting rid of check NOTEs
  leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
  pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
  nMoreExtreme=nGeEs=nLeEs=size=NULL
  leadingEdge=NULL
  .="damn notes"
  
  pval <- unlist(lapply(counts$pathway, function(c) min((1+sum(counts[c,]$leEs)) / (1 + sum(counts[c,]$leZero)),
                                                        (1+sum(counts[c,]$geEs)) / (1 + sum(counts[c,]$geZero)))))
  
  leZeroMean <- unlist(lapply(counts$pathway, function(d) sum(counts[d,]$leZeroSum) / sum(counts[d,]$leZero)))
  geZeroMean <- unlist(lapply(counts$pathway, function(e) sum(counts[e,]$geZeroSum) / sum(counts[e,]$geZero)))
  nLeEs <- unlist(lapply(counts$pathway, function(f) sum(counts[f,]$leEs)))
  nGeEs <- unlist(lapply(counts$pathway, function(g) sum(counts[g,]$geEs)))
  
  pvals <- data.frame(pval=pval, leZeroMean=leZeroMean, geZeroMean=geZeroMean, nLeEs=nLeEs, nGeEs=nGeEs)
  
  padj <- p.adjust(pvals$pval, method="fdr")
  ES <- pathwayScores
  NES <- ES / ifelse(ES > 0, pvals$geZeroMean, abs(pvals$leZeroMean))
  pvals$leZeroMean <- NULL
  pvals$geZeroMean <- NULL
  
  nMoreExtreme <- ifelse(ES > 0, pvals$nGeEs, pvals$nLeEs)
  pvals$nLeEs <- NULL
  pvals$nGeEs <- NULL
  
  size <- pathwaysSizes
  pathway <- names(pathwaysFiltered)
  
  leadingEdge <- sapply(leadingEdges, paste0, collapse = "; ")
  leadingEdge2 <- sapply(leadingEdge, function(x) strsplit(x, "; "))
  pathway.cpds <- sapply(pathwaysFiltered, attributes)
  
  matches <- mapply(intersect, leadingEdge2, pathway.cpds)
  
  leadingEdgeMatched <- sapply(matches, paste0, collapse = "; ")
  
  pvals.done <- cbind(pathway, pvals, padj, ES, NES, nMoreExtreme, size, leadingEdgeMatched)
  
  return(pvals.done)
}


#' MetaboAnalystR: A package for computating the notorious bar statistic.
#'
#' The MetaboAnalystR package provides a pipeline for metabolomics processing.
#' 
#' @section MetaboAnalystR functions:
#' The MetaboAnalystR functions ...
#'
#' @docType package
#' @name MetaboAnalystR
#' @useDynLib MetaboAnalystR, .registration=TRUE
NULL
#> NULL


######### ======= ------------- C function Bin: Batch Effect Module ----------- ======== #########

#' Internal C fucntion - C_imodwt_r
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.
C_imodwt_r <- function(y,z,N,j, L, ht, gt, XX){
  if (.on.public.web){ .C("imodwt", y, z, N, j, L, ht, gt, out=XX)$out} else{
    .C("imodwt", y, z, N, j, L, ht, gt, out=XX, PACKAGE = "MetaboAnalystR")$out}
}

#' Internal C fucntion - C_modwt_r
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.
C_modwt_r <- function(X,N,j, L, ht, gt,W, V){
  if (.on.public.web){ .C("modwt", X, N, as.integer(j), L,ht, gt, W = W, V = V)[7:8]} else {
    .C("modwt", X, N, as.integer(j), L, ht, gt, W = W, V = V, PACKAGE = "MetaboAnalystR")[7:8]}
}



############# ============ ------------- Bin bottom ----------- ============ ###########

#' Internal C fucntion - fmatch
#' @noRd
#' @references https://cran.r-project.org/web/packages/fastmatch/index.html
fmatch_r <- function(x, table, nomatch = NA_integer_, incomparables = NULL){
  .Call("fmatch", x, table, nomatch, incomparables, FALSE, PACKAGE = "MetaboAnalystR")
}
























CleanLipidNames <- function(qvec){
  
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.clean.lipid")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_lipid.Rc");    
    }
    return(my.clean.lipid(qvec));
  }else{
    return(my.clean.lipid(qvec));
  }
  
}

#' Perform detailed name match
#'@description Given a query, perform compound matching. 
#'@param mSetObj Input name of the created mSet Object.
#'@param q Input the query.
#'@export
#'
PerformDetailMatch <- function(mSetObj=NA, q){
  
  mSetObj <- .get.mSet(mSetObj);
  
  lipid = mSetObj$lipid.feats
  
  if(mSetObj$dataSet$q.type == "name"){
    PerformApproxMatch(mSetObj, q, lipid);
  }else{
    PerformMultiMatch(mSetObj, q, lipid);
  }
}

#' Perform multiple name matches
#'@description Given a query, performs compound name matching. 
#'@param mSetObj Input name of the created mSet Object.
#'@param q Input the query.
#'@export
#'
PerformMultiMatch <- function(mSetObj=NA, q, lipid){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
    cmpd.db <- .get.my.lib("lipid_compound_db.qs");
  }else if(anal.type == "utils"){
    cmpd.db <- .get.my.lib("master_compound_db.qs");
  }else{
    cmpd.db <- .get.my.lib("compound_db.qs");
  }
  
  matched.inx <- which(cmpd.db$kegg %in% q);
  if(length(matched.inx) > 0) {
    # record all the candidates,
    candidates <- cbind(matched.inx, cmpd.db$name[matched.inx]);
    mSetObj$dataSet$candidates <- candidates;
  }else{
    mSetObj$dataSet$candidates <- NULL;
  }
  return(.set.mSet(mSetObj));
}

#'Perform approximate compound matches
#'@description Given a query, perform approximate compound matching 
#'@param mSetObj Input the name of the created mSetObj.
#'@param q Input the q vector.
#'@export
#'
PerformApproxMatch <- function(mSetObj=NA, q, lipid){
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.approx.match")){ # public web on same user dir
      compiler::loadcmp("../../rscripts/metaboanalystr/_util_approx.Rc");    
    }
    return(my.approx.match(mSetObj, q, lipid));
  }else{
    return(my.approx.match(mSetObj, q, lipid));
  }
}

#'Set matched name based on user selection from all potential hits
#'@description Note: to change object in the enclosing enviroment, use "<<-"
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param query_nm Input the query name.
#'@param can_nm Input the candidate name.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

SetCandidate <- function(mSetObj=NA, query_nm, can_nm){
  
  mSetObj <- .get.mSet(mSetObj);
  
  lipid = mSetObj$lipid.feats
  
  query_inx <- which(mSetObj$name.map$query.vec == query_nm);
  
  can_mat <- mSetObj$dataSet$candidates;
  
  if(!is.null(can_mat)){
    
    if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
      cmpd.db <- .get.my.lib("lipid_compound_db.qs");
    }else if(anal.type == "utils"){
      cmpd.db <- .get.my.lib("master_compound_db.qs");
    }else{
      cmpd.db <- .get.my.lib("compound_db.qs");
    }
    
    can_inx <- which(can_mat[,2] == can_nm);
    
    if(can_inx <= nrow(can_mat)){
      can_inx <- which(cmpd.db$name == can_nm);
      hit <- cmpd.db[can_inx, ,drop=F];
      mSetObj$name.map$hit.inx[query_inx] <- can_inx;
      mSetObj$name.map$hit.values[query_inx] <- hit[,2];
      mSetObj$name.map$match.state[query_inx] <- 1;
      
      # re-generate the CSV file
      csv.res <- mSetObj$dataSet$map.table;
      if(ncol(csv.res) > 7){ # general utilities
        csv.res[query_inx, ]<-c(csv.res[query_inx, 1],
                                mSetObj$name.map$hit.values[query_inx],
                                hit$hmdb_id,
                                hit$pubchem_id,
                                hit$chebi_id,
                                hit$kegg_id,
                                hit$metlin_id,
                                hit$smiles,
                                1);
      }else{ # pathway analysis
        csv.res[query_inx, ]<-c(csv.res[query_inx, 1],
                                mSetObj$name.map$hit.values[query_inx],
                                hit$hmdb_id,
                                hit$pubchem_id,
                                hit$kegg_id,
                                hit$smiles,
                                1);
      }
      fast.write.csv(csv.res, file="name_map.csv", row.names=F);
      mSetObj$dataSet$map.table <- csv.res;
    }else{ #no match
      mSetObj$name.map$hit.inx[query_inx] <- 0;
      mSetObj$name.map$hit.values[query_inx] <- "";
      mSetObj$name.map$match.state[query_inx] <- 0;
      print("No name matches found.")
    }
  }
  
  if(.on.public.web){
    .set.mSet(mSetObj);
    return(query_inx);
  }else{
    return(.set.mSet(mSetObj));
  }
}

CrossReferencingAPI <- function(mSetObj=NA, inputType){
  
  mSetObj <- .get.mSet(mSetObj);
  
  toSend <- list(mSet = mSetObj, 
                 inputType = inputType,
                 analType = anal.type)
  
  load_httr()
  base <- api.base
  endpoint <- "/internal_mapcompounds"
  call <- paste(base, endpoint, sep="")
  print(call)
  
  saveRDS(toSend, "tosend.rds")
  request <- httr::POST(url = call, 
                        body = list(rds = upload_file("tosend.rds", "application/octet-stream")))
  
  # check if successful
  if(request$status_code != 200){
    AddErrMsg("Failed to connect to Xia Lab API Server!")
    return(0)
  }
  
  # now process return
  request <- httr::content(request, "raw")
  request <- unserialize(request)
  
  if(is.null(request$name.map)){
    AddErrMsg("Error! Compound name mapping via api.metaboanalyst.ca unsuccessful!")
    return(0)
  }else{
    mSetObj <- request
  }
  
  print("Compound name mapping via api.metaboanalyst.ca successful!")
  
  return(.set.mSet(mSetObj));
}

##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

#'Get all candidate compound names for a given index 
#'@description Returns 3 coloumns - inx, name, score
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

GetCandidateList <- function(mSetObj=NA, lipid){
  
  mSetObj <- .get.mSet(mSetObj);
  
  lipid = mSetObj$lipid.feats
  
  can_hits <- mSetObj$dataSet$candidates;
  
  if(is.null(can_hits)){
    can.mat <- matrix("", nrow=1, ncol= 6);
  }else{
    # construct the result table with cells wrapped in html tags
    # the unmatched will be highlighted in different background
    
    can.mat <- matrix("", nrow=nrow(can_hits)+1, ncol= 6);
    
    if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
      cmpd.db <- .get.my.lib("lipid_compound_db.qs");
    }else if(anal.type == "utils"){
      cmpd.db <- .get.my.lib("master_compound_db.qs");
    }else{
      cmpd.db <- .get.my.lib("compound_db.qs");
    }
    
    if(!lipid){
      # need to exclude lipids, to be consistent with approx matching part so that same index can be used to fetch db entries
      nonLipidInx <- cmpd.db$lipid == 0;
      cmpd.db <-cmpd.db[nonLipidInx,];
    }
    
    for (i in 1:nrow(mSetObj$dataSet$candidates)){
      hit.inx <- mSetObj$dataSet$candidates[i, 1];
      hit.name <- mSetObj$dataSet$candidates[i, 2];
      hit <- cmpd.db[hit.inx, ,drop=F];
      can.mat[i, ] <- c(hit.name,
                        paste(ifelse(hit$hmdb_id=="NA","", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")), sep=""),
                        paste(ifelse(hit$pubchem_id=="NA", "", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
                        paste(ifelse(hit$chebi_id=="NA","", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
                        paste(ifelse(hit$kegg_id=="NA","",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
                        paste(ifelse(hit$metlin_id=="NA","",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""));
    }
    # add "none" option
    can.mat[nrow(mSetObj$dataSet$candidates)+1,] <- c("None of the above", "", "", "", "", "");
  }
  # add the hit columns
  return.cols <- c(TRUE, mSetObj$return.cols);
  
  if(.on.public.web){
    return(as.vector(can.mat[,return.cols, drop=F]));
  }
  
  mSetObj$name.map$hits.candidate.list <- can.mat[,mSetObj$return.cols, drop=F]
  
  return(.set.mSet(mSetObj));
}

GetCanListRowNumber <- function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  if(is.null(mSetObj$dataSet$candidates)){
    return(1);
  }else{
    return(nrow(mSetObj$dataSet$candidates)+1); # include the "none" row
  }
}

GetQuery <- function(mSetObj=NA, inx){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$dataSet$cmpd[inx]);
}



#'Get mapping table
#'@description Return results from compound name mapping in a table
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export

GetMapTable <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  print(xtable::xtable(mSetObj$dataSet$map.table, caption="Result from Compound Name Mapping"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

#'Creates the mapping result table
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
CreateMappingResultTable <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  lipid <- mSetObj$lipid.feats
  
  if(lipid & anal.type == "msetqea"){
    qvec <- names(mSet$dataSet$url.var.nms);
  }else{
    qvec <- mSetObj$dataSet$cmpd;
  }
  
  if(is.null(qvec)){
    return();
  }
  # style for highlighted background for unmatched names
  pre.style<-NULL;
  post.style<-NULL;
  
  # style for no matches
  if(mSetObj$dataSet$q.type == "name"){
    no.prestyle<-"<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
    no.poststyle<-"</strong>";
  }else{
    no.prestyle<-"<strong style=\"background-color:red; font-size=125%; color=\"black\">";
    no.poststyle<-"</strong>";
  }
  
  hit.inx<-mSetObj$name.map$hit.inx;
  hit.values<-mSetObj$name.map$hit.values;
  match.state<-mSetObj$name.map$match.state;
  
  # construct the result table with cells wrapped in html tags
  # the unmatched will be highlighted in different background
  html.res <- matrix("", nrow=length(qvec), ncol=8);
  csv.res <- matrix("", nrow=length(qvec), ncol=9);
  colnames(csv.res) <- c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "SMILES", "Comment");
  
  if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
    cmpd.db <- .get.my.lib("lipid_compound_db.qs");
  }else if(anal.type == "utils"){
    cmpd.db <- .get.my.lib("master_compound_db.qs");
  }else{
    cmpd.db <- .get.my.lib("compound_db.qs");
  }
  
  for (i in 1:length(qvec)){
    if(match.state[i]==1){
      pre.style<-"";
      post.style="";
    }else{ # no matches
      pre.style<-no.prestyle;
      post.style<-no.poststyle;
    }
    hit <-cmpd.db[hit.inx[i], ,drop=F];
    html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                     paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$hmdb_id) || hit$hmdb_id=="" || hit$hmdb_id=="NA","-", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")),  sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$pubchem_id) || hit$pubchem_id=="" || hit$pubchem_id=="NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$chebi_id) || hit$chebi_id==""|| hit$chebi_id=="NA","-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$kegg_id) || hit$kegg_id==""|| hit$kegg_id=="NA","-",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$metlin_id) || hit$metlin_id==""|| hit$metlin_id=="NA","-",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""),
                     ifelse(match.state[i]!=1,"View",""));
    csv.res[i, ]<-c(qvec[i],
                    ifelse(match.state[i]==0, "NA", hit.values[i]),
                    ifelse(match.state[i]==0, "NA", hit$hmdb_id),
                    ifelse(match.state[i]==0, "NA", hit$pubchem_id),
                    ifelse(match.state[i]==0, "NA", hit$chebi_id),
                    ifelse(match.state[i]==0, "NA", hit$kegg_id),
                    ifelse(match.state[i]==0, "NA", hit$metlin_id),
                    ifelse(match.state[i]==0, "NA", hit$smiles),
                    match.state[i]);
  }
  # return only columns user selected
  
  # add query and match columns at the the beginning, and 'Detail' at the end
  return.cols <- c(TRUE, TRUE, mSetObj$return.cols, TRUE);
  html.res <- html.res[,return.cols, drop=F];
  csv.res <- csv.res[,return.cols, drop=F];
  
  # store the value for report
  mSetObj$dataSet$map.table <- csv.res;
  fast.write.csv(csv.res, file="name_map.csv", row.names=F);
  
  if(.on.public.web){
    .set.mSet(mSetObj);
    return(as.vector(html.res));
  }else{
    return(.set.mSet(mSetObj));
  }
}

GetHitsRowNumber<-function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(length(mSetObj$name.map$hit.inx));
}

GetPathNames<-function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$dataSet$path.res[,1]);
}

GetMatchedCompounds<-function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$dataSet$path.res[,2]);
}

GetIsLipids <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  is.lipid <- mSetObj$lipid.feats
  
  if(is.lipid){
    is.lipid <- "lipid"
  }else{
    is.lipid <- "met"
  }
  
  return(is.lipid)
}
