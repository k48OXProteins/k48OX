# accepts two tab-delimited gene list files, one tab-delimited gene universe file, and one CSV data table of characteristics
# The first column of each list should be the UNIPROT ID
# The second column of each list should be the ENSEMBL gene name
# requires: optparse  org.Sc.sgd.db

# usage: Rscript yPAD.R --list1 path/to/file --list2 path/to/file --univ path/to/file --data path/to/file --out output/dir

# 161108_yeast_data.csv is the default data table
# data/univ.txt is the default universe

library(optparse)

opts <- list(
  make_option('--list1',default = "data/POX.Upsig.H202.txt",type = 'character'),
  make_option("--list2",default="data/UBQ.Upsig.H202.txt",type = "character"),
  make_option("--univ",default="data/univ.txt",type = "character"),
  make_option("--data",default="161108_yeast_data.csv",type = "character"),
  make_option("--out",default="out",type = "character")
  )
# opts <- add_option(opts,"--list2",default="data/UBQ.Upsig.H2O2.txt",type = "character")
# opts <- add_option(opts,"--univ",default="data/univ.txt",type = "character")
# opts <- add_option(opts,"--data",default="161108_yeast_data.csv",type = "character")
opts <- parse_args(OptionParser(option_list = opts))
# if(is.na(opts$input)) stop()

# removes spaces from a filename
# could cause problems with reading filenames with spaces in them

pastelst <- function(lst,sep=' ') Reduce(function(x,y) paste(x,y,sep=sep), lst)
get.filename <- function(filename) pastelst(strsplit(filename,' ')[[1]],sep='_')

# concatenates path and filename into output and writes a data.frame to csv.
# any folders in the path that do not exist are created.
dir.csv <- function(x,filename, path = '~/Dropbox/vogel_christiano',...){
  path <- get.filename(path)
  if(!dir.exists(path)) dir.create(path,recursive = T)
  filename <- paste(path, filename, sep = '/')
  filename <- get.filename(filename)
  write.table(x,paste(filename,'csv',sep='.'), sep=',',...)
}

ensembl.to.uniprot <- function(t.half){
  require("org.Sc.sgd.db")
  genes <- select(org.Sc.sgd.db,row.names(t.half),'UNIPROT','ENSEMBL')
  ecounts <- table(genes$ENSEMBL)
  efilt <- names(ecounts[ecounts==1])
  pcounts <- table(genes$UNIPROT)
  pfilt <- names(pcounts[pcounts==1])
  new <- genes[genes$UNIPROT%in%pfilt&genes$ENSEMBL%in%efilt,]
  result <- t.half[new$ENSEMBL,]
  row.names(result) <- new$UNIPROT
  return(result)
}

# add #genes
# pull out genes that fall into categories
# add ls1 not in ls2
# sgd: published ds, curated data under downloads
# protein properties.tab
# Ingolia translation study 2009
# Arava
# TANGO: propensity to aggregate, disopred: predict disassembly
# change aa to relative
enrich.feature <- function(genes, univ,out=''){
    dat <- read.csv(opts$data,row.names = 1)
    if(!any(row.names(univ)%in%row.names(dat))) dat <- ensembl.to.uniprot(dat)
    dat <- dat[univ,]
    # convert raw aa freqs to relative aa freqs
    aas <- names(dat)[14:33]
    aafreq <- dat[,aas]/dat$PROTEIN.LENGTH
    dat[,aas] <- aafreq
    # features <- names(dat)[sapply(names(dat),function(x) 
    #   is.numeric(dat[genes,x]))]
    # result.tt <- as.data.frame(sapply(features,function(x) unlist(t.test(dat[univ,x],dat[genes,x]))))
    # result.wt <- as.data.frame(sapply(features,function(x) unlist(wilcox.test(dat[univ,x],dat[genes,x]))))
    result.tt <- featuretest(t.test,dat,univ,genes)
    result.wt <- featuretest(wilcox.test,dat,univ,genes)
    dir.csv(genes,'genes',out,row.names=F,col.names=F)
    dir.csv(result.tt,'t.test',out)
    dir.csv(result.wt,'wilcoxon.test',out)
}

featuretest <- function(fn,dat,x,y){
  features <- names(dat)[sapply(names(dat),function(i) 
    is.numeric(dat[y,i]))]
  result <- as.data.frame(
      t(sapply(
        features,function(i) fn(dat[x,i],dat[y,i]))))
  # result <- as.data.frame(sapply(tmp,unlist2))
  # row.names(result) <- row.names(tmp)
  result$FDR <- p.adjust(unlist2(result$p.value),'fdr')
  result <- t(as.data.frame(apply(result,1,unlist)))
  return(result)
}

lsdir <- function(lst){
  lst <- strsplit(lst,"/")[[1]]
  lst <- strsplit(lst[length(lst)],'.',T)[[1]]
  lst <- paste(lst[-length(lst)],collapse='.')
  return(lst)
}

ls1 <- read.table(opts$list1,fill = T,header = F)[,1]
ls2 <- read.table(opts$list2,fill = T,header = F)[,1]
univ <- read.table(opts$univ,fill = T,header = F)[,1]
enrich.feature(ls1,univ,paste(opts$out,lsdir(opts$list1),sep="/"))
enrich.feature(ls2,univ,paste(opts$out,lsdir(opts$list2),sep='/'))
enrich.feature(intersect(ls1,ls2),univ,opts$out)
enrich.feature(ls1[!ls1%in%ls2],univ,paste(opts$out,paste(lsdir(opts$list1),'exclusive',sep='.'),sep="/"))
enrich.feature(ls2[!ls2%in%ls1],univ,paste(opts$out,paste(lsdir(opts$list2),'exclusive',sep='.'),sep="/"))

# ox.sig <- read.table("~/grive/yPAD/TestFile/list_test_OXsigH2O2.lst",fill = T)[,1]
# ox.univ <- read.table("~/grive/yPAD/TestFile/list_test_OXuniverse.lst",fill = T)[,1]
# ub.sig <- read.table("~/grive/yPAD/TestFile/list_test_UBsigH2O2.lst",fill = T)[,1]
# ub.univ <- read.table("~/grive/yPAD/TestFile/list_test_UBuniverse.lst",fill = T)[,1]
# all.univ <- union(ox.univ,ub.univ)
# 
# enrich.feature(ox.sig,ox.univ,'~/grive/yPAD/Keira/ox')
# enrich.feature(ub.sig,ub.univ,'~/grive/yPAD/Keira/ub')
# enrich.feature(ox.sig,all.univ,'~/grive/yPAD/Keira/ox.all')
# enrich.feature(ub.sig,all.univ,'~/grive/yPAD/Keira/ub.all')
