pox <- read.csv("~/grive/yPAD/Keira/data/EBprot.POX_copy.csv",as.is = T,numerals = 'allow.loss')
ubq <- read.csv("~/grive/yPAD/Keira/data/EBprot.UBQ_copy.csv",as.is = T,numerals = 'allow.loss')
ebnames <- names(pox)
ebnames[7:length(ebnames)] <- c(sapply(ebnames[c(-1:-6,-grep("^X",ebnames))],function(i) rep(i,6)))
ebnames <- paste(ebnames,pox[1,],sep = "")
names(pox) <- ebnames
names(ubq) <- ebnames
pox <- pox[-1,]
ubq <- ubq[-1,]
row.names(pox) <- pox[,1]
row.names(ubq) <- ubq[,1]
lsts <- c("Upsig.H202","Upsig.WCE","Upsig.REC","Upsig.WCE_REC" )

poxlst <- lapply(lsts, function(i) pox[pox[,i]=="*",c("Comb.Protein","Comb.Gene")])
ubqlst <- lapply(lsts, function(i) ubq[ubq[,i]=="*",c("Comb.Protein","Comb.Gene")])
lapply(1:4,function(i) {
  write.table(poxlst[[i]],paste("~/grive/yPAD/Keira/data/POX",lsts[i],"txt",sep="."),row.names = F,col.names = F)
  write.table(ubqlst[[i]],paste("~/grive/yPAD/Keira/data/UBQ",lsts[i],"txt",sep="."),row.names = F,col.names = F)
})

poxuniv <- pox[,c(1,2)]
ubquniv <- ubq[,c(1,2)]
univ <- rbind(poxuniv,ubquniv[!ubquniv$Comb.Protein%in%poxuniv$Comb.Protein,])
write.table(univ,"~/grive/yPAD/Keira/data/univ.txt",row.names = F,col.names = F)
# write.table(poxuniv,paste("~/grive/yPAD/Keira/data/POX","univ","txt",sep="."),row.names = F,col.names = F)
# write.table(ubquniv,paste("~/grive/yPAD/Keira/data/UBQ","univ","txt",sep="."),row.names = F,col.names = F)