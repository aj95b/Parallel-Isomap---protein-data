calc_cor <- function(){
dist<-read.table("dist")
D<-read.table("D")
corcoef <- cor(dist,D)
write.table(corcoef,file="corcoef", append= FALSE,quote= FALSE,sep = " ",eol = "\n",na = "NA",dec = ".",row.names= FALSE,col.names= FALSE,qmethod= c("escape","double"))
}
calc_cor()
