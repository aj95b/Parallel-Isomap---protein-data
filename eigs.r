calc_eigs <- function(){
main_mat<-read.table("main_mat")
ev <- eigen(main_mat)
(val <- ev$values)
(vec <- ev$vectors)

write.table(val,file="val", append= FALSE,quote= FALSE,sep = " ",eol = "\n",na = "NA",dec = ".",row.names= FALSE,col.names= FALSE,qmethod= c("escape","double"))
write.table(vec,file="vec", append= FALSE,quote= FALSE,sep = " ",eol = "\n",na = "NA",dec = ".",row.names= FALSE,col.names= FALSE,qmethod= c("escape","double"))
}
calc_eigs()
