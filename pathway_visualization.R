library(org.Hs.eg.db)
library(pathview)
filename = 'EC24'#change filename here
reg_raw = read.table(filename,header = F,sep = '\t')
reg_2 = as.data.frame(reg_raw[,-1])
row.names(reg_2) = reg_raw[,1]
length(unique(reg_raw[,1]))
names(reg_2) = 'VALUE'


sym = row.names(reg_2)
ids = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
mvector = c()
for(i in 1:length(ids)){
  mvector = c(mvector,ids[[i]][1])
}

reg_3 = cbind(reg_2,mvector)
filter = reg_3$mvector != 'NA'
filter = filter %in% TRUE
reg_4 = reg_3[filter,]

reg_5 = as.data.frame(reg_4[,-2])
row.names(reg_5) = reg_4$mvector
reg_5 = as.matrix(reg_5)
mpathwayids = c('04620','04621','04622')
for(i in mpathwayids){
pv.out <- pathview(gene.data = reg_5[, 1], pathway.id = i, species = "hsa", out.suffix = paste(i,' ',filename), plot.col.key = F, kegg.native = T)
}
