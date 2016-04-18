library(limma)
setwd('/home/wenlei/')
targets = readTargets("GSE76599_target.txt")#change here
targets$FileName = targets$FileName
targets$FileName = gsub(' ','',targets$FileName)


#targets file describes the array names you want to analysze together
#example for targets file

#SampleNumber	FileName	Condition
#1	US83403541_252144010018_S01_GE1-v5_10_Apr08_1_2.txt	Brain
#2	US83403541_252144010020_S01_GE1-v5_10_Apr08_1_1.txt	Brain
#3	US83403541_252144010021_S01_GE1-v5_10_Apr08_1_1.txt	Lung
#4	US83403541_252144010021_S01_GE1-v5_10_Apr08_1_2.txt	Lung
#5	US83403541_252144010023_S01_GE1-v5_10_Apr08_1_1.txt	Liver
#6	US83403541_252144010023_S01_GE1-v5_10_Apr08_1_2.txt	Liver

gse_path = '/home/wenlei/microarray_dir/GSE76599'#change here
x = read.maimages(targets, path = gse_path, source = 'agilent', green.only = T)

#background correction, various parameters can be set
y <- backgroundCorrect(x, method="normexp", offset=16)

#Normalize between the arrays. Select a method appropriate for your data
y <- normalizeBetweenArrays(y, method="quantile")

# Use the avereps function to average genes with multiple probes.
y.ave <- avereps(y, ID=y$genes$ProbeName)#for some platform with no Genename, this should be ProbeName

#Build the design matrix for the linear modelling function.
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)


# Apply the intensity values to lmFit
fit <- lmFit(y.ave, design)

#contrast_GSE28166
#cmp_vector = c("VN1203_0H-Mock_0H","VN1203_3H-Mock_3H", "VN1203_7H-Mock_7H", "VN1203_12H-Mock_12H", "VN1203_18H-Mock_18H", "VN1203_24H-Mock_24H")

#contrast_GSE33267
#cmp_vector = c('SCL005_WT_0H-SCL005_mock_0H', 'SCL005_DORF6_0H-SCL005_mock_0H', 'SCL005_WT_3H-SCL005_mock_3H', 'SCL005_DORF6_3H-SCL005_mock_3H', 'SCL005_WT_7H-SCL005_mock_7H', 'SCL005_DORF6_7H-SCL005_mock_7H', 'SCL005_WT_12H-SCL005_mock_12H', 'SCL005_DORF6_12H-SCL005_mock_12H', 'SCL005_WT_24H-SCL005_mock_24H', 'SCL005_DORF6_24H-SCL005_mock_24H', 'SCL005_WT_30H-SCL005_mock_30H', 'SCL005_DORF6_30H-SCL005_mock_30H', 'SCL005_WT_36H-SCL005_mock_36H', 'SCL005_DORF6_36H-SCL005_mock_36H', 'SCL005_WT_48H-SCL005_mock_48H', 'SCL005_DORF6_48H-SCL005_mock_48H', 'SCL005_WT_54H-SCL005_mock_54H', 'SCL005_DORF6_54H-SCL005_mock_54H', 'SCL005_WT_60H-SCL005_mock_60H', 'SCL005_DORF6_60H-SCL005_mock_60H', 'SCL005_WT_72H-SCL005_mock_72H', 'SCL005_DORF6_72H-SCL005_mock_72H')

#contrast_GSE37571
#cmp_vector = c('ICL006_CA04_0h-ICL006_mock_0h', 'ICL006_CA04_3h-ICL006_mock_3h', 'ICL006_CA04_7h-ICL006_mock_7h', 'ICL006_CA04_12h-ICL006_mock_12h', 'ICL006_CA04_18h-ICL006_mock_18h', 'ICL006_CA04_24h-ICL006_mock_24h', 'ICL006_CA04_30h-ICL006_mock_30h', 'ICL006_CA04_36h-ICL006_mock_36h', 'ICL006_CA04_48h-ICL006_mock_48h')

#contrast_GSE69026
#cmp_vector = c('ICL102_AH1_0hr-ICL102_Mock_0hr', 'ICL102_FM_0hr-ICL102_Mock_0hr', 'ICL102_691_0hr-ICL102_Mock_0hr', 'ICL102_AH1_7hr-ICL102_Mock_7hr', 'ICL102_FM_7hr-ICL102_Mock_7hr', 'ICL102_691_7hr-ICL102_Mock_7hr', 'ICL102_AH1_12hr-ICL102_Mock_12hr', 'ICL102_FM_12hr-ICL102_Mock_12hr', 'ICL102_691_12hr-ICL102_Mock_12hr', 'ICL102_AH1_24hr-ICL102_Mock_24hr', 'ICL102_FM_24hr-ICL102_Mock_24hr', 'ICL102_691_24hr-ICL102_Mock_24hr')

#contrast_GSE76599
cmp_vector = c('ICL103_VN1203_0hr-ICL103_Mock_0hr', 'ICL103_627E_0hr-ICL103_Mock_0hr', 'ICL103_NS1_0hr-ICL103_Mock_0hr', 'ICL103_VN1203_7hr-ICL103_Mock_7hr', 'ICL103_627E_7hr-ICL103_Mock_7hr', 'ICL103_NS1_7hr-ICL103_Mock_7hr', 'ICL103_VN1203_12hr-ICL103_Mock_12hr', 'ICL103_627E_12hr-ICL103_Mock_12hr', 'ICL103_NS1_12hr-ICL103_Mock_12hr', 'ICL103_VN1203_24hr-ICL103_Mock_24hr', 'ICL103_627E_24hr-ICL103_Mock_24hr', 'ICL103_NS1_24hr-ICL103_Mock_24hr')

contrast.matrix <- makeContrasts(contrasts = cmp_vector, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

prefix = 'GSE76599'#change here for outputname prefix
setwd('/home/wenlei/tmp2/')

de_genes_wt <- c('')
de_genes_orf6 <- c('')

for(i in cmp_vector){

output <- topTable(fit2, adjust="BH", coef= i, genelist=y.ave$genes, number=Inf)

sig_genes <- output[output$adj.P.Val<0.05 & abs(output$logFC) > 1,]
sig_gene_names <- sig_genes$GeneName

if(grepl('WT',i)){
  de_genes_wt <- c(de_genes_wt,sig_gene_names)
  
}else{
  de_genes_orf6 <- c(de_genes_orf6,sig_gene_names)
  
}

up_sig_genes <- output[output$adj.P.Val<0.05 & output$logFC > 1,]
down_sig_genes <- output[output$adj.P.Val<0.05 & output$logFC < -1,]
write.table(down_sig_genes, file=paste(prefix,'downregulated',i,sep = '_'), sep="\t", quote=FALSE)
write.table(up_sig_genes, file=paste(prefix,'upregulated',i,sep = '_'), sep="\t", quote=FALSE)

print(i)
print(length(up_sig_genes$ProbeName))
#print(length(unique(up_sig_genes$GeneName)))
print('======')
print(length(down_sig_genes$ProbeName))
#print(length(unique(down_sig_genes$GeneName)))
}

print(length(unique(de_genes_orf6)))
print(length(unique(de_genes_wt)))
