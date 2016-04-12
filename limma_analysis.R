library(limma)
targets = readTargets("GSE28166_targets.txt")
targets$FileName = paste(targets$FileName,'.txt')
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

x = read.maimages(targets, path = '/home/wenlei/microarray_dir', source = 'agilent', green.only = T)

#background correction, various parameters can be set
y <- backgroundCorrect(x, method="normexp", offset=16)

#Normalize between the arrays. Select a method appropriate for your data
y <- normalizeBetweenArrays(y, method="quantile")

# Use the avereps function to average replicate spots.
y.ave <- avereps(y, ID=y$genes$ProbeName)

#Build the design matrix for the linear modelling function.
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)


# Apply the intensity values to lmFit
fit <- lmFit(y.ave, design)

contrast.matrix <- makeContrasts("Mock_0H-VN1203_0H","Mock_3H-VN1203_3H", "Mock_7H-VN1203_7H", "Mock_12H-VN1203_12H", "Mock_18H-VN1203_18H", "Mock_24H-VN1203_24H", levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, adjust="BH", coef="Mock_3H-VN1203_3H", genelist=y.ave$genes, number=Inf)

sig_genes <- output[output$adj.P.Val<0.05 & abs(output$logFC)>1,]
