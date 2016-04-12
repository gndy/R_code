library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
raw_data = read.table('expr_wgcna_EC_gt1_transformed',header = TRUE)
expr_data = as.data.frame(t(raw_data[,-1]))
names(expr_data) = raw_data$SAMPLE
enableWGCNAThreads()
powers = c(c(1:10),seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(expr_data,powerVector = powers, verbose = 5)
sizeGrWindow(9,5)
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")


net = blockwiseModules(expr_data, power = 8,
                       TOMType = "signed", minModuleSize = 40,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "mTOM_rnaseq", 
                       verbose = 3)


moduleColors = labels2colors(net$colors)
table(moduleColors)
modules = unique(moduleColors)

#resorted_color = read.table('color_sc_from_ec',sep = '\t')

sizeGrWindow(20, 9)
#merged_colors = resorted_color$V2
# Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
par(ps = 12)
#plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
plotDendroAndColors(net$dendrograms, moduleColors,
                                        
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#inmodule = (moduleColors == 'turquoise')
#modprobes = names(expr_data)[inmodule]
MEs0 = moduleEigengenes(expr_data, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

transcript_rate = as.data.frame(c(0.91,0.91,0.91,0.91,0.91,0.91,0.75,0.75,0.3,0.3))
names(transcript_rate) = 'rate'
nSamples = nrow(expr_data)
moduleTraitCor = cor(MEs, transcript_rate, use = 'p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);



TOM = TOMsimilarityFromExpr(expr_data, power = 18);
module = c("purple");
probes = names(expr_data)
inModule = (moduleColors==module)
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0.2
                            )



par(ps = 3)
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = 'rate',
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#geneModuleMembership
weight = transcript_rate
names(weight) = "weight"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(expr_data, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


#geneTraitSignificance
geneTraitSignificance = as.data.frame(cor(expr_data, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                   abs(geneTraitSignificance[moduleGenes, 1]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = "Gene significance for body weight",
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#filterGenes = abs(geneModuleMembership$MMblue)>0.9 & abs(geneTraitSignificance)>0.9

#write.csv(all_genes,file='all_module_gene_SC')
