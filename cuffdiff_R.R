library(cummeRbund)
sample_names = c('DEAB_SC0','DEAB_SC2','DEAB_SC4','DEAB_SC12','DEAB_SC24')
#sample_names = c('DEAB_EC0','DEAB_EC2','DEAB_EC4','DEAB_EC12','DEAB_EC24')
filename = '/home/wenlei/RNASeq/cuffdiff_out_SC/'
#filename = '/home/wenlei/RNASeq/cuffdiff_out_EC/'
baseline_label = 'DEAB_SC0' # to change
print(filename)
mindex = 2
print(sample_names[mindex])
output_name = sample_names[mindex]
cuff_data = readCufflinks(filename)
#gene.repfpkm = repFpkm(genes(cuff_data))
gene_diff_data = diffData(genes(cuff_data),baseline_label,sample_names[mindex])
#up_gene_data = subset(gene_diff_data, log2_fold_change>1&significant == 'yes')
#down_gene_data = subset(gene_diff_data, log2_fold_change< -1&significant == 'yes')
de_gene_data = subset(gene_diff_data, abs(log2_fold_change)>0.58&significant == 'yes')
#csDensity(genes(cuff_data))
#fpkmSCVPlot(genes(cuff_data))

#csScatter(genes(cuff_data),baseline_label,sample_names[mindex])
#print(nrow(up_gene_data)+nrow(down_gene_data))
#write.table(gene.repfpkm,paste('/home/wenlei/','SC','_replicates_fpkm_table',sep = ''),sep='\t')
#write.table(up_gene_data,paste('/home/wenlei/regulation/',output_name,'_upregulated',sep = ''),sep='\t')
#write.table(down_gene_data,paste('/home/wenlei/regulation/',output_name,'_downregulated',sep = ''),sep='\t')
write.table(de_gene_data,paste('/home/wenlei/',output_name,'_regulated',sep = ''),sep='\t')
