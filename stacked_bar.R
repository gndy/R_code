library(scales)
library(ggplot2)
library(reshape2)
theme_set(theme_gray(base_size = 25))
mapping_rate = read.table('mapping_rate', sep = ' ')

ec = mapping_rate[1,-1]
names(ec) = c('0_1','0_2','2_1','2_2','4_1','4_2','12_1','12_2','24_1','24_2')
rownames(ec) = 'host'
virus = 100-ec[1,]
rownames(virus) = 'mers'
dat = rbind(ec,virus)

datm <- melt(cbind(dat, ind = rownames(dat)), id.vars = c('ind'))

ggplot(datm,aes(x = variable, y = value,fill = ind, width = 0.5)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format())

