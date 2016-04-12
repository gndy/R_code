expr = read.table('Spink6.Laval.SNP.txt')
tt = expr$V1 == 2 
tc = expr$V1 == 1 
cc = expr$V1 == 0 

tt = expr$V2[tt]
tc = expr$V2[tc]
cc = expr$V2[cc]
length(expr$V1)
length(expr$V2)

x = expr$V1
y = expr$V2
#mylogit = glm(V1~V2,data = expr,family = poisson())
#summary(mylogit)

mcorr = cor.test(x,y)
pvalue_1 = mcorr$p.value

dev.new(width = 1, height = 6)
par(mar = c(5,1,3,1))
par(oma = c(2,5,2,5))
par(mfrow=c(1,3))
boxplot(tt,tc,cc,names = c('TT','TC','CC'), boxwex = 0.3, ylim = c(0.01,0.07), whisklty = 1)
points(1,mean(tt), pch = 18)
points(2,mean(tc), pch = 18)
points(3,mean(cc), pch = 18)

expr = read.table('Spink6.Groningen.SNP.txt')
tt = expr$V1 == 2 
tc = expr$V1 == 1 
cc = expr$V1 == 0 

tt = expr$V2[tt]
tc = expr$V2[tc]
cc = expr$V2[cc]

length(expr$V1)
length(expr$V2)

x = expr$V1
y = expr$V2
mcorr = cor.test(x,y)
pvalue_2 = mcorr$p.value

boxplot(tt,tc,cc,names = c('TT','TC','CC'), ylim = c(0.01,0.07), yaxt='n', boxwex = 0.3,whisklty = 1)
points(1,mean(tt), pch = 18)
points(2,mean(tc), pch = 18)
points(3,mean(cc), pch = 18)

expr = read.table('Spink6.UBC.SNP.txt')
tt = expr$V1 == 2 
tc = expr$V1 == 1 
cc = expr$V1 == 0 

tt = expr$V2[tt]
tc = expr$V2[tc]
cc = expr$V2[cc]

length(expr$V1)
length(expr$V2)

x = expr$V1
y = expr$V2
cor.test(x,y)
mcorr = cor.test(x,y)
pvalue_3 = mcorr$p.value

boxplot(tt,tc,cc,names = c('TT','TC','CC'), ylim = c(0.01,0.07), yaxt='n', boxwex = 0.3,whisklty = 1)
points(1,mean(tt), pch = 18)
points(2,mean(tc), pch = 18)
points(3,mean(cc), pch = 18)

par(ps = 15)
mtext('Normalized SPINK6 Expression', side = 2, outer = TRUE, line = 2)
mtext('rs1432689', side = 1, outer = TRUE, line = -1)
mtext(expression('P = 9.43 x 10'^'-12'), side = 3, outer = TRUE, line = -2)


library(metap)
pvalues = c(pvalue_1,pvalue_2,pvalue_3)
sumlog(pvalues)
