tt = read.table('TT',header = TRUE)
tc = read.table('TC',header = TRUE)
cc = read.table('CC',header = TRUE)
dev.new(width = 0.01, height = 6)


par(mar = c(5,3,3,1))
par(mfrow=c(1,2))
boxplot(tt$TT,tc$TC,cc$CC,names = c('TT','TC','CC'), ylim = c(-0.2,0.2), boxwex = 0.3,whisklty = 1)
points(1,mean(tt$TT),pch = 18)
points(2,mean(tc$TC), pch = 18)
points(3,mean(cc$CC), pch = 18)

boxplot(tt$TT,tc$TC,cc$CC,names = c('TT','TC','CC'), ylim = c(-0.2,0.2), boxwex = 0.3,whisklty = 1)

points(1,mean(tt$TT), pch = 18)
points(2,mean(tc$TC), pch = 18)
points(3,mean(cc$CC), pch = 18)

