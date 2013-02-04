#########################################################################################
#set constants
niter=1000000
options(scipen=niter*10)
npop=5
datafolder = 'BudowleMoretti/'
N = 1824085 #from jan 2012 update http://www.fbi.gov/about-us/lab/codis/ndis-statistics#California

#########################################################################################
#read in sib LR test files
runfolder = paste('sib_', format(niter), 'iter_simtheta0.010_calctheta0.000/', sep='')

#read in sim unrel LRs and CIs
unrelsibLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'fpr/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  unrelsibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim po LRs and CIs
posibLRs = array(0, dim=c(niter,npop,npop))
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'p-o/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  posibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim sib LRs and CIs
sibsibLRs = array(0, dim=c(niter,npop,npop))
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'sib/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  sibsibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim cousin LRs and CIs
cousinsibLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'cousin/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  cousinsibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim halfcousin LRs and CIs
halfcousinsibLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'halfcousin/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  halfcousinsibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim second cousin LRs and CIs
cousin2sibLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, '2cousin/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  cousin2sibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim half sib LRs and CIs
halfsibsibLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'halfsib/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  halfsibsibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim half sib LRs and CIs with subsampling in each population
sshalfsibsibLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'halfsib/', runfolder, 'LRs_iter', format(i), '_subsamp.res', sep='') 
  sshalfsibsibLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#########################################################################################
#read in p-o files
runfolder = paste('p-o_', format(niter), 'iter_simtheta0.010_calctheta0.000/', sep='')

#read in sim unrel LRs and CIs
unrelpoLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'fpr/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  unrelpoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim po LRs and CIs
popoLRs = array(0, dim=c(niter,npop,npop))
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'p-o/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  popoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim sib LRs and CIs
sibpoLRs = array(0, dim=c(niter,npop,npop))
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'sib/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  sibpoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim cousin LRs and CIs
cousinpoLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'cousin/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  cousinpoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim halfcousin LRs and CIs
halfcousinpoLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'halfcousin/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  halfcousinpoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim second cousin LRs and CIs
cousin2poLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, '2cousin/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  cousin2poLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim half sib LRs and CIs
halfsibpoLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'halfsib/', runfolder, 'LRs_iter', format(i), '.res', sep='') 
  halfsibpoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}

#read in sim half sib LRs and CIs with subsampling
sshalfsibpoLRs = array(0, dim=c(niter,npop,npop)) #LRs[niter, sim pop, assumed pop]
for (i in 0:(niter-1)) {
  filename = paste('results/autYsim/', datafolder, 'halfsib/', runfolder, 'LRs_iter', format(i), '_subsamp.res', sep='') 
  sshalfsibpoLRs[i+1,,] = array(scan(filename, quiet=TRUE), dim=c(5, 5))
}



#########################################################################################
#make npop x npop plots of lCIs for related and unrelateds for sibs
filename = paste('results/autYsim/', datafolder, 'sibComparativeDistributions.eps', sep='') 
postscript(filename, horizontal=FALSE, onefile=FALSE)#, paper='special', width=9, height=9)
par(mfrow=c(npop,npop), mar=c(1, 1, 0, 0), oma=c(3, 3, 4, 4))

popgroups = c('Vietnamese Am.', 'African Am.', 'European Am.', 'Latino Am.', 'Native Am.')
for(i in 1:5) {
  for(j in 1:5) {  
    plot(x=c(), y=c(), xlim=c(-20, 20), ylim=c(0,.14), xaxt='n', yaxt='n')
    axis(2, at=c(0, .5, .1), labels=FALSE)
    axis(1, at=c(-10, 0, 10, 20), labels=FALSE)
    if(j==1) { axis(2, tick=FALSE, line=FALSE, outer=FALSE, cex.axis=1, at=c(0, .1)) }
    if(j==1) { axis(2, tick=FALSE, line=FALSE, outer=TRUE, cex.axis=1, at=c(.06), labels='density') }
    if(i==5) { axis(1, tick=FALSE, line=FALSE, outer=FALSE, cex.axis=1, at=c(-10, 0, 10, 20)) }
    if(i==5) { axis(1, tick=FALSE, line=FALSE, outer=TRUE, cex.axis=1, at=c(5), labels=c(expression(paste('log(LR)')))) }

    if(j==5) { axis(4, tick=FALSE, line=FALSE, outer=FALSE, cex.axis=1, at=c(.07), labels=popgroups[i]) }
    if(i==1) { axis(3, tick=FALSE, line=FALSE, outer=FALSE, cex.axis=1, at=c(0), labels=popgroups[j]) }
    


    lines(density(sapply(log(posibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, lty='solid')
    lines(density(sapply(log(halfsibsibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, lty='dashed')
    lines(density(sapply(log(cousinsibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, lty='dotdash')
    #lines(density(sapply(log(halfcousinsibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, lty='dotdash')
    lines(density(sapply(log(cousin2sibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, lty='dotted')
    #lines(density(sapply(log(unrelsibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, lty='dotted')
    abline(v=0, lty='dotted')
    lines(density(sapply(log(sibsibLRs[,j,i]/N), function(x) { if(is.na(x)) -20 else x })), lwd=3, col='red')
  }
}
mtext("assumed population sample", side=4, cex=1.1, outer=TRUE, adj=.5, line=2.5)
mtext("true population sample", side=3, cex=1.1, outer=TRUE, adj=.5, line=2.5)

dev.off()



#########################################################################################
# calc power/false positives with CA method
unrelsibPRs = array(0, dim=c(npop))
posibPRs = array(0, dim=c(npop))
sibsibPRs = array(0, dim=c(npop))
cousinsibPRs = array(0, dim=c(npop))
halfcousinsibPRs = array(0, dim=c(npop))
cousin2sibPRs = array(0, dim=c(npop))
halfsibsibPRs = array(0, dim=c(npop))
sshalfsibsibPRs = array(0, dim=c(npop))
unrelpoPRs = array(0, dim=c(npop))
popoPRs = array(0, dim=c(npop))
sibpoPRs = array(0, dim=c(npop))
cousinpoPRs = array(0, dim=c(npop))
halfcousinpoPRs = array(0, dim=c(npop))
cousin2poPRs = array(0, dim=c(npop))
halfsibpoPRs = array(0, dim=c(npop))
sshalfsibpoPRs = array(0, dim=c(npop))

for(j in 1:5) {
  for(k in 1:niter) {
    if(sum(sibsibLRs[k,j,2:4]/N>.1)==3 && max(sibsibLRs[k,j,2:4]/N)>1.0) {
      sibsibPRs[j] = sibsibPRs[j] + 1
    }
    tmpLRs = sapply(unrelsibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      unrelsibPRs[j] = unrelsibPRs[j] + 1
    }
    tmpLRs = sapply(posibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      posibPRs[j] = posibPRs[j] + 1
    }
    tmpLRs = sapply(cousinsibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      cousinsibPRs[j] = cousinsibPRs[j] + 1
    }
    tmpLRs = sapply(halfcousinsibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      halfcousinsibPRs[j] = halfcousinsibPRs[j] + 1
    }
    tmpLRs = sapply(cousin2sibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      cousin2sibPRs[j] = cousin2sibPRs[j] + 1
    }
    tmpLRs = sapply(halfsibsibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      halfsibsibPRs[j] = halfsibsibPRs[j] + 1
    }
    #tmpLRs = sapply(sshalfsibsibLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    #if(sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
    #  sshalfsibsibPRs[j] = sshalfsibsibPRs[j] + 1
    #}

    if(sum(popoLRs[k,j,2:4]/N>.1)==3 && max(popoLRs[k,j,2:4]/N)>1.0) {
      popoPRs[j] = popoPRs[j] + 1
    }
    tmpLRs = sapply(unrelpoLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      unrelpoPRs[j] = unrelpoPRs[j] + 1
    }
    tmpLRs = sapply(sibpoLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      sibpoPRs[j] = sibpoPRs[j] + 1
    }
    tmpLRs = sapply(cousinpoLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      cousinpoPRs[j] = cousinpoPRs[j] + 1
    }
    tmpLRs = sapply(halfcousinpoLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      halfcousinpoPRs[j] = halfcousinpoPRs[j] + 1
    }
    tmpLRs = sapply(cousin2poLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      cousin2poPRs[j] = cousin2poPRs[j] + 1
    }
    tmpLRs = sapply(halfsibpoLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
      halfsibpoPRs[j] = halfsibpoPRs[j] + 1
    }
    #tmpLRs = sapply(sshalfsibpoLRs[k,j,2:4], function(x) { if(is.infinite(x)) { 0 } else { x }})
    #if(sum(is.na(tmpLRs))==0 && sum(tmpLRs/N>.1)==3 && max(tmpLRs/N)>1.0) {
    #  sshalfsibpoPRs[j] = sshalfsibpoPRs[j] + 1
    #}
  }
}


halfsibpoPRs
sshalfsibpoPRs

halfsibsibPRs
sshalfsibsibPRs


#########################################################################################
# plot positive rates over pops and true relationships

filename = paste('results/autYsim/', datafolder, 'relativeFPRs.eps', sep='') 
postscript(filename, horizontal=FALSE, onefile=FALSE, paper='special', width=9, height=8)
par(mfrow=c(1,2))

colors = rep(c('tomato', 'orange', 'orchid', 'steelblue', 'seagreen'), 5)
popgroups = c('Vietnamese Am.', 'African Am.', 'European Am.', 'Latino Am.', 'Native Am.')
pstyles = rep(seq(1,5,1), 5)

xs = c(rep(1,npop),rep(2,npop),rep(3,npop),rep(4,npop),rep(5,npop),rep(6,npop))

ys = c(cousin2poPRs, halfcousinpoPRs, cousinpoPRs, halfsibpoPRs, sibpoPRs, popoPRs)/niter
plot(x=xs, y=ys, col=colors, pch=pstyles, lwd=3, ylab='positive identification rate', xlab='true relationship', main='parent-offspring test', xaxt='n', ylim=c(0,1))
axis(1, at=seq(from=1,to=6,by=1), labels=c('\nsecond\ncousins', '\nhalf\ncousins', 'cousins', '\nhalf-\nsiblings', 'siblings', 'parent-\noffspring'), tick=FALSE, cex.axis=.6)
for (i in 1:npop) {
  mtext(popgroups[i], col=colors[i], line=-i, adj=0, cex=1, at=1.7)
  points(x=c(1.5), y=c(1-(i-1.4)*.033), col=colors[i], pch=pstyles[i], lwd=3)
}

ys = c(cousin2sibPRs, halfcousinsibPRs, cousinsibPRs, halfsibsibPRs, posibPRs, sibsibPRs)/niter
plot(x=xs, y=ys, col=colors, pch=pstyles, lwd=3, ylab='positive identification rate', xlab='true relationship', main='sibling test', xaxt='n', ylim=c(0,1))
axis(1, at=seq(from=1,to=6,by=1), labels=c('\nsecond\ncousins', '\nhalf-\ncousins', 'cousins', '\nhalf-\nsiblings', '\nparent-\noffspring', 'siblings'), tick=FALSE, cex.axis=.6)
for (i in 1:npop) {
  mtext(popgroups[i], col=colors[i], line=-i, adj=0, cex=1, at=1.7)
  points(x=c(1.5), y=c(1-(i-1.4)*.033), col=colors[i], pch=pstyles[i], lwd=3)
}

dev.off()

