#########################################################################################
#set constants
#niter=10000
#npop=5
#datafolder = 'BudowleMoretti/'

#########################################################################################
#read in sib files
#runfolder = paste('sib_100000iter_simtheta0.010_calctheta0.000/', sep='')

#read in true unrel LRs and CIs
#LRs = array(0, dim=c(niter,npop,npop,npop)) #LRs[niter, sim1 pop, sim2 pop, assumed pop]
#for (i in 1:niter) {
#  filename = paste('results/autYsimUnrelPopPair/', datafolder, 'fpr/', runfolder, 'LRs_iter', i-1, '.res', sep='') 
#  LRs[i,,,] = array(scan(filename), dim=c(5, 5, 5))
#}

#########################################################################################
#read in p-o files
#runfolder = paste('p-o_100000iter_simtheta0.010_calctheta0.000/', sep='')

#read in true unrel LRs and CIs
#poLRs = array(0, dim=c(niter,npop,npop,npop)) #LRs[niter, sim1 pop, sim2 pop, assumed pop]
#for (i in 1:niter) {
#  filename = paste('results/autYsimUnrelPopPair/', datafolder, 'fpr/', runfolder, 'LRs_iter', i-1, '.res', sep='') 
#  poLRs[i,,,] = array(scan(filename), dim=c(5, 5, 5))
#}


#########################################################################################
# calc unrel false positives with CA method
#sibFPRs = array(0, dim=c(npop, npop))
#poFPRs = array(0, dim=c(npop, npop))

#for(i in 1:5) {
#  for(j in 1:5) {
#    for(k in 1:niter) {
#      tmpLRs = sapply(LRs[k,i,j,2:4], function(x) { if(is.na(x) || is.infinite(x)) { 0 } else { x }})
#      if(sum(tmpLRs/niter>.1)==3 && max(tmpLRs/niter)>1.0) {
#        sibFPRs[i,j] = sibFPRs[i,j] + 1
#      }

#      tmpLRs = sapply(poLRs[k,i,j,2:4], function(x) { if(is.na(x) || is.infinite(x)) { 0 } else { x }})
#      if(sum(tmpLRs/niter>.1)==3 && max(tmpLRs/niter)>1.0) {
#        poFPRs[i,j] = poFPRs[i,j] + 1
#      }
#    }
#  }
#}



#########################################################################################
# read in fprs from autYsimUnrelPopPair runs

npop=5
#sibFPRs = t(array(scan('results/autYsimUnrelPopPair/BudowleMoretti/fpr/sib_100000000iter_simtheta0.010_calctheta0.000/FPRtab.res'), dim=c(npop,npop)))
#poFPRs = t(array(scan('results/autYsimUnrelPopPair/BudowleMoretti/fpr/p-o_100000000iter_simtheta0.010_calctheta0.000/FPRtab.res'), dim=c(npop,npop)))
sibFPRs = t(array(0, dim=c(npop,npop)))
poFPRs = t(array(0, dim=c(npop,npop)))

for(i in 1:2) {
  sibFPRs = sibFPRs + t(array(scan(paste('results/autYsimUnrelPopPair/BudowleMoretti/fpr/sib_100000000iter_simtheta0.010_calctheta0.000_FPRtab_rep',i,'.res',sep='')), dim=c(npop,npop)))
  poFPRs = poFPRs + t(array(scan(paste('results/autYsimUnrelPopPair/BudowleMoretti/fpr/p-o_100000000iter_simtheta0.010_calctheta0.000_FPRtab_rep',i,'.res',sep='')), dim=c(npop,npop)))
}

sibFPRs = sibFPRs/2
poFPRs = poFPRs/2

#########################################################################################
# calc fp breakdown by group

CAstate = c(.128, .067, .42, .373, .012)
CAprison = c(.006, .303, .268, .413, .009)

CAstate = c(.130, .062, .401, .376, .010) #from US census 2010  http://quickfacts.census.gov/qfd/states/06000.html
CAprison = c(.006, .290, .257, .396, .009) #from Yun's excel sheet

CAstate = c(.132789, .063329, .409602, .384065, .010215) #from US census 2010, sums to 1.0
CAprison = c(.006263, .302714, .268267, .413361, .009395) #from Yun's excel sheet, sums to 1.0


statepoppairprop = array(0, dim=c(npop,npop))
prisonpoppairprop = array(0, dim=c(npop,npop))
for(i in 1:5) { #found profile loop
  for(j in 1:5) { #database profile loop
    statepoppairprop[i,j] = CAstate[i]*CAstate[j]
    prisonpoppairprop[i,j] = CAprison[i]*CAprison[j]
  }
}

statesibFPRs = array(0, dim=c(npop))
prisonsibFPRs = array(0, dim=c(npop))
statepoFPRs = array(0, dim=c(npop))
prisonpoFPRs = array(0, dim=c(npop))
for(i in 1:npop) {
  statesibFPRs[i] = sum(c((statepoppairprop*sibFPRs)[i,], (statepoppairprop*sibFPRs)[,i])) - (statepoppairprop*sibFPRs)[i,i]
  prisonsibFPRs[i] = sum(c((prisonpoppairprop*sibFPRs)[i,], (prisonpoppairprop*sibFPRs)[,i])) - (prisonpoppairprop*sibFPRs)[i,i]
  statepoFPRs[i] = sum(c((statepoppairprop*poFPRs)[i,], (statepoppairprop*poFPRs)[,i])) - (statepoppairprop*poFPRs)[i,i]
  prisonpoFPRs[i] = sum(c((prisonpoppairprop*poFPRs)[i,], (prisonpoppairprop*poFPRs)[,i])) - (prisonpoppairprop*poFPRs)[i,i]
}

statesibFPRpercent = statesibFPRs/sum(statepoppairprop*sibFPRs)
prisonsibFPRpercent = prisonsibFPRs/sum(prisonpoppairprop*sibFPRs)
statepoFPRpercent = statepoFPRs/sum(statepoppairprop*poFPRs)
prisonpoFPRpercent = prisonpoFPRs/sum(prisonpoppairprop*poFPRs)


sum(statepoppairprop*sibFPRs)
sum(prisonpoppairprop*sibFPRs)
sum(statepoppairprop*poFPRs)
sum(prisonpoppairprop*poFPRs)




