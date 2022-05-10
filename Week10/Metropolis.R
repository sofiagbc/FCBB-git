set.seed(3)

#Generate data: create the observation
N=100
n=rpois(N,50)
#rpois(# observations,rate-average events per period)
# n is the number of trials in each observation. around 50
y1=rbinom(N,n,0.6)  #CpG
# rbinom(# observations, # trails/observation, probability of success )
# N is number of observation=number of martian genomes
# n is number of trials per observation=leangth of each reading
# 0.6 is prob of the sequence being a cpg-sequence if we are in island
# y1 is vector with random number of sucesses=number of sequences found in an island

y2=rbinom(N,n,0.1)  #non-CpG
# y1 is number of found sequences in a non-island
y=ifelse(runif(N)<0.2, y1, y2)
# uniform [0,1]
# runif(# observations)
# 20% of the times, you will be in island. lamba=0.2 is prob of being in island
# y is our observation: number of cpg-sequences found in 50 genomes
#===================================


## logposterior
logpost = function(theta) {
	lambda=theta[1]
	p1=theta[2]
	p2=theta[3]
	if(lambda>=1|p1>1|p2>=1|lambda<=0|p1<=0|p2<0|p1<p2) # probabilities and proportion is [0,1]. make sure they are valid
    # p1 has to be > than p2. higher prob in island
	  return(-999999)
  return(sum(log(lambda*p1^y*(1-p1)^(n-y) + (1-lambda)*p2^y*(1-p2)^(n-y))))
	# weight average of probability of your number of observations being island and non-island
	#LTP=(proportion of island)*(probability of your observation if you are in island)+(proportion of non-island)*(probability of your observation if you are in non-island)
	# take log of product to get a sum
}


#===================================

proposal = function(theta) {
## jumping distribution
## choosing the standard deviation of the jumping distribution is discretionary and will affect the acceptance ratio. The ideal acceptance ratio is 10-40%. 
New=theta+rnorm(3,0,0.5)*c(0.01,0.01,0.01)
return(New) 
}

#===================================

NREP = 3000
## starting values
lambda = 0.2
p1 = 0.6
p2 = 0.1

# new parameters far from correct
lambda = 0.4
p1 = 0.85
p2 = 0.4
mchain = data.frame(lambda=rep(NA,NREP),p1=rep(NA,NREP), p2=rep(NA,NREP))

#===================================

mchain[1,] = theta = c(lambda, p1, p2)
## keep track of acceptance rate
acc = 0; 
for(i in 2:NREP) {
  ## lambda
  thetaCandidate = proposal(theta)
  alpha = logpost(thetaCandidate)-logpost(theta)
  if( runif(1) <= exp(alpha) ) {
      acc = acc+1
      theta=thetaCandidate
  }
  ## update chain components
  mchain[i,] = theta
}

#===================================

accept.ratio=acc/NREP
print(paste('Acceptance Ratio:', accept.ratio, sep=' '))
par(mfrow=c(1,3))
plot(mchain[100:NREP,1],type="l",main="lambda",ylim=0:1)
plot(mchain[100:NREP,2],type="l",main="p1",ylim=0:1)
plot(mchain[100:NREP,3],type="l",main="p2",ylim=0:1)

# dont exclude 100 first
plot(mchain[1:NREP,1],type="l",main="lambda",ylim=0:1)
plot(mchain[1:NREP,2],type="l",main="p1",ylim=0:1)
plot(mchain[1:NREP,3],type="l",main="p2",ylim=0:1)

print(paste('Lambda Mean:', mean(mchain[100:NREP,1]), sep=' '))
print(paste('p1 Mean:', mean(mchain[100:NREP,2]), sep=' '))
print(paste('p2 Mean:', mean(mchain[100:NREP,3]), sep=' '))

