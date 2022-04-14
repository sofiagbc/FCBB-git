import numpy as np
import scipy 
import pandas as pd
import matplotlib
from scipy.stats import poisson, binom , norm, uniform 
import matplotlib.pyplot as plt
import argparse

np.random.seed(3)

N=100;
n=poisson.rvs(mu=50, size=N)
y1=binom.rvs(n=n, p=.6, size=N)
y2=binom.rvs(n=n, p=.1, size=N)
val=uniform.rvs(0,1, size = N) 
y = np.where(val<0.2, y1, y2)

def logpost(theta):
    [lam,p1,p2]=theta
    if(lam>=1 or p1>1 or p2>=1 or lam<=0 or p1<=0 or p2<0 or p1<p2):
      return -999999
    suma=np.sum(np.log(lam*p1**y*(1-p1)**(n-y) + (1-lam)*p2**y*(1-p2)**(n-y)))
    return suma

def proposal(theta):
    New=theta+norm.rvs(0,1,size=3)*[0.01,0.01,0.01]
    return New
    
NREP = 3000
## starting values
lam = 0.2
p1 = 0.6
p2 = 0.1

df={'lambda': np.repeat(np.nan,NREP), 'p1': np.repeat(np.nan,NREP),'p2': np.repeat(np.nan,NREP)}
mchain = pd.DataFrame(data=df)

mchain[0:1]=theta=[lam,p1,p2]
acc = 0

for i in range(1,NREP) :
  ## lambda
    thetaCandidate = proposal(theta)
    alpha= logpost(thetaCandidate)-logpost(theta)
    if( uniform.rvs(0,1) <= np.exp(alpha) ) :
        acc = acc+1
        theta=thetaCandidate
    ## update chain components
    mchain[i:i+1] = theta

accept_ratio=acc/NREP
print('Acceptance Ratio: %.3f' % (accept_ratio))
print('Lambda Mean: %.3f' % (np.mean(mchain['lambda'][100:NREP])))
print('p1 Mean: %.3f' % (np.mean(mchain['p1'][100:NREP])))
print('p2 Mean: %.3f' % (np.mean(mchain['p2'][100:NREP])))

plt.subplot(1, 3, 1)
plt.plot(mchain['lambda'][100:NREP])
plt.xlabel("Index")
plt.ylabel("mchain['lambda'][100:NREP]")
plt.ylim(0,1)
plt.title("lambda")

plt.subplot(1, 3, 2)
plt.plot(mchain['p1'][100:NREP])
plt.xlabel("Index")
plt.ylabel("mchain['p1'][100:NREP]")
plt.ylim(0,1)
plt.title("p1")

plt.subplot(1, 3, 3)
plt.plot(mchain['p2'][100:NREP])
plt.xlabel("Index")
plt.ylabel("mchain['p2'][100:NREP]")
plt.ylim(0,1)
plt.title("p2")

plt.tight_layout()

my_parser = argparse.ArgumentParser()

# add arguments
my_parser.add_argument('-o')
files = my_parser.parse_args()
plt.savefig(files.o)
