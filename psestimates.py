# -*- coding: utf-8 -*-
import time
from datetime import datetime
import itertools
import numpy as np
import statsmodels.api as sm

def getData(path):
	data=np.genfromtxt(path,names=True,delimiter=',')
	return data


def PSestimates(data,Xb,Kb,treat,tol):
	while len(Xb)>0:
		#Calculate first log-likelihood
		if len(Kb)==0: L1=-1000
		else: L1=sm.Logit(data[treat],sm.add_constant(data[Kb])).fit(disp=0).llf
		maxLike=0
		maxVar=''
		for x in range(0,len(Xb)):
			#Calculate seconf log-likelihood
			L2=sm.Logit(data[treat],sm.add_constant(data[Kb+Xb[x:x+1]])).fit(disp=0).llf
			#print aux.summary()
			#Likelihood-ratio test
			D=2*(L2-L1)
			if D>maxLike:
				maxLike=D
				maxVar=Xb[x]
		if maxLike>tol:
			#print "We added var: ", (maxVar), ", with likelihood ratio of", maxLike
			Kb.append(maxVar)
			Xb=list(set(Xb)-set(Kb))
		else: break
	return Kb


def main():
	data=getData('nswre.txt')
	print type(data)
	Xb=list( data.dtype.names)
	treat=Xb[0]
	#define Kb
	Kb=[]
	Kl=[]
	del Xb[0] #always the treatment has to be first
	Xb=list(set(Xb)-set(Kb))
	tol=1
	Kb=PSestimates(data,Xb,Kb,treat,tol)
	print Kb
	Xq=list(itertools.product(Kb, Kb))

	print Xq




if __name__ == '__main__':
	print str(datetime.now())
	t1=time.time()
	main()
	print "The run time was ", (time.time()-t1), " seconds."