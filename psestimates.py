# -*- coding: utf-8 -*-
import numpy as np
import statsmodels.api as sm
import time
import itertools
from datetime import datetime
from numpy.lib import recfunctions as rfn

# Input data
def getData(path):
	data = np.genfromtxt(path, names = True, delimiter=',')
	return data

# Main algorithm
def main(data,Xb,Kb,treatvar,tol):
	while len(Xb)>0:
		# Compute first log-likelihood
		if len(Kb)==0: L1=-1000
		else: L1=sm.Logit(data[treatvar],sm.add_constant(data[Kb])).fit(disp=0).llf
		maxLike=0
		maxVar=''
		for x in range(0,len(Xb)):
			# Compute second log-likelihood
			L2=sm.Logit(data[treatvar],sm.add_constant(data[Kb+Xb[x:x+1]])).fit(disp=0).llf
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

def genQuad(Kb,data):
	Xq=[]
	for x in range(1,len(Kb)):
		for y in range(0,x+1):
			if x==y:
				aux=set(data[Kb[x]])
				if len(aux)==2 and 0 in aux and 1 in aux:
					print Kb[x]," is a dummy"
				else:
					Xq.append(str(Kb[x]+'#'+Kb[y]))
					NewVar=data[Kb[x]]*data[Kb[y]]
					data=rfn.append_fields(data,names=str(Kb[x]+'#'+Kb[y]),data=NewVar,usemask=False)
			else:
				Xq.append(str(Kb[x]+'#'+Kb[y]))
				NewVar=data[Kb[x]]*data[Kb[y]]
				data=rfn.append_fields(data,names=str(Kb[x]+'#'+Kb[y]),data=NewVar,usemask=False)
	return (Xq,data)

def PSestimate(data,treatvar,totry=None,clin=1,cquad=2.71):
	# Get data:
	data = getData(data)
	# Define list of covariates to try (default is all):
	if totry is None:
		Xb = list(data.dtype.names)
		del Xb[Xb.index(treatvar)]
	else:
		Xb = list(totry)
	# Define Kb
	Kb = []
	Kl = []

	# First order terms:
	Xb=list(set(Xb)-set(Kb))
	tol=clin
	Kb=main(data,Xb,Kb,treatvar,tol)
	print "First order variables chosen: ", Kb

	# Second order terms:
	(Xq,data)=genQuad(Kb,data)
	tol=cquad
	Kb=main(data,Xq,Kb,treatvar,tol)
	print "All covariates chosen: ", Kb


if __name__ == '__main__':
	print str(datetime.now())
	t1=time.time()
	PSestimate(data='nswre.txt', treatvar='treat', totry=['black', 'ed', 'nodeg', 're74'])
	print "Run time was", (time.time()-t1), "seconds."