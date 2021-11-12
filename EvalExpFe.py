import numpy as np
from ctm4fitLFe_final import ctm4fitLFe
#Para testeos----
#from scipy.io import loadmat
#import matplotlib.pyplot as plt

def EvalExpFe(Data,CoeffNames,Stt,funf):
	x=Data.X

	for jj in np.arange(len(CoeffNames)):
		StrEval = "%s = %.10g" %(CoeffNames[jj],Stt[jj])
		exec(StrEval)

	f = eval(funf)

	return f

#Testeos
"""
JobCell=loadmat("Ejemplo_Batch_Cluster_python.mat",struct_as_record=False)
Job=JobCell['Job'][0,0]
Data = Job.Data[0,0]
CoeffNames=[]
for i in Job.Model[0,0].Global[0,0].CoeffNames[0,:]:
	CoeffNames.append(i[0])
#Para Stt se usó start values from Job.Model.Global.Parameters.Start
Stt=Job.Model[0,0].Global[0,0].Parameters[0,0].Start[:,0]
funf=Job.Model[0,0].Global[0,0].Function[0]

FF=EvalExpFe(Data, CoeffNames,Stt,funf)

x=Data.X[:,0]

plt.plot(x[x<1000],FF[0:299])
plt.show()

mask = np.logical_and(x > 1000, x < 2000)
plt.plot(x[mask],FF[300:599])
plt.show()

FF=np.append(FF,[0,0])
plt.plot(x[x>2000],FF[600:899])
plt.show()"""