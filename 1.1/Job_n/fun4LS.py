from ctm4fitLFe30 import ctm4fitLFe30
from Edge import Edge
import numpy as np

def fun4LS(Opt,x,Y,CoeffNames,Fix_in_cicle_var,funf,W):
	Stt=np.append(Opt,Fix_in_cicle_var)
	for jj in np.arange(len(CoeffNames)):
		StrEval = "%s = %.10g" %(CoeffNames[jj],Stt[jj])
		exec(StrEval)
	return np.sqrt(W)*(eval(funf) - Y)

def LS_Y_out(Opt,x,CoeffNames,Fix_in_cicle_var,funf):
	Stt=np.append(Opt,Fix_in_cicle_var)
	for jj in np.arange(len(CoeffNames)):
		StrEval = "%s = %.10g" %(CoeffNames[jj],Stt[jj])
		exec(StrEval)
	return eval(funf)