import numpy as np 
from scipy.stats import f
from scipy import linalg

def predint(level,p,dfe,rmse,ypred,jac,Rinv):
	E=jac@Rinv
	delta=np.sqrt(np.sum(E*E,axis=1))
	crit=np.sqrt(p*f.ppf(level,p,dfe))
	delta=delta*rmse*crit
	return [ypred-delta,ypred+delta]