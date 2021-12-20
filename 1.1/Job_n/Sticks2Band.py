import numpy as np

def Sticks2Band(x,E,f,MinL,MaxL,Split,Gap,Shape):

	B = np.zeros(x.size)

	l=0
	
	for i in E:

		if i <= Split:
			W = MinL
		elif i>Split and i <= Split+Gap:
			W = MinL + (MaxL-MinL)*(i-Split)/Gap
		elif i > Split + Gap:
			W = MaxL

		B =  B+f[l]*((1-Shape)/(np.pi*W*(1+((x-i)/W)**2))+Shape*(np.sqrt(2*np.log(2)))/(W*np.sqrt(2*np.pi))*np.exp(-((x-i))**2/(2*(W/(np.sqrt(2*np.log(2))))**2)))

		l=l+1

	return B
