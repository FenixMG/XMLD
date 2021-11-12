from scipy.io import loadmat
from scipy.io import savemat
import numpy as np
from datetime import datetime
from itertools import product 
from time import time
from scipy.optimize import least_squares
from EvalExpFe import EvalExpFe
from fun4LS import *
from scipy.stats import t
from scipy import linalg
from loadmat_dic import loadmat_dict
from predint import *
import os
from monitor import *

def BGridFe(Sf):
	current_dir=os.path.dirname(os.path.realpath(__file__))
	#[ ] Tal vez sea necesario poner todos los vectores columna en forma de vectores renglón.
	try:
		Fits = Sf['Fits']
		Skip=1
		Aux_Sf = loadmat_dict(Sf)
		Aux_Fits = Aux_Sf['Fits']
		try:
			L = len(Aux_Fits['Cycle'][:])
		except:
			Aux_Fits['Cycle']=[Aux_Fits['Cycle']]
			L=len(Aux_Fits['Cycle'][:])
		monitor_print('Resuming Job at Cycle '+str(L)+' and Attempted Fit '+str(Aux_Fits['Cycle'][L-1]['Counter']),new=True)
	except:
		Skip = 0


	if Skip==0:
		Job = Sf['Job']
	else:
		Job = Fits.Job

	Data = Job.Data
	Parameters = Job.Parameters
	Model = Job.Model
	funf = Model.Global.Function
	lb = Model.Global.Parameters.lb
	ub = Model.Global.Parameters.ub
	np_ = Model.Global.Parameters.np
	gr = Model.Global.Parameters.gr
	Ord = Model.Global.Parameters.ord

	try:
		Cluster=Job.Cluster
		Cluster=1
	except:
		Cluster=0

	Name = Job.Name

	if Skip==0:
		monitor_print(str(datetime.now())+' Starting Job '+Name,new=True)
	else:
		monitor_print(str(datetime.now())+' Starting Job '+Name)

	ex_= Parameters.ex_
	ex=ex_

	W_ = Parameters.W_

	X = Data.X
	Y = Data.Y

	NumFits = Parameters.NumFits

	if Skip == 0:
		Fits = {'fitobj':[],
				'Start':{'Coefficients':[]},
				'YFits':[],
				'Coefficients':[],
				'GOF':{'rsquare':[], 'adjrsquare':[],'RMSE':[],'SSE':[]},
				'CIntLowLine':[],
				'CIntUpLine':[],
				'CIntLowB':[],
				'CIntUpB':[],
				'Time':[],
				'Cycle':[],
				'Job':Job}
	else:
		Fits = Aux_Fits


	with np.errstate(divide='ignore'):
		pN = np.around((ub-lb)/np_+1)

	try:
		fx = Model.Global.Parameters.fx
	except:
		fx = np.zeros(gr.size)

	Coeff_1=[]

	if Skip==0:
		for i in Model.Global.CoeffNames:
			Coeff_1.append(i)
	else:
		sum_aux=gr.size-gr.sum()
		aux_coeff=0
		aux_coeff_2=0
		for i in Model.Global.CoeffNames:
			if gr[int(aux_coeff)]==0:
				Coeff_1.append(Model.Global.CoeffNames[int(aux_coeff_2)].strip())
				aux_coeff_2+=1
			else:
				Coeff_1.append(Model.Global.CoeffNames[int(sum_aux)].strip())
				sum_aux+=1
			aux_coeff+=1

	gub=ub
	glb=lb

	for kkkk in range(5):
		if Skip==1:
			if (kkkk+1)<L:
				continue

		try:
			aux0=Job.Cluster
			if kkkk+1 != aux0.Cycle:
				continue
		except:
			pass

		comb = 1

		Lower = []
		Upper = []
		Order = []
		Coeff_ = []
		Probl_ = []

		Gr = []
		Gr_ = []

		if kkkk>0:
			if (Skip==1 and (kkkk+1)==L) or Cluster==1:
				lb=Fits['Cycle'][kkkk]['Bounds']['Lower']
				ub=Fits['Cycle'][kkkk]['Bounds']['Upper']
			else:
				for jjjj in range(len(Coeff_1)):
					for iiii in range(len(Model.Global.CoeffNames)):
						if Coeff_1[jjjj]==Model.Global.CoeffNames[iiii].strip():
							Stdd = np.std([fila[iiii] for fila in Fits['Coefficients']])
							Ave = np.mean([fila[iiii] for fila in Fits['Coefficients']])
							if Stdd <= 5e-3:
								Stdd=abs(0.1*(gub[jjjj]-glb[jjjj]))
							if lb[jjjj]==0:
								pass
							else:
								if fx[jjjj]==0:
									if lb[jjjj]>0:
										if Ave-2*Stdd>0:
											lb[jjjj]=round((Ave-2*Stdd)*100)/100
										else:
											lb[jjjj]=0
									else:
										lb[jjjj]=round((Ave-2*Stdd)*100)/100
								else:
									AAA=round((Ave-2*Stdd)*100)/100
									if AAA>glb[jjjj]:
										lb[jjjj]=AAA
							if ub[jjjj]==100:
								pass
							else:
								if fx[jjjj]==0:
									if ub[jjjj]<0:
										if Ave+2*Stdd<0:
											ub[jjjj]=round((Ave+2*Stdd)*100)/100
										else:
											ub[jjjj]=0
									else:
										ub[jjjj]=round((Ave+2*Stdd)*100)/100
								else:
									AAA=round((Ave-2*Stdd)*100)/100
									if AAA<gub[jjjj]:
										ub[jjjj]=AAA

		for ii in range(len(lb)):
			if gr[ii]==1:
				Probl_.append(Coeff_1[ii])
				Gr_.append(np.linspace(lb[ii],ub[ii],int(pN[ii])).conj().transpose())
				Order.append(Ord[ii])
				comb = comb*(len(Gr_[-1]))
			else:
				Lower.append(lb[ii])
				Upper.append(ub[ii])
				Coeff_.append(Coeff_1[ii])
		aux1=0
		aux_prob=tuple(Probl_)
		for f in Order:
			Probl_[aux1] = aux_prob[f-1]
			aux1+=1

		for ii in np.arange(len(Gr_)):
			Gr.append(Gr_[Order[ii]-1])

		Coeff = Coeff_
		Probl = Probl_
		if Skip==1 and L ==kkkk+1:
			pass
		else:
			Fits['Cycle'].append({'Bounds':{'Lower':lb}})
			Fits['Cycle'][kkkk]['Bounds']['Upper']=ub
			if kkkk==0:
				Fits['Cycle'][kkkk]['Bounds']['Avg']=0
				Fits['Cycle'][kkkk]['Bounds']['Std']=0

		if Skip==0:
			Model.Global.CoeffNames=Coeff_+Probl_

		LimG=np.ones(20)
		for ii in range(len(Gr)):
			LimG[ii]=len(Gr[ii])

		cnt=0
		if Skip==1 and kkkk+1>L:
			Fits['Cycle'][kkkk]['Counter']=cnt

		for v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20 in product(range(int(LimG[0])),range(int(LimG[1])),range(int(LimG[2])),range(int(LimG[3])),range(int(LimG[4])),range(int(LimG[5])),range(int(LimG[6])),range(int(LimG[7])),range(int(LimG[8])),range(int(LimG[9])),range(int(LimG[10])),range(int(LimG[11])),range(int(LimG[12])),range(int(LimG[13])),range(int(LimG[14])),range(int(LimG[15])),range(int(LimG[16])),range(int(LimG[17])),range(int(LimG[18])),range(int(LimG[19]))):

			cnt+=1
			if Skip==1:
				if cnt<Fits['Cycle'][kkkk]['Counter']:
					continue

			if Cluster==1:
				if cnt<Job.Cluster.Min or cnt>Job.Cluster.Max:
					continue

			Fits['Cycle'][kkkk]['Counter']=cnt

			ProbVal = []

			for ii in range(len(Gr)):
			   ProbVal.append(eval("Gr[ii][v%g]" %(ii+1)))

			PValues = ProbVal

			mc = Job.Parameters.MonteCarlo
			StartPoint = np.zeros((len(lb),mc))
			SSE= np.zeros(mc)

			ti=time()
			Lower = np.array(Lower)
			Upper = np.array(Upper)

			for kk in range(mc):
				SP_ = Lower + np.random.rand(len(Lower))*(Upper-Lower)

				StartPoint[:,kk] = np.append(SP_,ProbVal)
				CoeffNames_for_fun = np.char.strip(Job.Model.Global.CoeffNames)

				try:
					Stt = np.append(SP_,ProbVal)
					f=EvalExpFe(Data,CoeffNames_for_fun,Stt,funf)
				except:
					pass

				Y = np.array(Y)
				W_ = np.array(W_)
				try:
					mask = np.logical_not(ex)

					SqrErr = W_[mask]*(f[mask]-Y[mask])**2
				except:
					SqrErr = W_ * (f - Y) ** 2

				SSErr = sum(SqrErr)
				SSE[kk]=SSErr

			tags = np.argsort(SSE)
			SSE.sort()

			StartPoint = StartPoint[:,tags]
			SSE_SP = SSE[0]

			Fits['GrindN']=cnt
			monitor_print('%s- Attempting Fit %g of %g (cycle %g of 5)...'%(datetime.now(),cnt,comb,kkkk+1))

			Nombres_coeff=Coeff+Probl 
			opt0=StartPoint[0:len(Lower),0]

			try:
				if Job.Parameters.MaxIter>0:
					Fit = least_squares(fun4LS,opt0,args=(X,Y,Nombres_coeff,PValues,funf,W_),method='trf',verbose=2,max_nfev=Parameters.MaxFunEvals,xtol=Parameters.TolX,ftol=Parameters.TolFun,bounds=(Lower,Upper))
				else:
					Fit = float("NaN")
			except:
				pass

			cntF=len(Fits['GOF']['SSE'])

			if cntF <= NumFits and (kkkk+1)<2:
				Fits['fitobj'].append(Fit)
				Fits['Time'].append(time()-ti)
				Fits['Start']['Coefficients'].append(StartPoint[:,0])

				if Job.Parameters.MaxIter>0:
					YFits=LS_Y_out(Fit.x,X,Nombres_coeff,PValues,funf)
					Fits['YFits'].append(YFits)
					vec_join_aux=np.append(Fit.x,PValues)
					Fits['Coefficients'].append(vec_join_aux)
					SSE=Fit.cost*2
					dfe=len(Fit.fun)-len(opt0)
					Fits['GOF']['SSE'].append(SSE)
					RMSE=np.sqrt(SSE/dfe)
					Fits['GOF']['RMSE'].append(RMSE)
					Fits['GOF']['rsquare'].append(1-(SSE/(sum(W_*(Y-np.mean(Y))**2))))
					Fits['GOF']['adjrsquare'].append(1-(SSE*(len(Fit.fun)-1))/((sum(W_*(Y-np.mean(Y))**2))*(len(Fit.fun)-len(opt0))))
					Level=Job.Parameters.CintL/100
					Fits['Level']=Level
					#No regresa valores parecidos a los de Matlab, sobre todo los últimos ya que te regresa dos órdenes de magnitud más chicas
					try:
						t_lim=t.ppf(1-(1-Level)/2,dfe)
						q,R =np.linalg.qr(Fit.jac)
						Rinv=linalg.inv(R)
						V=np.sum(Rinv**2,axis=1)*(SSE/dfe)
						db=t_lim*np.sqrt(V)
						lim_sup_coef=Fit.x+db
						lim_inf_coef=Fit.x-db
						try:
							pint=predint(Level,len(opt0),dfe,RMSE,YFits,Fit.jac,Rinv)
						except:
							pass
					except:
						lim_inf_coef=np.zeros(len(Fit.x))*float('NaN')
						lim_sup_coef=lim_inf_coef
						pint=[np.zeros(len(X))*float('NaN'),np.zeros(len(X))*float('NaN')]

					vec_join_aux=np.append(lim_inf_coef,PValues)
					Fits['CIntLowB'].append(vec_join_aux)
					vec_join_aux=np.append(lim_sup_coef,PValues)
					Fits['CIntUpB'].append(vec_join_aux)
					Fits['CIntLowLine'].append(pint[0][:])
					Fits['CIntUpLine'].append(pint[1][:])
				else:
					Stt=StartPoint[:,0]
					f=EvalExpFe(Data,CoeffNames_for_fun,Stt,funf)

					Fits['YFits'].append(f)
					mask = np.logical_not(ex)
					SqrErr = W_[mask]*(f[mask]-Y[mask])**2
					SSErr=sum(SqrErr)
					Fits['GOF']['SSE'].append(SSErr)
					Fits['Coefficients'].append(StartPoint[:,0])
					Fits['GOF']['rsquare'].append(float('NaN'))
					Fits['GOF']['adjsquare'].append(float('NaN'))
					Fits['GOF']['RMSE'].append(float('NaN'))
					Fits['CIntLowB'].append(np.zeros(len(StartPoint[:,0]))*float('NaN'))
					Fits['CIntUpB'].append(np.zeros(len(StartPoint[:,0]))*float('NaN'))
					Fits['CIntLowLine'].append(np.zeros(len(X))*float('NaN'))
					Fits['CIntUpLine'].append(np.zeros(len(X))*float('NaN'))

				Fits['Model']=Model

			else:
				SSE=Fits['GOF']['SSE'][:]
				tag = np.argsort(SSE)
				SSE.sort()
				tag=tag[len(SSE)-1]
				for ii in SSE:

					if Job.Parameters.MaxIter>0:
						if ii > Fit.cost*2:
							Fits['fitobj'][tag]=Fit
							Fits['Time'][tag]=time()-ti
							Fits['Start']['Coefficients'][tag][:]=StartPoint[:,0]
							YFits=LS_Y_out(Fit.x,X,Nombres_coeff,PValues,funf)
							Fits['YFits'][tag][:]=YFits
							Fits['Coefficients'][tag][:]=np.append(Fit.x,PValues)
							SSE4Fits=Fit.cost*2
							dfe=len(Fit.fun)-len(opt0)
							Fits['GOF']['rsquare'][tag]=1-(SSE4Fits/(sum(W_*(Y-np.mean(Y))**2)))
							Fits['GOF']['adjrsquare'][tag]=1-(SSE4Fits*(len(Fit.fun)-1))/((sum(W_*(Y-np.mean(Y))**2))*(dfe))
							RMSE=np.sqrt(SSE4Fits/dfe)
							Fits['GOF']['RMSE'][tag]=RMSE
							Fits['GOF']['SSE'][tag]=SSE4Fits
							Level=Job.Parameters.CintL/100
							try:
								t_lim = t.ppf(1 - (1 - Level) / 2, dfe)
								q, R = np.linalg.qr(Fit.jac)
								Rinv = linalg.inv(R)
								V = np.sum(Rinv ** 2, axis=1) * (SSE / dfe)
								db = t_lim * np.sqrt(V)
								lim_sup_coef = Fit.x + db
								lim_inf_coef = Fit.x - db
								try:
									pint = predint(Level, len(opt0), dfe, RMSE, YFits, Fit.jac, Rinv)
								except:
									pass
							except:
								lim_inf_coef = np.zeros(len(Fit.x)) * float('NaN')
								lim_sup_coef = lim_inf_coef
								pint = [np.zeros(len(X)) * float('NaN'), np.zeros(len(X)) * float('NaN')]


							Fits['CIntLowB'][tag][:] = np.append(lim_inf_coef,PValues)
							Fits['CIntUpB'][tag][:] = np.append(lim_sup_coef,PValues)
							Fits['CIntLowLine'][tag][:]=pint[0][:]
							Fits['CIntUpLine'][tag][:]=pint[1][:]

							break

						else:
							if ii>SSE_SP:
								Fits['fitobj'][tag]=Fit
								Fits['Time'][tag]=time()-ti
								Fits['Start']['Coefficients'][tag][:]=StartPoint[:,0]

								Stt=StartPoint[:,0]
								f=EvalExpFe(Data,CoeffNames_for_fun,Stt,funf)

								Fits['YFits'][tag][:]=f
								Fits['Coefficients'][tag][:]=StartPoint[:,0]
								Fits['GOF']['rsquare'][tag]=float('NaN')
								Fits['GOF']['adjsquare'][tag]=float('NaN')
								Fits['GOF']['RMSE'][tag]=float('NaN')
								Fits['GOF']['SSE'][tag]=SSErr

								Fits['CIntLowB'][tag][:]=np.zeros(len(StartPoint[:,0]))*float('NaN')
								Fits['CIntUpB'][tag][:]=np.zeros(len(StartPoint[:,0]))*float('NaN')
								Fits['CIntLowLine'][tag][:]=np.zeros(len(X))*float('NaN')
								Fits['CIntUpLine'][tag][:]=np.zeros(len(X))*float('NaN')

								break

			if not not Fits['Time']:
				tr=np.mean(Fits['Time'])*((comb-cnt+1)+(5-(kkkk+1))*comb)

				if tr<60:
					monitor_print('Estimated time to complete Job %s is %g seconds...'%(Name,round(tr)))
				elif tr<3600 and tr>=60:
					tr=tr/60
					monitor_print('Estimated time to complete Job %s is %g minutes...'%(Name,round(tr)))
				elif tr<86400 and tr>=3600:
					tr=tr/3600
					monitor_print('Estimated time to complete Job %s is %3.1f hours...'%(Name,tr))
				else:
					tr=tr/86400
					if tr<100:
						monitor_print('Estimated time to complete Job %s is %3.1f days...'%(Name,tr))
					else:
						monitor_print('Estimated time to complete Job %s is way too many days!...'%(Name))

			dictonary4save={'Fits':Fits,'Parameters':Parameters,'Model':Model,'Data':Data}

			savemat(current_dir+'/'+Name+'_FitData.mat',dictonary4save)
		
		try:

			FNameDef=Parameters.FileName
			StrFile=FNameDef+'_Cycle'+str(kkkk+1)

			dictonary4save={'Fits':Fits,'Parameters':Parameters,'Model':Model,'Data':Data}

			savemat(current_dir+'/'+StrFile+'_FitData.mat',dictonary4save)
		except:
			try:

				FNameDef = Job.Name
				StrFile = FNameDef + '_Cycle' + str(kkkk + 1)

				dictonary4save = {'Fits': Fits, 'Parameters': Parameters, 'Model': Model, 'Data': Data}

				savemat(current_dir+'/'+StrFile + '_FitData.mat', dictonary4save)
			except:
				monitor_print('Hubo un error al guardar el contenido del ciclo número '+str(kkkk+1))

	monitor_print(str(datetime.now())+'-Job'+Name+' is Complete')
