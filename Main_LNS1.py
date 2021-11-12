from BGridFe import BGridFe
from scipy.io import loadmat
import os

#Primera solución se obtiene el directorio actual
directorio_actual= os.path.dirname(os.path.realpath(__file__))

#No olvidar añadir .mat al final del nombre y entre comillas
#segunda solución se obtiene el path absoluto al archivo en automático
Batch_name=directorio_actual+"/Ejemplo_Batch_Cluster_3.mat"
JobCell=loadmat(Batch_name,struct_as_record=False,squeeze_me=True)
BGridFe(JobCell)