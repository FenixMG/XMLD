from BGridFe import BGridFe
from scipy.io import loadmat
import os

#Primera solución se obtiene el directorio actual
directorio_actual= os.path.dirname(os.path.realpath(__file__))

#No olvidar añadir .mat al final del nombre y entre comillas
#segunda solución se obtiene el path absoluto al archivo en automático
nombre_base="Fe3O4_L30_D3d_r1_Clu_"
l=0
for i in range(1000):
    l+=1
    try:
        try:
            Batch_name = directorio_actual + "/" + nombre_base + str(l) + "_FitData.mat"
            JobCell = loadmat(Batch_name, struct_as_record=False, squeeze_me=True)
            print('Se encontró un FitData.mat se reanudará el cáuculo desde donde se quedó')
            print('Batch cargado correctamente')
        except:
            Batch_name = directorio_actual + "/" + nombre_base + str(l) + ".mat"
            JobCell=loadmat(Batch_name,struct_as_record=False,squeeze_me=True)
            print('Ningún FitData.mat fue encontrado, se procede a calcular desde el inicio')
            print('Batch cargado correctamente')
        break
    except:
        pass
    if i+1==1000:
        print('después de 1000 ciclos el archivo no fue encontrado, revise que el nombre_base esté bien escrito, por ejemplo para Fe3O4_D3d_r3_Cluster_1.mat el nombre base es Fe3O4_D3d_r3_Cluster_')

BGridFe(JobCell)