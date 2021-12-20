import numpy as np
from scipy.special import erf

def Edge(X,I1,O1,G1,W1):
    E = I1*((1-G1)*(1/np.pi*np.arctan((X-O1)/W1)+1/2)+G1/2*(1+erf((X-O1)/(W1/(np.sqrt(2*np.log(2)))*np.sqrt(2)))))
    return E

#[*] ¿?Falta mandar a llamar a la librería en la función que llama a las funciones de FeL30
#[*] Faltan pruebas de ctmL30
#[ ] Probar que funcione con BGrid igual que matlab junto, ver en dónde y cómo se llama Edge y ver que se llame y funcione igual