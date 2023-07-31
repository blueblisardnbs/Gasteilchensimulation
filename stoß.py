import numpy as np
import scipy.spatial.distance as distance
import math
#from simulation import seitenlänge_Käfig

epsilon = 1
sigma = 1

def harter_Stoß(x, v, teilchen1, teilchen2, anzahl_dimensionen, m_1, m_2):
    
    v_1 = v[teilchen1, :]
    v_2 = v[teilchen2, :]
        
    v_rel = v_2 - v_1
    
    e_rel = (x[teilchen2,:] - x[teilchen1,:])
    e_rel /= np.linalg.norm(e_rel)
    
    v_rel_parallel_skalarprod = v_rel @ e_rel
    
    v_rel_parallel = v_rel_parallel_skalarprod * e_rel
    
    v_rel_neu = v_rel - 2 * v_rel_parallel
    
    v_cm = np.zeros(anzahl_dimensionen)
    
    for d in range(anzahl_dimensionen):
        v_cm[d] = (m_1*v_1[d] + m_2*v_2[d])/(m_1+m_2)
    
    v_1_neu = np.zeros(anzahl_dimensionen)
    v_2_neu = np.zeros(anzahl_dimensionen)
    
    for d in range(anzahl_dimensionen):
        v_1_neu[d] = v_cm[d] - (m_2/(m_1 + m_2) *v_rel_neu[d])
        v_2_neu[d] = v_cm[d] + (m_1/(m_1 + m_2) *v_rel_neu[d])
    
    return v_1_neu, v_2_neu

def wechselwirkungen(teilchen1, teilchen2, teilchen1_alt, teilchen2_alt, v1, v2, groesse_zeitschritt, seitenlänge_Käfig):
        
    r = teilchen2 - teilchen1

    e_r = r/np.linalg.norm(r)
    a = epsilon * (12 * pow(sigma, 12)/pow(np.linalg.norm(r), 13) - 6 * pow(sigma, 6)/pow(np.linalg.norm(r), 7)) * e_r #Masse immer 1
    Epot = epsilon * (pow(sigma, 12)/pow(np.linalg.norm(r),12) - pow(sigma,6)/pow(np.linalg.norm(r), 6))
    
    v2_neu = v2 + a * groesse_zeitschritt
    v1_neu = v1 - a * groesse_zeitschritt
    
    return v1_neu, v2_neu, Epot