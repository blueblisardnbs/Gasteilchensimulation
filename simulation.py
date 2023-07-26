import numpy as np
import math

groesse_zeitschritt = 0.001

anzahl_zeitschritte = 10000
anzahl_teilchen = 10
anzahl_dimensionen = 2

m_1 = 1
m_2 = 1

radius_teilchen = 0.005

def abstand(teilchen1, teilchen2):
    
    abstand_dimension = [0 for i in range(anzahl_dimensionen)]
    for d in range(anzahl_dimensionen):
        abstand_dimension[d] += orte[schritt-1, teilchen1, d]
        abstand_dimension[d] -= orte[schritt-1, teilchen2, d]
        
        abstand_dimension[d] = pow(abstand_dimension[d], 2)
    
    return math.sqrt(sum(abstand_dimension))
    
def stoß(x, v, teilchen1, teilchen2):
    
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
    
def zeitschritt(x, v, dt):

    v_neu = np.copy(v)
    for t in range(anzahl_teilchen):
        for d in range(anzahl_dimensionen):
            if x[t, d] > 1 or x[t, d] < 0:
                v_neu[t,d] *= -1

    for teilchen1 in range(anzahl_teilchen):
        for teilchen2 in range(anzahl_teilchen):
            if teilchen1 <= teilchen2:
                continue
            if abstand(teilchen1, teilchen2) < 2 * radius_teilchen:
                v_neu[teilchen1, :], v_neu[teilchen2, :] = stoß(x, v, teilchen1, teilchen2)
    
    x_neu = x + v_neu * dt

    return x_neu, v_neu

anfangspositionen = np.random.random((anzahl_teilchen, anzahl_dimensionen))
for position1 in range(len(anfangspositionen)):
    for position2 in range(len(anfangspositionen)):
        if position1 == position2:
            continue
        if np.linalg.norm(anfangspositionen[position1] - anfangspositionen[position2]) < 2*radius_teilchen:
            anfangspositionen[position1] = np.random.random(anzahl_dimensionen)
            position1 -= 1


anfangsgeschwindigkeiten = 2*np.random.random((anzahl_teilchen, anzahl_dimensionen))-1
orte = 1.8*np.zeros([anzahl_zeitschritte, anzahl_teilchen, anzahl_dimensionen])-0.9
geschwindigkeiten = np.zeros([anzahl_zeitschritte, anzahl_teilchen, anzahl_dimensionen])


orte[0] = anfangspositionen
geschwindigkeiten[0] = anfangsgeschwindigkeiten

zeiten = np.array([schritt * groesse_zeitschritt for schritt in range(anzahl_zeitschritte)])
for schritt in range(1, anzahl_zeitschritte):
    orte[schritt], geschwindigkeiten[schritt] = zeitschritt(orte[schritt-1], geschwindigkeiten[schritt-1], groesse_zeitschritt)
    print(f"{schritt+1}/{anzahl_zeitschritte}\t:\t{round(schritt/anzahl_zeitschritte/10,3)*1000}%")

header = "Teilchenpositionen"
for teilchen in range(anzahl_teilchen):
    np.savetxt("Ort_Teilchen" + str(teilchen + 1) + ".txt", np.column_stack([zeiten, orte[:, teilchen, :]]), header=header)
    np.savetxt("Geschwindigkeit_Teilchen" + str(teilchen + 1) + ".txt", np.column_stack([geschwindigkeiten[:, teilchen, :]]), header=header)
