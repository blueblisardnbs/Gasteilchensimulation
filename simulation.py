import numpy as np
import scipy.spatial.distance as distance

groesse_zeitschritt = 0.005

anzahl_zeitschritte = 10000
anzahl_teilchen = 5000
anzahl_dimensionen = 2

m_1 = 1
m_2 = 1

radius_teilchen = 0.05

seitenlänge_Käfig = 100

np.savetxt("Konfiguration.txt", (anzahl_teilchen, anzahl_zeitschritte))

def generate_random_initial_pos(num_particles, num_dimensions, radius_particles):

    box_left_wall  = [0.0,]*num_dimensions
    box_right_wall = [seitenlänge_Käfig,]*num_dimensions
    
    pos = np.zeros([num_particles, num_dimensions])
        
    md = 2.5 * radius_particles
    
    coord_arrays = [ ]
    for d in range(num_dimensions):
        l = box_left_wall[d] + md
        r = box_right_wall[d] - md
        num_sites = min( int((r - l) / md), 100)
        
        coord_arrays.append( np.linspace(l, r, num_sites) )
    
    coords = np.meshgrid(*coord_arrays)
    n = coords[0].size
    indices = np.random.choice(n, num_particles, replace=False)
    
    for d in range(num_dimensions):
        pos[:,d] = np.ravel(coords[d])[indices]
                
    return pos

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
            if x[t, d] > seitenlänge_Käfig or x[t, d] < 0:
                v_neu[t,d] *= -1
                
    distanzen = distance.pdist(x)
    q = distance.squareform(distanzen)

    t1, t2 = np.where(q < 2*radius_teilchen)

    for tt1, tt2 in zip(t1, t2):
        if tt1 < tt2:
            v_neu[tt1, :], v_neu[tt2, :] = stoß(x, v, tt1, tt2)
            
    x_neu = x + v_neu * dt

    return x_neu, v_neu

anfangspositionen = generate_random_initial_pos(anzahl_teilchen,anzahl_dimensionen, radius_teilchen)

distanzen = distance.pdist(anfangspositionen)


for position1 in range(len(anfangspositionen)):
    for position2 in range(len(anfangspositionen)):
        if position1 == position2:
            continue
        while np.linalg.norm(anfangspositionen[position1] - anfangspositionen[position2]) < 2*radius_teilchen:
            anfangspositionen[position1] = 0.9*np.random.random(anzahl_dimensionen)

anfangsgeschwindigkeiten = 2*np.random.random((anzahl_teilchen, anzahl_dimensionen))-1
orte = np.zeros([anzahl_zeitschritte, anzahl_teilchen, anzahl_dimensionen])
geschwindigkeiten = np.zeros([anzahl_zeitschritte, anzahl_teilchen, anzahl_dimensionen])


orte[0] = anfangspositionen
geschwindigkeiten[0] = anfangsgeschwindigkeiten

zeiten = np.array([schritt * groesse_zeitschritt for schritt in range(anzahl_zeitschritte)])
for schritt in range(1, anzahl_zeitschritte):
    orte[schritt], geschwindigkeiten[schritt] = zeitschritt(orte[schritt-1], geschwindigkeiten[schritt-1], groesse_zeitschritt)
    if (schritt+1)%1000 == 0:
        print(f"{schritt+1}/{anzahl_zeitschritte}\t:\t{round(schritt/anzahl_zeitschritte/10,3)*1000}%")

header = "Teilchenpositionen"
for teilchen in range(anzahl_teilchen):
    np.savetxt("Ort_Teilchen" + str(teilchen + 1) + ".txt", np.column_stack([zeiten, orte[:, teilchen, :]]), header=header)
    np.savetxt("Geschwindigkeit_Teilchen" + str(teilchen + 1) + ".txt", np.column_stack([geschwindigkeiten[:, teilchen, :]]), header=header)
    if teilchen%10 == 0: print(f"s\t:\t{round((teilchen/anzahl_teilchen)*100,1)}%")

print("\nFertig")
