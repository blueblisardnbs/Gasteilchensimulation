import numpy as np

def wand(teilchen_pos, teilchen_v, seitenlänge_Käfig, anzahl_dimensionen):
    kharm = 1e3
    E_pot_w = 0
    acc = np.zeros(anzahl_dimensionen)
    for d in range(anzahl_dimensionen):
        if teilchen_pos[d] < 0: 
            #print("bpoing")
            abstand_wand = abs(teilchen_pos[d])
            acc[d] = + kharm * abstand_wand
            E_pot_w = 1/2 * kharm * abstand_wand**2
        elif teilchen_pos[d] > seitenlänge_Käfig:
            abstand_wand = abs(seitenlänge_Käfig-teilchen_pos[d])
            acc[d] = - kharm * abstand_wand
            E_pot_w = 1/2 * kharm * abstand_wand**2
        
    return acc, E_pot_w