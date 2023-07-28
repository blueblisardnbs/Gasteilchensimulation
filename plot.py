import numpy as np
import matplotlib.pyplot as plt

anzahl_schritte = 10000
anzahl_teilchen = 1000

anzahl_teilchen, anzahl_schritte = [int(i) for i in np.loadtxt("Konfiguration.txt")]
#anzahl_teilchen = int(np.loadtxt("Konfiguration.txt")[0])

v_max = 5

orte_Teilchen = []
geschwindigkeiten_Teilchen = []

def ort_Anzeige():
    plt.figure()

    for teilchen in range(anzahl_teilchen):
        plt.plot([t[1] for t in orte_Teilchen[teilchen]], [t[2] for t in orte_Teilchen[teilchen]], label=f'Teilchen {teilchen+1}')

    plt.xlabel('Position x')
    plt.ylabel('Position y')
    plt.title('Bahnkurven der Teilchen')
    plt.show()

def energie_Anzeige():
    

    plt.figure()
    
    E_t = []
    
    for schritt in range((anzahl_schritte)):
        E_ges = 0

        for v_t in geschwindigkeiten_Teilchen:
            E_ges += np.linalg.norm(v_t[schritt])**2

        E_ges *= 0.5
        E_t.append(E_ges)

    plt.plot([schritt for schritt in range(anzahl_schritte)], E_t)
    
    plt.xlabel('Schritt')
    plt.ylabel('Energie E')
    plt.title('Gesamtenergie')
    plt.show()


def geschwindigkeit_Anzeige(schritt):

    plt.figure()
    
    geschwindigkeit_Intervalle = np.zeros(20)
    
    for intervall in range(20):
        for teilchen in range(len(geschwindigkeiten_Teilchen)):
            #print(np.linalg.norm(geschwindigkeiten_Teilchen[teilchen, schritt]))
            if np.linalg.norm(geschwindigkeiten_Teilchen[teilchen][schritt]) > v_max*(intervall/20-1/20) and np.linalg.norm(geschwindigkeiten_Teilchen[teilchen][schritt]) < v_max*intervall/20:
                geschwindigkeit_Intervalle[intervall] += 1
                
    
    #print(geschwindigkeit_Intervalle)
    
    plt.plot([round(i/10+1/10*v_max,1) for i in range(10)], [geschwindigkeit_Intervalle[i] for i in range(10)])
    
    #print(geschwindigkeit_Intervalle)
    
    plt.xlabel('Geschwindigkeit v')
    plt.ylabel('Häufigkeit n')
    plt.title('Geschwindigkeitsverteilung')
    plt.show()



for teilchen in range(anzahl_teilchen):
    orte_Teilchen.append(np.loadtxt(f"Ort_Teilchen{teilchen+1}.txt"))
    geschwindigkeiten_Teilchen.append(np.loadtxt(f"Geschwindigkeit_Teilchen{teilchen+1}.txt"))

#print(geschwindigkeiten_Teilchen)


geschwindigkeit_Anzeige(100)

ort_Anzeige()

energie_Anzeige()