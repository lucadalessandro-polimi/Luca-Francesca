import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_timing_data():
    # Cerca il file nella stessa cartella dello script
    file_path = os.path.join(os.path.dirname(__file__), 'timing_results.csv')


    # Leggi i dati dal file CSV
    df = pd.read_csv(file_path)
    
    # Estrai i dati
    num_points = df['NumPoints']
    times = df['TimeElapsed(ms)']  
    
    # Crea un grafico log-log
    plt.figure(figsize=(8, 6))
    plt.loglog(num_points, times, label='Tempo computazionale', marker='o')
    
    # Aggiungi curve di complessità teorica per confronto
    n = np.array(num_points)
    times = np.array(times)

    # Punto di riferimento per il riscalamento (il primo punto di times)
    n_min = n[0]
    time_min = times[0]

    # Scalare tutte le curve teoriche per farle partire da (n_min, time_min)
    scaling_factor_nlogn = time_min / (n_min * np.log(n_min))
    scaling_factor_n = time_min / n_min
    scaling_factor_n2 = time_min / (n_min ** 2)
    scaling_factor_n3 = time_min / (n_min ** 2*np.log(n_min))

    plt.figure(figsize=(8, 6))
    plt.loglog(n, scaling_factor_nlogn * n * np.log(n), label=r'$\mathcal{O}(n \log n)$', linestyle='dashed')
    plt.loglog(n, scaling_factor_n * n, label=r'$\mathcal{O}(n)$', linestyle='dashed')
    plt.loglog(n, scaling_factor_n2 * n**2, label=r'$\mathcal{O}(n^2)$', linestyle='dashed')
    plt.loglog(n, scaling_factor_n3 * n**2*np.log(n), label=r'$\mathcal{O}(n^2 \log n)$', linestyle='dashed')
    plt.loglog(n, times, 'b-', label='Times')
    

    # Dettagli grafico
    plt.title('Tempo Computazionale vs Numero di Punti')
    plt.xlabel('Numero di Punti')
    plt.ylabel('Tempo (millisecondi)')
    plt.legend()
    plt.grid(True)
    plt.savefig("timing_plot.png", dpi=300)


if __name__ == "__main__":
    plot_timing_data()

