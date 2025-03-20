import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_timing_data():
    # Cerca il file nella stessa cartella dello script
    file_path = os.path.join(os.path.dirname(__file__), 'timing_results.csv')
    
    # Controllo se il file esiste
    if not os.path.exists(file_path):
        print(f"❌ Il file {file_path} non esiste! Assicurati che sia stato generato correttamente.")
        return

    # Leggi i dati dal file CSV
    df = pd.read_csv(file_path)
    
    # Estrai i dati
    num_points = df['NumPoints']
    times = df['TimeElapsed(ms)']   # Convertiamo i millisecondi in secondi
    
    # Crea un grafico log-log
    plt.figure(figsize=(8, 6))
    plt.loglog(num_points, times, label='Tempo computazionale', marker='o')
    
    # Aggiungi curve di complessità teorica per confronto
    n = np.array(num_points)
    plt.loglog(n, n * np.log(n), label=r'$\mathcal{O}(n \log n)$', linestyle='dashed')
    plt.loglog(n, n, label=r'$\mathcal{O}(n)$', linestyle='dashed')
    plt.loglog(n, n**2, label=r'$\mathcal{O}(n^2)$', linestyle='dashed')

    # Dettagli grafico
    plt.title('Tempo Computazionale vs Numero di Punti')
    plt.xlabel('Numero di Punti')
    plt.ylabel('Tempo (secondi)')
    plt.legend()
    plt.grid(True)
    plt.savefig("timing_plot.png", dpi=300)


if __name__ == "__main__":
    plot_timing_data()

