import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

def plot_timing_data():
    folder = os.path.dirname(__file__)
    file_paths = glob.glob(os.path.join(folder, 'timing_*.csv'))

    if not file_paths:
        print("No timing_*.csv files found.")
        return

    for file_path in file_paths:
        df = pd.read_csv(file_path)

        # Estrai nome file per usarlo nel nome dell'immagine
        filename = os.path.splitext(os.path.basename(file_path))[0]
        suffix = filename.removeprefix("timing_")

        # Estrai dati
        num_points = df['NumPoints']
        times = df['TimeElapsed(ms)']

        # Conversione in array per i modelli teorici
        n = np.array(num_points)
        t = np.array(times)

        # Calcolo modelli teorici (normalizzati sul primo punto)
        n_min = n[0]
        time_min = t[0]

        scaling_factor_nlogn = time_min / (n_min * np.log(n_min))
        scaling_factor_n = time_min / n_min
        scaling_factor_n2 = time_min / (n_min ** 2)
        scaling_factor_n2logn = time_min / (n_min ** 2 * np.log(n_min))

        # Plot
        plt.figure(figsize=(8, 6))
        plt.loglog(n, t, '-', label='Measured time')
        plt.loglog(n, scaling_factor_n * n, '--', label=r'$\mathcal{O}(n)$')
        plt.loglog(n, scaling_factor_nlogn * n * np.log(n), '--', label=r'$\mathcal{O}(n \log n)$')
        plt.loglog(n, scaling_factor_n2 * n**2, '--', label=r'$\mathcal{O}(n^2)$')
        plt.loglog(n, scaling_factor_n2logn * n**2 * np.log(n), '--', label=r'$\mathcal{O}(n^2 \log n)$')

        plt.title(f'Computational cost: {suffix}')
        plt.xlabel('Number of Points')
        plt.ylabel('Time (ms)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        # Salvataggio del grafico con nome file specifico
        output_path = os.path.join(folder, f"{filename}.png")
        plt.savefig(output_path, dpi=300)
        plt.close()

        print(f"Saved plot for {filename} : {output_path}")


if __name__ == "__main__":
    plot_timing_data()

