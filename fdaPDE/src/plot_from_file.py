import numpy as np
import matplotlib.pyplot as plt

def plot_hexagon_from_file(filename="hexagon.txt"):
    try:
        # Carichiamo i dati dal file
        data = np.loadtxt(filename)

        # Aggiungiamo il primo punto alla fine per chiudere il poligono
        data = np.vstack([data, data[0]])

        # Disegniamo l'esagono
        plt.figure(figsize=(5, 5))
        plt.plot(data[:, 0], data[:, 1], 'bo-', label="Esagono")
        plt.fill(data[:, 0], data[:, 1], 'skyblue', alpha=0.3)
        plt.scatter(data[:-1, 0], data[:-1, 1], color='red', label="Vertici")

        # Etichettiamo i vertici
        for i, (x, y) in enumerate(data[:-1]):
            plt.text(x, y, str(i), fontsize=12, ha='right', va='bottom', color='black')

        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)
        plt.gca().set_aspect('equal')
        plt.legend()
        plt.title("Esagono generato da make_polygon")
        plt.grid(True)

        # Mostriamo la figura
        plt.savefig("hexagon_plot.png")  # Salva l'immagine in un file
        print("Grafico salvato come 'hexagon_plot.png'")


    except Exception as e:
        print(f"Errore nella lettura del file: {e}")

if __name__ == "__main__":
    plot_hexagon_from_file()
