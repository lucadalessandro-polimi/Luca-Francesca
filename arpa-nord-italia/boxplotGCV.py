import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from pathlib import Path

# =========================
# CONFIGURAZIONE
# =========================

files = {
    "Uniform mesh": "kFold_data/errors_mesh1.csv",
    "Data-driven mesh": "kFold_data/errors_mesh2.csv",
    "Residual-based mesh": "kFold_data/errors_mesh3.2.csv",
}

output_path = "4175_boxplot.png"
dpi = 300

method_colors = ["#74cb77", "#de77ae", "#6baed6"]

# =========================
# LETTURA DATI
# =========================

labels = list(files.keys())
data = []

for label in labels:
    df = pd.read_csv(files[label])
    values = pd.to_numeric(df.iloc[:, 0], errors="coerce").dropna().values
    values = np.sqrt(values)   # se i file contengono MSE
    data.append(values)

uniform = data[0]
datadriven = data[1]
residual = data[2]

# =========================
# TEST STATISTICI (paired)
# =========================

alpha = 0.05

p_dd = wilcoxon(uniform, datadriven, alternative='greater').pvalue
p_res = wilcoxon(uniform, residual, alternative='greater').pvalue


print(f"Uniform vs Data-driven p-value: {p_dd:.4f}")
print(f"Uniform vs Residual p-value: {p_res:.4f}")

# =========================
# BOXPLOT
# =========================

fig, ax = plt.subplots(figsize=(8, 7))

bp = ax.boxplot(
    data,
    patch_artist=True,
    widths=0.5,
    showfliers=True,
    medianprops=dict(color="black", linewidth=1.2),
)

for box, color in zip(bp["boxes"], method_colors):
    box.set_facecolor(color)

ax.set_ylim(3.5, 11)
ax.set_xticklabels(labels, fontsize=19)
ax.set_ylabel("RMSE", fontsize=19)
ax.set_title("Cross-validated RMSE under Different Mesh Discretizations", fontsize=22)

ax.grid(axis="y", linestyle="--", linewidth=0.7, alpha=0.6)

# =========================
# AGGIUNTA ASTERISCHI
# =========================

y_max = max([np.max(d) for d in data])
offset = 0.3

def add_star(x_pos, y_pos):
    ax.text(x_pos, y_pos, "*", ha="center", va="bottom", fontsize=22, color="black")

# box positions: 1, 2, 3
if p_dd < alpha and np.mean(datadriven) < np.mean(uniform):
    add_star(2, y_max + offset)

if p_res < alpha and np.mean(residual) < np.mean(uniform):
    add_star(3, y_max + offset)

plt.tight_layout()
plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
print(f"Figura salvata in: {Path(output_path).resolve()}")

plt.close()