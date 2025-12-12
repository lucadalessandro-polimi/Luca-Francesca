#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path



'''
CSV = Path("Meshes/Delaunay/timing_adaptivity.csv")

def _to_num(s): 
    return pd.to_numeric(s, errors="coerce")

# ====== (già avevi) stima dell'esponente α su log–log ======
def fit_alpha(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.count_nonzero(m) < 2: 
        return np.nan
    A = np.vstack([np.log(x[m]), np.ones(np.count_nonzero(m))]).T
    a, _ = np.linalg.lstsq(A, np.log(y[m]), rcond=None)[0]
    return float(a)

# ====== (già avevi) scaling LS per tracce guida ======
def ls_scaled(base, y):
    """
    Least–squares scale: find a so that y ≈ a * base (component-wise),
    returning a*base. Works if y is a vector OR a scalar (will broadcast).
    """
    base = np.asarray(base, float)
    y_arr = np.asarray(y, float)

    # allow scalar y: broadcast to base's shape
    if y_arr.ndim == 0:
        y_arr = np.full_like(base, float(y_arr))

    m = np.isfinite(base) & np.isfinite(y_arr) & (base > 0)
    if np.count_nonzero(m) < 2:
        return np.full_like(base, np.nan, dtype=float)

    a = np.dot(y_arr[m], base[m]) / np.dot(base[m], base[m])
    return a * base


# ===================== NUOVE FUNZIONI HELPER =====================

def _agg_median(df, by):
    """Aggrega per 'by' tenendo la mediana del Time_ms e ordina per le stesse colonne."""
    g = (df.groupby(by, dropna=False)["Time_ms"].median()
            .reset_index().sort_values(by))
    return g

def _geom_mid(x):
    x = np.asarray(x, float)
    return np.sqrt(x[0] * x[-1])          # centro geometrico (adatto per assi log)

def _log_interp(x, y, xq):
    lx, ly, lq = np.log(x), np.log(y), np.log(xq)
    return np.exp(np.interp(lq, lx, ly))  # interpolazione in log-log

def _complexity_guides(ax, x, y_ref, letter="n"):
    xg = np.asarray(x, float)
    xg = xg[np.isfinite(xg) & (xg > 1)]
    if xg.size < 2 or not np.isfinite(y_ref):
        return

    s      = xg / xg[0]
    yN     = y_ref * s
    yNlogN = y_ref * s * (np.log(xg)/np.log(xg[0]))
    yN2    = y_ref * s**2

    ax.loglog(xg, yN,     "--", alpha=0.35, linewidth=1.0, label="_nolegend_")
    ax.loglog(xg, yNlogN, "--", alpha=0.35, linewidth=1.0, label="_nolegend_")
    ax.loglog(xg, yN2,    "--", alpha=0.35, linewidth=1.0, label="_nolegend_")

    # >>> metti le etichette “verso destra” (85% del range in scala log) <<<
    frac = 0.965
    xt = np.exp(np.log(xg[0]) + frac*(np.log(xg[-1]) - np.log(xg[0])))

    # valore delle curve in xt (interpolazione log–log semplice)
    def y_at(xarr, yarr, xq):
        lx, ly, lq = np.log(xarr), np.log(yarr), np.log(xq)
        return np.exp(np.interp(lq, lx, ly))

    ax.text(xt, y_at(xg, yN,     xt) * 1.10,  fr"$\mathcal{{O}}({letter})$",            ha="center", va="center", fontsize=11, alpha=0.6)
    ax.text(xt, y_at(xg, yNlogN, xt) * 0.92,  fr"$\mathcal{{O}}({letter}\log {letter})$",ha="center", va="center", fontsize=11, alpha=0.6)
    ax.text(xt, y_at(xg, yN2,    xt) * 0.78,  fr"$\mathcal{{O}}({letter}^2)$",          ha="center", va="center", fontsize=11, alpha=0.6)


def _fit_and_overlay(ax, x, y, label_base):
    """Fit log–log (y ~ x^α), plottando dati e retta di fit tratteggiata. Ritorna α."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    a = fit_alpha(x, y)
    # dati
    ax.loglog(x, y, marker="o", linestyle="-", label=f"{label_base}  (α≈{a:.2f})")
    # linea di fit y = k * x^a
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.count_nonzero(m) >= 2 and np.isfinite(a):
        lx, ly = np.log(x[m]), np.log(y[m])
        k = np.exp(ly.mean() - a*lx.mean())
        xfit = np.linspace(x[m].min(), x[m].max(), 200)
        yfit = k * xfit**a
        ax.loglog(xfit, yfit, linestyle="--", linewidth=1.2, alpha=0.9)
    return a

# ============================== MAIN ==============================

def main():
    df = pd.read_csv(CSV)

    # nomi colonne attese
    if "TimeElapsed(ms)" in df.columns:
        df = df.rename(columns={"TimeElapsed(ms)": "Time_ms"})
    for col in ["DataN", "Time_ms"]:
        if col not in df.columns:
            raise KeyError(f"Missing column '{col}' in CSV")

    # opzionali
    df["NumPoints"] = _to_num(df.get("NumPoints", np.nan))
    df["AreaScale"] = _to_num(df.get("AreaScale", np.nan))

    # tipi + pulizia
    df["DataN"]   = _to_num(df["DataN"])
    df["Time_ms"] = _to_num(df["Time_ms"])
    df = df.dropna(subset=["DataN", "Time_ms"]).copy()

    # per sicurezza, sostituiamo 0 con eps per scala log dei plot
    df_plot = df.copy()
    df_plot["Time_ms"] = df_plot["Time_ms"].replace(0, 1e-9)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))

    # -------- LEFT: Time vs NumPoints (curve per DataN) --------
    if df_plot["NumPoints"].notna().any():
        alphas = []
        for data_n, g in df_plot.groupby("DataN"):
            g = g.dropna(subset=["NumPoints"])
            if g.empty: 
                continue
            gg = _agg_median(g, by=["NumPoints"])  # mediana per NumPoints
            x = gg["NumPoints"].to_numpy()
            y = gg["Time_ms"].to_numpy()
            a = _fit_and_overlay(axL, x, y, label_base=f"DataN={int(data_n)}")
            alphas.append((f"{int(data_n)}", a))

        # guide di complessità (scalate su un riferimento semplice)
        n_all = np.sort(df_plot["NumPoints"].dropna().unique())
        if n_all.size >= 2:
            n0 = n_all[0]
            y0 = df_plot.loc[df_plot["NumPoints"] == n0, "Time_ms"].median()
            letter = "n"
            _complexity_guides(axL, n_all, y0, letter)

        axL.set_xlabel("NumPoints (mesh nodes)")
        axL.set_ylabel("Time (ms)")
        axL.set_title("Time vs NumPoints (per DataN)")
        axL.grid(True, which="both", ls=":", alpha=0.5)
        axL.legend(fontsize=9, loc="upper left")

    else:
        axL.axis("off")
        axL.text(0.5, 0.5, "NumPoints not present in CSV", ha="center", va="center")

    # ---- RIGHT: Time vs DataN (curve per AreaScale, se presente) ----
    if df_plot["AreaScale"].notna().any():
        agg = (df_plot.groupby(["AreaScale", "DataN"], dropna=False)
                        .agg(Time_ms=("Time_ms", "median"))
                        .reset_index().sort_values(["AreaScale", "DataN"]))
        alphas = []
        for scale, g in agg.groupby("AreaScale", sort=True):
            x = g["DataN"].to_numpy()
            y = g["Time_ms"].to_numpy()
            a = _fit_and_overlay(axR, x, y, label_base=f"scale={scale:g}")
            alphas.append((f"{scale:g}", a))

        x_all = np.sort(agg["DataN"].unique())
        y_med = agg.groupby("DataN")["Time_ms"].median().reindex(x_all).to_numpy()
        y_ref = np.nanmedian(y_med)
        letter = "N"
        _complexity_guides(axR, x_all, y_ref, letter)

        axR.set_title("Time vs DataN (curves per AreaScale)")
        axR.set_xlabel("DataN (data size)")
        axR.set_ylabel("Time (ms)")
        axR.grid(True, which="both", ls=":", alpha=0.5)
        axR.legend(fontsize=9, loc="upper left")


    else:
        # fallback: curve per NumPoints se AreaScale manca
        agg = (df_plot.groupby(["NumPoints", "DataN"], dropna=False)
                        .agg(Time_ms=("Time_ms", "median"))
                        .reset_index().sort_values(["NumPoints", "DataN"]))
        alphas = []
        for npts, g in agg.groupby("NumPoints", sort=True):
            x = g["DataN"].to_numpy()
            y = g["Time_ms"].to_numpy()
            a = _fit_and_overlay(axR, x, y, label_base=f"NumPts={int(npts)}")
            alphas.append((f"{int(npts)}", a))

        x_all = np.sort(agg["DataN"].unique())
        y_med = agg.groupby("DataN")["Time_ms"].median().reindex(x_all).to_numpy()
        y_ref = np.nanmedian(y_med)
        _complexity_guides(axR, x_all, y_ref, "N")

        axR.set_title("Time vs DataN (fallback: per NumPoints)")
        axR.set_xlabel("DataN (data size)")
        axR.set_ylabel("Time (ms)")
        axR.grid(True, which="both", ls=":", alpha=0.5)
        axR.legend(fontsize=9, loc="upper left")


    plt.tight_layout()
    out = CSV.with_name(CSV.stem + "_complexity.png")
    plt.savefig(out, dpi=150)
    print("Saved:", out)

if __name__ == "__main__":
    main()
'''



CSV = Path("Meshes/Delaunay/complexity_test.csv")
df = pd.read_csv(CSV)
df.columns = df.columns.str.strip().str.replace(" ", "")

# conversione numerica
for col in ["Lmin", "Lmax", "tol", " n_data", " n_iter", "n_nodes", "time_ms"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# filtriamo
df = df.dropna(subset=["n_data", "time_ms"])
df = df[df["time_ms"] > 0]

# === 1️⃣ Boxplot: tempo per n_data ===
plt.figure(figsize=(7,5))
sns.boxplot(data=df, x="n_data", y="time_ms", color="lightsteelblue", fliersize=3)
plt.yscale("log")
plt.xlabel("Number of data points N")
plt.ylabel("Computation time (ms)")
plt.title("Computation time vs N (log scale)")
plt.grid(True, which="both", ls=":", alpha=0.4)
plt.tight_layout()
plt.savefig("Meshes/Delaunay/time_vs_ndata_boxplot.png", dpi=150)
plt.show()

# === 2️⃣ Fit log–log medio + curve guida ===
agg = df.groupby("n_data", as_index=False)["time_ms"].median()
x = agg["n_data"].to_numpy()
y = agg["time_ms"].to_numpy()

lx, ly = np.log(x), np.log(y)
a, b = np.polyfit(lx, ly, 1)  # fit log–log
alpha = a
k = np.exp(b)
print(f"Fitted complexity exponent α ≈ {alpha:.2f}")

x_fit = np.linspace(x.min(), x.max(), 100)
y_fit = k * x_fit ** alpha

plt.figure(figsize=(7,5))
plt.loglog(x, y, "o", label="Median data", color="royalblue")
plt.loglog(x_fit, y_fit, "k--", color="royalblue", lw=2, label=f"Fit: O(N^{alpha:.2f})")

# Curve guida normalizzate
x_ref = np.array([x.min(), x.max()])
y_ref = np.array([y.min(), y.max()])
s = x_fit / x_fit[0]
for exp, ls in zip([1, 1.0*np.log(x_fit)/np.log(x_fit[0]), 2], ["-", "--", ":"]):
    pass  # non serve, le curve guida le mettiamo sotto

# Curve teoriche
base = y_fit[0]
plt.loglog(x_fit, base * (x_fit/x_fit[0]), "--", color="gray", alpha=0.6, label="O(N)")
plt.loglog(x_fit, base * (x_fit/x_fit[0])*np.log(x_fit/x_fit[0]), "--", color="gray", alpha=0.6, label="O(N log N)")
plt.loglog(x_fit, base * (x_fit/x_fit[0])**2, "--", color="gray", alpha=0.6, label="O(N²)")

plt.xlabel("Number of data points N")
plt.ylabel("Computation time (ms)")
plt.title("Adaptivity cicle: computation time vs. # data points")
plt.grid(True, which="both", ls=":", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.savefig("Meshes/Delaunay/time_vs_ndata_fit.png", dpi=150)
plt.show()



#---------------------------------------------------------------------------

df = pd.read_csv("Meshes/Delaunay/complexity_test.csv")

# Cleaning and numeric conversion
for col in ["n_nodes", " n_iter", "time_ms"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=["n_nodes", "time_ms"])
df = df[(df["n_nodes"] > 0) & (df["time_ms"] > 0)]

# === SELECT X (n_cells or n_iter) ===
use_iters = False     # <- cambia a True se vuoi niter sull'asse x
x = df[" n_iter"] if use_iters else df["n_nodes"]
x_label = "Number of iterations" if use_iters else "Number of nodes"

y = df["time_ms"]

# === FIT legge di potenza y = C * x^alpha ===
logx, logy = np.log(x), np.log(y)
alpha, logC = np.polyfit(logx, logy, 1)
C = np.exp(logC)
print(f"Fitted complexity exponent α ≈ {alpha:.2f}")

# === Prepare line for fit and theoretical guides ===
x_fit = np.linspace(x.min(), x.max(), 200)
y_fit = C * x_fit**alpha

# Reference base for guide lines
ref_y = y.min()
ref_x = x.min()

# === PLOT ===
plt.figure(figsize=(8, 6))
plt.loglog(x, y, "o", alpha=0.7, color="royalblue", label="Data (samples)")
plt.loglog(x_fit, y_fit, "--", color="black", lw=2.0, label=f"Fit: O(n^{alpha:.2f})")

# Theoretical guide lines
plt.loglog(x_fit, ref_y * (x_fit / ref_x), "--", color="gray", alpha=0.6, label="O(n)")
plt.loglog(x_fit, ref_y * (x_fit / ref_x)**1.1, "--", color="gray", alpha=0.6, label="~O(n log n)")
plt.loglog(x_fit, ref_y * (x_fit / ref_x)**2, "--", color="gray", alpha=0.6, label="O(n²)")

plt.xlabel(x_label)
plt.ylabel("Computation time (ms)")
plt.title("Adaptivity cycle: computation time vs. # nodes")
plt.legend()
plt.grid(True, which="both", ls=":", alpha=0.4)
plt.savefig("Meshes/Delaunay/time_vs_nodes.png", dpi=150)
plt.tight_layout()
plt.show()