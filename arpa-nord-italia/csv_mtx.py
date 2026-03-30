import pandas as pd

def locs_csv_to_mtx(csv_path: str, mtx_path: str, x_col: str = "V1", y_col: str = "V2") -> None:
    """
    Convert a locs.csv (R write.csv style: "",X,Y with quoted strings) into MatrixMarket
    coordinate format compatible with your C++ read_locs_mtx (rows x 2, nnz=2*rows).
    """
    df = pd.read_csv(csv_path)

    if x_col not in df.columns or y_col not in df.columns:
        raise ValueError(f"Expected columns '{x_col}' and '{y_col}'. Found: {list(df.columns)}")

    # Convert to numeric (handles strings like " 7.6782")
    x = pd.to_numeric(df[x_col], errors="coerce")
    y = pd.to_numeric(df[y_col], errors="coerce")

    # Safety: drop rows with NaN coords (if any)
    mask = x.notna() & y.notna()
    x = x[mask].to_numpy()
    y = y[mask].to_numpy()

    n = len(x)
    nnz = 2 * n

    with open(mtx_path, "w", encoding="utf-8") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        # optional comment lines are allowed; your reader skips lines starting with '%'
        f.write("% generated from locs.csv (X,Y)\n")
        f.write(f"{n} 2 {nnz}\n")

        # Column 1: X values
        for i, xi in enumerate(x, start=1):   # MatrixMarket is 1-based
            f.write(f"{i} 1 {xi:.16g}\n")

        # Column 2: Y values
        for i, yi in enumerate(y, start=1):
            f.write(f"{i} 2 {yi:.16g}\n")

if __name__ == "__main__":
    locs_csv_to_mtx("data/2018.12.11_space_locs.csv", "../Meshes/Delaunay/2018.12.11_space_locs.mtx")

