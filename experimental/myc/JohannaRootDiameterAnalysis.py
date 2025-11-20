import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import openpyxl

df = pd.read_excel("Analysis_Root_Morphology_Johanna.xlsx")
# print(df)
result = df.loc[
    ~df["File.Name"].str.contains("Pick", case=False, na=False),
    df.columns.str.contains("Length.Diameter", case=False) + df.columns.str.contains("File.Name", case=False)
]
result = result.loc[
    ~result["File.Name"].str.contains("S", case=False, na=False)
]
# print(result)

# Werte der Root.Diameter-Spalten extrahieren
root_cols = result.columns[result.columns.str.contains("Length.Diameter")]

# x-Achse = Durchmesser-Bins
x = np.arange(len(root_cols))

# =========================================================
# 1. Plot – erste 5 Dateien speichern
# =========================================================
first5 = result.head(5)

plt.figure(figsize=(10, 5))
for idx, row in first5.iterrows():
    values = row[root_cols].astype(float).values
    plt.plot(x, values, alpha=0.7, label=row["File.Name"])

plt.xticks(x, root_cols, rotation=45, ha="right")
plt.title("RMC (erste 5 Files) – Root Diameter Distribution")
plt.xlabel("Root Diameter Klassen")
plt.ylabel("Root Length")
plt.legend()
plt.tight_layout()

plt.savefig("results/rmc_root_distribution.png", dpi=200)
plt.close()


# =========================================================
# 2. Plot – letzte 5 Dateien speichern
# =========================================================
last5 = result.tail(5)

plt.figure(figsize=(10, 5))
for idx, row in last5.iterrows():
    values = row[root_cols].astype(float).values
    plt.plot(x, values, alpha=0.7, label=row["File.Name"])

plt.xticks(x, root_cols, rotation=45, ha="right")
plt.title("WT (letzte 5 Files) – Root Diameter Distribution")
plt.xlabel("Root Diameter Klassen")
plt.ylabel("Root Length")
plt.legend()
plt.tight_layout()

plt.savefig("results/wt_root_distribution.png", dpi=200)
plt.close()


# =========================================================
# 3. Mittelwert + Standardabweichung berechnen pro Klasse
# =========================================================

# für RMC
rmc_mean = first5[root_cols].astype(float).mean()
rmc_std  = first5[root_cols].astype(float).std()

# für WT
wt_mean = last5[root_cols].astype(float).mean()
wt_std  = last5[root_cols].astype(float).std()

# alles in ein schönes DataFrame
summary = pd.DataFrame({
    "RMC_Mean": rmc_mean,
    "RMC_SD": rmc_std,
    "WT_Mean": wt_mean,
    "WT_SD": wt_std
})

print(summary)
