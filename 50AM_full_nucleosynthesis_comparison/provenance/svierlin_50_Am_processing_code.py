# Read-only extracted notebook code: /scratch/svierlin/50_Am_processing.ipynb
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

csv_files = glob.glob("incremental_xtimes/*.csv")
print(f"Found {len(csv_files)} CSV files.")

peak_temps = []
last_na24 = []

for i, fpath in enumerate(csv_files):
    try:
        df = pd.read_csv(fpath, usecols=['temp_GK', 'Na24'])
        peak_temps.append(df['temp_GK'].max())
        last_na24.append(df['Na24'].iloc[-1])
    except Exception as e:
        print(f"Error on {fpath}: {e}")

    if i % 10000 == 0:
        print(f"  {i}/{len(csv_files)} processed...")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(peak_temps, last_na24, s=1, alpha=0.3, color='steelblue')
ax.set_xlabel("Peak Temperature (GK)")
ax.set_ylabel("Final $^{24}$Na")
ax.set_title("Final $^{24}$Na vs. Peak Temperature")
ax.set_yscale('log')
plt.tight_layout()
plt.savefig("Na24_vs_peaktemp.png", dpi=150)
plt.show()
print("Saved to Na24_vs_peaktemp.png")

