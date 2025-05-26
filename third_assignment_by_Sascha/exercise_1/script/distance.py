import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

# Read the file as text to preserve date strings
distance_df = pd.read_csv(
    os.path.join(script_dir, "..", "data", "distance_earth_mars.txt"),
    sep=";",
    skiprows=1,
    names=["Date", "Distance"],
    dtype={"Date": str, "Distance": float},
    engine="python"
)

# Strip whitespace from Date column
distance_df["Date"] = distance_df["Date"].str.strip()

# Convert the 'Date' column to datetime objects, coerce errors to NaT
distance_df["Date"] = pd.to_datetime(distance_df["Date"], format="%d.%m.%Y", errors="coerce")

# Drop rows where Date could not be parsed
distance_df = distance_df.dropna(subset=["Date"])

plt.figure(figsize=(10, 5))
plt.scatter(distance_df["Date"], distance_df["Distance"], label="Distance", color="blue", marker="x",linewidths=0.8)
plt.xlabel("Date")
plt.ylabel("Distance (AU)")
plt.title("Distance between Earth and Mars Over Time")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()