import matplotlib.pyplot as plt
import numpy as np

# データをファイルから読み込む
file_path_2 = 'seismic_2.dat'
data_2 = []
file_path_3 = 'seismic_3.dat'
data_3 = []

with open(file_path_2, 'r') as file:
    for line in file:
        data_2.append(float(line.strip()))

with open(file_path_3, 'r') as file:
    for line in file:
        data_3.append(float(line.strip()))

fig, ax = plt.subplots(figsize=(7, 2.2))
ax.set_title("Seismic wave")
ax.set_xlabel("x")
ax.set_ylabel("y")

ax.plot(data_2, data_3, c="blue", marker=".", markersize= 1.8, linewidth=0)
ax.grid(True)
plt.show()