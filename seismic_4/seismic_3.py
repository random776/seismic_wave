# 走時曲線

import matplotlib.pyplot as plt
import numpy as np

# データをファイルから読み込む
file_path = 'seismic_curve.dat'
data = np.loadtxt(file_path)  # ファイルを数値データとして読み込む

# x, y を取得（z は今回は使わない）
x = data[:, 0]  # 1列目
y = data[:, 1]  # 2列目

# グラフの作成
fig, ax = plt.subplots(figsize=(7, 2.2))
ax.set_title("travel-time curve")
ax.set_xlim(0, 15.7)
ax.set_ylim(0, 6)
ax.set_xlabel("x")
ax.set_ylabel("t")

ax.plot(y, x, c="blue", marker=".", markersize=1.0, linewidth=0)
ax.grid(True)
plt.show()