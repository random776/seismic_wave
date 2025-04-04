import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
from scipy.interpolate import interp1d

# データファイルのパス
file_pat_x = 'seismic_x.dat'
file_pat_y = 'seismic_y.dat'
file_pat_t = 'seismic_t.dat'
file_speed = 'seismic_speed.dat'

# データの読み込み（エラーを避けるため `filling_values=np.nan` を指定）
data_x = np.genfromtxt(file_pat_x, delimiter=None, invalid_raise=False, filling_values=np.nan)
data_y = np.genfromtxt(file_pat_y, delimiter=None, invalid_raise=False, filling_values=np.nan)
data_t = np.genfromtxt(file_pat_t, delimiter=None, invalid_raise=False, filling_values=np.nan)

# wave_speed データの読み込み
data_speed = np.loadtxt(file_speed)

# y と wave_speed(y) の取得
y_speed = data_speed[:, 0]
speed = data_speed[:, 1]

# グリッド作成（カラーマップ用）
y_min, y_max = -10, 3
x_min, x_max = -5, 20
y_grid = np.linspace(y_min, y_max, 200)
x_grid = np.linspace(x_min, x_max, 200)
X, Y = np.meshgrid(x_grid, y_grid)

# wave_speed の補間（背景用）
wave_speed_interp = interp1d(y_speed, speed, bounds_error=False, fill_value="extrapolate")
speed_grid = wave_speed_interp(Y)  # Y に対応する wave_speed を計算

# 図の設定
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_title("Seismic wave propagation with wave speed")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.grid()

# 背景のカラーマップを表示
im = ax.imshow(speed_grid, extent=[x_min, x_max, y_min, y_max], origin="lower", cmap="jet", alpha=0.7, aspect="auto")

# カラーバーを追加
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Wave Speed")

# プロットの初期化
line, = ax.plot([], [], c="blue", marker=".", markersize=0.6, linewidth=0)

# アニメーションの更新関数
def update(frame):
    # NaN（データなし）を除外
    valid_x = data_x[frame][~np.isnan(data_x[frame])]
    valid_y = data_y[frame][~np.isnan(data_y[frame])]
    
    # 更新
    line.set_data(valid_x, valid_y)
    return line,

# アニメーションの作成
anim = animation.FuncAnimation(fig, update, frames=len(data_x), interval=100, blit=True)
rc("animation", html="jshtml")
plt.show()