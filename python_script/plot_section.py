import numpy as np
import matplotlib.pyplot as plt
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
data = np.loadtxt(os.path.join(script_dir, "section_dos.dat"), skiprows=1)

z_avg = data[:, 0]
energy = data[:, 1]
total = data[:, 2]

n_layers = 3
n_energy = len(data) // n_layers

z_avg_vals = z_avg[::n_energy]
energy_vals = energy[:n_energy]
total_2d = total.reshape(n_layers, n_energy).T

fig, ax = plt.subplots(figsize=(7, 5))
x_edges = np.arange(n_layers + 1)
y_edges = np.linspace(energy_vals.min(), energy_vals.max(), n_energy + 1)

mesh = ax.pcolormesh(x_edges, y_edges, total_2d, cmap="Reds", shading="flat")
cbar = fig.colorbar(mesh, ax=ax, label="Total DOS")

ax.set_xticks(0.5 + np.arange(n_layers))
ax.set_xticklabels([f"{z:.4f}" for z in z_avg_vals])
ax.set_xlabel("z_avg (Å)")
ax.set_ylabel("Energy (eV)")

plt.tight_layout()
plt.savefig(os.path.join(script_dir, "section_dos.png"), dpi=200)
print("  -> section_dos.png")
