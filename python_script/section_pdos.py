def add_tag_column(input_path, output_path):
    with open(input_path, "r") as f:
        lines = f.readlines()

    header_marker = "Position, move_x, move_y, move_z"
    section_marker = "Force"
    in_atom_section = False

    out_lines = []
    for line in lines:
        stripped = line.strip()

        if header_marker in stripped:
            in_atom_section = True
        elif stripped == section_marker or stripped.startswith(section_marker):
            in_atom_section = False

        if in_atom_section and header_marker not in stripped:
            if len(line.rstrip("\n").split()) >= 8:
                out_lines.append(line.rstrip("\n").rsplit(None, 1)[0] + "  1\n")
            else:
                out_lines.append(line.rstrip("\n") + "  1\n")
        else:
            out_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(out_lines)


def split_by_layers(input_path, layers, output_prefix="atom_layer", dos_input_path=None):
    with open(input_path, "r") as f:
        lines = f.readlines()

    header_marker = "Position, move_x, move_y, move_z"
    section_marker = "Force"
    in_atom_section = False

    atom_indices = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if header_marker in stripped:
            in_atom_section = True
            continue
        elif stripped == section_marker or stripped.startswith(section_marker):
            in_atom_section = False
            continue
        if in_atom_section:
            atom_indices.append(i)

    z_vals = [float(lines[idx].split()[3]) for idx in atom_indices]

    sorted_idx = sorted(range(len(z_vals)), key=lambda i: z_vals[i])
    n_atoms = len(atom_indices)
    base = n_atoms // layers
    remainder = n_atoms % layers

    layer_boundaries = []
    start = 0
    for i in range(layers):
        size = base + (1 if i < remainder else 0)
        end = start + size
        if end > n_atoms:
            end = n_atoms
        layer_boundaries.append((start, end))
        start = end

    layer_atom_sets = []
    for start, end in layer_boundaries:
        indices = set(sorted_idx[start:end])
        layer_atom_sets.append(indices)

    import os

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "tmp")

    for layer_idx in range(layers):
        out_lines = []
        in_atom_section = False
        atom_counter = 0
        for line in lines:
            stripped = line.strip()
            if header_marker in stripped:
                in_atom_section = True
                out_lines.append(line)
                continue
            elif stripped == section_marker or stripped.startswith(section_marker):
                in_atom_section = False
                out_lines.append(line)
                continue

            if in_atom_section:
                if atom_counter in layer_atom_sets[layer_idx]:
                    parts = line.rstrip("\n").rsplit(None, 1)
                    out_lines.append(parts[0] + "  1\n")
                else:
                    parts = line.rstrip("\n").rsplit(None, 1)
                    out_lines.append(parts[0] + "  0\n")
                atom_counter += 1
            else:
                out_lines.append(line)

        sub_dir = os.path.join(out_dir, f"PDOS{layer_idx + 1}")
        os.makedirs(sub_dir, exist_ok=True)
        output_path = os.path.join(sub_dir, "atom.config")
        with open(output_path, "w") as f:
            f.writelines(out_lines)
        print(f"  -> {output_path}")

        if dos_input_path is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            dos_input_path = os.path.join(script_dir, "DOS.input")
        if os.path.exists(dos_input_path):
            import shutil
            dst = os.path.join(sub_dir, "DOS.input")
            shutil.copy2(dos_input_path, dst)
            print(f"  -> {dst}")

        import glob

        for src in (
            [os.path.join(script_dir, "REPORT")] +
            glob.glob(os.path.join(script_dir, "*.UPF")) +
            [os.path.join(script_dir, "OUT.EIGEN")] +
            glob.glob(os.path.join(script_dir, "bpsiio*")) +
            [os.path.join(script_dir, f) for f in ("OUT.SYMM", "OUT.IND_EXT_KPT")]
        ):
            if os.path.exists(src):
                dst = os.path.join(sub_dir, os.path.basename(src))
                if not os.path.exists(dst):
                    import shutil
                    shutil.copy2(src, dst)
                    print(f"  -> {dst}")


def write_section_dos(atom_config_path, layers, tmp_dir, output_path):
    import numpy as np
    import os

    with open(atom_config_path, "r") as f:
        lines = f.readlines()

    lz = float(lines[4].split()[2])

    header_marker = "Position, move_x, move_y, move_z"
    section_marker = "Force"
    in_atom_section = False
    atom_indices = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if header_marker in stripped:
            in_atom_section = True
            continue
        elif stripped == section_marker or stripped.startswith(section_marker):
            in_atom_section = False
            continue
        if in_atom_section:
            atom_indices.append(i)

    z_vals = np.array([float(lines[idx].split()[3]) for idx in atom_indices]) * lz
    sorted_idx = np.argsort(z_vals)
    sorted_z = z_vals[sorted_idx]

    n_atoms = len(z_vals)
    base = n_atoms // layers
    remainder = n_atoms % layers

    z_avgs = []
    start = 0
    for i in range(layers):
        size = base + (1 if i < remainder else 0)
        end = start + size
        z_avgs.append(np.mean(sorted_z[start:end]))
        start = end

    all_data = []
    for i in range(layers):
        dos_path = os.path.join(tmp_dir, f"PDOS{i + 1}", "DOS.totalspin")
        data = np.loadtxt(dos_path, skiprows=1)
        energy = data[:, 0]
        total = data[:, 1]
        z_avg = z_avgs[i]
        for e, t in zip(energy, total):
            all_data.append((z_avg, e, t))

    with open(output_path, "w") as f:
        f.write(f"# {'z_avg':>15} {'energy':>15} {'total':>15}\n")
        for z_avg, e, t in all_data:
            f.write(f"  {z_avg:18.8E}  {e:18.8E}  {t:18.8E}\n")
    print(f"  -> {output_path}")


if __name__ == "__main__":
    import os
    import sys
    import argparse
    import numpy as np
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description="Process DOS layers and plot")
    parser.add_argument("--layers", type=int, default=3, help="number of layers to split (default: 3)")
    parser.add_argument("--cmax", type=float, default=None, help="colormap maximum value")
    parser.add_argument("--enerange", type=float, nargs=2, default=None, metavar=("EMIN", "EMAX"), help="energy range for plot, e.g. -20 0")
    parser.add_argument("--dpi", type=int, default=300, help="image DPI (default: 300)")
    parser.add_argument("--fermi", action="store_true", default=False, help="align energy to Fermi level from OUT.FERMI and draw dashed line")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)

    atom_config = os.path.join(script_dir, "atom.config")
    atom_tag_config = os.path.join(script_dir, "atom_tag.config")
    dos_input = os.path.join(script_dir, "DOS.input")

    add_tag_column(atom_config, atom_tag_config)

    split_by_layers(atom_tag_config, layers=args.layers, dos_input_path=dos_input)

    tmp_dir = os.path.join(script_dir, "tmp")
    for i in range(1, args.layers + 1):
        pdos_dir = os.path.join(tmp_dir, f"PDOS{i}")
        import subprocess
        subprocess.run(["plot_DOS_interp.x"], cwd=pdos_dir)

    dat_path = os.path.join(script_dir, "section_dos.dat")
    write_section_dos(atom_tag_config, layers=args.layers, tmp_dir=tmp_dir, output_path=dat_path)

    data = np.loadtxt(dat_path, skiprows=1)
    z_avg = data[:, 0]
    energy = data[:, 1]
    total = data[:, 2]

    if args.fermi:
        fermi_path = os.path.join(script_dir, "OUT.FERMI")
        with open(fermi_path) as f:
            efermi = float(f.read().split("=")[1].split()[0])
        energy -= efermi

    n_energy = len(data) // args.layers
    z_avg_vals = z_avg[::n_energy]
    energy_vals = energy[:n_energy]
    total_2d = total.reshape(args.layers, n_energy).T

    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["axes.titleweight"] = "bold"

    fig, ax = plt.subplots(figsize=(7, 5))
    x_edges = np.arange(args.layers + 1)
    y_edges = np.linspace(energy_vals.min(), energy_vals.max(), n_energy + 1)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("white_red", [(1, 1, 1), (1, 0, 0)])
    mesh = ax.pcolormesh(x_edges, y_edges, total_2d, cmap=cmap, shading="flat", vmax=args.cmax)
    fig.colorbar(mesh, ax=ax, label="Total DOS")
    ax.set_xticks(0.5 + np.arange(args.layers))
    ax.set_xticklabels([f"{z:.2f}" for z in z_avg_vals], rotation=30, ha="right")
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    ax.set_xlabel("z_avg (Å)")
    ax.set_ylabel("Energy (eV)")
    if args.enerange:
        ax.set_ylim(*args.enerange)
    if args.fermi:
        ax.axhline(y=0, color="black", linestyle="--", linewidth=1)
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, "section_dos.png"), dpi=args.dpi)
    print("  -> section_dos.png")
