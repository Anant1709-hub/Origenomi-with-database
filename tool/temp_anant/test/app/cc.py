import os
import re
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

# ---------- GC SKEW FUNCTION ----------
def fast_gc_skew(seq, win=1000):
    """Compute GC skew quickly using NumPy."""
    seq = np.frombuffer(seq.upper().encode(), dtype='S1')
    g = (seq == b'G').astype(int)
    c = (seq == b'C').astype(int)
    cum_g = np.cumsum(g)
    cum_c = np.cumsum(c)
    skew = np.empty(len(seq)//win)
    for i in range(len(skew)):
        start = i * win
        end = start + win
        g_win = cum_g[end - 1] - (cum_g[start - 1] if start > 0 else 0)
        c_win = cum_c[end - 1] - (cum_c[start - 1] if start > 0 else 0)
        denom = g_win + c_win
        skew[i] = (g_win - c_win) / denom if denom > 0 else 0
    return skew


# ---------- READ POSITIONS ----------
def read_positions(report):
    """Extract dnaA and oriC positions from the origami report."""
    dnaA, oriC = None, None
    with open(report) as f:
        for line in f:
            if "Chosen dnaA:" in line:
                next(f)
                coords = next(f)
                m = re.search(r"Coords:\s*(\d+)\s*-\s*(\d+)", coords)
                if m:
                    dnaA = int(m.group(1))
            if "Chosen oriC:" in line:
                next(f)
                coords = next(f)
                m = re.findall(r"(\d+)\s*-\s*(\d+)", coords)
                if m:
                    oriC = int(m[0][0])
    return dnaA, oriC


# ---------- MAIN CIRCOS FUNCTION ----------
def plot_circos(fna_path, report_path):
    record = next(SeqIO.parse(fna_path, "fasta"))
    seq = str(record.seq)
    dnaA, oriC = read_positions(report_path)
    genome_len = len(seq)

    # GC skew computation
    skew = fast_gc_skew(seq, win=1000)
    angles = np.linspace(0, 2*np.pi, len(skew), endpoint=False)

    # --- Setup plot ---
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)  # ensures 0.0 Mb is at top
    ax.set_axis_off()

    # --- Outer circle ---
    outer_r = 1.0
    ax.plot(np.linspace(0, 2*np.pi, 2000), [outer_r]*2000, color='navy', lw=4)

    # --- GC skew bars ---
    bar_r = 0.8  # inner GC radius
    bar_scale = 0.25
    ax.bar(angles[skew > 0], skew[skew > 0]*bar_scale, width=2*np.pi/len(skew),
           bottom=bar_r, color='lightblue', linewidth=0)
    ax.bar(angles[skew < 0], -skew[skew < 0]*bar_scale, width=2*np.pi/len(skew),
           bottom=bar_r - bar_scale*np.abs(skew[skew < 0]), color='orange', linewidth=0)

    
    # --- Tick logic (Always 28 ticks) ---
    num_ticks = 28
    tick_spacing = genome_len / num_ticks
    label_offset = 0.10
    for i in range(num_ticks):
        pos = tick_spacing * i
        theta = 2 * np.pi * (pos / genome_len)
        ax.plot([theta, theta], [outer_r, outer_r + 0.03], color='navy', lw=1.5)
        label = f"{pos / 1e6:.2f} Mb"
        ax.text(theta, outer_r + label_offset, label,
                ha='center', va='center', fontsize=9, rotation=0)
    # --- dnaA and oriC markers ---
    if dnaA:
        theta_dnaA = 2 * np.pi * (dnaA / genome_len)
        ax.plot([theta_dnaA], [outer_r], marker='s', color='blue', markersize=10, zorder=10)
    if oriC:
        theta_oriC = 2 * np.pi * (oriC / genome_len)
        ax.plot([theta_oriC], [outer_r], marker='o', color='red', markersize=10, zorder=10)

    # --- Legend (top-left corner) ---
    fig.text(0.05, 0.93, u'\u25A0 dnaA', color='blue', fontsize=13, fontweight='bold')
    fig.text(0.05, 0.89, u'\u25CF oriC', color='red', fontsize=13, fontweight='bold')

    # --- Title inside circle ---
    description = record.description.split(",")[0]
    parts = description.split(" ", 1)
    line1 = parts[0]
    line2 = parts[1] if len(parts) > 1 else ""
    ax.text(0, 0, f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha='center', va='center', fontsize=10, fontweight='bold')

    # --- Save output ---
    base = os.path.splitext(os.path.basename(fna_path))[0]
    plt.savefig(f"{base}_final_circos_v9.png", dpi=400, bbox_inches='tight')
    plt.savefig(f"{base}_final_circos_v9.svg", bbox_inches='tight')
    plt.close()
    print(f"✅ Saved {base}_final_circos_v9.[png/svg]")


# ---------- MAIN EXECUTION ----------
if __name__ == "__main__":
    fna = next((f for f in os.listdir() if f.endswith(".fna") and not f.startswith("origami_")), None)
    report = next((f for f in os.listdir() if f.endswith("_original_report.txt")), None)
    if not fna or not report:
        raise FileNotFoundError("FASTA or report file missing.")
    plot_circos(fna, report)

