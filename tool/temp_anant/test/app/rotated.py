import os
import re
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

# === FAST GC SKEW CALCULATION ===
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

# === READ MARKER POSITIONS ===
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

# === MAIN FUNCTION ===
def plot_circos(fna_path, report_path):
    record = next(SeqIO.parse(fna_path, "fasta"))
    seq = str(record.seq)
    dnaA, oriC = read_positions(report_path)
    genome_len = len(seq)

    if dnaA is None or oriC is None:
        raise ValueError("dnaA or oriC coordinates missing in report file.")

    # Decide which marker defines the origin (top)
    if oriC < dnaA:
        shift_bp = oriC  # oriC at 0.0 Mb
    else:
        shift_bp = dnaA  # dnaA at 0.0 Mb

    # === Compute GC skew ===
    skew = fast_gc_skew(seq, win=1000)
    genome_bins = len(skew)
    bin_size = genome_len / genome_bins

    # Rotate GC skew so dnaA/oriC is at 0.0 Mb
    shift_bins = int(shift_bp / bin_size)
    skew = np.roll(skew, -shift_bins)

    # Convert to polar coordinates
    angles = np.linspace(0, 2*np.pi, genome_bins, endpoint=False)

    # === PLOTTING ===
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)  # top = 0.0 Mb
    ax.set_axis_off()

    outer_r = 1.0
    inner_r = 0.8
    bar_scale = 0.25

    # === Outer circular genome ===
    ax.plot(np.linspace(0, 2*np.pi, 2000), [outer_r]*2000, color='navy', lw=4)

    # === GC Skew Curve ===
    ax.bar(angles[skew > 0], skew[skew > 0]*bar_scale, width=2*np.pi/len(skew),
           bottom=inner_r, color='lightblue', linewidth=0)
    ax.bar(angles[skew < 0], -skew[skew < 0]*bar_scale, width=2*np.pi/len(skew),
           bottom=inner_r - bar_scale*np.abs(skew[skew < 0]),
           color='orange', linewidth=0)

    # === FIXED TICKS (exact 0.00 Mb alignment) ===
    n_ticks = 28
    tick_spacing = genome_len / n_ticks
    tick_positions = [i * tick_spacing for i in range(n_ticks)]
    tick_len = 0.03
    label_offset = 0.08

    for i, pos in enumerate(tick_positions):
        # θ = π/2 for top; move clockwise by 2π/n per tick
        theta = 0 + (2 * np.pi * i / n_ticks)
        ax.plot([theta, theta], [outer_r, outer_r + tick_len],
                color='navy', lw=1.5, zorder=6)
        label = f"{(pos / 1e6):.2f} Mb"
        ax.text(theta, outer_r + tick_len + label_offset, label,
                ha='center', va='center', fontsize=9, rotation=0, zorder=6)

    # === Mark dnaA and oriC after rotation ===
    dnaA_rot = (dnaA - shift_bp) % genome_len
    oriC_rot = (oriC - shift_bp) % genome_len
    ax.plot([2*np.pi*(dnaA_rot/genome_len)], [outer_r],
            marker='s', color='blue', markersize=10, zorder=10)
    ax.plot([2*np.pi*(oriC_rot/genome_len)], [outer_r],
            marker='o', color='red', markersize=10, zorder=10)

    # === Title inside ===
    desc = record.description.split(",")[0]
    parts = desc.split(" ", 1)
    line1 = parts[0]
    line2 = parts[1] if len(parts) > 1 else ""
    ax.text(0, 0, f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha='center', va='center', fontsize=11, fontweight='bold')

    # === Legend ===
    fig.text(0.05, 0.93, u'\u25A0 dnaA', color='blue', fontsize=13, fontweight='bold')
    fig.text(0.05, 0.89, u'\u25CF oriC', color='red', fontsize=13, fontweight='bold')

    # === Save ===
    base = os.path.splitext(os.path.basename(fna_path))[0]
    plt.savefig(f"{base}_rotated_origin_aligned.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{base}_rotated_origin_aligned.svg", bbox_inches='tight')
    plt.close()
    print(f"✅ Saved {base}_rotated_origin_aligned.[png/svg] — perfectly aligned to 0.00 Mb at dnaA/oriC.")

# === MAIN ENTRY ===
if __name__ == "__main__":
    fna = next((f for f in os.listdir() if f.endswith('.fna') and not f.startswith('origami_')), None)
    report = next((f for f in os.listdir() if f.endswith('_original_report.txt')), None)
    if not fna or not report:
        raise FileNotFoundError("FASTA or report file missing.")
    plot_circos(fna, report)
