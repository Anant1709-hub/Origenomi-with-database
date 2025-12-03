#!/usr/bin/env python3
import os
import shutil
import subprocess
from Bio import SeqIO

from .circos_tools import plot_circos, plot_rotated_circos
from .dot import plot_dotplot
from .linear import run_linear_synteny_plot



# HELPERS
# def extract_fna_from_gbff(gbff_path: Path, output_fna: Path):
#     version_line = None
#     definition_line = None
#     sequence_lines = []
#     in_origin = False

#     with open(gbff_path, "r") as f:
#         for line in f:
#             if line.startswith("VERSION"):
#                 version_line = ">" + line.split("VERSION", 1)[1].strip()
#             elif line.startswith("DEFINITION"):
#                 definition_line = line.split("DEFINITION", 1)[1].strip()
#             elif line.startswith("ORIGIN"):
#                 in_origin = True
#             elif in_origin:
#                 if line.startswith("//"):
#                     break
#                 seq = "".join(filter(str.isalpha, line)).upper()
#                 sequence_lines.append(seq)

#     if not version_line or not definition_line or not sequence_lines:
#         raise ValueError(f"Invalid or incomplete GBFF file: {gbff_path}")

#     with open(output_fna, "w") as out:
#         out.write(version_line + " " + definition_line + "\n")
#         out.write("".join(sequence_lines) + "\n")

#     print(f"[INFO] Extracted .fna saved to {output_fna}")

def parse_report(report_path):
    """
    Parses the Origami original report and returns per-sequence dnaA/oriC status.
    Works for multiline blocks:

    >>seqid
    Chosen dnaA:
      ID: ...
      Coords: A - B
    """
    results = {}
    current_seq = None
    mode = None  # "dnaA" or "oriC" or None

    with open(report_path, "r") as f:
        for line in f:
            line = line.strip()

            # New sequence
            if line.startswith(">>"):
                current_seq = line[2:].split()[0]
                results[current_seq] = {
                    "dnaA_found": False,
                    "oriC_found": False,
                    "dnaA_coords": None,
                    "oriC_coords": None,
                }
                mode = None
                continue

            # Detect section headers
            if line.startswith("Chosen dnaA"):
                mode = "dnaA"
                continue

            if line.startswith("Chosen oriC"):
                mode = "oriC"
                continue

            # Detect coords inside the mode
            if line.startswith("Coords:") and mode:
                parts = line.replace("Coords:", "").strip().split("-")
                start, end = int(parts[0]), int(parts[1])

                if mode == "dnaA":
                    results[current_seq]["dnaA_found"] = True
                    results[current_seq]["dnaA_coords"] = (start, end)

                elif mode == "oriC":
                    results[current_seq]["oriC_found"] = True
                    results[current_seq]["oriC_coords"] = (start, end)

                continue

    return results



# PLOT GENERATION (SEQUENTIAL)
def generate_all_plots(fasta_output_path, report_path, rotated_fasta_path, work_dir, gcf_fasta_path):

    plot_dir = os.path.join(work_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # 1. Parse report to get per-sequence results
    report_info = parse_report(report_path)

    print("[LOG] Parsed report results:")
    for seq_id, data in report_info.items():
        print("  ", seq_id, data)

    # 2. Loop over sequences
    for seq_id, data in report_info.items():

        if not (data["dnaA_found"] and data["oriC_found"]):
            print(f"[LOG] {seq_id}: dnaA or oriC missing — skipping plots")
            continue

        print(f"[LOG] {seq_id}: Generating plots")

        # ===== ORIGINAL CIRCOS =====
        orig_plot = plot_circos(gcf_fasta_path, report_path)
        shutil.copy(orig_plot, os.path.join(plot_dir, f"{seq_id}_original_circos.png"))

        # ===== ROTATED CIRCOS =====
        rot_result = plot_rotated_circos(rotated_fasta_path, report_path)
        rot_png = rot_result[0] if isinstance(rot_result, (tuple, list)) else rot_result
        shutil.copy(rot_png, os.path.join(plot_dir, f"{seq_id}_rotated_circos.png"))

        # ===== DOTPLOT =====
        dot_png = plot_dotplot(
            os.path.basename(fasta_output_path),
            os.path.basename(rotated_fasta_path),
            os.path.basename(report_path),
            output_dir=work_dir,
        )
        if dot_png:
            shutil.copy(dot_png, os.path.join(plot_dir, f"{seq_id}_dotplot.png"))

        # ===== LINEAR SYNTENY =====
        synteny_dir = os.path.join(work_dir, "linear_synteny")
        os.makedirs(synteny_dir, exist_ok=True)

        synteny_png = run_linear_synteny_plot(synteny_dir)
        if synteny_png:
            shutil.copy(synteny_png, os.path.join(plot_dir, f"{seq_id}_linear_synteny.png"))

    print("[LOG] All multi-sequence plot generation complete.")
