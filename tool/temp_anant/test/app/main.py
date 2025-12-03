#!/usr/bin/env python3
import os
import shutil
import base64
import subprocess
from pathlib import Path
from fastapi import FastAPI, Form, File, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from Bio import SeqIO

from .circos_tools import plot_circos, plot_rotated_circos
from .dot import plot_dotplot
from .linear import run_linear_synteny_plot


app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)



def extract_fna_from_gbff(gbff_path: Path, output_fna: Path):
    version_line = None
    definition_line = None
    sequence_lines = []
    in_origin = False

    with open(gbff_path, "r") as f:
        for line in f:
            if line.startswith("VERSION"):
                version_line = ">" + line.split("VERSION", 1)[1].strip()
            elif line.startswith("DEFINITION"):
                definition_line = line.split("DEFINITION", 1)[1].strip()
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif in_origin:
                if line.startswith("//"):
                    break
                seq = "".join(filter(str.isalpha, line)).upper()
                sequence_lines.append(seq)

    if not version_line or not definition_line or not sequence_lines:
        raise ValueError(f"Invalid or incomplete GBFF file: {gbff_path}")

    with open(output_fna, "w") as out:
        out.write(version_line + " " + definition_line + "\n")
        out.write("".join(sequence_lines) + "\n")

    print(f"[INFO] Extracted .fna saved to {output_fna}")


def compute_genome_stats(fasta_file):
    def is_plasmid(header):
        h = header.lower()
        return "plasmid" in h or "plsm" in h or "pld" in h or "pls" in h

    seqs, current_seq, headers = [], [], []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    seqs.append("".join(current_seq))
                    current_seq = []
                headers.append(line.strip())
            else:
                current_seq.append(line.strip())
        if current_seq:
            seqs.append("".join(current_seq))

    total_size = sum(len(s) for s in seqs)
    gc = sum(s.count("G") + s.count("g") + s.count("C") + s.count("c") for s in seqs)
    genome_gc = round((gc / total_size) * 100, 2) if total_size > 0 else 0.0

    contigs = len(seqs)
    plasmids = sum(1 for h in headers if is_plasmid(h))
    chromosomes = contigs - plasmids

    return {
        "total_bp": total_size,
        "contigs": contigs,
        "plasmids": plasmids,
        "chromosomes": chromosomes,
        "gc_percent": genome_gc,
    }


def encode_file_to_base64(file_path):
    if not file_path:
        return None
    # file_path may already be a base64 string in some contexts; ensure path exists
    try:
        if isinstance(file_path, str) and file_path.startswith("data:image/"):
            return file_path
        if not os.path.exists(file_path):
            return None
        with open(file_path, "rb") as f:
            data = base64.b64encode(f.read()).decode("utf-8")
        return f"data:image/png;base64,{data}"
    except Exception as e:
        print(f"[WARN] Failed to encode {file_path}: {e}")
        return None


def run_origami_sync(fasta_input, work_dir):
    print(f"[INFO] Running Origami synchronously on {fasta_input}")
    result = subprocess.run(
        ["origami", "run", "-i", fasta_input.name, "--keep-temp","--keep-plot"],
        cwd=work_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        print(f"[ERROR] Origami failed: {result.stderr}")
        raise RuntimeError(result.stderr)

    print("[INFO] Origami finished successfully.")
    return result.stdout


def generate_all_plots(
    r,
    fasta_output_path,
    report_path,
    rotated_fasta_path,
    work_dir,
    gcf_fasta_path
):
    try:
        fasta_output_path = os.path.abspath(fasta_output_path)
        report_path = os.path.abspath(report_path)
        rotated_fasta_path = os.path.abspath(rotated_fasta_path) if rotated_fasta_path else None
        gcf_fasta_path = os.path.abspath(gcf_fasta_path) if gcf_fasta_path else None
        work_dir = os.path.abspath(work_dir)

        # ----- LOGGING -----
        print("[LOG] fasta_output_path:", fasta_output_path)
        print("[LOG] gcf_fasta_path:", gcf_fasta_path)
        print("[LOG] rotated_fasta_path:", rotated_fasta_path)
        print("[LOG] report_path:", report_path)
        print("[LOG] work_dir:", work_dir)

        dnaA_ok = r.get("dnaA_coords") != "Not found"
        oriC_ok = r.get("oriC_coords") != "Not found"

        # STOP ALL PLOTS IF dnaA OR oriC NOT FOUND
        if not (dnaA_ok and oriC_ok):
            print("[LOG] dnaA or oriC missing — skipping ALL plots")
            r["no_rotation_needed"] = True
            r["original_plots"] = []
            r["rotated_plots"] = []
            r["dot_plots"] = []
            r["pgv_plot"] = None
            # ensure rotation metadata fields exist but are None
            r["rotated_from"] = None
            r["rotated_shift_bp"] = None
            return r

        # --------------------------------------------------
        # ORIGINAL CIRCOS (GCF FASTA)
        # --------------------------------------------------
        print("[LOG] Running original circos with:", gcf_fasta_path)
        orig_plot_path = plot_circos(gcf_fasta_path, report_path)

        # --------------------------------------------------
        # ROTATED CIRCOS (Origami FASTA)
        # --------------------------------------------------
        print("[LOG] Running rotated circos with:", rotated_fasta_path)
        # plot_rotated_circos now returns (png_path, rotated_from, shift_bp)
        rot_result = plot_rotated_circos(rotated_fasta_path, report_path)

        # Support both old single-return and new triple-return signatures
        if isinstance(rot_result, tuple) or isinstance(rot_result, list):
            if len(rot_result) == 3:
                rot_plot_path, rotated_from, shift_bp = rot_result
            elif len(rot_result) == 1:
                rot_plot_path = rot_result[0]
                rotated_from = None
                shift_bp = None
            else:
                # unexpected shape
                rot_plot_path = rot_result[0]
                rotated_from = None
                shift_bp = None
        else:
            rot_plot_path = rot_result
            rotated_from = None
            shift_bp = None

        r["original_plots"] = [encode_file_to_base64(orig_plot_path)]
        r["rotated_plots"] = [encode_file_to_base64(rot_plot_path)]
        r["rotated_from"] = rotated_from
        r["rotated_shift_bp"] = shift_bp

        # --------------------------------------------------
        # NO ROTATED FASTA → STOP HERE
        # --------------------------------------------------
        if not rotated_fasta_path or not os.path.exists(rotated_fasta_path):
            print("[LOG] Rotated FASTA missing — skipping dot/synteny")
            r["dot_plots"] = []
            r["pgv_plot"] = None
            return r

        # --------------------------------------------------
        # DOTPLOT
        # --------------------------------------------------
        print("[LOG] Running dotplot")
        dotplot_path = plot_dotplot(
            os.path.basename(fasta_output_path),
            os.path.basename(rotated_fasta_path),
            os.path.basename(report_path),
            output_dir=work_dir,
        )

        if dotplot_path:
            print("[LOG] dotplot output:", dotplot_path)
            r["dot_plots"] = [encode_file_to_base64(dotplot_path)]

        # --------------------------------------------------
        # LINEAR SYNTENY
        # --------------------------------------------------
        print("[LOG] Running linear synteny plot")
        synteny_dir = os.path.join(work_dir, "linear_synteny")
        os.makedirs(synteny_dir, exist_ok=True)

        orig_copy = os.path.join(synteny_dir, os.path.basename(fasta_output_path))
        rot_copy = os.path.join(synteny_dir, os.path.basename(rotated_fasta_path))

        for src, dest in [(fasta_output_path, orig_copy), (rotated_fasta_path, rot_copy)]:
            if not os.path.exists(dest):
                try:
                    os.symlink(src, dest)
                except Exception:
                    shutil.copy(src, dest)

        pgv_png = run_linear_synteny_plot(synteny_dir,r)
        print("[LOG] PGV synteny output:", pgv_png)

        if pgv_png and os.path.exists(pgv_png):
            r["pgv_plot"] = encode_file_to_base64(pgv_png)
        else:
            r["pgv_plot"] = None

        # Cleanup synteny files
        for p in Path(synteny_dir).glob("pgv_out"):
            shutil.rmtree(p, ignore_errors=True)

        for p in Path(synteny_dir).glob("*.fna"):
            try:
                p.unlink()
            except:
                pass

    except Exception as e:
        print("[ERROR] generate_all_plots failed:", str(e))
        r["error"] = str(e)

    return r


# ============================================================
# MAIN ENDPOINT (SYNCHRONOUS)
# ============================================================
@app.post("/run-origami")
def run_origami(
    num: str = Form(...),
    refseq_id: str = Form(None),
    file: UploadFile = File(None),
):
    try:
        work_dir = Path(f"test_lazy_{num}")
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        print(f"Work Directory :{work_dir}")
        # -----------------------------
        # INPUT HANDLING
        # -----------------------------
        if file is not None:
            fasta_input = work_dir / file.filename
            with open(fasta_input, "wb") as f_out:
                f_out.write(file.file.read())
            refseq_id = Path(file.filename).stem

            if file.filename.endswith(".gbff"):
                output_fna = work_dir / f"{refseq_id}.fna"
                extract_fna_from_gbff(fasta_input, output_fna)
                fasta_input = output_fna

        elif refseq_id:
            common = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
            id1, id2, id3 = refseq_id[4:7], refseq_id[7:10], refseq_id[10:13]
            link = f"{common}{id1}/{id2}/{id3}/{refseq_id}/{refseq_id}_genomic.fna.gz"

            fasta_gz = work_dir / f"{refseq_id}_genomic.fna.gz"
            fasta_input = work_dir / f"{refseq_id}_genomic.fna"

            subprocess.run(["wget", "-q", "-O", str(fasta_gz), link])
            if not fasta_gz.exists():
                subprocess.run(["curl", "-L", "-o", str(fasta_gz), link], check=True)
            subprocess.run(["gunzip", str(fasta_gz)], check=True)

        else:
            return JSONResponse({"status": "error", "error": "No file or RefSeq ID provided."})

        # -----------------------------
        # RUN ORIGAMI (SYNC)
        # -----------------------------
        run_origami_sync(fasta_input, work_dir)

        # -----------------------------
        # FIND OUTPUTS
        # -----------------------------
        files = os.listdir(work_dir)
        report_file = next((f for f in files if f.endswith("_original_report.txt")), None)

        # GCF*.fna
        gcf_fasta = next((f for f in files if f.startswith("GCF") and f.endswith(".fna")), None)

        # *_original.fna
        fasta_output = next((f for f in files if f.endswith("_original.fna")), None)

        rotated_fasta = next(
            (f for f in files if f.startswith("origami_") and f.endswith(".fna") and not f.endswith("_original.fna")),
            None
        )

        if not report_file or not fasta_output or not gcf_fasta:
            print("[ERROR] Required files missing:", report_file, fasta_output, gcf_fasta)
            return JSONResponse({"status": "error", "error": "Required Origami or GCF files missing"})

        report_path = work_dir / report_file
        fasta_output_path = work_dir / fasta_output
        gcf_fasta_path = work_dir / gcf_fasta
        rotated_fasta_path = work_dir / rotated_fasta if rotated_fasta else None

        # ----- LOGGING -----
        print("[LOG] report_path:", report_path)
        print("[LOG] fasta_output_path (*_original.fna):", fasta_output_path)
        print("[LOG] gcf_fasta_path (GCF file):", gcf_fasta_path)
        print("[LOG] rotated_fasta_path:", rotated_fasta_path)

        # -----------------------------
        # STATS
        # -----------------------------
        stats = compute_genome_stats(fasta_input)
        stats["organism"] = refseq_id
        stats["tax_id"] = refseq_id.split("_")[1] if "_" in refseq_id else "uploaded"

        record = next(SeqIO.parse(fasta_input, "fasta"))
        seq = str(record.seq)
        orig_size = len(seq)

        record = next(SeqIO.parse(fasta_output_path, "fasta"))
        seq = str(record.seq)
        trimmed_size = len(seq)

        trimmed_portion = abs(orig_size - trimmed_size)

        stats["original_genome_size_bp"] = orig_size
        stats["trimmed_genome_size_bp"] = trimmed_size
        stats["trimmed_portion_bp"] = trimmed_portion

        # -----------------------------
        # PARSE REPLICONS
        # -----------------------------
        replicons = []
        current = None

        with open(report_path) as f:
            for line in f:
                line = line.strip()

                if line.startswith(">>"):
                    if current:
                        replicons.append(current)
                    parts = line[2:].strip().split(" ", 1)
                    current = {
                        "nc_id": parts[0],
                        "organism_name": parts[1].replace(", complete sequence", "") if len(parts) > 1 else "Unknown",
                        "dnaA_coords": "Not found",
                        "oriC_coords": "Not found",
                        "original_plots": [],
                        "rotated_plots": [],
                        "dot_plots": [],
                        "pgv_plot": None,
                        "no_rotation_needed": False,
                        # placeholders for rotation metadata
                        "rotated_from": None,
                        "rotated_shift_bp": None,
                    }

                elif "Chosen dnaA:" in line:
                    for _ in range(5):
                        next_line = next(f, "").strip()
                        if "Coords:" in next_line:
                            current["dnaA_coords"] = next_line.split("Coords:", 1)[1].strip()
                            break

                elif "Chosen oriC:" in line:
                    for _ in range(5):
                        next_line = next(f, "").strip()
                        if "Coords:" in next_line:
                            current["oriC_coords"] = next_line.split("Coords:", 1)[1].strip()
                            break

        if current:
            replicons.append(current)


        # PLOT GENERATION
        updated_replicons = []
        print(f"Work Directory:{work_dir}")
        for r in replicons:
            result = generate_all_plots(
                r,
                fasta_output_path,
                report_path,
                rotated_fasta_path,
                str(work_dir),
                gcf_fasta_path,
            )
            updated_replicons.append(result)

        # -----------------------------
        # FINAL RESPONSE
        # -----------------------------
        return JSONResponse({
            "status": "success",
            "refseq_id": refseq_id,
            "stats": stats,
            "replicons": updated_replicons,
            "download_report": f"/download/{work_dir.name}/{report_file}",
            "download_rotated_fasta": (
                f"/download/{work_dir.name}/{rotated_fasta}" if rotated_fasta else None
            ),
        })

    except Exception as e:
        print("[ERROR] run_origami failed:", str(e))
        return JSONResponse({"status": "failed", "error": str(e)})


# DOWNLOAD ENDPOINT
@app.get("/download/{folder}/{filename}")
def download_file(folder: str, filename: str):
    file_path = Path(folder) / filename
    if file_path.exists():
        return FileResponse(file_path, filename=filename)
    return JSONResponse({"error": "File not found"}, status_code=404)
