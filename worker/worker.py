import requests_unixsocket
import psycopg2
import os
import json
import sys

SOCK_PATH = "/tmp/origenomi/app.sock"
ENDPOINT = "run-origami"

DB_HOST = "postgres"
DB_NAME = "origenomi_database"
DB_USER = "anant"
DB_PASS = "8721"



# read input file from command line 
if len(sys.argv) < 2:
    print("Usage: python worker.py <input.fna>")
    sys.exit(1)

FNA_FILE_PATH = sys.argv[1]

if not os.path.exists(FNA_FILE_PATH):
    print(f"Error: File not found: {FNA_FILE_PATH}")
    sys.exit(1)

print(f"Processing file: {FNA_FILE_PATH}")


#  connect to UNIX socket API 
session = requests_unixsocket.Session()

url = f"http+unix://{SOCK_PATH.replace('/', '%2F')}/{ENDPOINT}"

print("Sending .fna file to:", url)

data = {"num": 1}  

with open(FNA_FILE_PATH, "rb") as f:
    files = {
        "file": (os.path.basename(FNA_FILE_PATH), f, "application/octet-stream")
    }
    response = session.post(url, data=data, files=files)

data = response.json()



status = data.get("status")
refseq_id = data.get("refseq_id")
stats = data.get("stats")
replicons = data.get("replicons", [])



conn = psycopg2.connect(
    host=DB_HOST,
    dbname=DB_NAME,
    user=DB_USER,
    password=DB_PASS
)
cur = conn.cursor()


cur.execute(
    """
    INSERT INTO refseq (refseq_id, status, stats)
    VALUES (%s, %s, %s)
    ON CONFLICT (refseq_id) DO NOTHING;
    """,
    (refseq_id, status, json.dumps(stats))
)


for r in replicons:
    cur.execute(
        """
        INSERT INTO replicon (
            refseq_id,
            nc_id,
            organism_name,
            dnaA_coords,
            oriC_coords,
            no_rotation_needed,
            rotated_from,
            rotated_shift_bp,
            original_plots,
            rotated_plots,
            dot_plots,
            pgv_plot
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
        """,
        (
            refseq_id,
            r.get("nc_id"),
            r.get("organism_name"),
            r.get("dnaA_coords"),
            r.get("oriC_coords"),
            r.get("no_rotation_needed"),
            r.get("rotated_from"),
            r.get("rotated_shift_bp"),
            json.dumps(r.get("original_plots")),
            json.dumps(r.get("rotated_plots")),
            json.dumps(r.get("dot_plots")),
            json.dumps(r.get("pgv_plot")),
        )
    )


conn.commit()
conn.close()

print(f"Stored results for {refseq_id} successfully.")
