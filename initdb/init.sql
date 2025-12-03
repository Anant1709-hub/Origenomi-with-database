CREATE TABLE IF NOT EXISTS refseq (
    refseq_id  TEXT PRIMARY KEY,
    status     TEXT NOT NULL,
    stats      JSONB NOT NULL,
    created_at TIMESTAMP DEFAULT NOW()
);

CREATE TABLE IF NOT EXISTS replicon (
    id SERIAL PRIMARY KEY,
    refseq_id TEXT REFERENCES refseq(refseq_id) ON DELETE CASCADE,
    nc_id TEXT NOT NULL,
    organism_name TEXT NOT NULL,
    dnaA_coords TEXT,
    oriC_coords TEXT,
    no_rotation_needed BOOLEAN,
    rotated_from TEXT,
    rotated_shift_bp INT,
    original_plots TEXT,
    rotated_plots TEXT,
    dot_plots TEXT,
    pgv_plot TEXT,
    created_at TIMESTAMP DEFAULT NOW()
);
