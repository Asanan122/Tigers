#!/usr/bin/env python3
import os
import shutil
from pathlib import Path

def setup_snpeff_database():
    root = Path(__file__).resolve().parent.parent
    snpeff_data = root / "tools" / "snpEff" / "data"
    db_dir = snpeff_data / "panthera_tigris"
    db_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(root / "data" / "reference" / "panthera_tigris.fasta", db_dir / "sequences.fa")
    shutil.copy2(root / "data" / "reference" / "panthera_tigris_scaffolds.gff", db_dir / "genes.gff")
    shutil.copy2(root / "data" / "reference" / "chromosomes.map", db_dir / "chromosomes.map")

    data_dir_abs = str(snpeff_data.resolve()) + "/"
    ref_abs = str((db_dir / "sequences.fa").resolve())
    gff_abs = str((db_dir / "genes.gff").resolve())
    map_abs = str((db_dir / "chromosomes.map").resolve())

    config_content = f"""# SnpEff configuration file for Panthera tigris
data.dir = {data_dir_abs}
panthera_tigris.genome : Panthera tigris
panthera_tigris.reference : {ref_abs}
panthera_tigris.genes.gff : {gff_abs}
panthera_tigris.chromosomes.map : {map_abs}
"""

    config_path = root / "tools" / "snpEff" / "snpEff.config"
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(config_content)

    print("SnpEff database setup complete!")
    print(f"Database directory: {db_dir}")
    print(f"Configuration file updated: {config_path}")

if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    setup_snpeff_database()
