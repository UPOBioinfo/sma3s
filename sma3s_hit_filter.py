#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
import re
import sqlite3
import sys
from pathlib import Path
from datetime import datetime


def open_maybe_gzip(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def detect_delimiter(path):
    with open_maybe_gzip(path, "rt") as handle:
        for line in handle:
            if line.strip():
                if "\t" in line:
                    return "\t"
                if "," in line:
                    return ","
                return "\t"
    return "\t"


def normalize_colname(name):
    return name.strip().lstrip("#").upper()


def build_column_map(fieldnames):
    return {normalize_colname(c): c for c in fieldnames}


def split_semicolon_field(value):
    if value is None:
        return []
    value = value.strip()
    if value == "":
        return []
    return [x.strip() for x in value.split(";")]


def parse_float(value):
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def is_perfect_hit(identity, qcov, tcov, min_identity, min_qcov, min_tcov):
    identity = parse_float(identity)
    qcov = parse_float(qcov)
    tcov = parse_float(tcov)

    if identity is None or qcov is None or tcov is None:
        return False

    return (
        identity >= min_identity
        and qcov >= min_qcov
        and tcov >= min_tcov
    )


def get_required_columns(fieldnames):
    colmap = build_column_map(fieldnames)

    required = [
        "ANNOTATION_HIT_SEQUENCE",
        "ANNOTATION_HIT_IDENTITY",
        "ANNOTATION_HIT_QCOV",
        "ANNOTATION_HIT_TCOV",
    ]

    missing = [c for c in required if c not in colmap]
    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
            + "\nDetected columns: "
            + ", ".join(fieldnames)
        )

    return colmap


def extract_gene_id(row, colmap):
    if "ID" in colmap:
        return row.get(colmap["ID"], "")
    return ""


def extract_gene_name(row, colmap):
    if "GENENAME" in colmap:
        return row.get(colmap["GENENAME"], "")
    return ""


def extract_description(row, colmap):
    if "DESCRIPTION" in colmap:
        return row.get(colmap["DESCRIPTION"], "")
    return ""


def normalize_species(os_lines):
    species = " ".join(os_lines).strip()
    species = " ".join(species.split())
    if species.endswith("."):
        species = species[:-1]
    return species if species else "NA"


def build_uniprot_species_index(
    uniprot_dat,
    species_db,
    overwrite=False,
    batch_size=50000,
    progress_every=100000,
):
    species_db = Path(species_db)

    if species_db.exists():
        if overwrite:
            species_db.unlink()
        else:
            raise FileExistsError(
                f"SQLite index already exists: {species_db}\n"
                f"Use --overwrite if you want to recreate it."
            )

    print(f"[INFO] Creating SQLite index: {species_db}", file=sys.stderr)
    print(f"[INFO] Reading UniProt .dat: {uniprot_dat}", file=sys.stderr)

    conn = sqlite3.connect(str(species_db))
    cur = conn.cursor()

    cur.execute("PRAGMA journal_mode=WAL;")
    cur.execute("PRAGMA synchronous=OFF;")
    cur.execute("PRAGMA temp_store=MEMORY;")
    cur.execute("PRAGMA cache_size=-200000;")

    cur.execute("""
        CREATE TABLE uniprot_species (
            uniprot_id TEXT PRIMARY KEY,
            species TEXT,
            taxid TEXT,
            entry_name TEXT,
            id_type TEXT
        ) WITHOUT ROWID;
    """)

    cur.execute("""
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        );
    """)

    cur.execute(
        "INSERT INTO metadata(key, value) VALUES (?, ?)",
        ("source_dat", str(uniprot_dat)),
    )
    cur.execute(
        "INSERT INTO metadata(key, value) VALUES (?, ?)",
        ("created_at", datetime.now().isoformat(timespec="seconds")),
    )

    conn.commit()

    entry_name = None
    accessions = []
    os_lines = []
    taxid = "NA"

    batch = []
    n_records = 0
    n_ids = 0

    taxid_re = re.compile(r"NCBI_TaxID=(\d+)")

    def flush_record():
        nonlocal batch, n_records, n_ids
        nonlocal entry_name, accessions, os_lines, taxid

        if not entry_name:
            return

        species = normalize_species(os_lines)

        rows = []

        rows.append((
            entry_name,
            species,
            taxid,
            entry_name,
            "ID",
        ))

        for ac in accessions:
            rows.append((
                ac,
                species,
                taxid,
                entry_name,
                "AC",
            ))

        batch.extend(rows)
        n_records += 1
        n_ids += len(rows)

        if len(batch) >= batch_size:
            cur.executemany(
                """
                INSERT OR REPLACE INTO uniprot_species
                (uniprot_id, species, taxid, entry_name, id_type)
                VALUES (?, ?, ?, ?, ?)
                """,
                batch,
            )
            conn.commit()
            batch = []

        if progress_every and n_records % progress_every == 0:
            print(
                f"[INFO] Parsed {n_records:,} UniProt records; "
                f"indexed {n_ids:,} IDs/accessions...",
                file=sys.stderr,
            )

    with open_maybe_gzip(uniprot_dat, "rt") as handle:
        for line in handle:
            if line.startswith("ID   "):
                entry_name = line[5:].split()[0].strip()
                accessions = []
                os_lines = []
                taxid = "NA"

            elif line.startswith("AC   "):
                ac_text = line[5:].strip()
                acs = [x.strip() for x in ac_text.split(";") if x.strip()]
                accessions.extend(acs)

            elif line.startswith("OS   "):
                os_lines.append(line[5:].strip())

            elif line.startswith("OX   "):
                match = taxid_re.search(line)
                if match:
                    taxid = match.group(1)

            elif line.startswith("//"):
                flush_record()
                entry_name = None
                accessions = []
                os_lines = []
                taxid = "NA"

    if batch:
        cur.executemany(
            """
            INSERT OR REPLACE INTO uniprot_species
            (uniprot_id, species, taxid, entry_name, id_type)
            VALUES (?, ?, ?, ?, ?)
            """,
            batch,
        )
        conn.commit()

    cur.execute(
        "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
        ("n_records", str(n_records)),
    )
    cur.execute(
        "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
        ("n_indexed_ids", str(n_ids)),
    )
    conn.commit()

    conn.close()

    print("[DONE] SQLite index created.", file=sys.stderr)
    print(f"[DONE] UniProt records parsed: {n_records:,}", file=sys.stderr)
    print(f"[DONE] IDs/accessions indexed: {n_ids:,}", file=sys.stderr)


def collect_perfect_hit_ids(
    input_tsv,
    delimiter,
    min_identity,
    min_qcov,
    min_tcov,
):
    perfect_hit_ids = set()

    total_genes = 0
    genes_with_perfect_hit = 0
    total_perfect_hits = 0
    rows_with_mismatched_hit_fields = 0

    with open_maybe_gzip(input_tsv, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("Input file has no header.")

        colmap = get_required_columns(reader.fieldnames)

        hit_col = colmap["ANNOTATION_HIT_SEQUENCE"]
        identity_col = colmap["ANNOTATION_HIT_IDENTITY"]
        qcov_col = colmap["ANNOTATION_HIT_QCOV"]
        tcov_col = colmap["ANNOTATION_HIT_TCOV"]

        for row in reader:
            total_genes += 1

            hits = split_semicolon_field(row.get(hit_col, ""))
            identities = split_semicolon_field(row.get(identity_col, ""))
            qcovs = split_semicolon_field(row.get(qcov_col, ""))
            tcovs = split_semicolon_field(row.get(tcov_col, ""))

            lengths = [len(hits), len(identities), len(qcovs), len(tcovs)]

            if len(set(lengths)) != 1:
                rows_with_mismatched_hit_fields += 1

            n = min(lengths)
            row_has_perfect_hit = False

            for i in range(n):
                if is_perfect_hit(
                    identities[i],
                    qcovs[i],
                    tcovs[i],
                    min_identity,
                    min_qcov,
                    min_tcov,
                ):
                    if hits[i]:
                        perfect_hit_ids.add(hits[i])
                        total_perfect_hits += 1
                        row_has_perfect_hit = True

            if row_has_perfect_hit:
                genes_with_perfect_hit += 1

    stats = {
        "total_genes": total_genes,
        "genes_with_perfect_hit": genes_with_perfect_hit,
        "total_perfect_hits": total_perfect_hits,
        "unique_perfect_hit_ids": len(perfect_hit_ids),
        "rows_with_mismatched_hit_fields": rows_with_mismatched_hit_fields,
    }

    return perfect_hit_ids, stats


def load_species_from_sqlite(species_db, wanted_ids, chunk_size=900):
    species_by_id = {}

    if not wanted_ids:
        return species_by_id

    conn = sqlite3.connect(str(species_db))
    cur = conn.cursor()

    wanted_ids = list(wanted_ids)

    for i in range(0, len(wanted_ids), chunk_size):
        chunk = wanted_ids[i:i + chunk_size]
        placeholders = ",".join(["?"] * len(chunk))

        query = f"""
            SELECT uniprot_id, species, taxid
            FROM uniprot_species
            WHERE uniprot_id IN ({placeholders})
        """

        cur.execute(query, chunk)

        for uniprot_id, species, taxid in cur.fetchall():
            species_by_id[uniprot_id] = {
                "species": species if species else "NA",
                "taxid": taxid if taxid else "NA",
            }

    conn.close()
    return species_by_id


def write_outputs(
    input_tsv,
    delimiter,
    output_prefix,
    species_by_id,
    min_identity,
    min_qcov,
    min_tcov,
):
    output_prefix = Path(output_prefix)

    genes_out = Path(str(output_prefix) + ".perfect_genes.tsv")
    hits_out = Path(str(output_prefix) + ".perfect_hits.long.tsv")
    summary_out = Path(str(output_prefix) + ".summary.tsv")

    total_genes = 0
    genes_with_perfect_hit = 0
    total_perfect_hits = 0
    rows_with_mismatched_hit_fields = 0
    perfect_hits_without_species = 0

    with open_maybe_gzip(input_tsv, "rt") as in_handle, \
            open(genes_out, "w", newline="") as genes_handle, \
            open(hits_out, "w", newline="") as hits_handle:

        reader = csv.DictReader(in_handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("Input file has no header.")

        colmap = get_required_columns(reader.fieldnames)

        hit_col = colmap["ANNOTATION_HIT_SEQUENCE"]
        identity_col = colmap["ANNOTATION_HIT_IDENTITY"]
        qcov_col = colmap["ANNOTATION_HIT_QCOV"]
        tcov_col = colmap["ANNOTATION_HIT_TCOV"]

        alnlen_col = colmap.get("ANNOTATION_HIT_ALNLEN")
        evalue_col = colmap.get("ANNOTATION_HIT_EVALUE")

        extra_gene_cols = [
            "N_PERFECT_HITS",
            "PERFECT_HIT_IDS",
            "PERFECT_HIT_SPECIES",
            "PERFECT_HIT_TAXIDS",
        ]

        genes_writer = csv.DictWriter(
            genes_handle,
            delimiter="\t",
            fieldnames=reader.fieldnames + extra_gene_cols,
            extrasaction="ignore",
        )
        genes_writer.writeheader()

        hits_fieldnames = [
            "ID",
            "GENENAME",
            "DESCRIPTION",
            "HIT_ID",
            "HIT_IDENTITY",
            "HIT_QCOV",
            "HIT_TCOV",
            "HIT_ALNLEN",
            "HIT_EVALUE",
            "HIT_SPECIES",
            "HIT_TAXID",
        ]

        hits_writer = csv.DictWriter(
            hits_handle,
            delimiter="\t",
            fieldnames=hits_fieldnames,
        )
        hits_writer.writeheader()

        for row in reader:
            total_genes += 1

            hits = split_semicolon_field(row.get(hit_col, ""))
            identities = split_semicolon_field(row.get(identity_col, ""))
            qcovs = split_semicolon_field(row.get(qcov_col, ""))
            tcovs = split_semicolon_field(row.get(tcov_col, ""))

            alnlens = split_semicolon_field(row.get(alnlen_col, "")) if alnlen_col else []
            evalues = split_semicolon_field(row.get(evalue_col, "")) if evalue_col else []

            lengths = [len(hits), len(identities), len(qcovs), len(tcovs)]

            if len(set(lengths)) != 1:
                rows_with_mismatched_hit_fields += 1

            n = min(lengths)

            perfect_hit_ids = []
            perfect_hit_species = []
            perfect_hit_taxids = []

            gene_id = extract_gene_id(row, colmap)
            gene_name = extract_gene_name(row, colmap)
            description = extract_description(row, colmap)

            for i in range(n):
                if is_perfect_hit(
                    identities[i],
                    qcovs[i],
                    tcovs[i],
                    min_identity,
                    min_qcov,
                    min_tcov,
                ):
                    hit_id = hits[i]

                    species_info = species_by_id.get(
                        hit_id,
                        {"species": "NA", "taxid": "NA"},
                    )

                    species = species_info["species"]
                    taxid = species_info["taxid"]

                    if species == "NA":
                        perfect_hits_without_species += 1

                    alnlen = alnlens[i] if i < len(alnlens) else ""
                    evalue = evalues[i] if i < len(evalues) else ""

                    perfect_hit_ids.append(hit_id)
                    perfect_hit_species.append(species)
                    perfect_hit_taxids.append(taxid)

                    total_perfect_hits += 1

                    hits_writer.writerow({
                        "ID": gene_id,
                        "GENENAME": gene_name,
                        "DESCRIPTION": description,
                        "HIT_ID": hit_id,
                        "HIT_IDENTITY": identities[i],
                        "HIT_QCOV": qcovs[i],
                        "HIT_TCOV": tcovs[i],
                        "HIT_ALNLEN": alnlen,
                        "HIT_EVALUE": evalue,
                        "HIT_SPECIES": species,
                        "HIT_TAXID": taxid,
                    })

            if perfect_hit_ids:
                genes_with_perfect_hit += 1

                unique_species = []
                seen_species = set()

                for sp in perfect_hit_species:
                    if sp not in seen_species:
                        unique_species.append(sp)
                        seen_species.add(sp)

                unique_taxids = []
                seen_taxids = set()

                for tx in perfect_hit_taxids:
                    if tx not in seen_taxids:
                        unique_taxids.append(tx)
                        seen_taxids.add(tx)

                row["N_PERFECT_HITS"] = str(len(perfect_hit_ids))
                row["PERFECT_HIT_IDS"] = ";".join(perfect_hit_ids)
                row["PERFECT_HIT_SPECIES"] = ";".join(unique_species)
                row["PERFECT_HIT_TAXIDS"] = ";".join(unique_taxids)

                genes_writer.writerow(row)

    with open(summary_out, "w") as handle:
        handle.write("METRIC\tVALUE\n")
        handle.write(f"total_genes\t{total_genes}\n")
        handle.write(f"genes_with_perfect_hit\t{genes_with_perfect_hit}\n")
        handle.write(f"total_perfect_hits\t{total_perfect_hits}\n")
        handle.write(f"perfect_hits_without_species\t{perfect_hits_without_species}\n")
        handle.write(f"rows_with_mismatched_hit_fields\t{rows_with_mismatched_hit_fields}\n")
        handle.write(f"min_identity\t{min_identity}\n")
        handle.write(f"min_qcov\t{min_qcov}\n")
        handle.write(f"min_tcov\t{min_tcov}\n")

    return genes_out, hits_out, summary_out


def run_filter(args):
    if args.sep == "auto":
        delimiter = detect_delimiter(args.input)
    elif args.sep == "tab":
        delimiter = "\t"
    else:
        delimiter = ","

    print("[INFO] Collecting perfect hit IDs from annotation table...", file=sys.stderr)

    perfect_hit_ids, first_pass_stats = collect_perfect_hit_ids(
        input_tsv=args.input,
        delimiter=delimiter,
        min_identity=args.min_identity,
        min_qcov=args.min_qcov,
        min_tcov=args.min_tcov,
    )

    print(
        f"[INFO] Total genes: {first_pass_stats['total_genes']:,}",
        file=sys.stderr,
    )
    print(
        f"[INFO] Genes with at least one perfect hit: "
        f"{first_pass_stats['genes_with_perfect_hit']:,}",
        file=sys.stderr,
    )
    print(
        f"[INFO] Unique perfect hit IDs: "
        f"{first_pass_stats['unique_perfect_hit_ids']:,}",
        file=sys.stderr,
    )

    species_by_id = {}

    if args.species_db:
        print(f"[INFO] Loading species from SQLite index: {args.species_db}", file=sys.stderr)

        species_by_id = load_species_from_sqlite(
            species_db=args.species_db,
            wanted_ids=perfect_hit_ids,
        )

        print(
            f"[INFO] Species found for "
            f"{len(species_by_id):,} / {len(perfect_hit_ids):,} unique perfect hits.",
            file=sys.stderr,
        )
    else:
        print(
            "[INFO] No --species-db provided. Species will be reported as NA.",
            file=sys.stderr,
        )

    print("[INFO] Writing output files...", file=sys.stderr)

    genes_out, hits_out, summary_out = write_outputs(
        input_tsv=args.input,
        delimiter=delimiter,
        output_prefix=args.output_prefix,
        species_by_id=species_by_id,
        min_identity=args.min_identity,
        min_qcov=args.min_qcov,
        min_tcov=args.min_tcov,
    )

    print("[DONE] Output files:", file=sys.stderr)
    print(f"  - {genes_out}", file=sys.stderr)
    print(f"  - {hits_out}", file=sys.stderr)
    print(f"  - {summary_out}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter genes with at least one annotation hit having "
            "identity, query coverage and target/hit coverage above configurable "
            "minimum percentages (100/100/100 by default). "
            "Optionally use a local SQLite index built from a UniProt .dat file "
            "to recover species names."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 sma3s_hit_filter.py filter -i annotation_hitmetrics.tsv -o perfect_hits\n"
            "  python3 sma3s_hit_filter.py filter -i annotation_hitmetrics.tsv -o hits_90_90_90 "
            "--min-identity 90 --min-qcov 90 --min-tcov 90\n"
            "  python3 sma3s_hit_filter.py filter --help"
        ),
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser(
        "build-index",
        help="Build a local SQLite index from a UniProt .dat or .dat.gz file.",
    )

    build_parser.add_argument(
        "--uniprot-dat",
        required=True,
        help="Input UniProt .dat or .dat.gz file.",
    )

    build_parser.add_argument(
        "--species-db",
        required=True,
        help="Output SQLite DB. Example: uniprot_bacteria_id_species.sqlite",
    )

    build_parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing SQLite DB.",
    )

    build_parser.add_argument(
        "--batch-size",
        type=int,
        default=50000,
        help="Number of rows inserted per SQLite batch. Default: 50000",
    )

    build_parser.add_argument(
        "--progress-every",
        type=int,
        default=100000,
        help="Report progress every N UniProt records. Default: 100000",
    )

    filter_parser = subparsers.add_parser(
        "filter",
        help="Filter hits using configurable identity and coverage thresholds.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 sma3s_hit_filter.py filter -i annotation_hitmetrics.tsv -o perfect_hits\n"
            "  python3 sma3s_hit_filter.py filter -i annotation_hitmetrics.tsv -o hits_90_90_90 "
            "--min-identity 90 --min-qcov 90 --min-tcov 90"
        ),
    )

    filter_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input annotation table in TSV format. .gz accepted.",
    )

    filter_parser.add_argument(
        "-o", "--output-prefix",
        default="perfect_uniprot_hits",
        help="Output prefix. Default: perfect_uniprot_hits",
    )

    filter_parser.add_argument(
        "--species-db",
        default=None,
        help=(
            "Optional SQLite species index generated with build-index. "
            "Example: uniprot_bacteria_id_species.sqlite"
        ),
    )

    filter_parser.add_argument(
        "--sep",
        default="auto",
        choices=["auto", "tab", "comma"],
        help="Input separator. Default: auto.",
    )

    filter_parser.add_argument(
        "--min-identity",
        type=float,
        default=100.0,
        help="Minimum identity percentage. Default: 100.0",
    )

    filter_parser.add_argument(
        "--min-qcov",
        type=float,
        default=100.0,
        help="Minimum query coverage percentage. Default: 100.0",
    )

    filter_parser.add_argument(
        "--min-tcov",
        type=float,
        default=100.0,
        help="Minimum target/hit coverage percentage. Default: 100.0",
    )

    args = parser.parse_args()

    if args.command == "build-index":
        build_uniprot_species_index(
            uniprot_dat=args.uniprot_dat,
            species_db=args.species_db,
            overwrite=args.overwrite,
            batch_size=args.batch_size,
            progress_every=args.progress_every,
        )

    elif args.command == "filter":
        run_filter(args)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

