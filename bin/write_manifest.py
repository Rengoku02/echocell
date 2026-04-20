#!/usr/bin/env python3
"""Write manifest.json summarizing a pipeline run.

Captures enough provenance to reproduce or audit the run: resolved parameters,
input/output SHA-256 hashes, env versions, seed/threads, reference versions.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


def sha256_of_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def file_summary(path: Path) -> dict:
    if not path.exists():
        return {"path": str(path), "exists": False}
    return {
        "path": str(path),
        "sha256": sha256_of_file(path),
        "size_bytes": path.stat().st_size,
    }


def git_sha(cwd: Path) -> str:
    if shutil.which("git") is None:
        return "unknown"
    for candidate in (cwd, cwd.parent):
        try:
            out = subprocess.run(
                ["git", "-C", str(candidate), "rev-parse", "HEAD"],
                capture_output=True, text=True, check=True, timeout=3,
            )
            return out.stdout.strip()
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            continue
    return "unknown"


def r_version(rscript_bin: str) -> str:
    try:
        out = subprocess.run(
            [rscript_bin, "-e", 'cat(paste(R.Version()$major, R.Version()$minor, sep="."))'],
            capture_output=True, text=True, check=True, timeout=10,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def conda_packages(env_prefix: str, names: list[str]) -> dict:
    try:
        out = subprocess.run(
            ["conda", "list", "--prefix", env_prefix, "--json"],
            capture_output=True, text=True, check=True, timeout=30,
        )
        pkgs = {p["name"]: p["version"] for p in json.loads(out.stdout)}
        return {n: pkgs.get(n, "absent") for n in names}
    except Exception:
        return {n: "unknown" for n in names}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True, help="Pipeline --outdir")
    p.add_argument("--input", required=True, help="Original input file path")
    p.add_argument("--params-json", required=True,
                   help="JSON file with resolved CLI parameters")
    p.add_argument("--outputs", required=True,
                   help="Comma-separated list of output file paths")
    p.add_argument("--conda-lock", default="", help="Path to conda.lock.txt (optional)")
    p.add_argument("--doublet-sha", default="unknown",
                   help="DoubletFinder commit SHA pinned in setup.sh")
    p.add_argument("--started-at", default="", help="ISO8601 pipeline start time")
    args = p.parse_args()

    outdir = Path(args.outdir).resolve()
    params = json.loads(Path(args.params_json).read_text())
    output_files = [Path(f).resolve() for f in args.outputs.split(",") if f.strip()]

    env_prefix = os.environ.get("CONDA_PREFIX", "")
    env_name = os.path.basename(env_prefix) if env_prefix else "unknown"

    pkg_versions = conda_packages(
        env_prefix,
        ["scanpy", "anndata", "harmonypy", "numpy", "scipy", "scikit-learn",
         "r-base", "r-seurat", "r-anndata", "bioconductor-singler",
         "bioconductor-celldex"],
    ) if env_prefix else {}

    manifest = {
        "pipeline_version": git_sha(Path(__file__).resolve().parent),
        "started_at": args.started_at or None,
        "finished_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "host": {
            "hostname": platform.node(),
            "os": f"{platform.system()} {platform.release()}",
            "arch": platform.machine(),
        },
        "environment": {
            "conda_env": env_name,
            "conda_prefix": env_prefix,
            "python_version": platform.python_version(),
            "r_version": r_version(os.path.join(env_prefix, "bin", "Rscript"))
                         if env_prefix else "unknown",
            "seed": params.get("seed"),
            "threads": params.get("threads"),
            "thread_env": {
                k: os.environ.get(k)
                for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                          "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                          "VECLIB_MAXIMUM_THREADS", "PYTHONHASHSEED")
            },
        },
        "parameters": params,
        "input": file_summary(Path(args.input).resolve()),
        "outputs": [file_summary(f) for f in output_files],
        "references": {
            "singler_ref": params.get("singler_ref"),
            "singler_labels": params.get("singler_labels"),
            "celldex_version": pkg_versions.get("bioconductor-celldex", "unknown"),
            "singler_version": pkg_versions.get("bioconductor-singler", "unknown"),
            "doubletfinder_sha": args.doublet_sha,
        },
        "packages": pkg_versions,
        "conda_lock": (
            file_summary(Path(args.conda_lock).resolve())
            if args.conda_lock else None
        ),
    }

    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    print(f"[write_manifest] wrote {manifest_path}")


if __name__ == "__main__":
    main()
