#!/usr/bin/env python3
"""
Check VASP single-point calculations:
  - Move finished runs to 'complete/'
  - Classify unfinished runs and prepare them for restart

A single-point run is considered finished when OUTCAR contains:
    "General timing and accounting informations for this job"
"""

import os
import shutil
from pathlib import Path


# ------------------------- status detection -------------------------

COMPLETION_MARKER = "General timing and accounting informations for this job"


def tail_text(path: Path, nbytes: int = 8192) -> str:
    """Return the last `nbytes` of a file as text (ignoring decoding errors)."""
    try:
        with path.open("rb") as f:
            try:
                f.seek(-nbytes, os.SEEK_END)
            except OSError:
                f.seek(0)
            return f.read().decode(errors="ignore")
    except OSError:
        return ""


def classify(folder: Path) -> str:
    """
    Return one of:
      'finished'   -> OUTCAR shows successful completion
      'not_started'-> no OUTCAR yet (never ran, or queued)
      'crashed'    -> OUTCAR exists but no completion marker (killed, OOM, walltime, error)
      'no_inputs'  -> required input files missing, can't restart
    """
    required = ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]
    if not all((folder / f).is_file() for f in required):
        return "no_inputs"

    outcar = folder / "OUTCAR"
    if not outcar.is_file():
        return "not_started"

    if COMPLETION_MARKER in tail_text(outcar):
        return "finished"

    return "crashed"


# ------------------------- restart preparation -------------------------

# Files from a previous (failed) run that should be removed before restart.
# For single-point we don't reuse WAVECAR/CHGCAR by default — safer to start
# fresh. If you DO want to reuse them (e.g. to save SCF iterations), remove
# them from this list.
STALE_FILES = [
    "OUTCAR", "OSZICAR", "vasprun.xml", "CONTCAR", "XDATCAR",
    "PCDAT", "REPORT", "vaspout.h5", "WAVECAR", "CHGCAR", "CHG",
    "DOSCAR", "EIGENVAL", "PROCAR", "LOCPOT", "ELFCAR", "WAVEDER",
    # stdout/stderr logs vary by scheduler; add yours if needed
]


def prepare_restart(folder: Path, dry_run: bool = False) -> None:
    """Remove stale output files so VASP can be rerun cleanly."""
    for name in STALE_FILES:
        p = folder / name
        if p.exists():
            if dry_run:
                print(f"  [dry-run] would remove {p.name}")
            else:
                p.unlink()


# ------------------------- main loop -------------------------

def main(base_dir: str = ".",
         complete_dir_name: str = "complete",
         start: int = 0, end: int = 1000,
         restart: bool = False,
         dry_run: bool = True) -> None:

    base = Path(base_dir).resolve()
    complete_dir = base / complete_dir_name
    complete_dir.mkdir(exist_ok=True)

    buckets = {"finished": [], "crashed": [], "not_started": [],
               "no_inputs": [], "missing": []}

    for i in range(start, end + 1):
        name = f"struc{i:05d}"          # adjust width if your naming differs
        folder = base / name

        if not folder.is_dir():
            buckets["missing"].append(name)
            continue

        status = classify(folder)
        buckets[status].append(name)

        if status == "finished":
            dest = complete_dir / name
            if dest.exists():
                print(f"  ! {dest.name} exists, skipping move")
                continue
            if dry_run:
                print(f"  [dry-run] would move {name} -> {complete_dir_name}/")
            else:
                shutil.move(str(folder), str(dest))
                print(f"  ✓ moved {name} -> {complete_dir_name}/")

        elif status == "crashed" and restart:
            print(f"  ↻ preparing {name} for restart")
            prepare_restart(folder, dry_run=dry_run)

    # ---- summary ----
    print("\n=== Summary ===")
    for k, v in buckets.items():
        print(f"  {k:12s}: {len(v)}")

    # write lists for crashed / not_started so you can resubmit them
    for k in ("crashed", "not_started", "no_inputs"):
        if buckets[k]:
            (base / f"{k}.txt").write_text("\n".join(buckets[k]))
            print(f"  -> wrote {k}.txt ({len(buckets[k])} entries)")


if __name__ == "__main__":
    # Recommended workflow:
    #   1) dry_run=True, restart=False  -> see what's finished/crashed
    #   2) dry_run=False, restart=False -> actually move finished ones
    #   3) dry_run=False, restart=True  -> clean crashed folders for resubmit
    #main(start=0, end=1000, restart=False, dry_run=True)
    #main(start=0, end=2000, restart=True, dry_run=False)
    main(start=0, end=2000, restart=False, dry_run=True)
