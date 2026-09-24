"""Run freemol's built binaries and capture their output."""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


class BinaryNotBuiltError(RuntimeError):
    """Raised when a freemol binary hasn't been built yet."""


def repo_root() -> Path:
    """The freemol git repo root (mcp/ is normally a direct child of it).

    Set FREEMOL_ROOT to override -- needed if freemol-mcp was installed
    non-editably (e.g. `pip install ./mcp` rather than `pip install -e
    ./mcp`), since then this file no longer lives inside the checkout it
    needs to find binaries in.
    """
    override = os.environ.get("FREEMOL_ROOT")
    if override:
        return Path(override).resolve()
    return Path(__file__).resolve().parents[3]


def freemol_dir() -> Path:
    return repo_root() / "Freemol"


def binary_path(name: str) -> Path:
    path = freemol_dir() / "bin" / f"{name}.exe"
    if not path.is_file():
        raise BinaryNotBuiltError(
            f"{path} does not exist. Build freemol first -- see README.md's "
            '"Build and test" section (cd Freemol && ./config/configure '
            "... && make freemol)."
        )
    return path


@dataclass
class RunResult:
    stdout: str
    stderr: str
    output_file: str
    returncode: int


def run(name: str, input_text: str, timeout: float = 60.0) -> RunResult:
    """Run Freemol/bin/<name>.exe -i <input> -o <output> in a scratch dir.

    Mirrors Freemol/tests/run_csmg_regression.sh: each call gets its own
    temp directory, since these programs can write scratch files (e.g.
    CSMG's CSMD.bufferInput) into the current directory.
    """
    exe = binary_path(name)
    with tempfile.TemporaryDirectory(prefix=f"freemol-mcp-{name}-") as workdir:
        work = Path(workdir)
        infile = work / "input.txt"
        outfile = work / "output.txt"
        infile.write_text(input_text)
        proc = subprocess.run(
            [str(exe), "-i", str(infile), "-o", str(outfile)],
            cwd=work,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        output_text = outfile.read_text() if outfile.exists() else ""
    return RunResult(
        stdout=proc.stdout,
        stderr=proc.stderr,
        output_file=output_text,
        returncode=proc.returncode,
    )
