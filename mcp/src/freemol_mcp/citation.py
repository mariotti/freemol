"""Build a citation pointing at the exact Fortran routine that produced a
tool's output: file:lines, freemol's own version, and a GitHub permalink
pinned to the commit this server is actually running from.
"""

from __future__ import annotations

import functools
import subprocess
from dataclasses import dataclass

from .binary import freemol_dir, repo_root

GITHUB_REPO = "mariotti/freemol"


@functools.lru_cache(maxsize=1)
def _commit_sha() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root(),
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


@functools.lru_cache(maxsize=1)
def freemol_version() -> str:
    """freemol's own version string, from config/printversion."""
    result = subprocess.run(
        ["sh", "./config/printversion", "long"],
        cwd=freemol_dir(),
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


@dataclass(frozen=True)
class Citation:
    routine: str
    file: str  # path relative to the repo root, e.g. Freemol/programs/.../X.F90
    line_start: int
    line_end: int
    note: str | None = None

    def as_dict(self) -> dict:
        sha = _commit_sha()
        lines = (
            f"L{self.line_start}"
            if self.line_start == self.line_end
            else f"L{self.line_start}-L{self.line_end}"
        )
        url = f"https://github.com/{GITHUB_REPO}/blob/{sha}/{self.file}#{lines}"
        result = {
            "routine": self.routine,
            "source": f"{self.file}:{self.line_start}-{self.line_end}",
            "freemol_version": freemol_version(),
            "commit": sha,
            "url": url,
        }
        if self.note:
            result["note"] = self.note
        return result
