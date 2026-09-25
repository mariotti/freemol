"""Generic tools for freemol's shared `[section-name]` input-file format
(see README.md's "Sectioned input files" section). Unlike the other
tools in this package, these don't call a built binary -- see
`freemol_mcp.sections`'s module docstring for why that's fine here.
"""

from __future__ import annotations

from pathlib import Path

from ..citation import Citation
from ..sections import SectionNotFoundError, list_sections, read_section

CITATION = Citation(
    routine="osec_set (section-header scan) + the section-body boundary "
    "convention used throughout, e.g. csmmod.F90's x-csmg-symop reader",
    file="Freemol/modules/osec.F90",
    line_start=191,
    line_end=296,
    note=(
        "This tool is a Python reimplementation of osec_set's scanning "
        "convention, not a subprocess call into a built binary -- there's "
        "no numerics here, just text scanning. One deliberate divergence: "
        "osec_set matches by prefix against the rightmost '[' on a header "
        "line (index(line,'[',.true.), an accidental quirk of the "
        "Fortran); this tool takes the text between the first '[' and its "
        "matching ']' instead -- equivalent for every fixture under "
        "Freemol/data/, none of which put a second '[' before the "
        "section name."
    ),
)


def _read_text(file_content: str | None, file_path: str | None) -> str:
    if (file_content is None) == (file_path is None):
        raise ValueError("pass exactly one of file_content or file_path")
    if file_content is not None:
        return file_content
    return Path(file_path).read_text()


def list_freemol_sections(
    file_content: str | None = None,
    file_path: str | None = None,
) -> dict:
    """List every section in a freemol-format input file, in file order.

    Pass exactly one of file_content (the file's text inline) or
    file_path (a path to read it from -- absolute, or relative to the
    server's working directory). Returns each section's name, its inline
    header parameters (e.g. "nrec=5 format=nscxyz" for `[molecule]
    nrec=5 format=nscxyz`), the 1-based line its header is on, and how
    many body lines follow it before the next section or end of file.
    """
    text = _read_text(file_content, file_path)
    sections = [
        {
            "name": s.name,
            "params": s.params,
            "header_line": s.header_line,
            "body_line_count": s.body_line_count,
        }
        for s in list_sections(text)
    ]
    return {"sections": sections, "citation": CITATION.as_dict()}


def read_freemol_section(
    section_name: str,
    file_content: str | None = None,
    file_path: str | None = None,
    occurrence: int = 1,
) -> dict:
    """Read one named section's raw body out of a freemol-format file.

    section_name is matched case-insensitively (e.g. "molecule" matches
    a `[molecule]` header). occurrence selects which match if the same
    section name appears more than once (1 = first; no fixture in this
    repo actually repeats a name, but the underlying format doesn't
    forbid it). Pass exactly one of file_content or file_path, as in
    list_freemol_sections.

    Returns the section's inline header parameters plus its body both as
    raw text (`body`) and as a list of lines (`body_lines`) -- the body
    is everything between this section's header and the next `[...]`
    header or end of file, unparsed (freemol's own programs each know
    how to read their own sections' row shapes; this tool doesn't guess
    at that).
    """
    text = _read_text(file_content, file_path)
    try:
        section = read_section(text, section_name, occurrence)
    except SectionNotFoundError as exc:
        return {
            "error": str(exc),
            "sections_present": exc.sections_present,
            "citation": CITATION.as_dict(),
        }
    return {
        "name": section.name,
        "params": section.params,
        "body": section.body,
        "body_lines": section.body_lines,
        "header_line": section.header_line,
        "body_start_line": section.body_start_line,
        "body_end_line": section.body_end_line,
        "citation": CITATION.as_dict(),
    }
