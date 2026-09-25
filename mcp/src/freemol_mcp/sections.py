"""Pure-Python reader for freemol's shared sectioned input-file format.

This mirrors the scanning convention `osec_set` implements
(`Freemol/modules/osec.F90:191-296`) and the section-body boundary every
per-section reader uses (e.g. `Freemol/programs/CSMG/csmmod.F90:770-787`
for `x-csmg-symop`): a line whose first character is `[` starts a new
section, named by the text between that `[` and its matching `]`, with
whatever trails the `]` on that line as the section's own inline
parameters; the section's body is every line after that up to the next
such header or end of file.

Unlike the coordinate-transformation tools, this does not call a built
Fortran binary -- there's no numerics here, just text scanning, so a
direct Python implementation doesn't touch the "reproducible double
precision" claim the rest of this package is careful about.

One deliberate divergence from `osec_set`: it prefix-matches against the
*rightmost* `[` on a header line (`index(line,'[',.true.)`, an accidental
quirk of the Fortran, not a documented format guarantee). This module
takes the text between the *first* `[` and its matching `]` instead --
simpler, and behaviorally identical for every fixture under
`Freemol/data/`, none of which put a second `[` before the section name.
"""

from __future__ import annotations

from dataclasses import dataclass


class SectionNotFoundError(ValueError):
    def __init__(self, name: str, occurrence: int, sections_present: list[str]):
        self.name = name
        self.occurrence = occurrence
        self.sections_present = sections_present
        super().__init__(
            f"section {name!r} (occurrence {occurrence}) not found; "
            f"sections present: {sections_present}"
        )


@dataclass(frozen=True)
class SectionHeader:
    name: str
    params: str
    header_line: int
    body_line_count: int


@dataclass(frozen=True)
class SectionBody:
    name: str
    params: str
    body: str
    body_lines: list[str]
    header_line: int
    body_start_line: int
    body_end_line: int


def _header_line_numbers(lines: list[str]) -> list[int]:
    return [i for i, line in enumerate(lines) if line.startswith("[")]


def _split_header(line: str) -> tuple[str, str]:
    """('[name] params' | '[name]params') -> (name, params)."""
    close = line.index("]")
    name = line[1:close]
    params = line[close + 1 :].strip()
    return name, params


def list_sections(text: str) -> list[SectionHeader]:
    lines = text.splitlines()
    header_idxs = _header_line_numbers(lines)
    headers = []
    for pos, idx in enumerate(header_idxs):
        name, params = _split_header(lines[idx])
        next_idx = header_idxs[pos + 1] if pos + 1 < len(header_idxs) else len(lines)
        body_line_count = next_idx - idx - 1
        headers.append(
            SectionHeader(
                name=name,
                params=params,
                header_line=idx + 1,
                body_line_count=body_line_count,
            )
        )
    return headers


def read_section(text: str, name: str, occurrence: int = 1) -> SectionBody:
    lines = text.splitlines()
    header_idxs = _header_line_numbers(lines)
    target = name.strip().lower()

    matches = []
    for pos, idx in enumerate(header_idxs):
        section_name, params = _split_header(lines[idx])
        if section_name.strip().lower() == target:
            next_idx = (
                header_idxs[pos + 1] if pos + 1 < len(header_idxs) else len(lines)
            )
            matches.append((idx, section_name, params, next_idx))

    if occurrence < 1 or occurrence > len(matches):
        sections_present = [_split_header(lines[idx])[0] for idx in header_idxs]
        raise SectionNotFoundError(name, occurrence, sections_present)

    idx, section_name, params, next_idx = matches[occurrence - 1]
    body_lines = lines[idx + 1 : next_idx]
    return SectionBody(
        name=section_name,
        params=params,
        body="\n".join(body_lines),
        body_lines=body_lines,
        header_line=idx + 1,
        body_start_line=idx + 2,
        body_end_line=next_idx,
    )
