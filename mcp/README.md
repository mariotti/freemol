# freemol-mcp

MCP tools wrapping freemol's coordinate-transformation programs. Every
tool response includes a `citation` field: the exact Fortran routine
(file + line range), freemol's own version, and a GitHub permalink
pinned to the commit the server is running from.

Currently wraps `XY4PolySphere` (`poly2cart`), `XY4Coord` and
`ch4sym2cart` -- see the root [README.md](../README.md)'s Programs
section for what each one does. `fit1Dpol` and `CSMG` aren't wrapped
yet.

This is a thin subprocess wrapper: it does not reimplement any of the
numerics in Python. It calls the same built binaries
(`Freemol/bin/*.exe`) that `Freemol/tests/run_smoke.sh` regression-tests
in CI.

Two more tools, `list_freemol_sections` and `read_freemol_section`, give
generic access to the `[section-name]` input-file format every program
above reads (see the root README's
["Sectioned input files"](../README.md#sectioned-input-files-an-old-ini-style-format-from-before-molden)
section). These are a genuine exception to the "no numerics
reimplemented" rule above -- there's no floating point involved, just
text scanning, so they're a direct Python mirror of `osec_set`'s
scanning convention rather than a subprocess call, and don't need
freemol built at all.

## Build freemol first

These tools need `Freemol/bin/{XY4PolySphere,XY4coord,ch4sym2cart}.exe`
to already exist (note the lowercase "c" in `XY4coord.exe` -- the program
directory is `XY4Coord`, but the built binary isn't). See the root
README's
["Build and test"](../README.md#build-and-test) section:

    cd Freemol
    ./config/configure m_generic_linux gfortran $PWD   # macOS: m_generic_osx
    make freemol

## Install and run

    cd mcp
    python3 -m venv .venv && source .venv/bin/activate
    pip install -e .
    freemol-mcp                # runs the server over stdio

## Add to an MCP client

Point the client at the `freemol-mcp` command (after `pip install -e .`
in an active virtualenv), or run it directly:

```json
{
  "mcpServers": {
    "freemol": {
      "command": "/absolute/path/to/freemol/mcp/.venv/bin/freemol-mcp"
    }
  }
}
```

## Tests

    pip install -e ".[dev]"
    pytest

Each test checks a tool against the exact numeric values already
validated in `Freemol/data/*/tests/*.inp` and `Freemol/tests/run_smoke.sh`
(e.g. `xy4polysphere_to_cartesian` with the tetrahedral-CH4 input must
return bonds=1.1000, angles=109.4712).
