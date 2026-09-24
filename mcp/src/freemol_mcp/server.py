"""freemol MCP server: coordinate-transformation tools wrapping freemol's
built Fortran binaries. Requires freemol to be built first -- see the
root README.md's "Build and test" section.
"""

from __future__ import annotations

from mcp.server.mcpserver import MCPServer

from .tools.ch4sym2cart import ch4sym2cart_apply_displacement
from .tools.xy4coord import xy4coord_apply_displacement
from .tools.xy4polysphere import xy4polysphere_to_cartesian

mcp = MCPServer(
    "freemol",
    instructions=(
        "Coordinate-transformation tools wrapping freemol's Fortran "
        "programs (built from Freemol/bin/*.exe). Every result includes "
        "a citation pointing at the exact routine, file:lines, freemol "
        "version and commit that produced it."
    ),
)

mcp.tool()(xy4polysphere_to_cartesian)
mcp.tool()(xy4coord_apply_displacement)
mcp.tool()(ch4sym2cart_apply_displacement)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
