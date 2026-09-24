import math

from freemol_mcp.tools.ch4sym2cart import ch4sym2cart_apply_displacement


def test_symmetric_stretch():
    # data/ch4sym2cart/tests/ch4_s1.inp's values.
    result = ch4sym2cart_apply_displacement(0.08, 0, 0, 0, 0, 0, 0, 0, 0)
    assert "error" not in result
    coords = result["cartesian"]
    assert coords["C"] == [0.0, 0.0, 0.0]
    assert coords["H1"] == [1.1, 0.0, 0.0]

    for label in ("H2", "H3", "H4"):
        x, y, z = coords[label]
        dist = math.sqrt(x * x + y * y + z * z)
        assert abs(dist - 1.1) < 1e-3

    h1, h2 = coords["H1"], coords["H2"]
    d1 = math.sqrt(sum(c * c for c in h1))
    d2 = math.sqrt(sum(c * c for c in h2))
    cos_angle = sum(a * b for a, b in zip(h1, h2)) / (d1 * d2)
    angle = math.degrees(math.acos(cos_angle))
    assert abs(angle - 109.47) < 0.1

    assert result["citation"]["url"].startswith(
        "https://github.com/mariotti/freemol/blob/"
    )
