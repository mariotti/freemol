from freemol_mcp.tools.xy4polysphere import xy4polysphere_to_cartesian


def test_tetrahedral_ch4():
    # data/XY4PolySphere/tests/ch4_poly.inp's values.
    result = xy4polysphere_to_cartesian(
        1.10, 1.10, 1.10, 1.10, 109.4712, 109.4712, 109.4712, 60.0, 60.0
    )
    assert "error" not in result
    assert result["bonds"] == [1.1, 1.1, 1.1, 1.1]
    for angle in result["angles_degrees"].values():
        assert angle == 109.4712
    assert result["metpot4"] == 107.6015
    assert result["metpot4"] > result["reference_geometry_metpot4"]
    assert len(result["cartesian"]) == 5
    assert result["cartesian"][0]["element"] == "C"
    assert result["citation"]["source"].endswith(
        "XY4PolySphere.F90:1580-1696"
    )
    assert result["citation"]["url"].startswith(
        "https://github.com/mariotti/freemol/blob/"
    )


def test_reference_molecule_only_affects_orientation():
    from freemol_mcp.molecule import Atom

    base = xy4polysphere_to_cartesian(
        1.10, 1.10, 1.10, 1.10, 109.4712, 109.4712, 109.4712, 60.0, 60.0
    )
    alt_molecule = (
        Atom("C", 0.0, 0.0, 0.0),
        Atom("H", 1.0, 0.0, 0.0),
        Atom("H", 0.0, 1.0, 0.0),
        Atom("H", 0.0, 0.0, 1.0),
        Atom("H", -1.0, -1.0, -1.0),
    )
    alt = xy4polysphere_to_cartesian(
        1.10,
        1.10,
        1.10,
        1.10,
        109.4712,
        109.4712,
        109.4712,
        60.0,
        60.0,
        reference_molecule=alt_molecule,
    )
    assert alt["bonds"] == base["bonds"]
    assert alt["angles_degrees"] == base["angles_degrees"]
    assert alt["metpot4"] == base["metpot4"]
