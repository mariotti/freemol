from freemol_mcp.tools.xy4coord import xy4coord_apply_displacement


def test_symmetric_stretch():
    # data/XY4Coord/tests/ch4_s1.inp's values.
    result = xy4coord_apply_displacement(0.08, 0, 0, 0, 0, 0, 0, 0, 0)
    assert "error" not in result
    assert result["bonds"] == [1.1, 1.1, 1.1, 1.1]
    for angle in result["angles_degrees"].values():
        assert angle == 109.4712
    assert len(result["cartesian"]) == 5
    assert result["cartesian"][0]["label"] == "X"
    assert result["citation"]["url"].startswith(
        "https://github.com/mariotti/freemol/blob/"
    )


def test_equilibrium():
    # data/XY4Coord/tests/equilibrium.inp's values.
    result = xy4coord_apply_displacement(0, 0, 0, 0, 0, 0, 0, 0, 0)
    assert "error" not in result
    assert result["bonds"] == [1.09, 1.09, 1.09, 1.09]


def test_known_issue_pure_angle_displacement_is_a_clean_error():
    # The documented Known Issue: this must not silently return numbers
    # the program itself flagged as wrong.
    result = xy4coord_apply_displacement(0, 0, 0, 0, 0.03, 0, 0, 0, 0)
    assert "error" in result
    assert "bonds" not in result
