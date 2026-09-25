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


def test_pure_angle_displacement_uses_redundancy_solve():
    # S2a alone leaves a second-order residual the direct Sr=0 input
    # can't satisfy on its own (PRECISION_NOTES.md section 3) -- fixed
    # via XY4Coord's own Sr redundancy solve, not a tolerance change.
    result = xy4coord_apply_displacement(0, 0, 0, 0, 0.03, 0, 0, 0, 0)
    assert "error" not in result
    assert result["bonds"] == [1.09, 1.09, 1.09, 1.09]
    assert result["requested_sr"] == 0.0
    assert result["sr_used"] is not None and result["sr_used"] != 0.0
    assert "note" not in result  # requested Sr (default 0) matches sr_used


def test_sr_alone_is_geometrically_impossible_flagged_not_hidden():
    # Raising all six angles together (Sr alone) is geometrically
    # impossible, not a solver bug -- see README.md's Known issues. The
    # redundancy solve (computed from s2a/s2b/s4x/s4y/s4z alone, all
    # zero here) finds equilibrium instead; this must not be returned as
    # if it satisfied the actually-requested Sr=0.03.
    result = xy4coord_apply_displacement(0, 0, 0, 0, 0, 0, 0, 0, 0, sr=0.03)
    assert "error" not in result  # the redundancy solve finds Sr=0 instead
    assert result["requested_sr"] == 0.03
    assert result["sr_used"] == 0.0
    assert "note" in result  # explicitly flags the mismatch, doesn't hide it
