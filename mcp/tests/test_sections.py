from freemol_mcp import binary
from freemol_mcp.tools.sections import list_freemol_sections, read_freemol_section

XY4COORD_FIXTURE = (
    binary.repo_root() / "Freemol" / "data" / "XY4Coord" / "tests" / "ch4_s1.inp"
)
CSMG_FIXTURE = binary.repo_root() / "Freemol" / "data" / "CSMG" / "tests" / "G_C1.mld"


def test_list_sections_xy4coord_fixture():
    result = list_freemol_sections(file_path=str(XY4COORD_FIXTURE))
    assert "error" not in result
    sections = result["sections"]
    assert [s["name"] for s in sections] == ["molecule", "x-xy4-symmcoord"]
    assert sections[0]["params"] == "nrec=5 format=nscxyz"
    assert sections[0]["header_line"] == 1
    assert sections[0]["body_line_count"] == 5
    assert sections[1]["body_line_count"] == 1
    assert result["citation"]["url"].startswith(
        "https://github.com/mariotti/freemol/blob/"
    )


def test_read_section_returns_exact_body_lines():
    result = read_freemol_section("molecule", file_path=str(XY4COORD_FIXTURE))
    assert "error" not in result
    assert result["params"] == "nrec=5 format=nscxyz"
    assert len(result["body_lines"]) == 5
    assert result["body_lines"][0].startswith("c  1 6")
    assert result["body_lines"][4].startswith("h  5 1")


def test_read_section_case_insensitive():
    result = read_freemol_section("MOLECULE", file_path=str(XY4COORD_FIXTURE))
    assert "error" not in result
    assert result["name"] == "molecule"


def test_read_section_not_found_lists_real_names():
    result = read_freemol_section("no-such-section", file_path=str(XY4COORD_FIXTURE))
    assert "error" in result
    assert result["sections_present"] == ["molecule", "x-xy4-symmcoord"]


def test_file_content_and_file_path_agree():
    text = XY4COORD_FIXTURE.read_text()
    by_content = list_freemol_sections(file_content=text)
    by_path = list_freemol_sections(file_path=str(XY4COORD_FIXTURE))
    assert by_content["sections"] == by_path["sections"]


def test_requires_exactly_one_of_file_content_or_file_path():
    try:
        list_freemol_sections()
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")
    try:
        list_freemol_sections(file_content="x", file_path="y")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_csmg_fixture_three_sections():
    # A shape with more sections than the XY4 fixtures, and a leading
    # comment line before the first header.
    result = list_freemol_sections(file_path=str(CSMG_FIXTURE))
    sections = result["sections"]
    assert [s["name"] for s in sections] == [
        "molecule",
        "x-csmg-cgauss",
        "x-csmg-symop",
    ]
    assert sections[0]["body_line_count"] == 4
    assert sections[1]["body_line_count"] == 4
    assert sections[2]["body_line_count"] == 1

    symop = read_freemol_section("x-csmg-symop", file_path=str(CSMG_FIXTURE))
    assert symop["body_lines"][0].strip().startswith("E")
