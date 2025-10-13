import re
import string
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gmx_top4py.parsing import read_top, write_top
from gmx_top4py.constants import AA3
from gmx_top4py.utils import get_gmx_dir


## test topology parser
def test_parser_doesnt_crash_on_example(arranged_tmp_path, caplog):
    """Example file urea.top
    from <https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>
    """
    top_dict = read_top(Path("urea.top"))
    assert isinstance(top_dict, dict)


def test_doubleparse_urea(arranged_tmp_path):
    """Parsing it's own output should return the same top on urea.top"""
    top_dict = read_top(Path("urea.top"))
    top2_path = Path("pytest_urea.top")
    write_top(top_dict, top2_path)
    top2_dict = read_top(top2_path)

    top3_path = Path("pytest_urea2.top")
    write_top(top2_dict, top3_path)
    top3_dict = read_top(top3_path)

    assert top2_dict == top3_dict


@pytest.mark.skipif(
    not get_gmx_dir("gmx"),
    reason="Command 'gmx' not found, can't test gmx dir parsing.",
)
def test_ff_includes_with_gmxdir(arranged_tmp_path):
    top_dict = read_top(Path("urea.top"))
    gmx_dir = get_gmx_dir("gmx")
    assert gmx_dir is not None, "gmx dir not found"
    tip3_dict = read_top(gmx_dir / "top" / "amber99.ff" / "tip3p.itp")

    assert top_dict["atomtypes"]
    assert top_dict["bondtypes"]
    assert top_dict["moleculetype_SOL"] == tip3_dict["moleculetype_SOL"]


import os


def test_ff_includes_with_ff_in_cwd(arranged_tmp_path):
    top_dict = read_top(Path("hexala.top"))
    l = os.listdir("amber99sb-star-ildnp.ff")
    ions_dict = read_top(Path("amber99sb-star-ildnp.ff/ions.itp"))
    assert top_dict["atomtypes"]
    assert top_dict["bondtypes"]
    for ion in [
        "IB+",
        "CA",
        "CL",
        "NA",
        "MG",
        "K",
        "RB",
        "CS",
        "LI",
        "ZN",
    ]:
        assert top_dict[f"moleculetype_{ion}"] == ions_dict[f"moleculetype_{ion}"]


# test whether topology parsing is invertible
allowed_text = st.text(
    string.ascii_letters + string.digits + "!\"$%&'()+,-./:<=>?@\\^_`{|}~", min_size=1
)


@given(
    # a list of lists that correspond to a sections of a top file
    sections=st.lists(
        st.lists(
            allowed_text.filter(
                lambda x: x not in ["ffdir", "define", "moleculetype_"]
            ),
            min_size=2,
            max_size=5,
        ),
        min_size=1,
        max_size=5,
    )
)
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=1000)
def test_parser_invertible(sections, arranged_tmp_path):
    # flatten list of lists of strings to list of strings with subsection headers
    # use first element of each section as header
    for s in sections:
        header = s[0]
        header = re.sub(r"\d", "x", header)
        s[0] = f"[ {header} ]"
    ls = [l for s in sections for l in s]

    top_path = Path("topol.top")
    top2_path = Path("topol2.top")
    with open(top_path, "w") as f:
        f.write("\n".join(ls))
    top_dict = read_top(top_path)
    write_top(top_dict, top2_path)
    top2_dict = read_top(top2_path)
    assert top_dict == top2_dict


@given(ls=st.lists(allowed_text))
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=500)
def test_parser_fails_without_sections(ls, arranged_tmp_path):
    p = Path("random_content.top")
    with open(p, "w") as f:
        f.writelines(ls)
    with pytest.raises(ValueError):
        read_top(p)


## test ff file parsing
def test_parse_aminoacids_read_top():
    aminoacids_path = (
        Path(__file__).parent
        / "test_files"
        / "assets"
        / "amber99sb-star-ildnp.ff"
        / "aminoacids.rtp"
    )
    aminoacids_dict = read_top(aminoacids_path, use_gmx_dir=False)
    for aminoacid in AA3:
        assert (
            entry := aminoacids_dict.get(aminoacid)
        ), f"Aminoacid {aminoacid} not in {aminoacids_path.name}"
        ref_subsections = ["atoms", "bonds", "impropers"]
        subsections = list(entry["subsections"].keys())

        assert all(
            x in subsections for x in ref_subsections
        ), f"Aminoacid {aminoacid} does not have the subsections {ref_subsections} but {subsections}"
        assert all(len(x) == 4 for x in entry["subsections"]["atoms"]["content"])
        assert all(len(x) == 2 for x in entry["subsections"]["bonds"]["content"])
        assert all(
            len(x) in [4, 7] for x in entry["subsections"]["impropers"]["content"]
        )


def test_parse_ffbonded_read_top():
    ffbonded_path = (
        Path(__file__).parent
        / "test_files"
        / "assets"
        / "amber99sb-star-ildnp.ff"
        / "ffbonded.itp"
    )
    ffbonded_dict = read_top(ffbonded_path)
    ref_sections = ["bondtypes", "angletypes", "dihedraltypes"]
    ffbonded_sections = list(ffbonded_dict.keys())

    assert all(
        x in ffbonded_sections for x in ref_sections
    ), f"Sections {ref_sections} should be in ffbonded sections: {ffbonded_sections}"
    assert all(
        len(x) == 5 for x in ffbonded_dict["bondtypes"]["content"]
    ), "Unexpected number of elements in bondtypes"
    assert all(
        len(x) == 6 for x in ffbonded_dict["angletypes"]["content"]
    ), "Unexpected number of elements in angletypes"
    assert all(
        len(x) == 8 for x in ffbonded_dict["dihedraltypes"]["content"]
    ), "Unexpected number of elements in dihedraltypes"
