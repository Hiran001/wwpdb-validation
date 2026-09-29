"""The README's code examples must actually run.

A usage guide with broken examples is worse than none, because it costs the
reader time before they discover the problem. These tests exercise every
example in the README against synthetic fixtures, so they run offline.
"""
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from wwpdb_validation import Report

pd = pytest.importorskip("pandas")


def _fixture() -> Report:
    rows = []
    for i in range(1, 61):
        rama = "OUTLIER" if i % 17 == 0 else "Favored"
        rows.append(
            f'<ModelledSubgroup chain="A" resnum="{i}" resname="ALA" model="1"'
            f' rscc="{0.95 - i * 0.004:.3f}" rsrz="{-0.5 + i * 0.03:.2f}"'
            f' EDIAm="{0.9 - i * 0.008:.3f}" owab="{10 + i * 1.1:.1f}"'
            f' phi="-{60 + i % 40}" psi="{120 - i % 50}" rama="{rama}"'
            f' rota="Favored" avgoccu="1.00"/>')
    rows.append('<ModelledSubgroup chain="W" resnum="501" resname="HOH"'
                ' model="1" rscc="0.80" owab="30.0" avgoccu="1.00"/>')
    xml = ('<?xml version="1.0"?><wwPDB-validation-information>'
           '<Entry pdbid="1tst" PDB-resolution="1.80" PDB-R="0.170"'
           ' PDB-Rfree="0.205" clashscore="3.1"'
           ' PDB-deposition-date="2020-01-15"'
           ' attemptedValidationSteps="molprobity,eds"/>'
           + "".join(rows) + "</wwPDB-validation-information>")
    return Report(xml.encode())


def test_entry_fields_example():
    rep = _fixture()
    assert rep.entry.pdb_id == "1TST"
    assert rep.entry.resolution == 1.80
    assert rep.entry.r_free == 0.205
    assert "eds" in rep.entry.attempted_steps
    assert rep.has_density


def test_residue_walk_example():
    rep = _fixture()
    found = [r for r in rep.residues(exclude_waters=True) if r.is_rama_outlier]
    assert found, "fixture should contain Ramachandran outliers"
    for r in found:
        assert r.rscc is not None and r.owab is not None


def test_dataframe_example():
    df = _fixture().to_dataframe(exclude_waters=True)
    assert len(df) == 60
    for col in ("rscc", "rsrz", "EDIAm", "owab", "phi", "psi", "rama"):
        assert col in df.columns, f"README references column {col!r}"


def test_correctable_recipe():
    df = _fixture().to_dataframe(exclude_waters=True)
    sel = df[(df.rscc < 0.7) & (df.EDIAm > 0.6)]
    assert sel is not None


def test_ramachandran_plot_recipe():
    df = _fixture().to_dataframe(exclude_waters=True).dropna(subset=["phi", "psi"])
    assert len(df) > 0
    colours = (df.rama == "OUTLIER").map({True: "red", False: "grey"})
    assert set(colours.unique()) <= {"red", "grey"}


def test_b_stratified_recipe():
    df = _fixture().to_dataframe(exclude_waters=True).dropna(subset=["owab", "rscc"])
    df["Bdecile"] = pd.qcut(df.owab, 10, labels=False, duplicates="drop")
    df["poor"] = df.rscc < 0.7
    out = df.groupby(["Bdecile", df.rama == "OUTLIER"]).poor.mean().unstack()
    assert out.shape[0] > 1, "stratification should produce several deciles"


def test_none_is_not_zero():
    """The README warns about this; assert the library behaves that way."""
    xml = ('<?xml version="1.0"?><wwPDB-validation-information>'
           '<Entry pdbid="1tst" attemptedValidationSteps="molprobity"/>'
           '<ModelledSubgroup chain="A" resnum="1" resname="ALA" model="1"/>'
           '</wwPDB-validation-information>')
    r = Report(xml.encode()).residues()[0]
    assert r.rscc is None and r.rscc != 0
