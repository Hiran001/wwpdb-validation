"""Every test here corresponds to a way of getting validation reports wrong.

All fixtures are synthetic XML, so the suite runs offline and in CI.
"""
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from wwpdb_validation import Report


def build(entry_attrs: str, residues: str) -> bytes:
    return (f'<?xml version="1.0"?><wwPDB-validation-information>'
            f'<Entry {entry_attrs}/>{residues}'
            f'</wwPDB-validation-information>').encode()


def test_entry_without_density_step_is_flagged():
    """The single most consequential field in the whole report.

    An entry whose report records no attempted density step carries geometry
    assessment but no RSCC, RSRZ or EDIAm. Code that retrieves reports and
    analyses whatever it finds will drop these silently and then report the
    number of reports retrieved as its sample size. That overstates the sample
    and does so non-randomly.
    """
    with_eds = Report(build('pdbid="1abc" attemptedValidationSteps="molprobity,eds"', ""))
    without = Report(build('pdbid="2abc" attemptedValidationSteps="molprobity,xtriage"', ""))
    assert with_eds.has_density is True
    assert without.has_density is False


def test_alternate_conformations_are_not_counted_twice():
    """A residue with altlocs appears as several rows sharing a residue number."""
    res = ('<ModelledSubgroup chain="A" resnum="10" resname="SER" altcode="A"'
           ' model="1" avgoccu="0.6" rscc="0.90"/>'
           '<ModelledSubgroup chain="A" resnum="10" resname="SER" altcode="B"'
           ' model="1" avgoccu="0.4" rscc="0.70"/>')
    r = Report(build('pdbid="1abc"', res))
    assert len(r.residues(collapse_altloc=False)) == 2
    collapsed = r.residues(collapse_altloc=True)
    assert len(collapsed) == 1
    # the higher-occupancy conformer is the one kept
    assert collapsed[0].rscc == 0.90


def test_insertion_codes_distinguish_residues():
    """resnum alone does not identify a residue. 52 and 52A are different."""
    res = ('<ModelledSubgroup chain="H" resnum="52" resname="GLY" model="1"/>'
           '<ModelledSubgroup chain="H" resnum="52" icode="A" resname="SER" model="1"/>')
    r = Report(build('pdbid="1abc"', res))
    got = r.residues()
    assert len(got) == 2, "insertion code was ignored, so two residues merged into one"
    assert {x.icode for x in got} == {"", "A"}


def test_multiple_models_multiply_the_sample_unless_filtered():
    res = "".join(
        f'<ModelledSubgroup chain="A" resnum="{i}" resname="ALA" model="{m}"/>'
        for m in (1, 2) for i in (1, 2, 3))
    r = Report(build('pdbid="1abc"', res))
    assert len(r.residues()) == 6
    assert len(r.residues(model=1)) == 3


def test_worst_clash_takes_the_largest_overlap():
    res = ('<ModelledSubgroup chain="A" resnum="5" resname="LEU" model="1">'
           '<clash clashmag="-0.42"/><clash clashmag="-0.91"/><clash clashmag="0.10"/>'
           '</ModelledSubgroup>')
    r = Report(build('pdbid="1abc"', res))
    x = r.residues()[0]
    assert x.n_clashes == 3
    assert x.worst_clash == pytest.approx(0.91)


def test_no_clashes_gives_zero_not_an_error():
    res = '<ModelledSubgroup chain="A" resnum="5" resname="LEU" model="1"/>'
    x = Report(build('pdbid="1abc"', res)).residues()[0]
    assert x.n_clashes == 0 and x.worst_clash == 0.0


def test_missing_measures_are_none_not_zero():
    """A residue with no density data must not read as a residue with bad density."""
    res = '<ModelledSubgroup chain="A" resnum="5" resname="LEU" model="1"/>'
    x = Report(build('pdbid="1abc"', res)).residues()[0]
    assert x.rscc is None and x.rsrz is None and x.ediam is None
    assert x.owab is None


def test_waters_can_be_excluded():
    res = ('<ModelledSubgroup chain="A" resnum="1" resname="ALA" model="1"/>'
           '<ModelledSubgroup chain="W" resnum="501" resname="HOH" model="1"/>')
    r = Report(build('pdbid="1abc"', res))
    assert len(r.residues()) == 2
    assert len(r.residues(exclude_waters=True)) == 1


def test_seldom_used_fields_are_exposed():
    """EDIAm, OPIA, phi and psi are present in real reports and rarely read."""
    res = ('<ModelledSubgroup chain="A" resnum="5" resname="LEU" model="1"'
           ' EDIAm="0.31" OPIA="55.0" phi="-120.5" psi="140.2" owab="63.2"'
           ' rama="OUTLIER" rota="OUTLIER" cis_peptide="1"/>')
    x = Report(build('pdbid="1abc"', res)).residues()[0]
    assert x.ediam == 0.31 and x.opia == 55.0
    assert x.phi == -120.5 and x.psi == 140.2
    assert x.is_rama_outlier and x.is_rota_outlier
    assert x.attrs["cis_peptide"] == "1"


def test_entry_attributes_are_all_available():
    r = Report(build('pdbid="1abc" PDB-resolution="2.10" PDB-Rfree="0.241"'
                     ' clashscore="4.25" PDB-deposition-date="2020-01-15"'
                     ' attemptedValidationSteps="molprobity,eds"', ""))
    assert r.entry.pdb_id == "1ABC"
    assert r.entry.resolution == 2.10
    assert r.entry.r_free == 0.241
    assert r.entry.clashscore == 4.25
    assert r.entry.deposition_date == "2020-01-15"


def test_residue_key_is_unique_per_residue():
    res = ('<ModelledSubgroup chain="A" resnum="1" resname="ALA" model="1"/>'
           '<ModelledSubgroup chain="B" resnum="1" resname="ALA" model="1"/>')
    keys = {x.key for x in Report(build('pdbid="1abc"', res)).residues()}
    assert len(keys) == 2, "chain is part of a residue's identity"
