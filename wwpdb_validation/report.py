"""Parse a wwPDB validation report into usable records.

The report is the authority on how a deposited structure was assessed, and it
carries far more than most tools read: 36 per-residue attributes, around 90
entry-level attributes, and eight kinds of outlier sub-element. This module
exposes all of it.

Four traps this handles, each of which silently corrupts a result rather than
raising:

1. **Entries with no density step.** If ``eds`` is absent from the entry's
   ``attemptedValidationSteps`` there are no real-space measures at all. Code
   that retrieves reports and analyses what it finds will silently drop those
   entries and report the number retrieved rather than the number contributing.
   ``Report.has_density`` makes the distinction explicit.

2. **Alternate conformations.** A residue with altlocs appears as several
   ``ModelledSubgroup`` rows sharing a residue number. Counting rows counts
   those residues twice. ``residues(collapse_altloc=True)`` keeps the highest
   occupancy copy.

3. **Insertion codes.** ``resnum`` alone does not identify a residue. The key
   is (model, chain, resnum, icode).

4. **Multiple models.** Reports for structures with several models repeat every
   residue per model. Filter on ``model`` or you will multiply your sample.
"""
from __future__ import annotations

import gzip
import pathlib
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any, Iterator

__all__ = ["Report", "Residue", "EntrySummary"]

# Sub-elements that hang off a residue, each a class of outlier.
_OUTLIER_TAGS = ("clash", "bond-outlier", "angle-outlier", "plane-outlier",
                 "symm-clash", "mog-bond-outlier", "mog-angle-outlier",
                 "mog-torsion-outlier", "mog-ring-outlier", "chiral-outlier")

_FLOAT_KEYS = {
    "rsr", "rsrz", "rscc", "EDIAm", "OPIA", "owab", "avgoccu", "phi", "psi",
    "mogul_bonds_rmsz", "mogul_angles_rmsz", "NatomsEDS", "RNAscore",
    "ligand_num_clashes", "ligand_num_symm_clashes", "num-H-reduce",
    "mogul_rmsz_numbonds", "mogul_rmsz_numangles", "resnum", "model",
}


def _num(v: Any) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


@dataclass
class Residue:
    """One modelled residue, with every attribute the report gives it."""
    chain: str = ""
    resnum: int | None = None
    resname: str = ""
    icode: str = ""
    altcode: str = ""
    model: int | None = None
    attrs: dict = field(default_factory=dict)
    outliers: dict = field(default_factory=dict)

    # --- the fields people actually reach for -------------------------
    @property
    def rscc(self): return _num(self.attrs.get("rscc"))
    @property
    def rsrz(self): return _num(self.attrs.get("rsrz"))
    @property
    def rsr(self): return _num(self.attrs.get("rsr"))
    @property
    def ediam(self): return _num(self.attrs.get("EDIAm"))
    @property
    def opia(self): return _num(self.attrs.get("OPIA"))
    @property
    def owab(self): return _num(self.attrs.get("owab"))
    @property
    def occupancy(self): return _num(self.attrs.get("avgoccu"))
    @property
    def phi(self): return _num(self.attrs.get("phi"))
    @property
    def psi(self): return _num(self.attrs.get("psi"))
    @property
    def rama(self) -> str: return self.attrs.get("rama", "")
    @property
    def rota(self) -> str: return self.attrs.get("rota", "")
    @property
    def is_rama_outlier(self) -> bool: return self.rama == "OUTLIER"
    @property
    def is_rota_outlier(self) -> bool: return self.rota == "OUTLIER"
    @property
    def n_clashes(self) -> int: return len(self.outliers.get("clash", []))

    @property
    def worst_clash(self) -> float:
        """Largest van der Waals overlap in Angstrom, 0.0 if none."""
        mags = [abs(_num(c.get("clashmag")) or 0.0)
                for c in self.outliers.get("clash", [])]
        return max(mags, default=0.0)

    @property
    def key(self) -> tuple:
        """Identity of this residue. resnum alone is NOT unique."""
        return (self.model, self.chain, self.resnum, self.icode)

    def as_dict(self) -> dict:
        d = {"chain": self.chain, "resnum": self.resnum, "resname": self.resname,
             "icode": self.icode, "altcode": self.altcode, "model": self.model,
             "n_clashes": self.n_clashes, "worst_clash": self.worst_clash}
        for k, v in self.attrs.items():
            d[k] = _num(v) if k in _FLOAT_KEYS else v
        for t in _OUTLIER_TAGS:
            if t in self.outliers:
                d[f"n_{t.replace('-', '_')}"] = len(self.outliers[t])
        return d


@dataclass
class EntrySummary:
    """Entry-level attributes, all of them, plus a few named conveniences."""
    attrs: dict = field(default_factory=dict)

    @property
    def pdb_id(self) -> str: return (self.attrs.get("pdbid") or "").upper()
    @property
    def resolution(self): return _num(self.attrs.get("PDB-resolution"))
    @property
    def r_work(self): return _num(self.attrs.get("PDB-R"))
    @property
    def r_free(self): return _num(self.attrs.get("PDB-Rfree"))
    @property
    def clashscore(self): return _num(self.attrs.get("clashscore"))
    @property
    def deposition_date(self) -> str: return self.attrs.get("PDB-deposition-date", "")
    @property
    def attempted_steps(self) -> list[str]:
        return [s for s in (self.attrs.get("attemptedValidationSteps") or "").split(",") if s]


class Report:
    """A parsed wwPDB validation report."""

    def __init__(self, source: str | pathlib.Path | bytes):
        data = self._read(source)
        self._root = ET.fromstring(data)
        e = self._root.find("Entry")
        self.entry = EntrySummary(dict(e.attrib) if e is not None else {})

    @staticmethod
    def _read(source) -> bytes:
        if isinstance(source, bytes):
            return gzip.decompress(source) if source[:2] == b"\x1f\x8b" else source
        p = pathlib.Path(source)
        raw = p.read_bytes()
        return gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw

    # --- the check that matters most ---------------------------------
    @property
    def has_density(self) -> bool:
        """True when an electron-density step was attempted.

        When this is False the report carries geometry assessment but NO
        real-space measures, so RSCC, RSRZ and EDIAm are absent for every
        residue. Roughly a third of X-ray entries are in this state. Filter on
        it explicitly, or your sample will be smaller than you think and not
        randomly so.
        """
        return "eds" in self.entry.attempted_steps

    def residues(self, collapse_altloc: bool = True,
                 model: int | None = None,
                 exclude_waters: bool = False) -> list[Residue]:
        """Every modelled residue in the report.

        collapse_altloc keeps one row per (model, chain, resnum, icode),
        choosing the highest occupancy, because alternate conformations appear
        as separate rows and naive counting double counts them.
        """
        out: list[Residue] = []
        for m in self._root.iter("ModelledSubgroup"):
            a = dict(m.attrib)
            if exclude_waters and (a.get("resname") or "").upper() in {"HOH", "DOD", "WAT"}:
                continue
            mdl = _num(a.get("model"))
            mdl = int(mdl) if mdl is not None else None
            if model is not None and mdl != model:
                continue
            rn = _num(a.get("resnum"))
            r = Residue(chain=a.get("chain", ""),
                        resnum=int(rn) if rn is not None else None,
                        resname=(a.get("resname") or "").upper(),
                        icode=(a.get("icode") or "").strip(),
                        altcode=(a.get("altcode") or "").strip(),
                        model=mdl, attrs=a)
            for tag in _OUTLIER_TAGS:
                found = m.findall(tag)
                if found:
                    r.outliers[tag] = [dict(x.attrib) for x in found]
            out.append(r)
        if not collapse_altloc:
            return out
        best: dict[tuple, Residue] = {}
        for r in out:
            cur = best.get(r.key)
            if cur is None or (r.occupancy or 0) > (cur.occupancy or 0):
                best[r.key] = r
        return list(best.values())

    def to_dataframe(self, **kw):
        """Residues as a pandas DataFrame. Requires pandas."""
        import pandas as pd
        return pd.DataFrame([r.as_dict() for r in self.residues(**kw)])

    def __repr__(self) -> str:
        return (f"<Report {self.entry.pdb_id} {self.entry.resolution} A, "
                f"density={'yes' if self.has_density else 'NO'}>")
