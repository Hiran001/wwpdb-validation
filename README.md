# wwpdb-validation

[![tests](https://github.com/Hiran001/wwpdb-validation/actions/workflows/tests.yml/badge.svg)](https://github.com/Hiran001/wwpdb-validation/actions/workflows/tests.yml)
[![python](https://img.shields.io/badge/python-3.9%2B-blue)](https://pypi.org/project/wwpdb-validation/)
[![licence](https://img.shields.io/badge/licence-MIT-green)](LICENSE)

**Read wwPDB validation reports completely, and correctly.**

Every structure in the Protein Data Bank ships with a validation report. It is
the authority on how that structure was assessed, and it contains far more than
most code reads: **36 per-residue attributes, around 90 entry-level attributes,
and eight classes of outlier sub-element.**

This library exposes all of it, and handles four traps that silently corrupt
results rather than raising an error.

```python
from wwpdb_validation import Report, fetch

rep = Report(fetch("1onh"))

if rep.has_density:                       # <- check this first, see below
    df = rep.to_dataframe(exclude_waters=True)
    print(df[df.rama == "OUTLIER"][["chain", "resnum", "rscc", "EDIAm", "owab"]])
```

## The four traps

### 1. A third of X-ray entries have no density data at all

If `eds` is absent from the entry's `attemptedValidationSteps`, the report
carries geometry assessment but **no** RSCC, RSRZ or EDIAm for any residue.

In a random sample of 1000 X-ray entries we found **342 (34.2%)** in this state.

Code that retrieves reports and analyses whatever it finds will drop those
entries silently, then report the number of reports *retrieved* as its sample
size. That overstates the sample and does so non-randomly.

```python
rep.has_density        # True only when a density step was attempted
```

### 2. Alternate conformations appear as separate rows

A residue modelled with altlocs produces several `ModelledSubgroup` rows sharing
one residue number. Counting rows counts that residue more than once.

```python
rep.residues()                          # one row per residue, highest occupancy
rep.residues(collapse_altloc=False)     # every row, if you want them
```

### 3. `resnum` does not identify a residue

Residue 52 and residue 52A are different residues. The identity is
**(model, chain, resnum, icode)**, which is what `Residue.key` returns.

### 4. Multi-model entries repeat every residue

Without filtering you will multiply your sample by the number of models.

```python
rep.residues(model=1)
```

## Fields nobody reads

These are present in real reports and are rarely used:

| field | what it is |
|---|---|
| `EDIAm` | density support for the *presence* of atoms, not model fit |
| `OPIA` | percentage of atoms with adequate density support |
| `phi`, `psi` | the actual backbone torsions, so you can plot Ramachandran space from reports alone, with no coordinates |
| `owab` | occupancy-weighted B, the covariate that confounds most density comparisons |
| `NatomsEDS` | atoms that entered the density calculation |
| `cis_peptide` | cis peptide flag |
| `RNApucker`, `RNAsuite`, `RNAscore` | nucleic acid validation |

`owab` deserves particular mention. Poor model-to-density agreement rises
steeply with B, so any comparison of flagged against unflagged residues that
does not stratify on it will overstate the association substantially.

## Install

```bash
pip install wwpdb-validation
```

Or from source:

```bash
git clone https://github.com/Hiran001/wwpdb-validation
cd wwpdb-validation && pip install -e .
```

`pandas` is optional and needed only for `to_dataframe`.

## Usage

**One entry**

```python
from wwpdb_validation import Report, fetch

rep = Report(fetch("1onh"))
print(rep.entry.pdb_id, rep.entry.resolution, rep.entry.r_free)

for r in rep.residues(exclude_waters=True):
    if r.is_rama_outlier:
        print(r.chain, r.resnum, r.resname, r.rscc, r.ediam, r.owab)
```

**Many entries, counting the exclusions instead of hiding them**

```python
from wwpdb_validation import fetch_many

usable = missing = no_density = 0
for pdb_id, rep in fetch_many(["1onh", "4lzt", "3nir"], require_density=True):
    if rep is None:
        missing += 1            # inspect why, do not just skip
        continue
    usable += 1
```

**Straight to a DataFrame**

```python
df = rep.to_dataframe(exclude_waters=True, model=1)
```

## Fetching

`fetch()` downloads from the public PDB archive to `~/.cache/wwpdb-validation`,
one request at a time with a short delay and a descriptive user agent. A cached
report is never re-downloaded. Nothing is uploaded, and the library makes no
network call unless you ask it to.

A missing report returns `None` rather than raising, because not every entry has
one and that is a normal outcome you should count.

## Tests

Eleven tests, each corresponding to one of the failure modes above, all running
offline against synthetic fixtures.

```bash
pytest tests -q
```

## Licence

MIT.
