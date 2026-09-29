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

## Contents

- [The four traps](#the-four-traps)
- [Fields nobody reads](#fields-nobody-reads)
- [Install](#install)
- [How to use it](#how-to-use-it)
- [How to read the output](#how-to-read-the-output)
- [Recipes](#recipes)
- [Fetching](#fetching)

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

## How to use it

### A first look at an entry

```python
from wwpdb_validation import Report, fetch

rep = Report(fetch("1onh"))
print(rep)
# <Report 1ONH 1.38 A, density=yes>
```

The repr tells you the two things to check before anything else: the resolution,
and whether density data exists at all. **If `density=NO`, stop.** That entry has
geometry assessment only, and any density field you read will be `None`.

```python
rep.entry.pdb_id            # '1ONH'
rep.entry.resolution        # 1.38
rep.entry.r_work            # 0.164
rep.entry.r_free            # 0.196
rep.entry.clashscore        # 3.12
rep.entry.deposition_date   # '2003-04-24'
rep.entry.attempted_steps   # ['molprobity', 'eds', 'xtriage', ...]
rep.entry.attrs             # all ~90 attributes, unmodified
```

### Walking the residues

```python
for r in rep.residues(exclude_waters=True):
    if r.is_rama_outlier:
        print(f"{r.chain} {r.resname}{r.resnum}{r.icode} "
              f"RSCC={r.rscc} RSRZ={r.rsrz} EDIAm={r.ediam} B={r.owab}")
```

`residues()` returns one row per residue by default, collapsing alternate
conformations to the highest-occupancy copy. Pass `collapse_altloc=False` if you
want every conformer, `model=1` to pick one model, `exclude_waters=True` to drop
solvent.

### Straight to a table

```python
df = rep.to_dataframe(exclude_waters=True)
df.shape                    # (368, 38)
```

Every attribute in the report becomes a column, numeric fields already
converted. Nothing is renamed, so a column name in the DataFrame is the
attribute name in the XML and you can look it up in the wwPDB schema.

## How to read the output

The report measures a structure along three separate axes. They answer different
questions and are easy to conflate.

### Axis 1: does the model agree with the map?

| field | range | reading |
|---|---|---|
| `rscc` | 0 to 1 | Real-space correlation between the deposited model's density and the experimental map. Higher is better. Well-ordered protein at good resolution typically sits above 0.9. Below about 0.7 is conventionally poor. **Not corrected for resolution**, so a 3.2 Å structure will look worse than a 1.2 Å one for the same quality of modelling. |
| `rsr` | ~0 upward | Real-space residual. Lower is better. Same caveat about resolution. |
| `rsrz` | Z-score | `rsr` expressed relative to structures at comparable resolution. **This is the resolution-corrected one**, and `rsrz > 2` is the wwPDB outlier criterion. Prefer it to `rscc` when comparing across resolutions. |

**What poor agreement does NOT tell you.** A low `rscc` means model and map
disagree. That happens either because the density is weak or absent, so the
atoms are not determined by the experiment, or because the density is perfectly
good and the residue is modelled wrongly. **These measures cannot separate the
two**, and the difference matters: the first case cannot be fixed by rebuilding,
the second can.

### Axis 2: is there density there at all?

| field | range | reading |
|---|---|---|
| `EDIAm` | 0 to ~1 | Median EDIA across the residue's atoms. Assesses support for the **presence** of atoms rather than the fit of a model. Below about 0.4 is commonly read as unsupported. |
| `OPIA` | 0 to 100 | Percentage of the residue's atoms with adequate density support. |
| `NatomsEDS` | count | How many atoms entered the density calculation. A small number here means the rest of the residue was not assessed. |

This is the axis that speaks to the question axis 1 cannot answer. If `rscc` is
low **and** `EDIAm` is low, the atoms are probably not determined by the data. If
`rscc` is low but `EDIAm` is high, there is density and the model is likely in
the wrong place, which is the correctable case.

### Axis 3: geometry

| field | values | reading |
|---|---|---|
| `rama` | `Favored`, `Allowed`, `OUTLIER` | Backbone conformation against the Ramachandran distribution. |
| `rota` | `Favored`, `Allowed`, `OUTLIER`, `""` | Side-chain rotamer. Empty for residues without rotamers. |
| `phi`, `psi` | degrees | The actual torsions, so you can plot Ramachandran space directly. |
| `n_clashes` | count | Atom-atom overlaps involving this residue. |
| `worst_clash` | Å | Largest overlap, as a positive number. MolProbity treats 0.4 Å and above as a clash. |
| `cis_peptide` | flag | Present when the preceding peptide bond is modelled cis. |

### The covariate that confounds all of it

| field | reading |
|---|---|
| `owab` | Occupancy-weighted average B factor for the residue. Mobility and disorder. |

**Read this whenever you compare groups.** Poor model-to-density agreement rises
steeply with B: in a sample of 646 X-ray entries, the rate among residues with no
geometry flag rose from 1.9% in the third B decile to 20.5% in the highest. Any
comparison of flagged against unflagged residues that does not stratify on `owab`
will attribute to geometry what is actually mobility. In that same sample,
standardising for B approximately halved every apparent effect.

### Entry-level percentiles

`rep.entry.attrs` contains many `absolute-percentile-*` and
`relative-percentile-*` fields. These rank the entry against the archive, either
overall or against structures at similar resolution. They describe **where the
entry sits among its peers**, not whether it is correct. A structure in the 95th
percentile for clashscore is unusual, not necessarily wrong.

### Values that are absent

Any field may be `None`, and `None` is not zero. A residue with `rscc is None`
was not assessed against density; it is not a residue that fits badly. Filtering
with `df.rscc < 0.7` silently drops those, which is usually what you want, but
count them so you know how many you dropped.

## Recipes

**Count the exclusions instead of hiding them**

```python
from wwpdb_validation import fetch_many

usable, no_report, no_density = [], [], []
for pdb_id, rep in fetch_many(ids):
    if rep is None:
        no_report.append(pdb_id)
    elif not rep.has_density:
        no_density.append(pdb_id)
    else:
        usable.append(rep)

print(f"{len(usable)} usable, {len(no_density)} without density, "
      f"{len(no_report)} with no report")
```

**Residues where the model disagrees with good density, the correctable case**

```python
df = rep.to_dataframe(exclude_waters=True)
correctable = df[(df.rscc < 0.7) & (df.EDIAm > 0.6)]
```

**Ramachandran plot from the report alone, no coordinates needed**

```python
import matplotlib.pyplot as plt
df = rep.to_dataframe(exclude_waters=True).dropna(subset=["phi", "psi"])
plt.scatter(df.phi, df.psi, s=4,
            c=(df.rama == "OUTLIER").map({True: "red", False: "grey"}))
plt.xlabel("phi (degrees)"); plt.ylabel("psi (degrees)")
```

**Compare a flag class against the rest, stratified by B**

```python
import pandas as pd
df = rep.to_dataframe(exclude_waters=True).dropna(subset=["owab", "rscc"])
df["Bdecile"] = pd.qcut(df.owab, 10, labels=False, duplicates="drop")
df["poor"] = df.rscc < 0.7
print(df.groupby(["Bdecile", df.rama == "OUTLIER"]).poor.mean().unstack())
```

Comparing the two columns **within** a decile is the comparison worth making.
Comparing the overall means is the comparison that misleads.

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
