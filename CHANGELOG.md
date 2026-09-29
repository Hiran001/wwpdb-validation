# Changelog

## 0.1.0

First release.

- Parse wwPDB validation reports from a path, a gzipped path, or bytes.
- Expose all per-residue attributes and all entry-level attributes, rather than
  a chosen subset.
- `Report.has_density` distinguishes entries where an electron-density step was
  attempted from those where it was not. In a random sample of 1000 X-ray
  entries, 342 (34.2%) had no density step and therefore no RSCC, RSRZ or
  EDIAm for any residue.
- Collapse alternate conformations by occupancy, filter by model, exclude
  waters, and identify residues by (model, chain, resnum, icode) rather than by
  residue number alone.
- Outlier sub-elements parsed for ten categories, including clashes with their
  magnitudes.
- `to_dataframe()` for pandas users.
- `fetch()` and `fetch_many()` retrieve reports from the public PDB archive with
  a local cache, one request at a time. A missing report returns None rather
  than raising, because that is a normal outcome worth counting.
- Eleven offline tests, each corresponding to a specific failure mode.
