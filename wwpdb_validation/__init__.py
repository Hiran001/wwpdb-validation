"""Read wwPDB validation reports completely and correctly.

    from wwpdb_validation import Report, fetch

    rep = Report(fetch("1onh"))
    if rep.has_density:
        df = rep.to_dataframe(exclude_waters=True)
        print(df[df.rama == "OUTLIER"][["chain", "resnum", "rscc", "EDIAm", "owab"]])
"""
from .report import Report, Residue, EntrySummary
from .archive import fetch, fetch_many, default_cache

__version__ = "0.1.0"
__all__ = ["Report", "Residue", "EntrySummary", "fetch", "fetch_many",
           "default_cache", "__version__"]
