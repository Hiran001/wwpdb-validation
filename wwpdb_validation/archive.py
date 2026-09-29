"""Retrieve validation reports from the PDB archive, with a local cache.

Reports are public files. This fetches them politely: one request at a time,
a small delay between requests, a descriptive user agent, and a cache so a
re-run costs nothing. Nothing is uploaded.
"""
from __future__ import annotations

import pathlib
import time
import urllib.error
import urllib.request
from typing import Iterable, Iterator

from .report import Report

__all__ = ["fetch", "fetch_many", "default_cache"]

_BASE = "https://files.rcsb.org/pub/pdb/validation_reports"
_UA = {"User-Agent": "wwpdb-validation (https://github.com/Hiran001/wwpdb-validation)"}


def default_cache() -> pathlib.Path:
    p = pathlib.Path.home() / ".cache" / "wwpdb-validation"
    p.mkdir(parents=True, exist_ok=True)
    return p


def fetch(pdb_id: str, cache: pathlib.Path | None = None,
          delay: float = 0.15, timeout: int = 60) -> pathlib.Path | None:
    """Download one report, or return the cached copy. None if unavailable.

    A missing report is a normal outcome, not an error: not every entry has
    one. Callers should count the Nones rather than assume success.
    """
    pid = pdb_id.lower().strip()
    if len(pid) != 4:
        raise ValueError(f"expected a 4-character PDB id, got {pdb_id!r}")
    cache = cache or default_cache()
    dest = cache / f"{pid}.xml.gz"
    if dest.exists() and dest.stat().st_size > 2000:
        return dest
    url = f"{_BASE}/{pid[1:3]}/{pid}/{pid}_validation.xml.gz"
    try:
        data = urllib.request.urlopen(
            urllib.request.Request(url, headers=_UA), timeout=timeout).read()
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return None
    if len(data) < 2000:
        return None
    dest.write_bytes(data)
    time.sleep(delay)
    return dest


def fetch_many(pdb_ids: Iterable[str], cache: pathlib.Path | None = None,
               require_density: bool = False,
               delay: float = 0.15) -> Iterator[tuple[str, Report | None]]:
    """Yield (pdb_id, Report or None) for each id, in order.

    With require_density=True, entries whose report records no attempted
    density step yield None. That exclusion is then visible in the caller's
    own loop rather than hidden inside a filter, which is the mistake this
    library exists to prevent.
    """
    for pid in pdb_ids:
        path = fetch(pid, cache=cache, delay=delay)
        if path is None:
            yield pid, None
            continue
        try:
            rep = Report(path)
        except Exception:
            yield pid, None
            continue
        if require_density and not rep.has_density:
            yield pid, None
            continue
        yield pid, rep
