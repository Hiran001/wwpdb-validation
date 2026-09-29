# Contributing

Bug reports and pull requests are welcome.

## Reporting a parsing problem

The most useful bug report names a **PDB ID** and what the library returned
against what the report actually contains. Reports are public, so anyone can
reproduce it:

```python
from wwpdb_validation import Report, fetch
rep = Report(fetch("XXXX"))
```

If a field is parsed wrongly, quoting the relevant XML attribute helps a lot.

## Pull requests

- Every behavioural change needs a test. The suite is offline and uses
  synthetic XML fixtures, so tests stay fast and deterministic.
- Each existing test corresponds to a specific way of reading these reports
  incorrectly. If you fix a new one, add a test that names it in the docstring,
  so the reason survives longer than the memory of the person who found it.
- Run `pytest tests -q` before opening the PR. CI runs Python 3.9, 3.11 and 3.13.

## Scope

This library parses validation reports and nothing else. It does not fetch
coordinates, compute validation metrics, or judge structures. Keeping the scope
narrow is deliberate: it makes the library easy to trust and easy to depend on.

Requests to add analysis, thresholds or recommendations are likely to be
declined, kindly, as belonging in a layer above this one.
