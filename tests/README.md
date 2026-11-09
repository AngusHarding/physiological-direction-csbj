# Tests

Unit tests for the preparation stage (`../1_prepare_nhanes/`), one file per
subject (for example `test_cbc.py` for the blood count). Some test a script's
rules on small made-up examples; others check the prepared files, for example
the cohort size or that every participant appears once.

Run them from the top folder of this repository, after `1_prepare_nhanes/` has made its outputs:

```
python -B -m unittest discover -s tests -p "test_*.py"
```

107 tests, a few seconds. They only read files.
