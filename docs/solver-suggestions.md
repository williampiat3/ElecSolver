# Solver suggestions

- For **small or moderately sized systems**, `scipy.sparse.linalg.spsolve` is effective.
- For **large-scale temporal problems**, consider **MUMPS** through `python-mumps`.

MUMPS is more efficient when only the second member changes during time-stepping.

!!! tip

    See `tests.test_temporal_system` for an example of using `python-mumps` to solve the resulting system efficiently.
