"""Relax pandas' minimum SQLAlchemy version check for the pinned Mikado container.

The Mikado biocontainer bundles pandas>=2.2 with SQLAlchemy 1.4.x. Pandas
enforces a hard floor of SQLAlchemy>=2.0.0 (see
`pandas.compat._optional.VERSIONS`) before it will treat a SQLAlchemy Engine
as a proper `Connectable`; when the check fails it silently falls back to its
legacy sqlite3-only code path and crashes on a real Engine object with
`AttributeError: 'Engine' object has no attribute 'cursor'`
(`Mikado serialise` hits this while loading TransDecoder ORFs via
`pd.read_sql_table`). SQLAlchemy 1.4's Engine API is otherwise fully
compatible with pandas' SQLAlchemy-aware code path, so lower the floor to
match what's actually installed instead of upgrading dependencies inside the
pinned container. Staged into the task workdir as `sitecustomize.py` and
picked up automatically by the Python interpreter via `PYTHONPATH`.
"""
from pandas.compat._optional import VERSIONS

VERSIONS["sqlalchemy"] = "1.4.0"
