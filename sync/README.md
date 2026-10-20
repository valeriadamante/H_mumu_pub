# Run 3 H→μμ synchronisation tools

`standalone/` contains the comparison tools; they need only numpy, awkward, uproot and matplotlib (no analysis framework).

- `standalone/sync_compare.py`: yield tables, event-overlap tables and plots (bar charts, ratio to a reference group, overlap matrices) for all groups.
- `standalone/compare_run3_sync_groups.py`: detailed pair-by-pair comparison (distributions, migrations, per-pair SUMMARY.md).
- `standalone/HOWTO.txt`: step-by-step instructions to recreate the tables and plots from the common ntuples.

`python/` holds the scripts that produce the TAMU sync ntuples and the older helper scripts.
