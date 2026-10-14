# H→μμ sync: TAMU tuple production and group comparison

Entry point for the sync work (scripts and docs under `sync/`).

## Compare groups (standalone: numpy, awkward, uproot, matplotlib)
- Docs: `sync/RUN_ALL_COMPARISONS.md` (commands), `sync/GROUPS_COMPARISON.md` (how it works, how to add a group).
- Script: `sync/python/compare_run3_sync_groups.py`, reads `/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/<GROUP>/`.
- Quick start: `python3 sync/python/compare_run3_sync_groups.py --list`, then `--group1 TAMU --group2 PKU --output <dir>`.

## Produce the TAMU sync tuples (needs the repository, HTCondor)
- Docs: `sync/GGFVBF_SYNC_20261015.md` (preselection, submit, merge), `sync/TAMU_REFERENCE_AUDIT_20261013.md`, `sync/LESSONS_LEARNED_PISA.md`.
- Scripts: `special_sync_skim.py` (job), `special_sync_submit.py` (Condor), `special_sync_merge.py` (merge with entries check), config `config/Run3_2024/selections_sync_ggfvbf_pisa_20261013.yaml`.
- Current state of the work: `sync/GROUPS_SYNC_STATUS_20261013.md`.
- Caveat: plain `hadd` can silently drop job files with different branch/leaf types; `special_sync_merge.py` unifies them and checks entries.
