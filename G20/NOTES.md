# NOTES.md (G20)

Group G20 of the 20-group geomagnetic precursor dataset — see `../NOTES.md` for the shared IAGA-2002 data format and station table, and `../docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2025-12-24 ML6.1 (CWA; USGS Mww6.0), Taitung Beinan, extremely shallow (11.9km) + 2025-12-27 M7.0 (CWA; USGS Mww6.6 — one of the largest CWA/USGS magnitude gaps documented in this dataset, 0.4 units), offshore ESE Yilan (anchor) — CWA press statements called it Taiwan's largest earthquake since the 1999 921 earthquake and the 2024-04-03 Hualien earthquake (G10); the single largest-magnitude event in this entire 20-group registry.

## Contents (verified against disk 2026-08-09)

- 1,428 `.sec` files, 2025-09-22 ~ 2026-01-18.
- 12 stations, uniform 119 files each: `cnu, csg, hcn, kma, lnu, lyn, mtu, ncg, twu, xcg, yhg, zbn`. No gaps within the group.
- `ttn` (Beinan) entirely absent, consistent with the permanent post-2024-12-19 outage documented in the root NOTES.md (same as G12/G13).
- Full 12-station vector (X/Y/Z) pool, same as G19 — `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are available for this group.
