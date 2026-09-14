# CLAUDE.md (G23)

Group G23 of the 23-group geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../geomag_precursor/docs/candidate_fetch_ranges_from_GDMScatalog.md` for the candidate research that identified this date range (superseded by `../geomag_precursor/scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-20). New standalone group added 2026-08-20 alongside G21 and G22, from a user-supplied CWA GDMS regional catalog export (`GDMScatalog.json`) that surfaced M≥6 events the original per-event web search had missed.

## Event(s)

Merged group, mirroring the `G2_G3`/`G6_G7_G8` precedent — two events whose fetch windows overlapped each other:

- 2020-06-14 ML6.09, CWA, coordinates 24.2632N/122.4350E, depth 55.55km. Not the anchor.
- 2020-07-26 ML6.24 (anchor), CWA, coordinates 24.2552N/122.4215E, depth 53.59km. Higher magnitude than 06-14, so this is the group's anchor.

Neither event has a confirmed CWA place name found via web search (both too small to have dedicated news coverage indexed) — coordinates put both offshore NE Taiwan, close to each other; treat any place-name description as coordinate-derived, not an official CWA name. Combined baseline/aftermath window (2020-03-13~2020-08-17) didn't overlap G4's window (starts 2020-09-08).

## Contents (verified against disk 2026-08-20)

- 1,264 `.sec.gz` files, 2020-03-13 ~ 2020-08-17.
- 8 stations, uniform 158 files each: `csg, hcn, kmn, lyn, sme, ttn, twu, yhg`. No gaps within the group.
- **Notably smaller station pool than G21/G22 (8 vs. 11) despite being chronologically modern** — has `csg` (installed by this era, unlike G14–G16/G21/G22) and still uses the pre-transition codes `kmn`/`sme` (this window is entirely before G4's Nov–Dec 2020 `kmn→kma`/`sme→cnu` transition), but is missing `cnu, kma, lnu, mtu, ncg, pta, xcg, zbn` that the modern G4+ network otherwise has by this vintage. **Reason unconfirmed** — could be GDMS's regional-catalog batch export not including every station slot CWA had by 2020, or a real narrower deployment at this specific time; flagging as an open question rather than a settled explanation.
- **Has usable vector (X/Y/Z) stations**, unlike G21/G22 and G14–G18 — `ulf_analysis.py` ran successfully (`data/interim/G23/ulf_daily.csv`, `ulf_near_far_index.csv`, `ulf_spectrogram_{kmn,yhg}_eq_window.json`), and `H`/`Z`/`ULF_pc3` methods are all available. G23 was added to the `ULF_GROUPS` tuple (duplicated across 5 `geomag_precursor/scripts/*.py` files) 2026-08-20 alongside this group's registration.
