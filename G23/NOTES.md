# NOTES.md (G23)

**Update 2026-09-22:** this folder now holds two separate analysis groups, **G23** (2020-07-26 anchor) and **G24** (2020-06-14 anchor), which share this raw-data folder (`events.py::Group.folder`) but have their own `data/interim/G23/`, `data/interim/G24/` and their own near/far station pools chosen from each epicenter (in practice nearly identical here -- the two epicenters are only ~1.6km apart, far smaller than any inter-station distance). The folder was merged only because the two events' fetch windows overlapped (mechanical rule), not because they belong to one sequence -- confirmed no foreshock/mainshock/aftershock relationship. The text below predates the split and describes the folder's data as a whole.

Combined-folder description (historical): covering event sequences G23 and G24

Group covering event sequences G23 and G24 of the (now 24-group) geomagnetic precursor dataset — see `../NOTES.md` for the shared IAGA-2002 data format and station table, and `../docs/candidate_fetch_ranges_from_GDMScatalog.md` for the candidate research that identified this date range (superseded by `../scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-20). New standalone group added 2026-08-20 alongside G21 and G22, from a user-supplied CWA GDMS regional catalog export (`GDMScatalog.json`) that surfaced M≥6 events the original per-event web search had missed.

## Event(s)

- G23: 2020-07-26 ML6.24 (anchor), CWA, coordinates 24.2552N/122.4215E, depth 53.59km. Higher magnitude than 2020-06-14, so this was already the merged group's anchor -- unchanged by the 2026-09-22 split.
- G24: 2020-06-14 ML6.09 (anchor of its own group since the 2026-09-22 split), CWA, coordinates 24.2632N/122.4350E, depth 55.55km. Previously the non-anchor event of the merged group.

Neither event has a confirmed CWA place name found via web search (both too small to have dedicated news coverage indexed) — coordinates put both offshore NE Taiwan, ~1.6km apart; treat any place-name description as coordinate-derived, not an official CWA name. Combined baseline/aftermath window (2020-03-13~2020-08-17) didn't overlap G4's window (starts 2020-09-08).

## Contents (verified against disk 2026-08-20)

- 1,264 `.sec.gz` files, 2020-03-13 ~ 2020-08-17 (combined baseline/aftermath window; shared unchanged by the group split, since raw data is not moved/duplicated).
- 8 stations, uniform 158 files each: `csg, hcn, kmn, lyn, sme, ttn, twu, yhg`. No gaps within the group.
- **Notably smaller station pool than G21/G22 (8 vs. 11) despite being chronologically modern** — has `csg` (installed by this era, unlike G14–G16/G21/G22) and still uses the pre-transition codes `kmn`/`sme` (this window is entirely before G4's Nov–Dec 2020 `kmn→kma`/`sme→cnu` transition), but is missing `cnu, kma, lnu, mtu, ncg, pta, xcg, zbn` that the modern G4+ network otherwise has by this vintage. **Reason unconfirmed** — could be GDMS's regional-catalog batch export not including every station slot CWA had by 2020, or a real narrower deployment at this specific time; flagging as an open question rather than a settled explanation.
- **Has usable vector (X/Y/Z) stations**, unlike G21/G22 and G14–G18 — `ulf_analysis.py` ran successfully (`data/interim/{G23,G24}/ulf_daily.csv`, `ulf_near_far_index.csv`, `ulf_spectrogram_{kmn,yhg}_eq_window.json`), and `H`/`Z`/`ULF_pc3` methods are all available. G23 was added to the `ULF_GROUPS` tuple (duplicated across 5 `geomag_precursor/scripts/*.py` files) 2026-08-20 alongside this group's registration; G24 was added to the same tuple 2026-09-22 when the group split.
