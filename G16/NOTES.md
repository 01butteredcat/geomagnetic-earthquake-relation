# CLAUDE.md (G16)

Group G16 of the 20-group geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2013-10-31 M6.3 (USGS Mww6.3, 46km SSW of Hualien City; anchor) — CWA-attributed figures found via web search conflict (some report ML6.4/depth15.0km near Wanrong, others ML6.3/depth19.5km) and could not be resolved to a final CWA catalog value, so USGS's figure is used as the primary/consistent one here, flagged low-confidence pending direct CWA bulletin confirmation (see `events.py`).

## Contents (verified against disk 2026-08-09)

- 1,160 `.sec` files, 2013-07-30 ~ 2013-11-22.
- 10 stations, uniform 116 files each: `hcn, hln, kmn, lyn, ncg, slg, ttn, twu, yhg, yli`. No gaps within the group.
- Same station-pool caveats as G15: `csg` not yet installed (absent from G14–G16), `pta` already retired.
- Zero usable vector (X/Y/Z) stations — scalar `F`-only; `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped for this group.
