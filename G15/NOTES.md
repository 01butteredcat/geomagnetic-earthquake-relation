# NOTES.md (G15)

Group G15 of the 20-group geomagnetic precursor dataset — see `../NOTES.md` for the shared IAGA-2002 data format and station table, and `../docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2013-06-02 ML6.5 (CWA; USGS Mww6.2 — another large ML/Mw gap in this dataset), Nantou Puli/Yuchi (anchor). CWA's own decimal epicenter was only reported as a relative bearing and could not be located via web search; USGS coordinates are used as a substitute (see `events.py`).

## Contents (verified against disk 2026-08-09)

- 1,160 `.sec` files, 2013-03-01 ~ 2013-06-24.
- 10 stations, uniform 116 files each: `hcn, hln, kmn, lyn, ncg, slg, ttn, twu, yhg, yli`. No gaps within the group.
- `csg` not yet installed (absent from G14–G16); `pta` already retired here (it's present in G14, G21, and G22, but not in this group).
- Zero usable vector (X/Y/Z) stations — scalar `F`-only; `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped for this group.
