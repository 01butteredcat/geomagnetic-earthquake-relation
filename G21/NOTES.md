# CLAUDE.md (G21)

Group G21 of the 23-group geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../geomag_precursor/docs/candidate_fetch_ranges_from_GDMScatalog.md` for the candidate research that identified this date range (superseded by `../geomag_precursor/scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-20). New standalone group added 2026-08-20 alongside G22 and G23, from a user-supplied CWA GDMS regional catalog export (`GDMScatalog.json`) that surfaced M≥6 events the original per-event web search had missed.

## Event(s)

2010-11-21 ML6.14 (anchor), CWA, coordinates 23.8525N/121.6857E, depth 46.87km. No confirmed CWA place name found via web search (event too small/old to have dedicated news coverage indexed) — the coordinates put it offshore, roughly east of Hualien; treat any place-name description as coordinate-derived, not an official CWA name. This event's baseline/aftermath window (2010-08-20~2010-12-13) didn't overlap G14's window (ends 2010-03-26) or G15's (starts 2013-03-01), so it was registered as a new standalone group rather than merged into either.

## Contents (verified against disk 2026-08-20)

- 1,276 `.sec.gz` files, 2010-08-20 ~ 2010-12-13.
- 11 stations, uniform 116 files each: `hcn, hln, kmn, lyn, ncg, pta, slg, ttn, twu, yhg, yli`. No gaps within the group.
- Uses the era's retired station codes (`hln, kmn, slg, yli` — no confirmed successor mapping, see root CLAUDE.md); `csg` not yet installed (same era as G14–G16). `pta` (Majja) is present here — confirmed via file header (`Station Name Majja`) — correcting the earlier (now-fixed) root CLAUDE.md claim that `pta` only appeared in G14; it's also present in G22.
- Zero usable vector (X/Y/Z) stations — scalar `F`-only, same situation as G14–G18; `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped for this group.
