# CLAUDE.md (G14)

Group G14 of the 20-group geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../geomag_precursor/docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../geomag_precursor/scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2009-07-14 Mw6.3 + 2009-10-04 Mw6.1 (both offshore Hualien, USGS-sourced — CWA reports not located via web search) + 2009-12-19 ML6.9 (anchor, CWA, offshore Hualien) + 2010-03-04 ML6.4 Kaohsiung Jiaxian. The Jiaxian event is an inland collision-zone earthquake mixed into an otherwise offshore-Hualien group purely because its GDMS fetch window overlaps the anchor's on the calendar — see `events.py`'s note on this event for the full magnitude-tie/anchor-selection rationale.

## Contents (verified against disk 2026-08-09)

- 3,839 `.sec` files, 2009-04-12 ~ 2010-03-26.
- 11 stations, uniform 349 files each: `hcn, hln, kmn, lyn, ncg, pta, slg, ttn, twu, yhg, yli`. No gaps within the group.
- Uses the era's retired station codes (`hln, kmn, slg, yli` — no confirmed successor mapping, see root CLAUDE.md); `csg` not yet installed (absent from G14–G16). `pta` (Majja) is a historical scalar-only station slot retired before G1–G13/G19–G23's windows — also confirmed present in **G21** and **G22** (added 2026-08-20), so it's not unique to this group as earlier documented here.
- Zero usable vector (X/Y/Z) stations — scalar `F`-only, per root CLAUDE.md's "G14–G18 predate the vector-station network" note; `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped for this group.
