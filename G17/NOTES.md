# NOTES.md (G17)

Group G17 of the 20-group geomagnetic precursor dataset — see `../NOTES.md` for the shared IAGA-2002 data format and station table, and `../docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2014-12-11 M6.1 (USGS Mww6.1; CWA report not located via web search), very deep (256km) offshore NE Taiwan, roughly "offshore Yilan" + 2015-02-14 ML6.3 (CWA, anchor), offshore Taitung — CWA's own decimal coordinates not located, USGS's Mww6.2 coordinates used as a substitute (depths reasonably consistent between sources, see `events.py`).

## Contents (verified against disk 2026-08-09)

- 1,991 `.sec` files, 2014-09-09 ~ 2015-03-08.
- 11 stations, uniform 181 files each: `csg, hcn, hln, kmn, lyn, ncg, slg, ttn, twu, yhg, yli`. No gaps within the group.
- `csg` reappears here — installed sometime between G16's end (2013-11) and this group's start (2014-09), per root NOTES.md's station-code notes. Retired codes `hln, kmn, slg, yli` still in use.
- Zero usable vector (X/Y/Z) stations — scalar `F`-only; `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped for this group.
