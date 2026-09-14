# CLAUDE.md (G18)

Group G18 of the 20-group geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../geomag_precursor/docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../geomag_precursor/scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2016-02-06 ML6.6 (CWA, revised upward from an initial ML6.4 rapid report; anchor), Kaohsiung Meinong, inland collision-zone event — caused 117 deaths (mostly the Weiguan Jinlong Building collapse in Tainan), Taiwan's deadliest earthquake since 1999 at the time + 2016-05-31 M6.4 (USGS Mww6.4; CWA report not located), very deep (246km) offshore NE Taiwan, roughly "offshore Yilan".

## Contents (verified against disk 2026-08-09)

- 2,772 `.sec` files, 2015-11-05 ~ 2016-06-22.
- 12 stations, uniform 231 files each: `csg, hcn, hln, kmn, lyn, msi, ncg, slg, ttn, twu, yhg, yli`. No gaps within the group.
- `msi` (retired code, no confirmed successor — see root CLAUDE.md) appears here alongside `yhg`, consistent with the root doc's note that they coexisted in G14–G18's data.
- Zero usable vector (X/Y/Z) stations — scalar `F`-only; `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped for this group. This is the last chronological group in the dataset without a vector pool — `G19` onward has full vector coverage.
