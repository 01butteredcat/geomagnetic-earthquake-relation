# CLAUDE.md (G11)

Group G11 of the 13-event geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../docs/13_groups_fetch_ranges.md` for the full event/range rationale.

## Event(s)

2025-01-21 M6.4, Chiayi Dapu. **This is the only inland, non-subduction-zone event** in the 13-group dataset — the other 12 groups cluster around the Yilan-Hualien offshore subduction/plate-boundary zone. May warrant separate/stratified statistical treatment.

## Contents (verified against disk 2026-08-04)

- 1,452 `.sec` files, 2024-10-20 ~ 2025-02-12.
- 12 stations with full coverage, 116 files each: `cnu, csg, hcn, kma, lnu, lyn, mtu, ncg, twu, xcg, yhg, zbn`.
- **`ttn` (Beinan) gap**: only 60 files, 2024-10-20 ~ 2024-12-18 — stops well short of the group's end date. File mtimes show `ttn` was fetched in the same batch as every other station, so this isn't a partial/failed re-fetch — the source data itself stops there. `ttn` is also entirely absent from G12 and G13 (see those groups' notes and root CLAUDE.md). Working hypothesis: `ttn` went offline/was decommissioned around 2024-12-19; recommend confirming against the CWA GDMS portal before relying on `ttn` for cross-group analysis involving this or later groups.
