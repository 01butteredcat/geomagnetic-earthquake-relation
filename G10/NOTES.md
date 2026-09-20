# CLAUDE.md (G10)

Group G10 of the 13-event geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../docs/13_groups_fetch_ranges.md` for the full event/range rationale. This group is the original single-event (2024/4/3 M7.2) dataset that `geomag_precursor/` was first built around.

## Event(s)

2024-04-03 M7.2 mainshock (3 records same day) + 2024-04-23 (2 records) + 2024-05-10, all Hualien-area.

## Contents (verified against disk 2026-08-04, post-merge)

- 1,989 `.sec` files, 2024-01-01 ~ 2024-06-01.
- 13 stations, uniform 153 files each: `cnu, csg, hcn, kma, lnu, lyn, mtu, ncg, ttn, twu, xcg, yhg, zbn`. No gaps within the group.
- This is the first group to use `cnu` (replacing the retired code `sme`, last seen in G9).
- Range was extended from the original 2024-01-01~2024-04-25 fetch to the full 2024-01-01~2024-06-01 recommended in `13_groups_fetch_ranges.md`, to also cover the 04-23 aftershocks and the 05-10 M6.0.
- `geomag_precursor/scripts/common.py`'s `KNOWN_OUTAGE_WINDOWS` documents several genuine short data-outage windows within this range (network-wide brief outages, and per-station gaps for `twu`, `zbn`, `kma`, `yhg`) — distinct from the structural `88888.00`/`99999.00` sentinels, see that file before treating any dip as anomaly signal.
