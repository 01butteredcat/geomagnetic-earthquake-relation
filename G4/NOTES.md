# CLAUDE.md (G4)

Group G4 of the 13-event geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../docs/13_groups_fetch_ranges.md` for the full event/range rationale.

## Event(s)

2020-12-10 M6.7, offshore Yilan.

## Contents (verified against disk 2026-08-04)

- 1,313 `.sec` files, 2020-09-08 ~ 2021-01-01.
- Per-station counts are **uneven** — this window sits during the network's code-transition/expansion period (old codes phasing out, new codes coming online mid-window), not a merge/download problem:

  | Station | Files | Range |
  |---|---|---|
  | csg, hcn, lyn, sme, ttn, twu, yhg | 114 each | 2020-09-08 ~ 2021-01-01 (full) |
  | mtu | 111 | 2020-09-11 ~ 2021-01-01 |
  | kmn | 76 | 2020-09-08 ~ 2020-11-22 |
  | ncg | 81 | 2020-10-11 ~ 2021-01-01 |
  | xcg | 83 | 2020-10-09 ~ 2021-01-01 |
  | zbn | 84 | 2020-10-10 ~ 2021-01-01 |
  | lnu | 61 | 2020-11-01 ~ 2021-01-01 |
  | kma | 19 | 2020-12-11 ~ 2021-01-01 |

  `kmn` stops 2020-11-22 while `kma` starts 2020-12-11 — a confirmed code replacement (19-day gap). `lnu`, `xcg`, `zbn` all come online mid-window as **new** station slots (not confirmed renames of any earlier retired code — see root CLAUDE.md); `ncg` already existed in G1/G2_G3 but has a partial-coverage gap here.
