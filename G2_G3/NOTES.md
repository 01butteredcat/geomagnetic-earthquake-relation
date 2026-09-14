# CLAUDE.md (G2_G3)

Combined group covering event sequences G2 and G3 of the 13-event geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../geomag_precursor/docs/13_groups_fetch_ranges.md` for the full event/range rationale.

## Event(s)

- G2: 2019-04-18 M6.1, Hualien Xiulin.
- G3: 2019-08-08 M6.0, offshore Yilan.

## Contents (verified against disk 2026-08-04)

- 2,652 `.sec` files, 2019-01-15 ~ 2019-10-30.
- Per-station counts are **uneven** — this reflects real station code transitions during 2019, not a merge/download problem:

  | Station | Files | Range |
  |---|---|---|
  | csg, hcn, kmn, lyn, ttn, twu, yhg | 289 each | 2019-01-15 ~ 2019-10-30 (full) |
  | ncg | 227 | 2019-01-15 ~ 2019-08-29 |
  | yli | 249 | 2019-01-15 ~ 2019-09-20 |
  | hln | 85 | 2019-01-15 ~ 2019-04-09 |
  | slg | 45 | 2019-01-15 ~ 2019-02-28 |
  | sme | 23 | 2019-10-08 ~ 2019-10-30 |

  `hln` and `slg` both stop mid-window with no confirmed successor code in this dataset; `sme` starts late (2019-10-08) and continues through G9 before being replaced by `cnu` starting at G10 — see root CLAUDE.md's station-code table for what is and isn't a confirmed 1:1 transition.
