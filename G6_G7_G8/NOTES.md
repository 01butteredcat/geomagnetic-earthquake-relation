# NOTES.md (G6_G7_G8)

**Update 2026-09-20:** this folder now holds three separate analysis groups, **G6** (2021-10-24), **G7** (2022-01-03) and **G8** (2022-03-23 anchor, plus its 03-23b aftershock and the 2022-05-09 event), which share this raw-data folder (`events.py::Group.folder`) but have their own `data/interim/<group>/` and their own near/far station pools chosen from each epicenter. The folder was merged only because the events were fetched together for convenience. The text below predates the split and describes the folder's data as a whole.

Combined-folder description (historical): covering event sequences G6, G7, and G8

Combined group covering event sequences G6, G7, and G8 of the 13-event geomagnetic precursor dataset — see `../NOTES.md` for the shared IAGA-2002 data format and station table, and `../docs/13_groups_fetch_ranges.md` for the full event/range rationale.

## Event(s)

- G6: 2021-10-24 M6.5, Yilan City.
- G7: 2022-01-03 M6.0, offshore Yilan.
- G8: 2022-03-23 M6.6, offshore Hualien.

## Contents (verified against disk 2026-08-04)

- 3,432 `.sec` files, 2021-07-23 ~ 2022-04-14.
- 13 stations, uniform 264 files each: `csg, hcn, kma, lnu, lyn, mtu, ncg, sme, ttn, twu, xcg, yhg, zbn`. No gaps within the group.
