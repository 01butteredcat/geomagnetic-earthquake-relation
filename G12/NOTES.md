# CLAUDE.md (G12)

Group G12 of the 13-event geomagnetic precursor dataset — see `../CLAUDE.md` for the shared IAGA-2002 data format and station table, and `../docs/13_groups_fetch_ranges.md` for the full event/range rationale.

## Event(s)

Anchor 2025-06-11 ML6.42 offshore Hualien (since 2026-08-20; the data window was fetched around the earlier anchor, 2025-08-27 ML6.05 offshore Yilan). See `../scripts/events.py` for the full event list (incl. the M5 backfill).

## Contents (verified against disk 2026-08-04)

- 1,392 `.sec` files, 2025-05-26 ~ 2025-09-18.
- 12 stations, uniform 116 files each: `cnu, csg, hcn, kma, lnu, lyn, mtu, ncg, twu, xcg, yhg, zbn`.
- **`ttn` (Beinan) is entirely absent** from this group (12 stations instead of the usual 13), consistent with it going offline in late 2024 — see G11's notes and root CLAUDE.md for the full timeline.

## Known limitation: too little data before the anchor (accepted 2026-09-28)

The folder starts 2025-05-26 -- ~93 days before the old 2025-08-27 anchor, but only **16 days before the current 2025-06-11 anchor**. The daily index's 28-day trailing baseline (`compute_indices.TRAILING_WINDOW_DAYS`) therefore has fewer than 5 clean days (`MIN_CLEAN_POINTS`) for the days just before the anchor, so their H/Z index is NaN, and `verify_pipeline.py`'s `baseline_window_excludes_storms` check fails for G12. This is a data limitation, not a bug: no window length fixes it; only fetching data from ~2025-04 onward would. Accepted as-is -- G12's daily-scale pre-anchor result should be read as "no baseline", not as "no anomaly". The ULF Pc3 pre-event rank test cannot run for G12 either (1 non-storm day in its 30-day pre-event window; `rank_p` is None), so G12 does not enter the primary Pc3 test.

