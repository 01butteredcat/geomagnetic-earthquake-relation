# NOTES.md (G19)

Group G19 of the 20-group geomagnetic precursor dataset — see `../NOTES.md` for the shared IAGA-2002 data format and station table, and `../docs/candidate_groups_G14_G20.md` for the candidate research that identified this date range (superseded by `../scripts/events.py`'s finalized event/coordinate details, researched fresh after the raw data was fetched 2026-08-07).

## Event(s)

2024-08-16 M6.3 (CWA; USGS Mww6.1 coordinates substituted, coord confidence low), SSE Hualien (anchor). Occurs amid the ongoing G10 (2024-04-03) aftershock sequence's broader activity but falls outside G10's fetch window (which ends 2024-06-01); CWA/press treat it as a distinct notable event rather than merely an aftershock label.

## Contents (verified against disk 2026-08-09)

- 1,461 `.sec` files, 2024-05-15 ~ 2024-09-07.
- 13 stations: `cnu, csg, hcn, kma, lnu, lyn, mtu, ncg, ttn, twu, xcg, yhg, zbn`, all uniform at 116 files **except `ttn` (69/116)**.
- **`ttn` (Beinan) stops after 2024-07-22** — a previously-undocumented gap, discovered while writing this file, that leaves `ttn` unavailable for this group's 2024-08-16 anchor event and its ~47-day aftermath. This is a **separate, temporary outage** from the permanent post-2024-12-19 `ttn` gap documented in the root NOTES.md: G11's data shows `ttn` back online 2024-10-20 through 2024-12-18, so `ttn` clearly came back up between this group's end (2024-09-07) and G11's start — it is not evidence the permanent outage started earlier. (This also corrects the root NOTES.md's previous "`ttn` full coverage ... through G19" framing, updated alongside this file.)
- First group in the newly-added G14–G20 batch with a full 13-station vector (X/Y/Z) pool, matching G4-onward's modern network — `ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are available for this group (subject to the `ttn` caveat above, since `ttn` is scalar-only/`F` anyway and not part of the vector pool).
