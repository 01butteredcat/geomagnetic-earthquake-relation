"""Earthquake event / group registry for the 23-group multi-event
geomagnetic precursor pipeline (originally 13 groups; see the dated entries
below for how it grew to 23).

Precise epicenter coordinates/depth/magnitude were compiled from CWA
(Central Weather Administration) earthquake bulletins, cross-checked against
USGS where CWA's own decimal coordinates could not be located via web
search -- see each event's `coord_source`/`coord_confidence`/`note` fields.
This is real research data, not placeholders: where precision is genuinely
lower (CWA coordinates unavailable, using USGS as a substitute), that is
recorded explicitly rather than silently treated as equally precise.

Each group's `anchor` event is the single largest-magnitude event in that
group -- the one test point used for cross-group statistics (§7 of the
G1-G13 plan), to avoid treating a foreshock/mainshock/aftershock sequence's
members as pseudo-independent samples.

2026-09-23: 68 non-anchor events added to G11 (30), G12 (13), G13 (12), G19 (2),
and G20 (11) from a user-supplied CWA GDMS regional magnitude-report export
(`GDMScatalog.txt`, container-root-relative, M>=5.0, 2024-09-01~2026-07-31) --
this is a plain magnitude-report dump (date/time/lat/lon/depth/ML/nstn/.../
quality columns), a different export than `GDMScatalog.json`'s used elsewhere
in this file's history. Cross-checked programmatically against every existing
event in the five affected groups' `folder_events()` (+-6h / +-0.3 magnitude
tolerance, the same rule `fetch_earthquake_catalog.py::flag_known_events()`
uses) -- 7 catalog rows matched already-registered events (including all 4 of
this window's pre-existing anchors) and were skipped; the remaining 68 were
new. One of those, 2025-08-07 ML6.32 in G12's window, is a real M6+ event this
registry had not previously recorded at all (not flagged in any prior
candidate doc) but does not exceed G12's 2025-06-11 ML6.42 anchor, so no
anchor reassignment is triggered anywhere in this batch. 15 further catalog
rows fall in calendar gaps between these groups' fetch windows (no raw
geomagnetic data covers them) and were left out entirely -- not registered
here, and not a "candidate" list either since this note isn't itself a
candidate-groups doc; see `docs/candidate_events_gdms_2024_2026.md`.

2026-08-20 (yet later same day): 3 new standalone groups (G21, G22, G23) were
added, plus one non-anchor event each to G11, G17, and G6_G7_G8 -- all from
raw .sec data the user fetched from GDMS covering the candidate fetch ranges
proposed in `docs/candidate_fetch_ranges_from_GDMScatalog.md` (itself derived
from `GDMScatalog.json`'s Table 2). G21/G22/G23 are single-baseline-window
groups that didn't overlap any existing group's fetch range; G23 merges two
events (2020-06-14, 2020-07-26) whose windows overlapped each other, per the
usual mechanical merge rule. None of the three appended events (G11's
2025-04-08, G17's 2015-03-23, G6_G7_G8's 2022-05-09) exceed their group's
existing anchor magnitude, so no further anchor reassignment is triggered.

2026-08-20 (later same day): 11 additional non-anchor events, plus one anchor
reassignment (G12), were added from the same user-supplied CWA GDMS regional
catalog export (`GDMScatalog.json`) -- these are events the catalog contains
that this registry had never recorded at all (as opposed to the earlier same-day
update below, which only replaced coordinate/magnitude precision on events
already known). See `docs/candidate_groups_from_GDMScatalog.md` for the full
cross-check. G12's anchor moved from 2025-08-27 (ML6.05) to the newly-added
2025-06-11 (ML6.42), since the mechanical "largest magnitude in group" rule
(see below) applies to it too, mirroring the G17/G18 precedent.
**`data/interim/G12/` was regenerated under the new anchor 2026-08-21**
(confirmed via `data/interim/all_groups_run_summary.json`'s recorded
`anchor_date`/`anchor_magnitude`, which already show 2025-06-11/ML6.42) --
this is no longer stale, superseding an earlier version of this note that
said otherwise.

2026-08-20: 12 previously `coord_confidence="low"` (USGS-substitute) events
were upgraded to CWA-sourced coordinates/magnitudes using a user-supplied CWA
GDMS regional catalog export (see each event's note for per-event detail).
This changed G17's and G18's anchor events (2015-02-14 -> 2014-12-11 for
G17; 2016-02-06 -> 2016-05-31 for G18), since the anchor rule above is
applied mechanically to whichever event has the highest confirmed magnitude.
**`data/interim/G17/` and `data/interim/G18/` were regenerated under their
new anchors 2026-08-21** (confirmed via `all_groups_run_summary.json`'s
recorded `anchor_date`/`anchor_magnitude` for both groups) -- no longer
stale, superseding an earlier version of this note that said otherwise (see
geomag_precursor/README.md).
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta


@dataclass(frozen=True)
class Event:
    date: str  # Taiwan local calendar date commonly used to reference this event, YYYY-MM-DD
    time_local: str  # Taiwan local time (UTC+8), "YYYY-MM-DD HH:MM:SS"
    lat: float
    lon: float
    depth_km: float
    magnitude: float
    magnitude_type: str  # as reported by magnitude_source, e.g. "ML", "M"
    magnitude_source: str  # "CWA" (this dataset's primary convention) or other
    coord_source: str  # "CWA" or "USGS"
    coord_confidence: str  # "high" (CWA's own coordinates) or "low" (CWA coords not found, using USGS)
    anchor: bool = False
    note: str = ""

    @property
    def time_utc(self) -> str:
        dt = datetime.strptime(self.time_local, "%Y-%m-%d %H:%M:%S") - timedelta(hours=8)
        return dt.strftime("%Y-%m-%d %H:%M:%S")


@dataclass(frozen=True)
class Group:
    group_id: str
    folder: str  # subdirectory name under this project's root
    events: tuple[Event, ...]

    @property
    def anchor_event(self) -> Event:
        anchors = [e for e in self.events if e.anchor]
        assert len(anchors) == 1, f"{self.group_id} must have exactly one anchor event, got {len(anchors)}"
        return anchors[0]


GROUPS: dict[str, Group] = {
    "G1": Group("G1", "G1", (
        Event("2018-02-04", "2018-02-04 21:56:41", 24.15, 121.74, 10.6, 5.8, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="USGS reports this same event as Mww 6.1 (us1000cfn6) -- an unusually large ML/Mw gap "
                   "documented in seismological literature. CWA ML5.8 used here for consistency with the "
                   "rest of this registry, which prefers CWA throughout."),
        Event("2018-02-06", "2018-02-06 23:50:00", 24.10, 121.73, 6.3, 6.2, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Hualien. Mw 6.4 (USGS) / Mj 6.7 (JMA)."),
    )),
    # G2 and G3 share one raw-data folder (`G2_G3/`) purely because they were fetched together for
    # convenience -- they are independent events (112 days and ~45km apart), so each is its own
    # group with its own anchor. Split 2026-09-20 (previously one merged "G2_G3" group whose only
    # anchor was G2, leaving G3 out of every cross-group test). See folder_events() below for how
    # analyses that must exclude "every real event in this data" still see both.
    "G2": Group("G2", "G2_G3", (
        Event("2019-04-18", "2019-04-18 13:01:07", 24.06, 121.54, 18.8, 6.3, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Hualien Xiulin (花蓮秀林). Mw 6.1 (USGS)."),
    )),
    "G3": Group("G3", "G2_G3", (
        Event("2019-08-08", "2019-08-08 05:28:04", 24.44, 121.91, 24.2, 6.2, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Offshore Yilan. Mw 5.8 (USGS). Anchor of its own group since the 2026-09-20 split "
                   "of G2_G3 (previously the non-anchor second event of the merged group)."),
    )),
    "G4": Group("G4", "G4", (
        Event("2020-12-10", "2020-12-10 21:19:58", 24.74, 122.03, 76.8, 6.7, "M", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Offshore Yilan, deep subduction-zone event. Mww 6.1/depth 71km (USGS), M6.3/depth 86km "
                   "(JMA) -- depth estimate varies substantially by agency for this event; CWA's 76.8km "
                   "used as primary."),
    )),
    "G5": Group("G5", "G5", (
        Event("2021-02-07", "2021-02-07 01:36:03", 24.6632, 122.6068, 111.27, 6.21, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered independent "
                   "event falling within G5's existing fetch window. Magnitude checked against this "
                   "group's anchor: GDMScatalog.json's own record for the anchor is ML6.26, so this "
                   "ML6.21 event does not exceed it and no anchor reassignment is triggered (see the "
                   "anchor event's note -- CORRECTED 2026-09-14: the registry's anchor magnitude field "
                   "itself has now been updated to match, resolving a mismatch where this note claimed "
                   "6.26 in prose but the anchor was still registered as 6.2). See "
                   "docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2021-04-18", "2021-04-18 22:14:37", 23.8592, 121.48, 14.42, 6.26, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Hualien Shoufeng (花蓮壽豐). Mww 5.8 (USGS). CORRECTED 2026-09-14: magnitude/"
                   "coordinates/depth updated from 6.2/23.86N,121.48E/14.4km to 6.26/23.8592N,121.48E/"
                   "14.42km, confirmed against GDMScatalog.json's own record (2021-04-18 14:14:37.80 UTC, "
                   "quality B, 99 stations, exact origin-time match to the second) -- this resolves a "
                   "prior inconsistency where the 2021-02-07 event's note already asserted the anchor's "
                   "'true' magnitude was ML6.26, but the anchor's own registered magnitude field here "
                   "still said 6.2 (i.e. literally less than 2021-02-07's own registered 6.21), breaking "
                   "this registry's mechanical 'anchor = largest magnitude in group' rule on paper even "
                   "though the prose note claimed otherwise. Previous 23.86/121.48/14.4 values were "
                   "themselves just lower-precision roundings of the same event, not a different source. "
                   "NOTE: `data/interim/G5/` was regenerated 2026-09-19 with this corrected anchor -- "
                   "`run_pipeline.sh --group G5` was rerun (verify_pipeline.py: 7 pass/0 fail/1 "
                   "inconclusive, consistent with other groups) and confirmed the near/far station "
                   "ranking is unchanged (the coordinate shift is tiny, 5th/4th decimal degree, tens of "
                   "meters). `data/interim/all_groups_run_summary.json` was also regenerated by a "
                   "full run_all_groups.sh batch the same day and now records the corrected "
                   "ML6.26 anchor for G5."),
    )),
    # G6, G7, G8 share one raw-data folder (`G6_G7_G8/`) for fetch convenience only (71 and 79 days
    # apart, different epicenters); split into three groups 2026-09-20. G8 keeps its own aftershock
    # sequence (2022-03-23b, 2022-05-09) as non-anchor events.
    "G6": Group("G6", "G6_G7_G8", (
        Event("2021-10-24", "2021-10-24 13:11:34", 24.53, 121.78, 65.6, 6.5, "M", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Located in Nan'ao Township (南澳鄉), Yilan -- NOT Yilan City itself; 'Yilan City' in an "
                   "earlier internal doc was a location-name error. Mww 6.2 (USGS)."),
    )),
    "G7": Group("G7", "G6_G7_G8", (
        Event("2022-01-03", "2022-01-03 17:46:37", 24.0203, 122.1710, 22.35, 6.06, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations, exact origin-time match to the second) -- "
                   "supersedes the previous USGS Mww6.2 (us7000g8n3) substitute coordinates. Offshore "
                   "Yilan/Hualien."),
    )),
    "G8": Group("G8", "G6_G7_G8", (
        Event("2022-03-23", "2022-03-23 01:41:39", 23.40, 121.61, 25.7, 6.7, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Offshore Hualien. CWA revised this event after an initial rapid report of ML6.6/30.6km; "
                   "this entry uses the final/revised CWA catalog value (ML6.7, 25.7km). Mww 6.7 (USGS)."),
        Event("2022-03-23b", "2022-03-23 04:29:59", 23.4190, 121.4300, 22.57, 6.04, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- aftershock ~2h48m after the 2022-03-23 "
                   "anchor mainshock (01:41 local), same local calendar date hence the 'b' suffix (the "
                   "anchor itself keeps its plain date). A separate, smaller quality-C catalog record at "
                   "17:43:25 UTC (~107s after the anchor's own 17:41:38 UTC origin time) is likely a "
                   "duplicate/reprocessed solution of the anchor itself, not a distinct event, and is "
                   "intentionally not registered here. See docs/candidate_groups_from_GDMScatalog.md "
                   "Table 1."),
        Event("2022-05-09", "2022-05-09 14:23:03", 23.9702, 122.5257, 16.76, 6.27, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from raw .sec data the user fetched from GDMS covering the extended "
                   "fetch window proposed in docs/candidate_fetch_ranges_from_GDMScatalog.md (GDMScatalog.json "
                   "quality B, 99 stations) -- this event's baseline window overlapped G6_G7_G8's existing "
                   "fetch window, so it was registered as an extension of this group rather than a new "
                   "standalone one. ML6.27 does not exceed the 2022-03-23 anchor's ML6.7, so no anchor "
                   "reassignment is triggered."),
    )),
    "G9": Group("G9", "G9", (
        Event("2022-09-17", "2022-09-17 21:41:19", 23.08, 121.16, 8.6, 6.6, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Foreshock near Guanshan, Taitung. Correct CWA magnitude is ML6.6 -- an earlier internal "
                   "doc mistakenly cited 'M6.4', which traces back to South Korea KMA's figure, not CWA's."),
        Event("2022-09-18a", "2022-09-18 13:19:19", 23.1305, 121.1817, 12.13, 6.15, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered foreshock ~1h25m "
                   "before the 2022-09-18 anchor mainshock (14:44 local); shares the anchor's local "
                   "calendar date hence the 'a' suffix even though it precedes the anchor chronologically "
                   "(the anchor keeps its plain, un-suffixed date). See "
                   "docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2022-09-18", "2022-09-18 14:44:00", 23.14, 121.20, 7.8, 6.8, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Chishang/Guanshan mainshock. Mw 6.9 (USGS)."),
        Event("2022-09-19", "2022-09-19 10:07:45", 23.4410, 121.2995, 13.38, 6.02, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered aftershock the "
                   "day after the 2022-09-18 anchor mainshock. See "
                   "docs/candidate_groups_from_GDMScatalog.md Table 1."),
    )),
    "G10": Group("G10", "G10", (
        Event("2024-04-03", "2024-04-03 07:58:09", 23.8757, 121.5735, 19.72, 7.19, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Hualien mainshock. UPDATED 2026-09-21 from the user-supplied CWA GDMS regional catalog "
                   "export (GDMScatalog.json, 2024-04-02 23:58:09.94 UTC, quality B, 99 stations) -- "
                   "supersedes the original single-event pipeline's hardcoded ML7.2 / 23.819N,121.562E / "
                   "15.5km / 07:58:11, which were CWA's pre-revision values. CWA revised this event from "
                   "ML7.2 (depth 22.5km) to ML7.1 on 2025-02-01 after a second manual relocation; CWA's "
                   "scweb page now shows Mag 7.1, 23.88N/121.57E, 19.7km, which is the catalog's 7.19 / "
                   "23.8757N,121.5735E / 19.72km at published precision (7.19 -> 7.1 suggests truncation "
                   "rather than rounding; not confirmed). The catalog's two-decimal ML7.19 is used here, "
                   "matching how every other catalog-sourced event in this registry is recorded. CWA's own "
                   "2025-12-27 press statement (re: the G20 2025-12-27 Yilan earthquake) still quoted the "
                   "original 7.2 figure. NOTE: this changes the anchor's coordinates (~6km from the old "
                   "ones), which the near/far station ranking is keyed to -- G10's analysis outputs have "
                   "NOT been re-run against the new values yet. Sources: https://scweb.cwa.gov.tw/zh-tw/"
                   "earthquake/Parameters/2024040307580971019 ; zh.wikipedia.org \"2024年花蓮地震\"."),
        Event("2024-04-23a", "2024-04-23 02:26:52", 23.7400, 121.6127, 9.20, 6.16, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations, exact origin-time match to the second) -- "
                   "supersedes the previous USGS Mww6.1 (us6000mt0r) substitute coordinates."),
        Event("2024-04-23b", "2024-04-23 02:32:49", 23.85, 121.54, 5.5, 6.3, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Shoufeng Township. Mww 6.1 (USGS)."),
        Event("2024-04-23c", "2024-04-23 08:04:05", 23.8352, 121.5715, 11.68, 6.14, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered event on the "
                   "same local calendar date as 2024-04-23a/b, several hours later; 'c' continues that "
                   "lettering. See docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2024-04-27a", "2024-04-27 02:21:24", 24.1828, 121.6748, 35.54, 6.31, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered event, 28 "
                   "minutes before the next entry ('b'). See docs/candidate_groups_from_GDMScatalog.md "
                   "Table 1."),
        Event("2024-04-27b", "2024-04-27 02:49:29", 24.2505, 121.6865, 29.43, 6.00, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered event, 28 "
                   "minutes after the previous entry ('a'). See docs/candidate_groups_from_GDMScatalog.md "
                   "Table 1."),
        Event("2024-05-06", "2024-05-06 17:45:32", 23.7658, 121.5623, 27.71, 6.05, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered independent "
                   "event falling within G10's existing fetch window. See "
                   "docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2024-05-10", "2024-05-10 15:45:00", 24.2315, 121.8342, 7.72, 6.01, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality C, 99 stations, origin time matches to within 18s) -- "
                   "supersedes both the previous ML5.8 magnitude and the USGS Mww5.8 (us6000mxpi) "
                   "substitute coordinates; the catalog's ML6.01 is notably higher than the ML5.8 this "
                   "registry previously cited from a CWA report."),
    )),
    "G11": Group("G11", "G11", (
        Event("2025-01-21", "2025-01-21 00:17:26", 23.22, 120.55, 15.8, 6.4, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Chiayi Dapu / Tainan Nanxi (officially epicentered in Tainan's Nanxi District, adjacent "
                   "to Chiayi's Dapu Township where max intensity was felt). The only inland, "
                   "non-subduction-zone event of the original 13 groups (G1-G13) -- no longer the only "
                   "one registry-wide once G14-G20/G21-G23 were added: G14's 2010-03-04 Jiaxian event and "
                   "G18's 2016-02-06 Meinong event are also inland/collision-zone (see their notes and "
                   "parent CLAUDE.md's 'Inland, non-subduction-zone events' section). Mw 6.0 (USGS)."),
        Event("2025-04-08", "2025-04-08 23:26:35", 24.6573, 123.1025, 113.89, 6.15, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from raw .sec data the user fetched from GDMS covering the extended "
                   "fetch window proposed in docs/candidate_fetch_ranges_from_GDMScatalog.md (GDMScatalog.json "
                   "quality B, 99 stations) -- this event's baseline window overlapped G11's existing fetch "
                   "window, so it was registered as an extension of this group rather than a new standalone "
                   "one. ML6.15 does not exceed the 2025-01-21 anchor's ML6.4, so no anchor reassignment is "
                   "triggered."),
        Event("2024-10-24", "2024-10-24 23:59:32", 24.2405, 122.345, 69.61, 5.02, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-10-27", "2024-10-27 18:21:45", 24.0125, 121.6213, 31.47, 5.45, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-01", "2024-11-01 00:18:19", 23.5625, 121.539, 33.14, 5.50, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-02", "2024-11-02 14:27:30", 21.1612, 121.1728, 64.58, 5.20, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-07", "2024-11-07 01:19:39", 23.5643, 121.5512, 32.52, 5.56, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-14", "2024-11-14 10:38:16", 23.9023, 121.6213, 33.23, 5.46, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-22", "2024-11-22 20:40:17", 23.1897, 120.2015, 13.99, 5.53, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-27", "2024-11-27 23:32:24", 24.0757, 121.6687, 40.33, 5.36, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-11-30", "2024-11-30 17:07:03", 23.817, 123.4173, 83.77, 5.58, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-12-10", "2024-12-10 21:56:47", 24.8033, 122.4835, 96.24, 5.12, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-12-26", "2024-12-26 16:08:51", 23.962, 121.708, 26.19, 5.31, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-12-30", "2024-12-30 03:51:36", 23.5383, 120.6845, 15.08, 5.20, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-01-12", "2025-01-12 04:27:09", 23.7973, 121.4622, 11.34, 5.01, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-01-21b", "2025-01-21 00:26:25", 23.182, 120.5273, 12.73, 5.00, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'b' (following the G10 2024-04-23a/b/c "
                   "precedent) since this shares its calendar date with the 2025-01-21 anchor above -- "
                   "`seismometer_comparison.py`/`coseismic_step_analysis.py` key per-event output by "
                   "`(group_id, event.date)`, which would otherwise collide."),
        Event("2025-01-21c", "2025-01-21 01:42:31", 23.1662, 120.5722, 14.22, 5.13, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'c' for the same reason as "
                   "2025-01-21b's note above."),
        Event("2025-01-22", "2025-01-22 02:48:19", 21.1258, 121.1365, 50.71, 5.70, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-01-24", "2025-01-24 19:18:42", 23.1588, 120.5103, 15.7, 5.41, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-01-25a", "2025-01-25 06:01:57", 23.166, 120.5023, 15.22, 5.19, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'a' (G10 2024-04-23a/b/c precedent) since "
                   "this shares its calendar date with the 2025-01-25b event below -- see 2025-01-21b's note for "
                   "why."),
        Event("2025-01-25b", "2025-01-25 19:49:17", 23.2572, 120.4997, 9.2, 5.79, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'b' for the same reason as 2025-01-25a's "
                   "note above."),
        Event("2025-01-26a", "2025-01-26 00:10:50", 22.8205, 120.6782, 18.56, 5.43, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'a' (G10 2024-04-23a/b/c precedent) since "
                   "this shares its calendar date with the 2025-01-26b event below -- see 2025-01-21b's note for "
                   "why."),
        Event("2025-01-26b", "2025-01-26 07:38:53", 23.1613, 120.5143, 14.4, 5.84, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'b' for the same reason as 2025-01-26a's "
                   "note above."),
        Event("2025-01-27", "2025-01-27 08:18:50", 22.7688, 121.0303, 6.66, 5.34, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-01-30", "2025-01-30 10:11:54", 23.2347, 120.5832, 12.09, 5.68, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-02-08", "2025-02-08 00:46:58", 23.2773, 120.589, 11.41, 5.30, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-03-05", "2025-03-05 21:27:39", 23.1493, 120.497, 15.95, 5.25, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-03-13", "2025-03-13 13:09:39", 23.1102, 121.4093, 18.55, 5.72, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-03-25", "2025-03-25 09:40:48", 23.9265, 121.6245, 34.8, 5.23, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-04-01", "2025-04-01 19:26:02", 22.6255, 122.5623, 53.34, 5.10, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-04-09", "2025-04-09 09:53:25", 24.5877, 121.8153, 68.71, 5.82, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-04-23", "2025-04-23 08:02:54", 24.0708, 122.5285, 53.46, 5.08, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G11's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
    )),
    "G12": Group("G12", "G12", (
        Event("2025-06-11", "2025-06-11 19:00:29", 23.4257, 121.5323, 34.61, 6.42, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered independent "
                   "event falling within G12's existing fetch window. ML6.42 is higher than "
                   "2025-08-27's ML6.05, so `anchor` moves here from 2025-08-27 (see that event's note) "
                   "-- this registry's anchor rule is mechanical/magnitude-based (see module docstring), "
                   "mirroring the G17/G18 precedent. NOTE: `data/interim/G12/` was regenerated under "
                   "this (2025-06-11) anchor on 2026-08-21 (confirmed via all_groups_run_summary.json) "
                   "-- no longer stale, superseding an earlier version of this note that said otherwise. "
                   "See docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2025-08-27", "2025-08-27 21:11:00", 24.8665, 121.8983, 108.88, 6.05, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations, exact origin-time match to the second) -- "
                   "supersedes the previous USGS Mww5.3 (us7000qra8, depth 114.91km) substitute "
                   "coordinates. Offshore Yilan, deep event. ML6.05 is now confirmed lower than "
                   "2025-06-11's ML6.42, so this event is **no longer G12's anchor** (was anchor=True "
                   "before 2026-08-20 later same-day update; see 2025-06-11's note)."),
        Event("2025-05-30", "2025-05-30 08:02:24", 24.8692, 122.8572, 130.19, 5.36, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-06-12", "2025-06-12 00:01:04", 23.4043, 121.5322, 34.79, 5.11, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-06-20", "2025-06-20 04:18:39", 21.9103, 119.9202, 44.13, 5.20, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-06-24", "2025-06-24 02:00:57", 23.7833, 121.3857, 23.58, 5.06, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-07-07", "2025-07-07 20:01:28", 24.558, 121.7683, 62.05, 5.08, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-07-08", "2025-07-08 03:57:06", 23.9393, 121.4747, 19.47, 5.06, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-07-15", "2025-07-15 21:58:17", 21.6963, 120.9835, 38.51, 5.18, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-08-07", "2025-08-07 15:45:04", 24.5565, 122.9643, 108.29, 6.32, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-08-10", "2025-08-10 14:03:36", 23.7597, 121.5453, 21.93, 5.26, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-08-17", "2025-08-17 06:51:31", 21.2348, 119.9673, 76.25, 5.33, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-08-19", "2025-08-19 09:26:08", 23.9883, 121.7167, 25.77, 5.09, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-08-21", "2025-08-21 16:37:46", 23.2565, 120.5535, 15.44, 5.14, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-08-22", "2025-08-22 14:06:15", 23.1532, 120.5287, 16.92, 5.49, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G12's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
    )),
    "G13": Group("G13", "G13", (
        Event("2026-05-01", "2026-05-01 20:39:55", 24.93, 122.08, 98.3, 6.1, "M", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="CWA official report EQ115037. NE offshore Yilan, deep event. Mww 5.8 (USGS)."),
        Event("2026-02-24", "2026-02-24 12:37:15", 24.6567, 121.9048, 63.99, 5.68, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-03-03", "2026-03-03 19:52:48", 23.2078, 120.4723, 9.22, 5.09, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-03-12", "2026-03-12 20:14:13", 23.7757, 121.551, 23.39, 5.78, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-03-15", "2026-03-15 16:14:56", 24.3528, 121.999, 29.61, 5.59, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-03-20", "2026-03-20 21:32:42", 23.8347, 121.6112, 36.53, 5.20, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-04-05a", "2026-04-05 01:14:58", 24.0437, 121.6313, 25.43, 5.74, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'a' (G10 2024-04-23a/b/c precedent) since "
                   "this shares its calendar date with the 2026-04-05b event below -- see 2025-01-21b's note "
                   "(G11) for why."),
        Event("2026-04-05b", "2026-04-05 23:34:03", 24.821, 122.9833, 133.56, 5.26, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'b' for the same reason as 2026-04-05a's "
                   "note above."),
        Event("2026-05-04", "2026-05-04 00:18:55", 24.8417, 123.1363, 135.92, 5.69, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-05-09", "2026-05-09 03:26:39", 24.6535, 121.9998, 70.83, 5.00, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-05-12", "2026-05-12 14:53:33", 23.2922, 121.4407, 26.38, 5.61, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-05-13", "2026-05-13 18:43:30", 24.0278, 121.6147, 23.24, 5.10, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-05-17", "2026-05-17 08:46:14", 24.0065, 120.9998, 16.12, 5.16, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G13's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
    )),
    # --- G14-G20: added 2026-08-08. Candidate events were pre-researched (relative-bearing
    # locations + USGS-sourced approximate magnitudes only) in docs/candidate_groups_G14_G20.md;
    # precise CWA/USGS coordinates/depths/magnitudes below were researched fresh via WebSearch/
    # WebFetch against CWA bulletin pages, USGS FDSN, and zh-Wikipedia articles citing CWA, after
    # the user manually fetched G14-G20's raw .tgz data from GDMS on 2026-08-07.
    "G14": Group("G14", "G14", (
        Event("2009-07-14", "2009-07-14 02:05:01", 24.0228, 122.2193, 18.08, 6.00, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality C, 74 stations, origin time matches to within 1s) -- "
                   "supersedes the previous USGS Mwc6.3 (usp000gz7c) substitute magnitude/coordinates. "
                   "61km E of Hualien City per USGS."),
        Event("2009-10-04", "2009-10-04 01:36:06", 23.6475, 121.5790, 29.15, 6.09, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 71 stations, exact origin-time match to the second) -- "
                   "supersedes the previous USGS Mww6.1 (usp000h2b0) substitute magnitude/coordinates. "
                   "41km SSW of Hualien City per USGS."),
        Event("2009-11-05", "2009-11-05 17:32:57", 23.7890, 120.7187, 24.08, 6.15, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 75 stations) -- previously unregistered independent "
                   "event falling within G14's existing fetch window; not documented in "
                   "candidate_groups_G14_G20.md, only surfaced by the full regional-catalog cross-check. "
                   "WEB-VERIFIED 2026-09-14: this is a SECOND inland, non-subduction-zone event in G14 "
                   "(alongside 2010-03-04 Jiaxian below) -- epicenter is near Lugu Township, Nantou "
                   "County (南投縣鹿谷鄉), central Taiwan's inland collision/fold-thrust belt, ~90km+ "
                   "inland from the offshore-Hualien events elsewhere in this group (confirmed via USGS "
                   "usp000h3s7, \"1km SSW of Lugu, Taiwan\", 2009-11-05 09:32:56 UTC = 17:32 Taiwan local "
                   "-- same calendar date in both, and zh-Wikipedia's Taiwan earthquake list entry for "
                   "2009-11-05 17:32 Nantou, ML6.2/23.79N,120.72E/24.1km, both matching this record "
                   "closely). Commonly referenced simply as the 2009年11月5日南投(鹿谷)地震 -- no distinct "
                   "named designation like '921'. See docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2009-12-19", "2009-12-19 21:02:00", 23.79, 121.66, 43.8, 6.9, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="CWA ML6.9, 21.4km SE of Hualien City offshore, depth 43.8km (per zh-Wikipedia citing "
                   "the CWA bulletin). USGS reports this same event as Mwc6.4 (usp000h56g) -- a large "
                   "ML/Mw gap similar to the one already documented for G1's 2018-02-04 event. This "
                   "registry's mechanical 'largest CWA magnitude in group' anchor rule resolves what "
                   "docs/candidate_groups_G14_G20.md described as an unresolved tie between two "
                   "USGS-scale M6.4 events (this one and 2010-03-04 Jiaxian): using CWA's own ML6.9 "
                   "figure, this event is unambiguously G14's largest, not a tie. G14 also mixes "
                   "tectonic settings: this event plus 2009-07-14/2009-10-04 are Hualien-offshore "
                   "subduction-zone events, while 2010-03-04 (Jiaxian) AND 2009-11-05 (Nantou Lugu, "
                   "web-verified 2026-09-14 -- see that event's own note) are BOTH inland collision-zone "
                   "events -- comparable to G11's 'inland, non-subduction' caveat, but G14 has two such "
                   "events, not one. Near/far station ranking is keyed to this (offshore Hualien) anchor, "
                   "so proximity classification for the 2010-03-04/2009-11-05 inland events specifically "
                   "is not analytically meaningful (only the anchor event is tested per-group currently, "
                   "so low-priority unless a future analysis starts testing non-anchor events "
                   "individually). Pre-fetch data-availability confidence per candidate_groups_G14_G20.md: "
                   "unconfirmed; confirmed present in G14.tgz as of 2026-08-07."),
        Event("2010-03-04", "2010-03-04 08:18:51", 22.97, 120.71, 22.6, 6.4, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Kaohsiung Jiaxian (高雄甲仙), inland collision-zone earthquake -- NOT the G14 anchor "
                   "despite being tied on magnitude with 2009-12-19 in USGS-scale terms (both ~M6.4); "
                   "CWA's ML6.9 for the 2009-12-19 event makes that one the group's actual largest "
                   "magnitude event (see its note). This event caused real, documented damage (96 "
                   "injuries per zh-Wikipedia across Tainan/Chiayi/Kaohsiung; ~20 tilted buildings, 340 "
                   "damaged school facilities, 8 damaged heritage sites, HSR service disruption incl. one "
                   "derailment) -- more societal impact than the anchor event, but this registry's "
                   "magnitude-based anchor rule is kept mechanical/consistent across all groups rather "
                   "than substituting a damage-based judgment call. Mww 6.3 (USGS, usp000h8nd, "
                   "22.918,120.795, depth 21km)."),
    )),
    "G15": Group("G15", "G15", (
        Event("2013-03-27", "2013-03-27 10:03:19", 23.9022, 121.0527, 19.43, 6.24, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations) -- previously unregistered independent "
                   "event falling within G15's existing fetch window. See "
                   "docs/candidate_groups_from_GDMScatalog.md Table 1."),
        Event("2013-06-02", "2013-06-02 13:43:03", 23.8615, 120.9742, 14.54, 6.48, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Nantou Puli / Yuchi (南投埔里/魚池). Updated 2026-08-20 from the user-supplied CWA "
                   "GDMS regional catalog export (GDMScatalog.json, quality B, 99 stations, exact "
                   "origin-time match to the second) -- ML6.48 closely confirms the ML6.5 previously cited "
                   "from zh-Wikipedia citing the CWA bulletin, and supersedes the USGS Mww6.2 (usb000hbrt) "
                   "substitute coordinates that were used because CWA's own decimal epicenter had only "
                   "been reported as a relative bearing. Pre-fetch data-availability confidence per "
                   "candidate_groups_G14_G20.md: unconfirmed; confirmed present in G15.tgz as of "
                   "2026-08-07."),
    )),
    "G16": Group("G16", "G16", (
        Event("2013-10-31", "2013-10-31 20:02:09", 23.5662, 121.3485, 14.98, 6.42, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="46km SSW of Hualien City. Updated 2026-08-20 from the user-supplied CWA GDMS regional "
                   "catalog export (GDMScatalog.json, quality B, 99 stations, exact origin-time match to "
                   "the second) -- this resolves the previously-unresolved CWA ML6.4-vs-ML6.3 ambiguity "
                   "found via web search (ML6.42 rounds to the ML6.4 report) and supersedes the USGS "
                   "Mww6.3 (usc000ksdy) substitute coordinates used pending this direct GDMS confirmation. "
                   "Pre-fetch data-availability confidence per candidate_groups_G14_G20.md: unconfirmed; "
                   "confirmed present in G16.tgz as of 2026-08-07."),
    )),
    "G17": Group("G17", "G17", (
        Event("2014-12-11", "2014-12-11 05:03:39", 25.4492, 122.6070, 268.62, 6.70, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality C, 99 stations, exact origin-time match to the second) -- "
                   "supersedes the previous USGS Mww6.1 (usc000t5bu) substitute magnitude/coordinates. "
                   "Very deep event offshore NE Taiwan (77km NE of Jiufen per USGS) -- deep subduction-zone "
                   "event, roughly 'offshore Yilan' per candidate doc's characterization. ML6.70 makes "
                   "this G17's largest-magnitude event, so `anchor` moved here from 2015-02-14 (see that "
                   "event's note) -- this registry's anchor rule is mechanical/magnitude-based (see module "
                   "docstring), so the anchor tracks whichever event has the highest confirmed CWA "
                   "magnitude. NOTE: `data/interim/G17/` was regenerated under this (2014-12-11) anchor "
                   "on 2026-08-21 (confirmed via all_groups_run_summary.json) -- no longer stale, "
                   "superseding an earlier version of this note that said otherwise."),
        Event("2015-02-14", "2015-02-14 04:06:32", 22.6582, 121.3967, 27.78, 6.28, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="27.4km E-SE of Taitung County government, offshore. Updated 2026-08-20 from the "
                   "user-supplied CWA GDMS regional catalog export (GDMScatalog.json, quality C, 99 "
                   "stations, exact origin-time match to the second) -- supersedes the previous USGS "
                   "Mww6.2 (usb000tp6y) substitute coordinates. ML6.28 is now confirmed lower than "
                   "2014-12-11's ML6.70, so this event is **no longer G17's anchor** (was anchor=True "
                   "before 2026-08-20; see 2014-12-11's note). Pre-fetch data-availability confidence per "
                   "candidate_groups_G14_G20.md: fairly confident; confirmed present in G17.tgz as of "
                   "2026-08-07."),
        Event("2015-03-23", "2015-03-23 18:13:51", 23.7257, 121.6707, 38.40, 6.19, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-08-20 from raw .sec data the user fetched from GDMS covering the extended "
                   "fetch window proposed in docs/candidate_fetch_ranges_from_GDMScatalog.md (GDMScatalog.json "
                   "quality B, 99 stations) -- this event's baseline window overlapped G17's existing fetch "
                   "window, so it was registered as an extension of this group rather than a new standalone "
                   "one. ML6.19 does not exceed the 2014-12-11 anchor's ML6.70, so no anchor reassignment is "
                   "triggered."),
    )),
    "G18": Group("G18", "G18", (
        Event("2016-02-06", "2016-02-06 03:57:27", 22.93, 120.54, 16.7, 6.6, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Kaohsiung Meinong (高雄美濃). CWA revised this event's magnitude upward from an "
                   "initial ML6.4 rapid report to a final ML6.6 (max intensity also revised, from Meinong "
                   "district itself to Tainan Xinhua at intensity 7) -- this entry uses the final/revised "
                   "CWA catalog value, mirroring the G6_G7_G8 2022-03-23 event's revision precedent "
                   "already in this registry. Depth reported as 14.6km or 16.7km depending on source; "
                   "16.7km (NCREE technical report) used here. Caused 117 deaths (mostly the Weiguan "
                   "Jinlong Building collapse in Tainan) -- Taiwan's deadliest earthquake since the 1999 "
                   "921 earthquake at the time. Mww 6.4 (USGS, us20004y6h). Pre-fetch data-availability "
                   "confidence per candidate_groups_G14_G20.md: confident (a published paper confirms "
                   "≥10 CWA proton magnetometers operating that year); confirmed present in G18.tgz "
                   "as of 2026-08-07. Updated 2026-08-20: ML6.6 is now confirmed lower than 2016-05-31's "
                   "ML6.91 (per the user-supplied CWA GDMS regional catalog export) -- this event is **no "
                   "longer G18's anchor** (was anchor=True before 2026-08-20; see 2016-05-31's note)."),
        Event("2016-05-31", "2016-05-31 13:23:47", 25.4920, 122.6827, 256.89, 6.91, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality C, 99 stations, exact origin-time match to the second) -- "
                   "supersedes the previous USGS Mww6.4 (us20005zay) substitute magnitude/coordinates. "
                   "Very deep event offshore NE Taiwan, roughly 'offshore Yilan' per candidate doc's "
                   "characterization. ML6.91 makes this G18's largest-magnitude event, so `anchor` moved "
                   "here from 2016-02-06 (see that event's note) -- this registry's anchor rule is "
                   "mechanical/magnitude-based (see module docstring). NOTE: `data/interim/G18/` was "
                   "regenerated under this (2016-05-31) anchor on 2026-08-21 (confirmed via "
                   "all_groups_run_summary.json) -- no longer stale, superseding an earlier version of "
                   "this note that said otherwise."),
    )),
    "G19": Group("G19", "G19", (
        Event("2024-08-16", "2024-08-16 07:35:53", 23.7828, 121.7050, 19.36, 6.37, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Updated 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations, origin time matches to within 2s) -- "
                   "supersedes the previous USGS Mww6.1 (us7000n7b8, 30km SSE of Hualien City) substitute "
                   "coordinates/depth; also differs from the 9.7km depth previously cited from Focus "
                   "Taiwan's citation of a preliminary CWA report (19.36km here is the catalog's finalized "
                   "depth). Offshore SSE Hualien; occurs amid the ongoing 2024-04-03 (G10) aftershock "
                   "sequence's broader activity but is outside G10's fetch window (which ends "
                   "2024-06-01), and CWA/press treat it as a distinct notable event rather than merely an "
                   "aftershock label. Pre-fetch data-availability confidence per "
                   "candidate_groups_G14_G20.md: high (same era/network as G1-G13); confirmed present in "
                   "G19.tgz as of 2026-08-07."),
        Event("2024-09-02", "2024-09-02 16:26:25", 23.9255, 121.6402, 32.68, 5.63, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G19's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2024-09-07", "2024-09-07 13:16:50", 23.9572, 121.6922, 29.49, 5.35, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G19's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
    )),
    "G20": Group("G20", "G20", (
        Event("2025-12-24", "2025-12-24 17:47:06", 22.85, 121.15, 11.9, 6.1, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Taitung Beinan (台東卑南), 10.1km N of Taitung County government, extremely shallow. "
                   "CWA M6.1. Mww 6.0 (USGS, us7000rkjm)."),
        Event("2025-12-27", "2025-12-27 23:05:56", 24.66, 122.00, 67.7, 7.0, "M", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Offshore ESE Yilan. CWA's finalized report puts this at M7.0 (NOT M6.6 as "
                   "docs/candidate_groups_G14_G20.md's USGS-sourced candidate table characterized it). "
                   "CORRECTED 2026-09-14 (web-verified against CWA's own 2025-12-27 press statement, "
                   "https://news.ltn.com.tw/news/life/breakingnews/5291755 , quoting CWA directly: "
                   "\"921大地震規模達7.3、去年0403大地震規模達7.2，這起地震規模達7.0\"): this is NOT the "
                   "largest-magnitude event in this registry, as a prior version of this note incorrectly "
                   "claimed -- CWA itself explicitly ranks it third, behind the 1999 921 earthquake "
                   "(M7.3) and the 2024-04-03 Hualien earthquake (G10, M7.2 as originally reported and "
                   "still CWA's own press-comparison figure; see G10's note on CWA's later 2025-02-01 "
                   "revision to ML7.1). G10's 2024-04-03 mainshock remains the largest-magnitude event "
                   "in this 23-group registry either way (7.2 or 7.1, both > this event's 7.0). "
                   "Mww 6.6 (USGS, us7000rl2n) -- one of the largest CWA/USGS magnitude gaps documented "
                   "in this registry (0.4 magnitude units, using CWA's press-comparison M7.0 figure). "
                   "Depth also varies by source/report stage: 67.7km (CWA finalized report, scweb detail "
                   "page #2025156, https://scweb.cwa.gov.tw/zh-tw/earthquake/details/2025156) vs 72.8km "
                   "(an earlier CWA rapid-report figure circulated in press); 67.7km used here as the "
                   "more authoritative finalized value. Pre-fetch data-availability confidence per "
                   "candidate_groups_G14_G20.md: high (same era/network as G1-G13); confirmed present "
                   "in G20.tgz as of 2026-08-07."),
        Event("2025-09-29", "2025-09-29 05:23:50", 24.328, 123.2852, 80.82, 5.46, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-10-08a", "2025-10-08 07:52:13", 23.99, 121.5502, 10.53, 5.18, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'a' (G10 2024-04-23a/b/c precedent) since "
                   "this shares its calendar date with the 2025-10-08b event below -- see 2025-01-21b's note "
                   "(G11) for why."),
        Event("2025-10-08b", "2025-10-08 11:33:54", 23.5133, 122.0388, 37.31, 5.20, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered. `date` suffixed 'b' for the same reason as 2025-10-08a's "
                   "note above."),
        Event("2025-10-18", "2025-10-18 10:04:13", 24.2967, 122.083, 27.48, 5.44, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-10-31", "2025-10-31 08:25:07", 23.3465, 122.9392, 62.49, 5.11, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-11-10", "2025-11-10 13:00:48", 23.184, 120.6643, 9.07, 5.40, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-11-21", "2025-11-21 07:36:07", 25.232, 124.9578, 168.88, 5.90, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-12-08", "2025-12-08 19:24:57", 23.859, 121.6222, 24.52, 5.80, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2025-12-18", "2025-12-18 19:32:50", 24.1158, 121.7225, 31.0, 5.18, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-01-12", "2026-01-12 21:31:49", 24.6893, 122.0058, 70.54, 5.37, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
        Event("2026-01-18", "2026-01-18 16:47:46", 24.323, 121.8107, 7.67, 5.04, "ML", "CWA",
              coord_source="CWA", coord_confidence="high",
              note="Added 2026-09-23 from a user-supplied CWA GDMS regional magnitude-report export "
                   "(GDMScatalog.txt, covers 2024-09-01~2026-07-31) -- previously unregistered independent event "
                   "falling within G20's existing fetch window. Does not exceed this group's anchor magnitude, so no "
                   "anchor reassignment is triggered."),
    )),
    # --- G21-G23: added 2026-08-20. New standalone candidate groups proposed in
    # docs/candidate_fetch_ranges_from_GDMScatalog.md (itself derived from GDMScatalog.json's Table 2),
    # whose baseline/aftermath fetch windows did not overlap any of G1-G20's existing windows. The user
    # fetched each group's raw .sec data from GDMS on 2026-08-20, matching the proposed windows exactly.
    "G21": Group("G21", "G21", (
        Event("2010-11-21", "2010-11-21 20:31:45", 23.8525, 121.6857, 46.87, 6.14, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 76 stations). New standalone group -- this event's "
                   "baseline/aftermath window (2010-08-20~2010-12-13) did not overlap G14's window (ends "
                   "2010-03-26) or G15's (starts 2013-03-01). Raw .sec data fetched by the user from GDMS "
                   "2026-08-20, matching that proposed window exactly (1276 files)."),
    )),
    "G22": Group("G22", "G22", (
        Event("2012-06-10", "2012-06-10 05:00:18", 24.4590, 122.3068, 69.88, 6.62, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations). New standalone group -- this event's "
                   "baseline/aftermath window (2012-03-09~2012-07-02) did not overlap G21's window (ends "
                   "2010-12-13) or G15's (starts 2013-03-01). Raw .sec data fetched by the user from GDMS "
                   "2026-08-20, matching that proposed window exactly (1276 files)."),
    )),
    # G23 and G24 share one raw-data folder (`G23/`) purely because they were fetched together for
    # convenience -- confirmed no foreshock/mainshock/aftershock relationship (42 days, ~1.6km apart),
    # so each is its own group with its own anchor. Split 2026-09-22, mirroring the G2_G3/G6_G7_G8
    # precedent. Unlike G2/G3 (each has a confirmed CWA place name), NEITHER event here does -- both
    # are coordinate-derived only; do not add a place name to either note. See folder_events() below
    # for how analyses that must exclude "every real event in this data" still see both.
    "G23": Group("G23", "G23", (
        Event("2020-07-26", "2020-07-26 20:52:29", 24.2552, 122.4215, 53.59, 6.24, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations). ML6.24 is higher than 2020-06-14's ML6.09, "
                   "so this was already the merged group's anchor -- unchanged by the 2026-09-22 split. "
                   "Merged group's combined baseline/aftermath window (2020-03-13~2020-08-17) did not "
                   "overlap G4's window (starts 2020-09-08). Raw .sec data fetched by the user from GDMS "
                   "2026-08-20, matching that proposed window exactly (1264 files, shared with G24's "
                   "raw-data folder)."),
    )),
    "G24": Group("G24", "G23", (
        Event("2020-06-14", "2020-06-14 04:18:59", 24.2632, 122.4350, 55.55, 6.09, "ML", "CWA",
              coord_source="CWA", coord_confidence="high", anchor=True,
              note="Added 2026-08-20 from the user-supplied CWA GDMS regional catalog export "
                   "(GDMScatalog.json, quality B, 99 stations). Anchor of its own group since the "
                   "2026-09-22 split of the merged G23 group -- previously the non-anchor event: its own "
                   "+22-day window (through 2020-07-06) overlapped 2020-07-26's -93-day baseline (from "
                   "2020-04-24), so per the usual mechanical merge rule these two events were registered "
                   "as one group rather than two standalone ones, mirroring the G2_G3/G6_G7_G8 precedent. "
                   "Confirmed the two are not foreshock/mainshock/aftershock of the same sequence -- "
                   "1.6km/42 days apart is a coincidental window overlap only. Shares G23's raw-data "
                   "folder and combined fetch window (2020-03-13~2020-08-17, 1264 files)."),
    )),
}

ALL_GROUP_IDS: tuple[str, ...] = tuple(GROUPS.keys())


def get_group(group_id: str) -> Group:
    try:
        return GROUPS[group_id]
    except KeyError:
        raise KeyError(f"unknown group_id {group_id!r}; valid: {sorted(GROUPS)}") from None


def sibling_group_ids(group_id: str) -> tuple[str, ...]:
    """Every group (including `group_id` itself) whose raw data lives in the same
    `folder`. Groups normally have a folder of their own, so this is just
    `(group_id,)`; G2/G3 share `G2_G3/`, G6/G7/G8 share `G6_G7_G8/` (both split
    2026-09-20), and G23/G24 share `G23/` (split 2026-09-22) -- each from what used
    to be one merged group."""
    folder = get_group(group_id).folder
    return tuple(gid for gid, g in GROUPS.items() if g.folder == folder)


def folder_events(group_id: str) -> tuple[Event, ...]:
    """Every registered event in `group_id`'s raw-data folder, i.e. its own events plus
    its sibling groups'. Use this wherever the question is "which real earthquakes are in
    this data?" (excluding them from null/random reference draws, flagging catalog rows as
    already known, capping a search window at the nearest other event) rather than "which
    events belong to this group's analysis?" -- for the latter use `get_group(id).events`.
    Identical to `get_group(id).events` for every group that has a folder to itself."""
    return tuple(ev for gid in sibling_group_ids(group_id) for ev in GROUPS[gid].events)


def assign_group_for_time(group_id: str, when) -> str:
    """Which of `group_id`'s folder-sharing groups an extra (non-registered) earthquake at
    time `when` (UTC datetime/Timestamp, naive or tz-aware) should be attributed to: the one
    whose anchor is nearest in time (ties go to the earlier anchor). Sibling groups query the
    same date range, so without this rule the same catalog event would be counted once per
    sibling."""
    when = when.replace(tzinfo=None) if getattr(when, "tzinfo", None) is not None else when
    def key(gid: str):
        a = datetime.strptime(GROUPS[gid].anchor_event.time_utc, "%Y-%m-%d %H:%M:%S")
        return (abs((when - a).total_seconds()), a)
    return min(sibling_group_ids(group_id), key=key)


if __name__ == "__main__":
    for gid, g in GROUPS.items():
        a = g.anchor_event
        print(f"{gid:10s} anchor={a.date} {a.magnitude_type}{a.magnitude} "
              f"({a.lat:.3f},{a.lon:.3f}) depth={a.depth_km}km confidence={a.coord_confidence} "
              f"n_events={len(g.events)}")
