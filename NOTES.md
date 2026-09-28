# CLAUDE.md — geomag_precursor

This file provides guidance to Claude Code when working inside `geomag_precursor/`. This is a **self-contained project**: both the raw `G1`...`G24` geomagnetic data folders and `seismometer/` (seismometer/accelerometer cross-check data) live directly inside this directory (moved here 2026-09-14 — `Gx` first, `seismometer/` shortly after — so this repo is fully self-contained and standalone).

## What this project is

A Python analysis pipeline testing whether Taiwan's CWA (Central Weather Administration) geomagnetic network shows statistically detectable earthquake-precursor or coseismic signals, across a 24-group/137-event multi-earthquake dataset (69 M≥6, 68 M5) (2009–2026, Taiwan M≥6.0 events). Two independent lines of analysis exist:

- **Daily-scale precursor screening** (`build_daily_features.py` → `compute_indices.py`/`ulf_analysis.py` → `cross_group_analysis.py`/`superposed_epoch_analysis.py`): looks for candidate anomalies in the days-to-weeks before each mainshock. Results: `data/interim/cross_group_summary.md`, `output/geomag_precursor_validation_report.html`.
- **Coseismic (at-origin-second) analysis** (`coseismic_step_analysis.py` → `coseismic_stacking_analysis.py` → `seismometer_comparison.py`): looks for a step/spike right at each earthquake's origin second in the raw 1Hz data, then asks whether that's a real field change or shaking-induced instrument noise. This is the newer, still-active line of work — see "Coseismic pipeline" below.

## Data Provenance

- All earthquake event parameters (time, location, magnitude) come from the CWA catalog / CWA GDMS. Never describe them as coming from web searches, news, or press releases in reports or docs.
- The event registry (`scripts/events.py`) must match the CWA catalog. When editing magnitudes or event lists, cross-check against the catalog and report any discrepancies.

## Git & Publishing

- Only `geomag_precursor/` is published to the public GitHub remote. Never push the container/parent repo (`precursor/`) to it.
- Use the repo-local git email (already configured). Do not change the global git config.
- Before committing, make sure no absolute local machine paths (e.g. `C:\Users\...`, `/home/...`) appear in code or docs. Use `PROJECT_DIR` in `scripts/common.py`, resolved relative to the repo.

## Project Facts (keep current)

- There are 24 station/event groups (G1–G24, after the G23/G24 split). Update this count whenever groups are split or merged. Pipeline-wide tasks such as storm cancellation or reruns apply to ALL groups unless stated otherwise.
- Full rerun: `bash scripts/run_all_groups.sh`, then run verify (`scripts/verify_pipeline.py`). Report which verify failures existed before the rerun and which are new.
- Tests: `.venv/bin/python3 -m pytest tests` (regression tests for past statistical fixes plus property tests that no null distribution sees the tested window/event). Run them after touching any null/baseline/test code.
- Magnitude tiers for validation: M≥6.0 / M≥5.5 / M≥5.0.
- Analyses are either "daily" or "second-scale (coseismic)". State which line a task belongs to before starting.

## The raw data (G1..G24)

69 CWA-catalog M≥6.0 earthquake records from 2009-07-14 through 2026-05-01 (49 original + 20 found missing by a reverse catalog check on 2026-09-27, 11 of them distant/deep offshore events that can never be an anchor; plus 68 M5 events backfilled 2026-09-23), grouped by proximity in time/location into **24 independent event sequences**, since several records are foreshock/mainshock/aftershock of the same sequence and treating them as independent samples would be pseudo-replication. Each group `Gx` is a data-fetch window of "~93 days before the (first) mainshock as baseline, ~22 days after the last event in the sequence." G1–G13 (2018–2026) were the original batch; G14–G20 (added 2026-08-08) extend the dataset backward to 2009 plus two more recent 2024/2025 events, crossing the lower edge of the 20–30-group threshold commonly cited in the literature for a statistically meaningful precursor test (still not sufficient on its own — see `data/interim/cross_group_summary.md`'s own caveats, e.g. no independent quiet-period control); G21–G23 (added 2026-08-20, alongside non-anchor events appended to G6_G7_G8, G11, G12, and G17's existing windows) came from a user-supplied CWA GDMS regional catalog export (`GDMScatalog.json`, not included in this repo) that surfaced further M≥6 events missing from the original event list.

Full event/date-range rationale and background: `docs/13_groups_fetch_ranges.md` (G1–G13), `docs/candidate_groups_G14_G20.md` (G14–G20's candidate research), and `docs/candidate_fetch_ranges_from_GDMScatalog.md` (G21–G23's candidate research) — all superseded by `scripts/events.py` once a group's data was actually fetched and registered. Group-specific details (actual file counts, stations, known gaps) live in each `Gx/CLAUDE.md`. Precise per-event epicenter coordinates/depth/magnitude, all from the CWA catalog / CWA GDMS (every event has `coord_source="CWA"`; USGS/JMA magnitudes appear only as reference values in some events' `note`), are in `scripts/events.py`.

### Layout

All 24 groups live in 20 flat directories of `.sec` files (no nested subfolders). **Three directories hold more than one group**, because those events were fetched together for convenience, not because they belong to one sequence: `G2_G3/` holds groups **G2** and **G3**, `G6_G7_G8/` holds **G6**, **G7** and **G8**, `G23/` holds groups **G23** and **G24**. Until 2026-09-20 each of the first two was registered as a single merged group (`G2_G3`, `G6_G7_G8`) with a single anchor, which left G3, G6 and G7 out of every cross-group test; `G23` was similarly a single merged group (its only anchor being the 2020-07-26 event) until the 2026-09-22 split that gave 2020-06-14 its own group, `G24`. All three are now separate groups (own anchor, own `data/interim/<group>/`, own near/far station pools chosen from their own epicenter) that share a raw-data `folder` (`events.py::Group.folder`). Anything that needs "every real event in this data" -- null/random reference draws, catalog `is_known_event` flagging, search-window caps -- uses `events.py::folder_events()` rather than `Group.events`; extra (unregistered) catalog events are attributed to one sibling by `assign_group_for_time()` so they are not counted once per group. Table rows below are per directory:

| Group | Date range (actual) | Event(s) |
|---|---|---|
| G1 | 2017-11-03 ~ 2018-02-28 | 2018-02-04 ML5.8 (CWA; Mww6.1 USGS) + 02-06 ML6.2 (anchor) Hualien |
| G2_G3 (groups G2, G3) | 2019-01-15 ~ 2019-10-30 | G2: 2019-04-18 ML6.3 (anchor) Hualien Xiulin; G3: 2019-08-08 ML6.2 offshore Yilan |
| G4 | 2020-09-08 ~ 2021-01-01 | 2020-12-10 M6.7 offshore Yilan |
| G5 | 2021-01-15 ~ 2021-05-10 | 2021-02-07 ML6.21 (non-anchor, added 2026-08-20) + 04-18 ML6.26 (anchor, corrected 2026-09-14 from a stale ML6.2 -- see events.py's note) Hualien Shoufeng |
| G6_G7_G8 (groups G6, G7, G8) | 2021-07-23 ~ 2022-05-31 | G6: 2021-10-24 M6.5 Nan'ao Township, Yilan; G7: 2022-01-03 ML6.06 offshore Yilan; G8: 2022-03-23 ML6.7 (CWA-revised from an initial ML6.6 rapid report, anchor) offshore Hualien + 03-23b ML6.04 (non-anchor, added 2026-08-20, aftershock ~2h48m after the anchor) + 05-09 ML6.27 (non-anchor, added 2026-08-20; ~24.0N,122.5E, coordinate-derived direction, no confirmed CWA place name) |
| G9 | 2022-06-17 ~ 2022-10-10 | 2022-09-17 ML6.6 foreshock + 09-18a ML6.15 (non-anchor, added 2026-08-20, foreshock ~1h25m before the mainshock) + 09-18 ML6.8 (anchor) Chishang/Guanshan mainshock + 09-19 ML6.02 (non-anchor, added 2026-08-20, aftershock the next day) |
| G10 | 2024-01-01 ~ 2024-06-01 | 2024-04-03 M7.2 mainshock (anchor) + 04-23a ML6.16 + 04-23b ML6.3 (Shoufeng Township) + 04-23c ML6.14 (non-anchor, added 2026-08-20) + 04-27a ML6.31 + 04-27b ML6.0 (both non-anchor, added 2026-08-20) + 05-06 ML6.05 (non-anchor, added 2026-08-20) + 05-10 ML6.01 (corrected 2026-08-20 from a stale ML5.8) Hualien -- 8 events total |
| G11 | 2024-10-20 ~ 2025-04-30 | 2025-01-21 M6.4 Chiayi Dapu (anchor) + 04-08 ML6.15 (non-anchor, added 2026-08-20; ~24.7N,123.1E, far NE offshore near the Ryukyu arc, coordinate-derived direction, no confirmed CWA place name) |
| G12 | 2025-05-26 ~ 2025-09-18 | 2025-06-11 ML6.42 offshore Hualien (anchor, added 2026-08-20, supersedes 2025-08-27 as largest-magnitude event in window) + 08-27 ML6.05 offshore Yilan |
| G13 | 2026-01-28 ~ 2026-05-23 | 2026-05-01 M6.1 NE offshore Yilan |
| G14 | 2009-04-12 ~ 2010-03-26 | 2009-07-14 ML6.0 (CWA; Mwc6.3 USGS) + 10-04 ML6.09 (CWA; Mww6.1 USGS) offshore Hualien + 11-05 ML6.15 (non-anchor, added 2026-08-20; Nantou Lugu (南投鹿谷), inland, web-verified 2026-09-14) + 12-19 ML6.9 (anchor) offshore Hualien + 2010-03-04 ML6.4 Kaohsiung Jiaxian (inland) -- G14 has **two** inland events, see note below |
| G15 | 2013-03-01 ~ 2013-06-24 | 2013-03-27 ML6.24 (non-anchor, added 2026-08-20; 23.9022N,121.0527E, coordinate-derived direction, no confirmed CWA place name) + 06-02 ML6.48 (anchor) Nantou Puli/Yuchi |
| G16 | 2013-07-30 ~ 2013-11-22 | 2013-10-31 ML6.42 (CWA, resolves a prior ML6.3-vs-6.4 ambiguity) SSW Hualien |
| G17 | 2014-09-09 ~ 2015-04-14 | 2014-12-11 ML6.7 (CWA, corrected 2026-08-20 from a stale USGS-substitute Mww6.1, anchor) offshore Yilan + 2015-02-14 ML6.28 offshore Taitung + 03-23 ML6.19 (non-anchor, added 2026-08-20; ~23.7N,121.7E, offshore Hualien, coordinate-derived direction, no confirmed CWA place name) |
| G18 | 2015-11-05 ~ 2016-06-22 | 2016-02-06 ML6.6 Kaohsiung Meinong + 2016-05-31 ML6.91 (CWA, corrected 2026-08-20 from a stale USGS-substitute Mww6.4, anchor) offshore Yilan |
| G19 | 2024-05-15 ~ 2024-09-07 | 2024-08-16 ML6.37 SSE Hualien |
| G20 | 2025-09-22 ~ 2026-01-18 | 2025-12-24 ML6.1 Taitung Beinan + 2025-12-27 M7.0 (CWA; Mww6.6 USGS) offshore ESE Yilan |
| G21 | 2010-08-20 ~ 2010-12-13 | 2010-11-21 ML6.14 (~23.9N,121.7E, offshore Hualien, coordinate-derived direction, no confirmed CWA place name) |
| G22 | 2012-03-09 ~ 2012-07-02 | 2012-06-10 ML6.62 (~24.5N,122.3E, far offshore NE Taiwan, coordinate-derived direction, no confirmed CWA place name) |
| G23 (groups G23, G24) | 2020-03-13 ~ 2020-08-17 | G23: 2020-07-26 ML6.24 (anchor), offshore NE Taiwan, coordinate-derived direction, no confirmed CWA place name; G24: 2020-06-14 ML6.09 (anchor of its own group since the 2026-09-22 split), same area, coordinate-derived direction, no confirmed CWA place name |

**Inland, non-subduction-zone events**: G11 (Chiayi Dapu, 2025), G14's 2009-11-05 (Nantou Lugu, added/web-verified 2026-08-20/21) and 2010-03-04 (Jiaxian) events, and G18's 2016-02-06 Meinong event are inland collision-zone earthquakes — 4 events across 3 groups (G11, G14, G18); the other 21 groups cluster around the Yilan-Hualien offshore subduction/plate-boundary zone (this now includes G21–G24, whose anchors are all offshore NE Taiwan, not inland). Statistically these may need separate/stratified treatment. Note G14 mixes **two** inland events (Nantou Lugu, Jiaxian) into an otherwise offshore-Hualien group purely because their GDMS fetch windows overlap on the calendar — see `scripts/events.py`'s notes on those two events for details.

### Seismometer/accelerometer data — not geomagnetic, don't confuse with `Gx`

`seismometer/` (moved here from the parent directory 2026-09-14, alongside `Gx`) holds 29 `GXX_MMDD` folders (e.g. `G10_0403`, `G9_0918`), each containing both its SAC PoleZero instrument-response files *and* its matching `GXX_MMDD_w.mseed` miniSEED waveform file together. **These are a different data type from everything else on this page**: instrument-response files + waveforms from nearby/co-located *seismometers/accelerometers* (not the CWA geomagnetic network), fetched manually to cross-check whether coseismic geomagnetic anomalies are real field changes or shaking-induced instrument noise. See "Seismometer comparison data" further below for the full data-quirks writeup. None of the "Data format"/"Station codes" sections below apply to these folders — they're not IAGA-2002 files at all.

### Data format

Each `.sec` file is an **IAGA-2002** format ASCII text file (CRLF line endings) containing one calendar day of 1-second geomagnetic field measurements from a single station.

- **Header** (12 fixed lines, pipe-terminated): source (Taiwan CWA), station name, geodetic lat/lon, elevation, reported/sensor-orientation channels (`XYZF`), sampling rate (1 second), data type (`Definitive`).
- **Column header**: `DATE TIME DOY <STA>X <STA>Y <STA>Z <STA>F`
- **Data rows**: one per second, 86400 per file (86413 lines including header), e.g.:
  ```
  2024-01-01 00:00:00.000 001     36390.46  -3023.81  26379.67  88888.00
  ```
  X/Y/Z are the magnetic field components in nT. For all stations except `ttn`, F is a placeholder (`88888.00`, IAGA-2002's "not reported" value) — real data is in X/Y/Z. **`ttn` (Beinan) is scalar-only**: F holds the real reading, X/Y/Z are always the placeholder. `99999.00` is a distinct sentinel meaning genuine data outage (as opposed to `88888.00`'s "channel not reported"). See `scripts/parser.py` and `common.py` for the reference implementation of this convention.

File naming: `<station><YYYYMMDD>dsec.sec`, e.g. `cnu20240101dsec.sec` = Chinan station, 2024-01-01. `:Zone.Identifier` sidecars are Windows browser download markers with no data value — ignore them when iterating over `.sec` files (glob `*.sec`, not all files).

### Station codes

The network's station codes changed over the years (older groups use retired codes; from ~G10 onward the current set is used):

| Code | Station   | Lat    | Lon     | Notes |
|------|-----------|--------|---------|---|
| cnu  | Chinan    | 23.957 | 120.928 | first appears G10; replaces retired code `sme` in the 13-station slot count (`sme` present G2_G3–G9 and, chronologically earlier, G23/G24; `cnu` present from G10 on, never overlapping). Also present G19–G20 |
| csg  | Chihshang | 23.111 | 121.226 | present in every group, G1–G13 and G17–G20; **absent from G14–G16, G21, G22** (2009–2013-era data, before this station slot's installation) — the slot appears to have been installed sometime between G16's end (2013-11) and G17's start (2014-09). Present in **G23/G24** (2020 data, chronologically after installation despite the high group numbers) |
| hcn  | Hengchun  | 21.940 | 120.814 | present in every group, G1–G24 |
| kma  | Kinmen    | 24.443 | 118.353 | replaces `kmn` — contiguous transition inside G4's window (`kmn` last file 2020-11-22, `kma` first file 2020-12-11). Also present G19–G20 |
| lnu  | Lanyu     | 22.037 | 121.558 | first appears G4 (2020-11-01) as a new station slot; not a confirmed rename of any earlier code (retired code `hln` stopped ~19 months earlier, in G2_G3). Also present G19–G20 |
| lyn  | Liyutan   | 24.346 | 120.780 | present in every group, G1–G24 |
| mtu  | Matsu     | 26.169 | 119.923 | first appears G4. Also present G19–G20 |
| ncg  | Neicheng  | 24.718 | 121.683 | present in every group, G1–G20 (partial in G1–G4 windows, see per-group notes) and G21–G22; **absent from G23/G24** (which have only an 8-station pool, see `G23/CLAUDE.md`) |
| pta  | Majja     | 22.703 | 120.653 | scalar-only (F channel), retired before G1–G13/G19–G23's windows; appears in **G14, G21, and G22** (2009–2012-era data) — a historical station slot not documented anywhere in G1–G13, first discovered when G14 extended the dataset back to 2009, later also confirmed present in G21/G22 (added 2026-08-20). Absent from G15–G18 and G23/G24 |
| ttn  | Beinan    | 22.818 | 121.080 | scalar-only (F channel); present G10, G14–G18, G19 (partial, through 2024-07-22 only — see gap note below), and **G21–G24** (full coverage, G23/G24 sharing the same raw files); absent G12–G13 and G20 (post-2024-12 permanent gap) |
| twu  | Wanqiu    | 23.185 | 120.529 | present in every group, G1–G24 |
| xcg  | Xincheng  | 24.038 | 121.609 | first appears G4 (2020-10-09) as a new station slot, not a rename — coexists with `sme` through G5–G9. Also present G19–G20 |
| yhg  | Yeheng    | 24.670 | 121.376 | present in every group, G1–G24 — **not** a rename of `yli`/`msi` (those coexisted with `yhg` in G1, and again in G14–G18's 2009–2016 data) |
| zbn  | Zhiben    | 22.739 | 121.064 | first appears G4. Also present G19–G20 |

Retired codes seen only in early-network-era data and never in the modern (G10+ vintage) network: `hln`, `kmn`, `msi`, `pta`, `slg`, `sme`, `yli`. This is a **chronological**, not group-number, distinction — the newly-added G14–G18 and G21–G23 carry group numbers overlapping or higher than G10–G13/G19–G20 but their raw data is chronologically older, so they legitimately contain these retired codes too (`hln`/`kmn`/`msi`/`slg`/`yli` all appear across G14–G18; `pta` in G14, G21, G22; `kmn`/`sme` also appear in **G23/G24**, whose 2020-03~08 data predates that year's Nov–Dec `kmn→kma`/`sme→cnu` transition inside G4's window — G23/G24 are "modern-era" chronologically (2020) but still pre-transition, distinct from G14–G18's genuinely early 2009–2016 vintage). Of these, only `kmn→kma` and `sme→cnu` are well-evidenced 1:1 transitions (contiguous or exact slot replacement); `hln`, `msi`, `pta`, `slg`, `yli` retire without a confirmed successor code in this dataset — don't assume they map onto `lnu`/`xcg`/`yhg`. Uneven per-station file counts within a group's date range (e.g. in G1–G4, G14–G18) generally reflect this real historical station churn, not merge or download damage.

### Known data issue: ttn (Beinan) gap from late 2024 onward

- **G10** (2024-01-01~06-01): `ttn` present, full coverage (153/153 files).
- **G11** (2024-10-20~2025-02-12): `ttn` present only through **2024-12-18** (60/116 files) — every other station in the group has full coverage through 2025-02-12, and file mtimes show all stations were fetched in the same batch, so this isn't a slow/retried fetch, `ttn` data simply stops there.
- **G12** and **G13** (2025-05-26 onward): `ttn` is **absent entirely** (12 stations instead of 13).

Working hypothesis: the `ttn` station went offline/was decommissioned around **2024-12-19**, rather than this being a fetch-script bug — recommend confirming against the CWA GDMS portal (`gdmsn.cwb.gov.tw`, see `docs/13_groups_fetch_ranges.md`) before relying on this in the cross-group statistical analysis.

`ttn` is present with full coverage in the newly-added historical groups **G14** (2009–2010) through **G18** (2015–2016) and, consistent with the same permanent gap, **absent from G20** (2025-09~2026-01). **G19** (2024-05-15~2024-09-07) is a separate case: `ttn` stops after **2024-07-22** (69/116 files), leaving it unavailable for G19's 2024-08-16 anchor event and its ~47-day aftermath. This is a *distinct, temporary* outage, not an early start of the permanent one — G11's data shows `ttn` back online 2024-10-20 through 2024-12-18, i.e. after G19's window and before the permanent gap. See `G19/CLAUDE.md` for details.

### Known data issue: G14–G18 predate the vector-station network

Empirically checking each new group's station pool (`common.py`'s XYZ/F pool discovery) found **G14, G15, G16, G17, and G18 all have zero usable vector (X/Y/Z) stations** — every station in these groups' `.sec` files reports only the scalar `F` channel as real data (`Reported`/`Sensor Orientation` header fields say `F`, not `XYZF`), the same situation as `G1`/`G2_G3`. This isn't a parsing bug: it reflects that Taiwan's CWA geomagnetic network was largely proton-magnetometer-based (scalar-only) before a vector-magnetometer upgrade that appears to have been substantially complete by the time of **G19** (2024) and **G20** (2025-2026), both of which have full 12-13-station vector pools like G4 onward. Concretely:

- **G14–G18** (2009–2016 data): scalar-only, `F`-method analysis only (`ULF_pc3`/`ULF_pc4`/`H`/`Z` methods are skipped — see `common.py`'s `MIN_STATIONS_FOR_METHOD` gating).
- **G19–G20** (2024-2025 data): full vector pools, same analysis methods available as G4–G13.

This means the dataset's "vector-sufficient" tier (usable for the professor-suggested ULF/ Pc3-Pc4 polarization method) is **G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G19, G20, G23, G24** (14 groups; G6/G7/G8 counted as one until the 2026-09-20 split, G23/G24 counted as one until the 2026-09-22 split) — see `scripts/*.py`'s `ULF_GROUPS` tuple (duplicated across 5 files; keep them in sync if this set changes again). `G23`/`G24` (2020 data) joined this tier 2026-08-20 despite their high group numbers, since the raw data is chronologically modern (post-vector-upgrade) even though it was fetched after G14–G20; `G21`/`G22` (2010/2012 data, also added 2026-08-20) remain scalar-only like G14–G18.

### Working with this data

- Files are large (~6 MB/day/station); when writing a parser, stream/parse a file at a time rather than loading a whole group into memory.
- Skip the 12-line header (or detect the `DATE       TIME` column line) before parsing data rows.
- No per-file checksums or manifest exist; verify row counts (should be 86400 data rows/file) as a basic integrity check.
- `build_daily_features.py` ingests via a **non-recursive** glob (`common.list_day_refs`, `*dsec.sec*` matching both plain `.sec` and gzip-compressed `.sec.gz`) against a single flat directory — this is why all `Gx` groups are kept flat rather than nested.

### `.tgz` batch downloads (GDMS batch-download feature)

A `Gx` folder can also hold one or more `.tgz` archives downloaded from GDMS's batch-download feature, directly alongside (or instead of) loose `.sec`/`.sec.gz` files — the pipeline reads `.sec` members straight out of a `.tgz` via Python's `tarfile` module, no manual `tar xzf` extraction needed. Relevant API: `common.py`'s `list_day_refs`/`resolve_day_ref`, returning `parser.DayFileRef` — a small reference that works uniformly whether the underlying data is a loose file or a `.tgz` member.

- **A `Gx` folder can have multiple `.tgz` archives.** GDMS batch exports appear to be split into several `.tgz` per request, and **archive filenames do not reliably reflect their date-range contents** (confirmed empirically: an archive named `20240403.tgz` actually contained files from 2024-05-08 through 2024-06-01). Never infer date coverage from a `.tgz`'s filename — the pipeline always enumerates the actual members inside.
- **Precedence when the same station+date exists in more than one source**: loose `.sec`/`.sec.gz` always wins over `.tgz` content (this is the expected steady state, not warned about); between multiple `.tgz` archives, the first one encountered in filename-sorted order wins, and `list_day_refs` prints a warning to stderr so overlapping batch downloads are noticeable.
- **`<name>.tgz.idx.json` sidecars**: the first `list_day_refs`/`resolve_day_ref` call against a `.tgz` builds and caches a member index next to the archive (gzip has no random access, so a full index scan of a large archive can take tens of seconds — cached indefinitely, keyed on the archive's size+mtime). Safe to delete to force a rescan, e.g. after replacing a `.tgz` with a different file of the same name.
- Because `.tgz` reads reopen/scan the whole archive per lookup unless done sequentially, `build_daily_features.py`'s batch pipeline opens each `.tgz` exactly once and extracts all wanted members in one sequential pass in the main process, handing already-extracted text to worker processes — not one `tarfile.open()` per file per worker, which measured ~600x slower on a real GDMS archive.

## Environment setup

```bash
cd geomag_precursor
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

`requirements.txt` (added 2026-08-13) is this project's first formal dependency record — before that, packages were installed ad hoc into `.venv` with nothing tracking what/why. `obspy` is the newest addition, added specifically for `seismometer_comparison.py` to read miniSEED waveforms and SAC PoleZero instrument-response files (see "Seismometer comparison data" below) — it's the project's first dependency that isn't numpy/pandas/scipy-adjacent. If you're running any script under `scripts/`, always invoke it via `.venv/bin/python3`, not a bare `python3`.

## Scripts pipeline (daily-scale, per group)

Orchestrated by `scripts/run_pipeline.sh --group <G1|...|G24> [--full-report]` (single group) / `scripts/run_all_groups.sh` (all 24 groups in 20 raw-data folders, steps 1–6 only):

1. `build_daily_features.py` — parses raw `.sec`/`.sec.gz`/`.tgz` day files → `data/interim/<group>/daily_features.csv` + `minute_series_<station>.parquet`.
2. `timezone_check.py` — confirms the `.sec` TIME column's timezone via diurnal-signal inference (dataset-wide constant `common.DATA_TIMEZONE = "UTC"`, not re-verified per group once established).
3. `fetch_space_weather.py` — Dst/Kp storm-day flags (cached-skip only if `storm_days.csv` exists AND `storm_days_summary.json` says `confidence: high`; a partial/failed earlier fetch is retried; a missing Dst month is now recorded in `dst_missing_months` and downgrades the confidence to `medium`).
4. `compute_indices.py` — daily H/Z/F near-far regression anomaly index + candidate-date flagging (MAD z-score, threshold 2.5).
5. `ulf_analysis.py` — Pc3/Pc4 polarization (XYZ-capable groups only). Since 2026-09-28 the near/far index normalizes each station's nightly Z/H by its own non-storm median before taking the pool median (the old raw-median version, kept as `*_raw`, shifted whenever a station joined or left the pool -- G11's zbn, Z/H 15–45, stopped 2024-12-23 and dropped the near index ~20× at the start of G11's pre-event window). `--from-daily` rebuilds only `ulf_near_far_index.csv` from the existing `ulf_daily.csv`; `--self-test` checks the dropout case. **Primary Pc3 rank test after the fix: p = 0.049** (was 0.016; one per family 0.048; without the discovery case G10 0.055; leave-one-out 0.020–0.103; placebo block just before each real window 0.96). G10's 2024-03-30 dip goes from z −4.64 to −3.29 and is no longer the series minimum (04-14, after the mainshock, is). SEA M≥5.5 day −4 no longer leaves the null band. See `docs/statistical_methods.md` "2026-09-28 修正".
6. `verify_pipeline.py` — hard-gates everything downstream on all checks passing.
7–8. (`--full-report` only, G10-specific narrative) `prepare_report_data.py` + `build_artifact.py` → `output/geomag_precursor_report.html`.

Cross-group validation on top of that (`scripts/run_validation_pipeline.sh`, assumes steps 1–6 already ran for every `ULF_GROUPS` member): `fetch_earthquake_catalog.py` (×2 magnitude thresholds) → `surrogate_test.py --all` → `superposed_epoch_analysis.py` (×2) → `backtest_rule.py` (×2) → `prepare_validation_report_data.py` + `build_validation_report.py` → `output/geomag_precursor_validation_report.html`.

**Interactive pipeline flowchart** (every step's summary, I/O, statistical method and math formula): `output/pipeline_flowchart.html`, built from the single content file `docs/flowchart/pipeline_flow.yaml`. After adding/renaming/moving anything under `scripts/`, run `.venv/bin/python3 docs/flowchart/check_flowchart.py` (verifies cited file:line/symbols, constants, and method IDs vs `docs/statistical_methods.md`), then `build_flowchart.py`. See `docs/flowchart/README.md`.

## Coseismic pipeline (second-scale, per event)

Not wired into the shell orchestrators above — run directly:

- `coseismic_step_analysis.py [--all]` — per-event (all 137, not just anchors) step/spike detector at the origin second, ±180s search window (`SCAN_HALF_SEC`, widened from ±120s 2026-08-16), null distribution from random reference times excluding every real event in that group. Output: `data/interim/coseismic_step_analysis/`.
- `coseismic_stacking_analysis.py [--all]` — cross-event superposed-epoch stacking of the above statistics (bootstrap CI90 over events + null band over random per-event reference times, mirroring `superposed_epoch_analysis.py`'s day-scale design at 1-second resolution instead). Output: `data/interim/coseismic_stacking_analysis/`.
- `seismometer_comparison.py [--all]` — compares the geomagnetic anomaly's timing against independent seismometer/accelerometer waveform data (see below) to help separate "instrument shaken by real ground motion" from "the magnetic field itself changed", using the same ±180s window (`SEARCH_HALF_SEC`) on its own independent computation. Output: `data/interim/seismometer_comparison/`. Since 2026-09-25 it also writes PGA (counts / PZ `SENSITIVITY`, max horizontal, own-event window only), a geomagnetic 1Hz noise ratio, a `geomag_dropout` flag, and `alignment_verdict_gated` (the verdict only where the event's own step30 p from `coseismic_step_analysis.py` is < 0.05 -- which is just 7 of 93 events -- G10 2024-04-03's xcg step30 p sits right on the boundary (0.044 with the per-key null streams of 2026-09-26; 0.0500/0.0505 in earlier runs whose shared stream reshuffled whenever any event changed) -- so the ungated verdicts are mostly timing noise).
- `coseismic_dose_response.py [--all]` — the label-free replacement for reading the verdicts: Spearman of log PGA vs. geomagnetic anomaly size, p from within-group permutation. Result (2026-09-25, numbers updated 2026-09-28): the 1Hz noise ratio rises with PGA in every subset (all ρ=0.33 p=0.009; M≥6 ρ=0.50 p=0.016; seismic station ≤2 km ρ=0.89 p=0.010), and 6 of the 7 significant anomalies sit at PGA ≥100 gal -- consistent with the magnetometers being shaken -- but a genuine field change coupled to ground motion would also scale with PGA, so this test alone can't tell the two apart; only the co-located seismometer check (n=11) and the G10 waveform point to shaking specifically. 7/93 significant events is itself at chance level (4.7 expected at α=0.05, P(≥7) ≈ 0.2). Output: `data/interim/coseismic_dose_response/`.
- `coseismic_onset_moveout.py [--all] [--min-mag 6] [--group ID]` — near- vs. far-station onset timing (added 2026-09-25), needs no seismometer data: for every event × every station (not just the nearest 2), the first second in −60..+180 s where the 1 Hz first difference exceeds 4σ (σ = MAD of −660..−60 s) in ≥3 of 5 s; each station's false-trigger rate and null onset distribution come from 200 random reference times in the same loaded data. External disturbances would start everywhere at once; shaking noise follows the seismic waves. **Result (2026-09-25, numbers updated 2026-09-28): the onsets follow the seismic waves** (shaking, or a field change carried by ground motion -- the test can't tell which). Since 2026-09-28 all stations of an event share the same random null reference times and the primary p is **event-level** (one shared reference time per event per simulation): the old per-trigger Poisson-binomial treated an event's stations as independent, but they fire together (null variance 1.4–1.7× the independent value; in the self-test's network-wide-disturbance H0 the old p is falsely significant 17.5 % of the time at 5 %, the event-level one 5.0 %). M≥6 (69 events): 59 stations in their arrival window [R/6 − 5, R/3.5 + 30] s vs. 15.4 expected at a random time (null SD 4.5, p<0.0001; clean stations 24 vs 1.45). The event-level count doesn't condition on triggering, so more triggering at the origin also raises it -- M<6 gets event-level p=0.02 while the timing-only test (Poisson-binomial with variance × dispersion ratio) is not significant (z=0.8, p=0.21; M≥6 z=7.8, clean 6.8). Timing-only numbers for M≥6: 59 of 120 onsets vs. 21.8 expected (clean stations 24/32 vs 6.5), pooled within-event Theil–Sen slope 0.22 s/km (CI90 0.12–0.31, p=0.0005); trigger rate falls from 92 % at <30 km to 9 % beyond 250 km. Median residual vs. S arrival is −0.8 s over the original 49 events (47/88 in window vs 16.1) but −12.6 s overall. **Prior-event shaking is not the cause** (sensitivity added 2026-09-28): events with another registered event ≤900 s before them (`PRIOR_EVENT_EXCLUSION_SEC`; G8 2022-03-23c 106 s, G10 2024-04-23b 356 s, G10 2024-04-03b 797 s, G11 2025-01-21b 539 s) have their baseline and search window inside that shaking, so `prior_event_shaking` flags them and `m6_no_prior_shaking`/`m5_no_prior_shaking` in `summary.json` repeat every statistic without them (main `m6`/`m5` unchanged, kept because the main test was fixed first) -- M≥6 gives 55 in window vs 14.3 at a random time (event-level p<0.0001; timing-only 55/112 vs 19.9), clean stations 24 vs 1.27, slope 0.20 s/km (CI90 0.10–0.28, p=0.001), and the S residual is still −13.5 s. The negative median comes from far and noisy stations instead: among M≥6 non-flagged onsets it is −0.5 s at R<150 km vs −33.9 s at R≥150 km, and +2.2 s for clean stations vs −37.1 s for the rest (clean = own null false-trigger rate < 5 %, 31 onsets) -- a noise trigger lands anywhere in −60..+180 s, which is early relative to a late far-field S arrival. Only registered events can be checked (M5 registered only from 2024-09, the GDMS json is M≥6 only), so smaller aftershocks, e.g. before G10 2024-04-03c (2 h after the M7.2), can't be ruled out. No clean station triggers network-wide before the origin. Against the co-located seismometer (stations with false rate <0.2, n=11) the geomagnetic onset is a median 0.6 s after the seismic onset. `xcg` triggers at ~80 % of random times, so its onsets carry almost no timing information. M<6 timing is not significant (25/95 vs 21.3). Output: `data/interim/coseismic_onset_moveout/` (`onsets.csv`, `moveout_per_event.csv`, `summary.json`, `travel_time.png`); `--group`/`--min-mag` runs go to a `subset_*` subfolder.

Both `coseismic_step_analysis.py` and `seismometer_comparison.py` cap that ±180s window per-event (`_effective_half_sec`, a local copy in each script) at half the gap to the nearest *other* real event in the same group, since G10's 2024-04-23a/b are only 357s apart — without the cap, one event's search would reach into the other's real anomaly. Every other group's events are hours-to-years apart, so this only ever narrows that one pair (to ~178s). 2026-08-19 investigated widening the shared window further, to 300s/360s, specifically to check whether G20's 2025-12-24 `obs_lag=-178s` (2s from the 180s boundary) was itself a truncation artifact — it wasn't (unchanged at 300s), but the wider window changed several *other*, non-boundary-pinned events' reported peaks (a wider search finds a larger max-of-N by chance even under pure noise), so the widening was reverted; the per-event cap was kept as a genuine fix. See plan `artifact-wobbly-kettle.md` for the full investigation.

Geomagnetic storms are **flagged, not excluded** on this line (added 2026-09-24): `coseismic_step_analysis.py::_storm_status` looks up each event's UTC date in its group's `storm_days.csv` and writes `is_storm_day` / `is_storm_onset` (`None` if the event falls outside the fetched space-weather range) into the per-event JSON, `summary.csv`, and `all_events_run_summary.json` (which also gets a storm vs. quiet `storm_sensitivity` split). `coseismic_stacking_analysis.py` re-stacks with flagged events dropped (own rng, so main combos are unchanged) → `stack_storm_sensitivity.csv`. Why not exclude: the 1hr detrend + same-window null already absorb a storm-raised noise floor, and ~half of all events (70/137 storm-or-recovery, 43 onset as of 2026-09-27; 59/117 and 41 onset since the 2026-09-25 re-fetch -- G6/G7/G8's and G17's cached storm_days.csv had been built before their folders were extended, leaving G8 2022-05-09 and G17 2015-03-23 unflagged; `fetch_space_weather.py --check-cache` now also requires the cache to cover the folder's date range) fall on flagged days under the daily pipeline's definition. The remaining risk is a single SSC/substorm transient inside the ±180s window, which is what the flag lets you check. **Result (2026-09-24 run): storms dilute, they don't create.** Fraction of tests with p<0.05 is 6.7% on storm-or-recovery days vs 12.5% on quiet days (137 events, 2026-09-27; 5.3% vs 13.1% over the earlier 117) (12.2% before the re-fetch); dropping storm days *strengthens* the near-station H stack at lag +10s (137 events: step30 p 0.43→0.037, step90 0.74→0.013, spike 0.175→0.042, n 96→44; none of these 48 stacking p-values is multiplicity-corrected; before the 2026-09-26 correction of four minute-truncated origin times the peak sat at +25s) while far-station/F combos stay null — consistent with storm noise masking a shaking-timed signal, not faking one.

All three take `--self-test` (synthetic/known-file sanity checks, run automatically before touching real data unless `--self-test` is passed alone) and `--group <ID>` (repeatable, restricts to a subset).

Findings from this line of analysis were written up as a published Artifact report (not stored in this repo — see the report itself for the link) rather than a `docs/` file; the two persistence-vs-recovery flagship cases are G9's 2022-09-18 mainshock (`csg` station, persistent non-recovering offset) and G10's 2024-04-03 mainshock (`xcg` station, noise burst that recovers).

## Seismometer comparison data (fetched separately, not by this pipeline)

`seismometer_comparison.py` reads data that lives under `seismometer/` (moved into this project 2026-09-14, shortly after `Gx`; resolved via `SEISMIC_ROOT = common.GX_DATA_ROOT / "seismometer"`), fetched manually by the user rather than by any script here. As of the 2026-08-17 reorganization, each event gets one `GXX_MMDD/` folder (e.g. `seismometer/G10_0403/`) holding *both*:

- SAC PoleZero instrument-response files, `SAC_PZs_TW_<STA>_<CHAN>_<LOC>_<start>_<end>`.
- The matching miniSEED waveform, `GXX_MMDD_w.mseed` (100Hz, `event_utc-60s` to `event_utc+600s`).

(Before 2026-08-17 the mseed files sat loose directly under this project's parent directory, one level up from their PZ folder; before 2026-08-16 only the 16 anchor events had been fetched at all.)

**Coverage is 111/137 events** (the 20 M≥6 events registered 2026-09-27 have no seismometer data fetched yet; updated 2026-09-25: the remaining 80 fetchable events -- 14 M>=6, 66 M5 -- came in one batch, laid out as `G<NN>_<MMDD><suffix>/` and auto-registered by `seismometer_comparison.py::_register_convention_dirs`; the 6 left are the pre-2012 permanent gap below, status `no_data_pre_2012`). History: coverage was 30/49 as of 2026-09-20 — the 11 non-anchor events (foreshocks/aftershocks within the originally-covered 16 groups) plus a previously-missing G9 foreshock were fetched and wired into `SEISMIC_DATA_DIRS` 2026-08-16 (27 events); G22's 1 event and G23's 1 event + G24's 1 event (split 2026-09-22 from a single merged "G23" with 2 events; mseed fetched 2026-09-19, PoleZero files added 2026-09-20) were wired in 2026-09-20. Of the 19 uncovered events, **6 are a permanent, structural gap**: G14's 5 (its earliest event, 2009-07-14, predates the seismic data source's ~2012 cutoff) plus G21's 1 (2010-11-21, same reason — confirmed 2026-09-19 that the data source's fetchable range only starts 2012-01-01). The remaining **13** (spread across G5, G6_G7_G8, G9, G10, G11, G12, G15, G17) are simply not yet fetched. The G22/G23 PZ sets are missing a few stations that appear in the mseed (G22: CHK/ELD loc 11 and HEN; G23: HEN and SSH) — those traces get `no_pz_epoch`, which doesn't affect the three events' results (the picked stations are ILA for G22, NSK for G23). `seismometer_comparison.py`'s `coverage_summary.json` output enumerates all 49 events with an honest status rather than silently implying full coverage.

**The folder-name ↔ mseed-filename mapping is hardcoded** in `seismometer_comparison.py::SEISMIC_DATA_DIRS` (keyed by `(group_id, event.date)`, mseed path resolved as `SEISMIC_ROOT / pz_dir / mseed_filename` since the reorg), not derived by pattern — most groups follow `G0N_MMDD` for both the PZ folder and `G0N_MMDD_w.mseed`, but **G9 doesn't**: its PZ folder is `G9_0918` (no leading zero) while its mseed file is `G09_0918_w.mseed` (leading zero). `G10`'s two 2024-04-23 events also share one PZ folder (`G10_0422`) but each keep their own mseed file inside it. Don't assume the pattern holds for any future additions; check `SEISMIC_DATA_DIRS` directly.

Other confirmed quirks handled in `seismometer_comparison.py` (see its module docstring for detail): mseed traces are not uniformly the full ±60s/+600s window (some stations are short triggered-accelerograph recordings, flagged `triggered_short_trace`); some `(station, channel)` pairs have duplicate traces in the same file (handled via `Stream.merge()`); PZ epoch selection matches on the trace's own `location` code, not a guess.

## Version control

This project got a fresh `git init` on 2026-09-14 (its earlier 2026-08-14 history was dropped when `G1`..`G24` and, shortly after, `seismometer/` were moved in from the parent directory and it was re-established as a self-contained standalone repo — this is the only repo among the moved material that's published). `.gitignore` excludes `.venv/`, `__pycache__/`, `data/` (regenerable analysis intermediates/outputs, ~850MB, reproducible from `scripts/` + the raw data with the project's fixed seed `20260805` — not worth version-controlling), each `Gx/*` raw day-file/batch-archive (~18GB, too big for git; each `Gx/CLAUDE.md` is explicitly kept via a `!` negation), and all of `seismometer/` (~384MB, all binary, no docs to keep). Tracked: `scripts/`, `docs/`, `output/` (the built HTML reports), `requirements.txt`, each `Gx/CLAUDE.md`, this file, `README.md`.

## Documentation

- Math in markdown uses `$...$` inline and `$$...$$` display. Never use `\(...\)` or `\[...\]`.
- Documentation, reports and figure labels are for Chinese grad-school application materials. Use Chinese labels on figures unless told otherwise.
- `docs/statistical_methods.md` must stay in sync with the code. When a formula or function changes, update its formula and line references.
- Files under `docs/` planning folders are historical records. Do not edit them unless asked.

## Workflow

- At the end of each work block: commit with clear split commits, update the progress notes file, and update memory with pending items (e.g. open verify failures).
- Rendered outputs (HTML slides, flowcharts, `.docx`, charts) must be visually verified before reporting them done. Render the BUILT output, not the template.
- For LICENSE files (CC-BY-4.0), do not fetch the license text from the web because it triggers content-filter errors. Ask the user to add it, or write a short header linking to the official URL.
