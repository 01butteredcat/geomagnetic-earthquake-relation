# geomag_precursor

用台灣中央氣象署（CWA）地磁觀測網的 1 秒解析度資料，統計檢定「地震前有沒有可偵測的地磁前兆」與「地震發生瞬間，地磁場有沒有真實變化（同震效應）」這兩個問題。資料涵蓋 2009–2026 年間 49 起 CWA M≥6.0 地震（依時間/地點分成 24 個獨立事件序列，避免同一序列的前震/主震/餘震被當成獨立樣本重複計算；詳細事件/日期範圍考證見本 repo 的 `CLAUDE.md`）。

這份文件是給**不熟悉這個 repo 的外部使用者**看的操作手冊：如何設置環境、去哪裡拿資料、怎麼跑整條分析、有哪些參數可以調整。如果你是 Claude Code（或其他 AI 助理）在這個 repo 裡工作，請改讀 `CLAUDE.md`（給 AI 助理的操作指南，假設了較多背景知識）。

## 這個專案在做什麼：兩條獨立的分析線

- **日尺度前兆篩選線**：找主震前「數天到數週」內，地磁場有沒有異常，可能是地震前兆訊號。核心方法是「近站－遠站」相減（濾掉太陽磁暴等全球共同背景訊號）+ MAD z-score 異常日標記，另外對向量測站足夠的組別跑 ULF Pc3/Pc4 極化分析。
- **同震（秒級）分析線**：找地震發生那一刻附近（±180 秒內），地磁場有沒有階躍/尖峰訊號，並試圖分辨這是磁場真的變了，還是測站被地動搖晃出的儀器雜訊。方法上更年輕、還在活躍推進，目前還沒有 shell script 自動串接。

兩條線各自獨立、互不依賴，可以只跑其中一條。兩者的結論目前都是：**沒有找到具統計顯著性的前兆/同震訊號**——但同震線的兩個案例研究（G9、G10 主震）呈現出跟獨立地震儀資料吻合的、值得繼續深究的訊號樣貌。詳見各自報告（第 7 節）。

## 環境設置

```bash
cd geomag_precursor
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

- 沒有 `setup.py`/`pyproject.toml`，也沒有明文寫死 Python 版本下限；開發與測試環境是 **Python 3.12.3**，建議至少 3.10+。
- `requirements.txt` 釘死四個套件版本：`numpy==2.5.1`、`pandas==3.0.3`、`scipy==1.18.0`、`obspy==1.5.0`（`obspy` 只有同震線的 `seismometer_comparison.py` 需要，用來讀 miniSEED 波形與 SAC PoleZero 儀器響應檔）。
- **一律用 `.venv/bin/python3` 執行 `scripts/` 底下的腳本，不要用系統的裸 `python3`**——三支 shell orchestrator 也是這樣寫的（優先找 `../.venv/bin/python3`，找不到才 fallback 到系統 `python3`）。

## 資料需求（怎麼拿到原始資料）

跑得動這條 pipeline 之前，必須先準備好原始資料。地磁資料（`G1`~`G24`）的資料夾本身在這個 repo 裡（跟著 `scripts/`、`docs/` 一起被 git 追蹤），但**實際的原始資料檔不在 git 裡**（`.gitignore` 排除，只留每組的 `CLAUDE.md`）——資料量約 18GB，遠超 GitHub 建議的 repo 大小。

| 資料 | 放置路徑 | 從哪裡拿 | 誰負責準備 |
|---|---|---|---|
| 地磁 1 秒解析度資料（IAGA-2002 格式 `.sec`，也接受 `.sec.gz`/`.tgz`） | `G1/`、`G2_G3/`、…、`G23/`（本 repo 底下，每個原始資料夾對應一或多個 group_id，命名 `<station><YYYYMMDD>dsec.sec`；2026-09-14 整併進本目錄，讓這個 repo 自成一個完整自包的專案） | CWA GDMS 地磁資料下載系統（`gdms.cwb.gov.tw`），需申請帳號登入 | **使用者手動下載**，repo 完全不含 |
| 地震儀/加速度計波形（SAC PoleZero 響應檔 + miniSEED） | `seismometer/<GXX_MMDD>/`（本 repo 底下，如 `seismometer/G10_0403/`；2026-09-14 跟著 `Gx` 一起整併進本目錄） | 使用者自行取得的地震儀網路資料 | **使用者手動下載**，且下載後要同步更新 `seismometer_comparison.py::SEISMIC_DATA_DIRS` 這個手動維護的對照表（folder 名稱 ↔ mseed 檔名，兩者並非永遠一致，例如 G9 的 folder 是 `G9_0918`、mseed 檔是 `G09_0918_w.mseed`） |
| 空間天氣指數（Dst/Kp） | `data/external/<group>/{kp,dst}.csv`（在 repo 內，但被 `.gitignore` 排除） | `kp.gfz.de`、`wdc.kugi.kyoto-u.ac.jp` | `fetch_space_weather.py` **自動連網抓取**；只有 Kp＋Dst 完整（摘要 `confidence: high`）的快取才會被沿用，抓取不完整會自動重試（要強制重抓可刪除該組的 `storm_days.csv`） |
| 擴充地震目錄 | `data/external/extended_catalog_m5.{0,5}.csv` | USGS FDSN Event API | `fetch_earthquake_catalog.py` **自動連網抓取** |

只有地磁資料與地震儀資料需要使用者自己張羅；後兩項腳本會自己處理（前提是要有對外網路連線）。

## 專案結構

```
geomag_precursor/
├── README.md              # 本文件
├── CLAUDE.md               # 給 AI 助理看的操作指南（背景知識、已知陷阱更詳細）
├── requirements.txt
├── .venv/                  # 自行建立，.gitignore 排除
├── scripts/                # 全部分析程式碼（見下）
├── docs/                   # 研究規劃筆記（資料抓取範圍表、候選事件評估），非操作手冊
├── data/
│   ├── external/            # 自動抓取的空間天氣、地震目錄快取
│   └── interim/             # 各腳本的中繼與最終輸出（見「輸出產物」一節）
│                             # 整個 data/ 都被 .gitignore 排除（~850MB，可從原始資料+固定 seed 重現）
└── output/                 # 3 份自含式 HTML 報告，這個資料夾「有」被 git 追蹤
```

## 如何執行

### 日尺度前兆篩選線

單一組別跑一次（6 步：解析 → 時區驗證 → 抓磁暴日 → 算異常指數 → ULF 極化 → 驗證關卡）：

```bash
cd scripts
./run_pipeline.sh --group G10
./run_pipeline.sh --group G10 --full-report   # 多跑 2 步，另外產生敘事型 HTML 報告（只有 G10 的手寫敘事文字，其他組別跑 --full-report 沒有意義）
```

全部 24 個 `group_id`（放在 20 個原始資料夾裡：G2、G3 共用 `G2_G3/`，G6、G7、G8 共用 `G6_G7_G8/`，G23、G24 共用 `G23/`；G2_G3/G6_G7_G8 於 2026-09-20、G23/G24 於 2026-09-22 各自拆成獨立的組別，每組有自己的 anchor 與 `data/interim/<group>/`）各跑一輪：

```bash
./run_all_groups.sh
```

單組失敗（`verify_pipeline.py` 檢查沒過）不會中止整批，最後彙整成 `data/interim/all_groups_run_summary.json` 並印出每組的驗證結果摘要。

跨組彙整（假設 `run_all_groups.sh` 已跑完）：

```bash
./cross_group_analysis.py   # 沒有命令列參數；輸出 data/interim/cross_group_summary.{json,md}
```

跨組正式驗證方法論（假設 `ULF_GROUPS` 這 14 組已跑完前 6 步）：

```bash
./run_validation_pipeline.sh   # 抓擴充地震目錄 → 洗牌檢定 → 疊加曆元分析 → 規則回測 → output/geomag_precursor_validation_report.html
```

**日尺度線各腳本的命令列參數：**

| 腳本 | 參數 | 說明 |
|---|---|---|
| `build_daily_features.py` / `timezone_check.py` / `fetch_space_weather.py` / `compute_indices.py` / `ulf_analysis.py` / `verify_pipeline.py` / `prepare_report_data.py` / `build_artifact.py` | `--group <ID>`（必填） | 單一 group_id，如 `G10`、`G3` |
| `fetch_earthquake_catalog.py` | `--min-mag <float>`（預設 5.5）、`--output <path>`、`--groups <ID...>`（預設全部 `ULF_GROUPS`） | 抓 USGS 擴充地震目錄 |
| `surrogate_test.py` | `--group <ID>`、`--all` | 洗牌顯著性檢定 |
| `superposed_epoch_analysis.py` / `backtest_rule.py` | `--catalog <path>`（必填）、`--label <str>`（必填，如 `m5.5`） | 疊加曆元分析／規則回測 |
| `cross_group_analysis.py` / `prepare_validation_report_data.py` / `build_validation_report.py` | （無參數） | 直接執行即可 |

### 同震（秒級）分析線

**這 4 支腳本目前沒有 shell script 自動串接**——要照順序個別直接執行：

```bash
cd scripts
.venv/bin/python3 coseismic_step_analysis.py --all        # 1. 49 起事件逐一 step/spike 偵測
.venv/bin/python3 coseismic_stacking_analysis.py --all    # 2. 跨事件疊加（不依賴步驟 1 的輸出檔，是獨立重算）
.venv/bin/python3 seismometer_comparison.py --all         # 3. 比對地震儀資料（依賴 SEISMIC_DATA_DIRS 手動維護表 + seismometer/ 底下的資料）
.venv/bin/python3 coseismic_joint_analysis.py --all        # 4. 聯合統計檢定（依賴步驟 3 輸出的 comparison_summary.csv）
```

**依賴關係**：步驟 2（疊加）獨立於步驟 1，重算一次自己的統計量，不吃步驟 1 的快取；步驟 3 依賴 `seismometer/` 底下的原始波形資料是否齊全；步驟 4 直接讀步驟 3 輸出的 `data/interim/seismometer_comparison/comparison_summary.csv`，所以**步驟 3 一定要先成功跑過**才能跑步驟 4。

**同震線各腳本的命令列參數：**

| 腳本 | 參數 | 說明 |
|---|---|---|
| `coseismic_step_analysis.py` | `--all`（跑全部 49 事件，預設只跑 G10 錨定事件）、`--no-injection`（跳過正對照組注入測試）、`--self-test` | **沒有 `--group`**，要嘛全部跑、要嘛只跑 G10 |
| `coseismic_stacking_analysis.py` | `--all`、`--group <ID>`（可重複指定）、`--self-test` | |
| `seismometer_comparison.py` | `--all`、`--group <ID>`（可重複指定）、`--geomag-station <站碼>`（手動排查單一測站用）、`--self-test` | |
| `coseismic_joint_analysis.py` | `--all`、`--self-test` | 沒有 `--group` |

四支都支援 `--self-test`（合成資料/已知檔案的健檢），且**除非單獨傳 `--self-test`，否則正式分析前會自動先跑一次 self-test，失敗就中止**——第一次在新環境執行前，建議先手動跑一次 `--self-test` 確認環境正常：

```bash
.venv/bin/python3 coseismic_step_analysis.py --self-test
```

## 關鍵參數（可調整）

所有參數都是各腳本檔案裡模組層級的大寫常數，改完直接重跑對應腳本即可生效，不需要額外設定檔或環境變數。

### 日尺度線

| 參數 | 定義檔案 | 目前值 | 用途 |
|---|---|---|---|
| `MIN_STATIONS_FOR_METHOD` | `common.py` | 5 | 一個測站池（F-only 或 XYZ）至少要幾站，該方法才會在這一組執行，否則整組跳過 |
| `N_NEAR_STATIONS` / `N_FAR_STATIONS` | `common.py` | 3 / 2 | 近站／遠站相減法各取幾座測站 |
| `OUTAGE_PCT_MISSING_THRESHOLD` | `common.py` | 0.05 | 站/日缺測比例超過此值，該天不列入「乾淨基線日」 |
| `SPIKE_THRESHOLD_NT` | `build_daily_features.py` | 300.0 | 秒級尖峰雜訊濾除門檻（nT） |
| `NIGHT_HOURS_UTC` | `compute_indices.py`／`ulf_analysis.py`（各自定義一份） | `{17,18,19}`（本地 01:00–03:59） | 只用夜間（磁場最安靜）時段算異常指數 |
| `MIN_NIGHT_MINUTES` | `compute_indices.py`（`verify_pipeline.py` 匯入同一值） | 90（夜間 180 分鐘的 50%） | 某站某夜有效分鐘數低於此值，該站當夜不貢獻指數；`verify_pipeline` 的 outage 判定也用同一標準（2026-09-20 新增，原本整天 `pct_missing>5%` 與實際使用的夜間窗口不一致） |
| `TRAILING_WINDOW_DAYS` | `compute_indices.py` | 21 | 滾動基線窗長（天） |
| `CANDIDATE_Z_THRESHOLD` | `compute_indices.py` | 2.5 | MAD z-score 超過此值標記為候選異常日 |
| `ZOOM_WINDOW_DAYS` | `ulf_analysis.py` | 12 | 錨定事件前後放大檢視的天數 |
| `KP_STORM_THRESHOLD` / `DST_STORM_THRESHOLD` | `fetch_space_weather.py` | 5.0 / −30.0 nT | 磁暴日判定門檻 |
| `WINDOWS_DAYS` | `cross_group_analysis.py`／`backtest_rule.py` | `[7, 14, 30]` | 震前重現率/回測用的窗口長度（天） |
| `N_SURROGATES` | `surrogate_test.py` | 2000 | 洗牌檢定抽樣次數 |
| `N_BOOTSTRAP` / `N_NULL` | `superposed_epoch_analysis.py` | 2000 / 1000 | 疊加曆元分析的 bootstrap CI 抽樣次數／null 抽樣次數 |
| `WINDOW_BEFORE_DAYS` / `WINDOW_AFTER_DAYS` | `superposed_epoch_analysis.py` | 30 / 10 | 疊加曆元窗口（震前/震後天數） |
| `FIXED_RULE_THRESHOLD` | `stat_utils.py` | −4.1 | 舊報告引用的固定判定門檻（z-score） |
| `ULF_GROUPS` | `fetch_earthquake_catalog.py`／`surrogate_test.py`／`superposed_epoch_analysis.py`／`backtest_rule.py`（**5 個檔案各自定義一份，同一組值**，另含 `prepare_validation_report_data.py`） | `("G4","G5","G6","G7","G8","G9","G10","G11","G12","G13","G19","G20","G23","G24")` | 哪些組別的向量站夠格跑跨組驗證。**改這個集合要記得同步改全部 5 個檔案**（`G23`/`G24` 的原始資料為 2026-08-20 補上，2020 年資料雖群組編號高但時間上晚於向量站升級，屬向量站齊全世代；G23/G24 於 2026-09-22 從合併的單一 `G23` 拆分而來） |

### 同震（秒級）線

| 參數 | 定義檔案 | 目前值 | 用途 |
|---|---|---|---|
| `SCAN_HALF_SEC` | `coseismic_step_analysis.py` | 180（秒，逐事件上限值，見下方提醒） | 在報震秒前後搜尋 step/spike 極值的目標半窗 |
| `STEP_WINDOWS_SEC` | `coseismic_step_analysis.py` | `(10, 30, 90)` | 三種 step 偵測統計量的窗長（秒） |
| `N_NULL` | `coseismic_step_analysis.py` | 2000 | null 分布抽樣次數 |
| `EXCLUSION_BUFFER_SEC` | `coseismic_step_analysis.py` | 600 | null 參考時間點必須離任何真實地震事件至少多遠（秒） |
| `MISSING_FRACTION_THRESHOLD` | `coseismic_step_analysis.py` | 0.10 | 掃描窗內 NaN 比例超過此值，該次檢定作廢 |
| `INJECTION_AMPLITUDES_SIGMA` | `coseismic_step_analysis.py` | `(1,2,3,5,10)` | 正對照組合成訊號注入的振幅倍率（相對 local null MAD） |
| `STACK_HALF_SEC` | `coseismic_stacking_analysis.py` | 300（秒） | 跨事件疊加的半窗，故意比 `SCAN_HALF_SEC` 寬，讓疊加後的形狀看得完整 |
| `N_BOOTSTRAP` / `N_NULL_STACK` | `coseismic_stacking_analysis.py` | 2000 / 1000 | 疊加分析的 bootstrap CI／null band 抽樣次數 |
| `MIN_OFF_EVENT_SAMPLES` | `coseismic_stacking_analysis.py`／`seismometer_comparison.py` | 200 | 離峰時段樣本數低於此值，中位數/MAD 估計太不穩定，判定該次無效 |
| `GEOMAG_HALF_SEC` | `seismometer_comparison.py` | 240（秒，目標值） | 地磁側 z-score step30 剖面的擷取窗 |
| `SEARCH_HALF_SEC` | `seismometer_comparison.py` | 180（秒，逐事件上限值） | 地磁側峰值搜尋子窗 |
| `PERSISTENCE_Z_THRESHOLD` | `seismometer_comparison.py` | 2.0 | 判定「震後異常持續未恢復」的 z-score 門檻 |
| `PERSISTENCE_TAIL_FRACTION` | `seismometer_comparison.py` | 0.3 | 用波形尾段多少比例判定持續性 |
| `ALIGNMENT_TOLERANCE_SEC` | `seismometer_comparison.py` | 5 | 判定地磁異常時間點是否「對齊」地動窗口的容忍秒數 |
| `N_PERM` | `coseismic_joint_analysis.py` | 2000 | 兩臂標籤置換檢定的排列次數 |

**三則重要提醒：**

1. **`SCAN_HALF_SEC`（`coseismic_step_analysis.py`）與 `SEARCH_HALF_SEC`（`seismometer_comparison.py`）是兩份各自獨立寫死的同名概念參數，不是共用變數——改一個不會自動同步另一個**，兩處都要改才會一致。
2. 這兩個參數其實是**逐事件的「目標值」，不是硬性套用的統一值**：兩支腳本都有一個 `_effective_half_sec()` 函式，會把單一事件的實際半窗，依「跟同組最近另一起事件的時間間隔」自動收窄（例如 G10 的 2024-04-23a/b 只相隔 357 秒，實際半窗會被收到約 178 秒），避免一起事件的搜尋窗吃到另一起事件自己的訊號。曾經試過把兩個目標值從 180 秒進一步放寬到 300/360 秒，結果發現搜尋窗越寬、單純靠雜訊也越容易找到更大的極值（look-elsewhere 效應），讓好幾起事件的結果變得不穩定，因此最終定案維持在 180/240 秒——**調整這兩個值之前，务必意識到這個副作用**。
3. 全專案共用同一組固定亂數種子 **`SEED = 20260805`**（`coseismic_step_analysis.py`、`surrogate_test.py`、`superposed_epoch_analysis.py`、`coseismic_joint_analysis.py` 都是這個值），確保重跑結果可重現；如果要測試「換一個種子會不會改變結論」，記得同步全部改。

## 輸出產物與如何解讀

- `data/interim/<group>/`（24 個資料夾，每個 group_id 一個）——各組日尺度線的中繼與最終產物：`daily_features.csv`、`minute_series_<station>.parquet`、`storm_days.csv`、`ulf_near_far_index.csv`、`local_anomaly_index.csv`、`candidate_windows.json`、`verification_report.json`。
- `data/interim/coseismic_step_analysis/`、`coseismic_stacking_analysis/`、`seismometer_comparison/`、`coseismic_joint_analysis/`——同震線四支腳本的輸出，各自有 `*_summary.csv`（總覽表）與逐事件/逐組合的 JSON 詳細檔。
- `data/interim/cross_group_summary.{json,md}`、`all_groups_run_summary.json`、`validation_report_data.json`——跨組彙整與驗證報告用的中繼資料。
- `output/`（**這個資料夾有被 git 追蹤**，其餘 `data/` 都沒有）：
  - `geomag_precursor_report.html` / `geomag_precursor_report_G10.html` — 日尺度線 G10 單事件敘事報告（手寫敘事文字只對 G10 有意義，`--full-report` 只建議對 G10 跑）
  - `geomag_precursor_validation_report.html` — 跨組驗證方法論報告
- **同震線目前沒有對應的 `output/` 報告產生腳本**：它的成果是手動整理成外部發布的網頁報告，不在這個 repo 裡（`data/interim/coseismic_*` 底下的 CSV/JSON 是完整可重現的原始數字/圖表資料，任何人都能拿這些資料自己重建報告）。

## 已知限制與踩雷提醒

- **G14 的 5 起事件永遠不會有地震儀比對資料**：地動資料源只回溯到約 2012 年，G14 最早的事件是 2009 年，結構性缺口，不是還沒抓而已。
- **測站有純量／向量世代分野**：G14–G18（2009–2016 年資料）與 G21、G22（2010、2012 年資料）全數是純量站（只有總磁場 F），完全跑不動 H/Z 篩選法與 ULF 極化分析；G4 之後（含 G19、G20）以及 G23、G24（2020 年資料，雖群組編號高但時間上已晚於向量站升級）才是向量站齊全的世代——這正是 `ULF_GROUPS` 只列這 14 組的原因（2026-09-20 `G6_G7_G8` 拆成 G6、G7、G8 後是 13 組，2026-09-22 `G23` 拆成 G23、G24 後變 14 組）。
- **G21 永遠沒有地震儀比對資料**：2010-11-21 的事件早於地動資料源的回溯起點（2012-01-01），與 G14 同屬結構性缺口。G22、G23 已於 2026-09-20 補齊 mseed 與 PoleZero 並納入 `SEISMIC_DATA_DIRS`（`seismometer/` 現有 29 個 `GXX_MMDD` 資料夾），G23 於 2026-09-22 拆出 G24 後兩者共用同一批已補齊的資料，目前 49 個事件中 30 個可做比對。
- **同資料夾的兄弟組並非完全獨立**：G2／G3（共用 `G2_G3/`）與 G6／G7／G8（共用 `G6_G7_G8/`）的原始資料相同，虛無抽樣、`is_known_event` 旗標與窗口上限已改用整個資料夾的所有事件（`events.py::folder_events()`），回測也把兄弟事件的震前窗口排除在命中／誤報／分母之外；但 `cross_group_analysis.py` 的 `sliding_baseline_rate` 基線窗口仍可能含到兄弟事件附近的日期，這一點尚未處理。也因此跨組檢定裡這 5 組的結果並非五個獨立樣本。
- **日尺度候選日對測站池很敏感**：G2 與 G3 用同一批原始資料，但近站各自依震央挑選（G2 為 hln,slg,sme；G3 為 ncg,hln,lyn），候選日從 5 天變成 13 天，14 天窗口的 F 法命中因此由 4/9 變 5/10（p 0.452→0.298，仍不顯著）。解讀單組結果時要記得這一點。
- **`verify_pipeline` 的 `storm_cancellation_test` 目前失敗的組別（2026-09-22，含 bootstrap CI 改版與 G23/G24 拆分後的最新結果）**：G1、G2、G3、G11、G12、G13、G14、G23、G24。
  - 判定邏輯：彙總該組所有磁暴日（Kp≥5 或 Dst≤−30）算出 `median|local_anomaly_index| / median|near_index|` 點估計後，另外對這批磁暴日做 90% bootstrap 信賴區間（2000 次重抽樣、固定種子 `20260805`），**pass 現在要求 CI 上界 < 0.6**，比單看點估計更嚴格（`scripts/verify_pipeline.py::check_storm_cancellation`，2026-09-22 改版）。磁暴日少於 5 天判 inconclusive（G4）；少於 15 天會多一則 `low_bootstrap_power` 警告但不影響判定（不改變 pass/fail，只是提醒 CI 精確度較粗）；磁暴訊號中位數低於 0.5 z 判 inconclusive。擬合天數少於 30 天只示警，不影響判定。
  - 這次改版是**非對稱收緊**：CI 上界必然 ≥ 點估計，所以原本點估計就 fail 的組別（G2、G3、G11）依然 fail；G12（點估計 0.51）、G13（0.40）、G14（0.48）原本點估計壓線 pass，但因磁暴日樣本數少、CI 偏寬（上界分別 0.70、0.78、0.72）而翻成 fail；G1（點估計 0.38，僅 5 個磁暴日）CI 上界達 0.72 同樣翻成 fail。
  - G23（2020-07-26 anchor）與 G24（2020-06-14 anchor，2026-09-22 從合併的單一 G23 拆分而來）目前**點估計與 CI 完全相同**（1.03，CI 上界 1.09）：兩者拆分後近／遠測站池剛好一致（震央只差 1.6km，遠小於任何站間距離），且這個 check 讀的 `local_anomaly_index.csv` 完全不依賴 anchor 日期，只依賴測站池，所以兩組數值必然一致，並非計算錯誤。遠站與近站相關性極低、回歸斜率約 0.12，這兩組的候選日都應視為磁暴訊號可能未被濾乾淨。
  - `baseline_window_excludes_storms`（2026-09-20 記錄的失敗組別：G6、G11、G19，尚未在本次改版後重新核對）：事前基線窗口內乾淨天數不足，是真實的低基線（補入 Dst 後 G11 的磁暴／恢復期占範圍約七成），只在文件註明。（附註：G23 拆分後這項檢查對 G23 為 pass、對 G24 為 inconclusive——G24 的震前窗口剛好沒有任何一天真的碰到磁暴日，屬於「這次沒機會測到失效情境」的 inconclusive，不是失敗。）
- **Kp/Dst 抓取失敗不再被快取沿用**：`fetch_space_weather.py` 會保留已抓到的月份並在 `storm_days_summary.json` 記錄 `dst_missing_months`，缺漏時 `confidence` 降為 `medium`；`run_pipeline.sh` 只在 `high` 信心時才沿用快取。G1／G2／G3／G11 原本缺 Dst 卻標 high（2026-09-20 已補抓，候選日 G1 6→8、G2 5→0、G3 13→5、G11 H 25→3／Z 31→6）。
- **`ttn`（卑南）測站有已知的資料缺口**：2024 年 12 月起疑似永久停站，G12、G13、G20 完全沒有這一站；G11、G19 各有一段較短的缺測期。細節見本 repo 的 `CLAUDE.md`。
- **`seismometer_comparison.py::SEISMIC_DATA_DIRS` 是手動維護的對照表**，新增地震儀資料要手動同步更新這個表，folder 名稱與 mseed 檔名的對應規則並非永遠一致（已知例外：G9）。
- 空間天氣抓取（`fetch_space_weather.py`）已存在的 `storm_days.csv` 會直接跳過重抓，要強制更新需手動刪除該檔案。

## 延伸閱讀

- `CLAUDE.md` — 給 AI 助理的操作指南，背景知識與已知陷阱記錄得更詳細（例如同震線的完整方法論演進過程，也包含原始地磁 `.sec` 資料本身的格式、測站代碼歷史、各組資料涵蓋範圍）。
- `docs/13_groups_fetch_ranges.md` — G1–G13 資料抓取時間範圍規劃與理由。
- `docs/candidate_groups_G14_G20.md` — G14–G20 候選地震評估紀錄。
- `docs/candidate_groups_from_GDMScatalog.md`、`docs/candidate_fetch_ranges_from_GDMScatalog.md` — G21–G23 候選事件評估與抓取範圍規劃（源自使用者提供的 CWA GDMS 區域目錄 `GDMScatalog.json`）。
