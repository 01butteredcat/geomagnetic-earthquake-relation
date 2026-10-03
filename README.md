# geomag_precursor

這個專案是我用中央氣象署（CWA）地磁觀測網的 1 秒資料，去看地震前地磁會不會有異常，然後地震發生那一刻地磁有沒有真的跳一下。資料是 2009 到 2026 年的 137 起地震，69 起是 M6 以上、68 起是 M5 級，因為很多其實是同一串的前震、主震、餘震，所以我把它們分成 24 組，免得同一串被當成很多個獨立樣本重複算。每起地震的詳細資料都在 `scripts/events.py`，地震參數全部來自 CWA 目錄。

這份 README 是寫給第一次看到這個 repo 的人，講怎麼裝、資料去哪裡拿、怎麼跑。比較細的背景和踩過的坑我都記在 `NOTES.md`。

## 我在看的兩件事

一個是日尺度的前兆，就是主震前幾天到幾週，震央附近的地磁有沒有怪怪的。做法是拿近震央的測站減掉遠的測站，把磁暴這種全台一起動的訊號扣掉，再用 MAD z-score 標出異常的日子，向量站夠多的組別另外跑 ULF Pc3/Pc4 的 Z/H 極化比。結果是大部分指標都不顯著，Pc3 的 Z/H 在震前 30 天有偏低，合併起來 p = 0.049，但拿掉最早發現的 G10 之後就變 0.055，而且方法是看過資料才定下來的，所以我只把它當探索性的結果。

另一個是同震，也就是地震發生前後 ±180 秒內地磁有沒有跳。結果是有跳，但跳的時間是跟著震波走的，離震央近的站先跳、遠的站後跳，雜訊大小也跟著地動加速度變大，最直接的解釋就是磁力儀被震到了，不太像磁場本身真的改變。

兩條線互相獨立，只想跑其中一條也可以。

## 環境

```bash
cd geomag_precursor
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

我是用 Python 3.12.3 開發的，3.10 以上應該都可以。套件只有 numpy、pandas、scipy 跟 obspy，版本都釘在 `requirements.txt`，obspy 只有同震線讀地震儀波形會用到。`scripts/` 底下的東西請一律用 `.venv/bin/python3` 跑，不要用系統的 `python3`，shell script 也是這樣寫的。

測試的話跑 `.venv/bin/python3 -m pytest tests`，裡面主要是之前修過的統計 bug 的回歸測試。

## 資料要自己準備

原始資料很大（地磁大概 18GB），所以沒有放進 git，`G1`～`G23` 這些資料夾在 repo 裡只留各自的 `NOTES.md`。要自己準備的有兩種：

| 資料 | 放哪裡 | 去哪裡拿 |
|---|---|---|
| 地磁 1 秒資料（IAGA-2002 的 `.sec`，`.sec.gz` 或 GDMS 批次下載的 `.tgz` 也可以直接丟） | `G1/`、`G2_G3/`、…、`G23/`，檔名像 `cnu20240101dsec.sec` | CWA GDMS（`gdms.cwb.gov.tw`），要先註冊帳號 |
| 地震儀波形和儀器響應（miniSEED + SAC PoleZero） | `seismometer/<GXX_MMDD>/`，例如 `seismometer/G10_0403/` | 也是 CWA GDMS，我是自己註冊帳號下載的 |

地震儀資料抓下來之後，如果資料夾名稱不是照 `G<NN>_<MMDD>` 的規則，要去 `seismometer_comparison.py` 的 `SEISMIC_DATA_DIRS` 手動加一筆對照，像 G9 的資料夾是 `G9_0918`，但 mseed 檔卻是 `G09_0918_w.mseed`。

剩下兩種腳本會自己上網抓，有網路就好：

- 空間天氣指數 Kp、Dst：`fetch_space_weather.py` 會抓到 `data/external/<group>/`，抓不完整下次會自己重抓，想強制重抓就把那組的 `storm_days.csv` 刪掉。
- 擴充地震目錄：`fetch_earthquake_catalog.py` 讀的是 repo 上一層的兩個 CWA GDMS 匯出檔，`../GDMScatalog_2009-2024.txt` 和 `../GDMScatalog.txt`，這兩個要自己從 GDMS 匯出，repo 裡沒有。時間超出這兩個檔的組別會直接報錯，真的要用 USGS 補才加 `--allow-usgs`。

## 資料夾

```
geomag_precursor/
├── README.md
├── NOTES.md          開發筆記
├── requirements.txt
├── scripts/          所有分析程式
├── tests/
├── docs/             統計方法說明、當初挑事件和抓資料範圍的紀錄
├── G1/ … G23/        原始地磁資料（git 只留 NOTES.md）
├── seismometer/      地震儀資料（不進 git）
├── data/             中間結果和輸出，不進 git，可以從原始資料重跑出來
└── output/           HTML 報告，這個有進 git
```

## 怎麼跑

### 日尺度前兆線

單跑一組：

```bash
cd scripts
./run_pipeline.sh --group G10
```

它會依序做解析、時區檢查、抓磁暴日、算異常指數、ULF 極化，最後跑 `verify_pipeline.py` 檢查。加 `--full-report` 會多產生一份敘事型的 HTML 報告，不過那段敘事是我專門寫給 G10 的，其他組跑了也沒什麼意義。

全部 24 組一起跑，然後做跨組彙整：

```bash
./run_all_groups.sh
./cross_group_analysis.py
```

某一組沒過檢查不會卡住整批，最後會整理在 `data/interim/all_groups_run_summary.json`。

全部跑完之後，可以再跑正式的跨組驗證（抓目錄、洗牌檢定、疊加曆元分析、規則回測，最後出報告）：

```bash
./run_validation_pipeline.sh
```

### 同震線

這條沒有 shell script 串起來，要照順序自己跑：

```bash
cd scripts
../.venv/bin/python3 coseismic_step_analysis.py --all       # 每起地震找 step/spike
../.venv/bin/python3 coseismic_stacking_analysis.py --all   # 跨事件疊加，自己重算，不吃上一步的結果
../.venv/bin/python3 seismometer_comparison.py --all        # 跟地震儀比對
../.venv/bin/python3 coseismic_joint_analysis.py --all      # 聯合檢定，要讀上一步的 comparison_summary.csv
```

第三步一定要先跑成功，第四步才跑得動。另外還有兩支後來加的：`coseismic_dose_response.py` 看地磁異常跟地動加速度的關係，`coseismic_onset_moveout.py` 比近站和遠站的起始時間，後者不需要地震儀資料。

這幾支都有 `--self-test`，正式跑之前也會自己先跑一次，沒過就停。在新環境第一次跑的話，建議先手動跑一次確定環境沒問題：

```bash
../.venv/bin/python3 coseismic_step_analysis.py --self-test
```

## 常調的參數

參數都是寫在腳本裡的大寫常數，改完重跑那支腳本就好，沒有設定檔。比較常動到的是這幾個：

| 參數 | 在哪裡 | 現在的值 | 做什麼 |
|---|---|---|---|
| `N_NEAR_STATIONS` / `N_FAR_STATIONS` | `common.py` | 3 / 2 | 近站、遠站各取幾站 |
| `TRAILING_WINDOW_DAYS` | `compute_indices.py` | 28 | 滾動基線的天數 |
| `CANDIDATE_Z_THRESHOLD` | `compute_indices.py` | 2.5 | z 超過多少算異常日 |
| `NIGHT_HOURS_UTC` | `compute_indices.py`、`ulf_analysis.py` | 17–19 UTC | 只用台灣半夜 1 點到 4 點最安靜的時段 |
| `KP_STORM_THRESHOLD` / `DST_STORM_THRESHOLD` | `fetch_space_weather.py` | 5 / −30 nT | 磁暴日的門檻 |
| `ULF_GROUPS` | 5 支腳本各一份 | 14 組 | 哪些組的向量站夠多，可以跑跨組驗證 |
| `SCAN_HALF_SEC` | `coseismic_step_analysis.py` | 180 秒 | 同震搜尋的半窗 |
| `SEARCH_HALF_SEC` | `seismometer_comparison.py` | 180 秒 | 同上，比對地震儀那支用的 |
| `SEED` | 幾支有亂數的腳本 | 20260805 | 固定亂數種子，結果才重現得出來 |

有幾個地方要小心。`ULF_GROUPS` 是在 `fetch_earthquake_catalog.py`、`surrogate_test.py`、`superposed_epoch_analysis.py`、`backtest_rule.py`、`prepare_validation_report_data.py` 各寫一份，改的話五個都要改。`SCAN_HALF_SEC` 和 `SEARCH_HALF_SEC` 雖然是同一個概念，但其實是兩個獨立的常數，也是要兩邊一起改。然後我有試過把搜尋窗放寬到 300 秒以上，結果窗口越寬，光靠雜訊就越容易找到更大的值，好幾起事件的結果反而變不穩，所以最後還是維持 180 秒。

## 跑完的東西在哪

日尺度線每一組的結果在 `data/interim/<group>/`，像 `daily_features.csv`、`local_anomaly_index.csv`、`candidate_windows.json`、`verification_report.json`。同震線每支腳本各有一個 `data/interim/<腳本名稱>/` 資料夾，裡面有總表 `*_summary.csv` 和每起事件的 JSON。跨組的彙整是 `data/interim/cross_group_summary.{json,md}`。

`output/` 裡有三份 HTML：G10 的單事件報告、跨組驗證報告，還有一份可以點的流程圖 `pipeline_flowchart.html`，每一步的輸入輸出和用到的公式都在裡面。同震線的結果我是另外整理成報告，沒有放在這個 repo。

## 會踩到的雷

- 2016 年以前的資料（G14～G18、G21、G22）都是只有總磁場 F 的純量站，H、Z 和 ULF 極化完全跑不了，只有 G4 之後跟 G19、G20、G23、G24 這些向量站齊全的組才行，所以 `ULF_GROUPS` 才只有 14 組。
- G14 和 G21 永遠不會有地震儀資料，因為地震儀資料只回溯到 2012 年，這兩組的地震都比那更早。目前 137 起裡面有 111 起可以跟地震儀比對。
- 放在同一個資料夾的組其實用的是同一批原始資料，像 G2/G3、G6/G7/G8、G23/G24，所以跨組檢定裡它們不算完全獨立的樣本。
- 日尺度的候選日對測站選哪幾站很敏感，G2 和 G3 資料一樣，只是近站不同，候選日數量就差很多。
- `verify_pipeline.py` 有幾組的檢查一直沒過，大多是磁暴扣不乾淨或震前的乾淨基線天數不夠，哪幾組、為什麼都記在 `NOTES.md`。
- 卑南站 `ttn` 從 2024 年底就沒資料了，G12、G13、G20 都沒有這一站。

## 其他文件

- `NOTES.md`：開發筆記，原始資料格式、測站代碼的歷史、各組資料範圍、同震線方法怎麼一路改過來的都在這裡。
- `docs/statistical_methods.md`：每個統計方法的公式和對應的程式碼位置。
- `docs/13_groups_fetch_ranges.md`、`docs/candidate_groups_G14_G20.md`、`docs/candidate_groups_from_GDMScatalog.md`：當初怎麼挑地震、抓哪段資料的紀錄。
