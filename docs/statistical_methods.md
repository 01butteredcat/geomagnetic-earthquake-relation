# 地磁資料統計分析方法總覽

本文彙整 `scripts/` pipeline 中實際用在地磁（及對照用地震）資料上的統計/訊號處理方法，依分析階段分組列出：**用途**、**數學式**、**優點**、**缺點/限制**、**程式碼位置**。每個 file:line 對應目前程式碼，可直接查證。

整體設計理念（見 `scripts/stat_utils.py` 開頭註解）：地磁日序列樣本數少（單一群組僅 ~93 天基線）、非常態、具自相關，因此本專案**捨棄傳統母數統計（t 檢定、Pearson 相關、PCA 等）**，改以「穩健統計量（中位數/MAD、Theil–Sen）＋重抽樣式顯著性檢定（bootstrap／替代資料／置換檢定）」為主軸，見文末對照小節。

---

## A. 資料前處理

### A1. 去尖峰濾波（Despiking，一階差分跳點偵測）
- **用途**：在計算日摘要統計前，先剔除 1Hz X/Y/Z/F 中的儀器/傳輸尖峰雜訊。
- **程式碼**：`scripts/build_daily_features.py:37-47`
- **數學式**：

  設原始序列為 $x_t$，一階差分 $d_t = |x_t - x_{t-1}|$。若 $d_t >$ `SPIKE_THRESHOLD_NT`（300 nT），則將 $x_t$ 標記為缺值（NaN）：

  $$
  x_t' = \begin{cases} \text{NaN} & d_t > 300\text{ nT} \\ x_t & \text{otherwise} \end{cases}
  $$

- **優點**：計算極輕量（O(n)）、對單點尖峰異常非常敏感、不需假設資料分布。
- **缺點**：固定絕對閾值（300 nT）沒有隨測站/季節雜訊水準自適應；連續多點的「階梯」異常不會被偵測（差分本身可能小於閾值）；可能誤刪真實的快速磁暴變化。

### A2. 去趨勢（1 小時置中滾動平均）
- **用途**：在做 ULF 頻帶濾波（見 C6）或 coseismic 階躍偵測（見 D15）前，先移除長週期背景趨勢（如日變 Sq 曲線），避免其能量洩漏進目標頻帶。
- **程式碼**：`scripts/ulf_analysis.py:52-54`；`scripts/coseismic_step_analysis.py:219-223`（獨立複製的同一公式）
- **數學式**：

  $$
  \bar{x}_t = \frac{1}{3601}\sum_{i=t-1800}^{t+1800} x_i, \qquad x_t^{\text{detrend}} = x_t - \bar{x}_t
  $$

  （置中、視窗寬 3601 秒 ≈ 1 小時）
- **優點**：實作簡單、線性相位（不會造成時間位移）、足以壓制日變等長週期成分。
- **缺點**：視窗邊界（資料頭尾）平滑效果較弱（`min_periods=1`）；對比真正的高通濾波器，滾動平均的頻率響應較差（旁瓣大），並非嚴謹的頻域設計。

---

## B. 穩健標準化與異常偵測（precursor screening 主幹）

### B1. 中位數/MAD 穩健 z-score
- **用途**：本專案幾乎所有異常判定的共同基礎統計量——將夜間平均 H/Z/F 值、ULF near−far 差值序列等標準化，用來偵測偏離「正常」範圍的日子。
- **程式碼**：`scripts/stat_utils.py:27-37`（全序列版，各腳本共用）；`scripts/compute_indices.py:74-89`（21 天滾動窗版，`TRAILING_WINDOW_DAYS=21`、`MIN_CLEAN_POINTS=5`）
- **數學式**：

  $$
  \text{med} = \operatorname{median}(x), \qquad
  \text{MAD} = 1.4826 \times \operatorname{median}\bigl(|x_i - \text{med}|\bigr)
  $$
  $$
  z_i = \frac{x_i - \text{med}}{\text{MAD}}
  $$

  （常數 1.4826 使 MAD 在真常態分布下與標準差同尺度）
- **優點**：對離群值/厚尾分布穩健（中位數與 MAD 的崩潰點遠高於平均值/標準差）；不需假設常態性；計算快速、易於逐日滾動更新。
- **缺點**：樣本數過少（如滾動窗內乾淨天數 < 5）時 MAD 估計本身不穩定；MAD 接近 0（序列非常平坦）時 z 值會爆炸或無定義；仍是「單變量、逐點」統計量，未考慮序列的自相關結構（因此後續才需要替代資料法補強顯著性推論，見 D）。

### B2. Theil–Sen 穩健回歸（近-遠測站共模訊號校正）
- **用途**：地磁日變化中很大一部分是全網共有的太陽風/磁暴訊號（共模雜訊），而非局部異常。此法用「遠測站指標」回歸「近測站指標」，殘差即為扣除共模訊號後的「局部異常指標」（local anomaly index）。
- **程式碼**：`scripts/compute_indices.py:92-165`，呼叫 `scipy.stats.theilslopes`
- **數學式**：

  設 $y_t=$ 近測站群 MAD-z 中位數（near_index）、$x_t=$ 遠測站群 MAD-z 中位數（far_index），僅用「乾淨日」擬合：

  $$
  \hat{\beta} = \operatorname{median}_{i<j}\left(\frac{y_j-y_i}{x_j-x_i}\right), \qquad \hat{\alpha} = \operatorname{median}(y_t - \hat{\beta} x_t)
  $$
  $$
  \text{local\_anomaly}_t = y_t - (\hat{\beta} x_t + \hat{\alpha})
  $$

  候選異常日旗標：$|\text{local\_anomaly}_t| > 2.5 \times \text{MAD}(\text{residual})$
- **優點**：Theil–Sen 斜率估計對離群值穩健（崩潰點約 29%，遠高於最小二乘法的 0%）；不需誤差常態假設；能有效扣除全網共模的磁暴/日變訊號，凸顯局部異常。
- **缺點**：假設 near/far 關係為線性、且斜率在整個視窗內恆定，若局部異常本身很大也會被線性擬合部分「吸收」而低估；`n_fit_days < 10` 時退化為零斜率（不校正），可能失真；仍只是雙變量迴歸，未考慮測站間的空間相關結構。

---

## C. 訊號處理與頻譜分析（ULF Pc3/Pc4 極化）

### C1. Butterworth 零相位帶通濾波
- **用途**：從 1Hz X/Y/Z 中萃取 Pc3（10–45 秒週期）與 Pc4（45–150 秒週期）地磁脈動頻帶，這是文獻中常見的「震前 ULF 磁擾動」候選頻段。
- **程式碼**：`scripts/ulf_analysis.py:46-49`，`scipy.signal.butter` + `filtfilt`
- **數學式**：4 階 Butterworth 帶通濾波器，正規化截止頻率為

  $$
  f_{\text{low}} = \frac{1}{150\text{s}},\quad f_{\text{high}} = \frac{1}{45\text{s}} \ (\text{Pc4}); \qquad
  f_{\text{low}} = \frac{1}{45\text{s}},\quad f_{\text{high}} = \frac{1}{10\text{s}} \ (\text{Pc3})
  $$

  再以 `filtfilt`（正向+反向濾波）消除相位延遲，等效於 8 階零相位濾波器。
- **優點**：通帶內響應平坦（Butterworth 特性）；`filtfilt` 保證零相位失真，不影響事件時間對齊；頻帶邊界對應標準 IAGA ULF 波段命名，方便與文獻比較。
- **缺點**：資料中若有缺值（gap）會破壞 `filtfilt` 的連續性假設，程式選擇直接跳過該日而非插值修補，可能造成資料選擇性缺失；4 階濾波器過渡帶仍有一定寬度，頻帶邊緣能量會互相滲漏。

### C2. Hilbert 轉換包絡線與 Z/H 極化比
- **用途**：計算每個頻帶在「本地夜間」（UTC 17–19 時）的訊號包絡能量，並取垂直分量對水平分量的比值（Z/H ratio）——地磁極化分析中常見的震前異常指標，理論基礎是地下電性結構變化會改變波的極化特性。
- **程式碼**：`scripts/ulf_analysis.py:57-87`，`scipy.signal.hilbert`
- **數學式**：

  $$
  \text{env}_c(t) = |\mathcal{H}(x_c(t))| \quad (c \in \{X,Y,Z\})，\qquad
  \text{env}_H(t) = \sqrt{\text{env}_X(t)^2+\text{env}_Y(t)^2}
  $$
  $$
  \text{RMS}_H = \sqrt{\frac{1}{N}\sum_{t\in\text{night}}\text{env}_H(t)^2}, \quad
  \text{RMS}_Z = \sqrt{\frac{1}{N}\sum_{t\in\text{night}}\text{env}_Z(t)^2}
  $$
  $$
  \text{ZH ratio} = \frac{\text{RMS}_Z}{\text{RMS}_H}
  $$

  最終候選前兆序列為近−遠測站群中位數 ZH ratio 之差：`diff_zh = median(near_zh) − median(far_zh)`（`ulf_analysis.py:153-170`）。
- **優點**：Hilbert 包絡是萃取瞬時振幅的標準方法，不需額外視窗化假設；限定於本地夜間可避開人為電磁雜訊（工業用電、交通）干擾最嚴重的時段；near−far 差分設計再次扣除共模效應。
- **缺點**：僅取夜間 3 小時樣本，樣本量偏少，單日估計變異大；Z/H ratio 對水平分量接近零時數值不穩定（雖設有 `1e-9` 保護）；物理機制（地下電性/壓磁效應）本身在文獻中仍有爭議，此為觀測型指標而非嚴格物理模型。

### C3. 短時傅立葉頻譜圖（STFT Spectrogram）
- **用途**：於事件前後 ±12 天窗口內，視覺化 X 分量在 ULF 頻段的功率隨時間變化，作為診斷/報告用圖，輔助判讀是否有頻譜異常。
- **程式碼**：`scripts/ulf_analysis.py:118-150`，`scipy.signal.spectrogram`（`nperseg=3600, noverlap=1800`）
- **數學式**：對每個長度 3600 秒的分段做 FFT 取功率譜密度 $S_{xx}(f,t)$，圖示為 $\log_{10}(S_{xx}+10^{-6})$（避免 log(0)）。
- **優點**：直觀呈現時頻演化，適合人工判讀是否有異常頻譜特徵；`filtfilt` 之外的獨立驗證管道。
- **缺點**：純視覺化/診斷用途，未產生量化統計檢定；時間解析度與頻率解析度受限於分段長度（Heisenberg 測不準原理的工程版本，1 小時分段對應約 0.28 mHz 頻率解析度）。

---

## D. 假設檢定與顯著性推論（以重抽樣為主）

### D1. Shapiro–Wilk 常態性檢定
- **用途**：檢驗 ULF near−far 差分序列是否服從常態分布，用來論證後續為何改用重抽樣式檢定而非傳統母數檢定。
- **程式碼**：`scripts/surrogate_test.py:31,68`，`scipy.stats.shapiro`
- **數學式**：

  $$
  W = \frac{\left(\sum_{i=1}^n a_i x_{(i)}\right)^2}{\sum_{i=1}^n (x_i-\bar{x})^2}
  $$

  （$x_{(i)}$ 為排序統計量，$a_i$ 為常態分位數推導出的係數）
- **優點**：小樣本下檢定力（power）優於 K-S 檢定等其他常態性檢定；廣泛被接受為標準做法。
- **缺點**：僅檢驗「是否顯著偏離常態」，不告訴你偏離的方向/程度；樣本數很小或很大時分別可能低估/過度偵測偏態；此處只作為**動機說明**，並非分析結論本身。

### D2. 傅立葉相位隨機化替代資料法（FFT Phase Randomization Surrogate）
- **用途**：產生與觀測序列**振幅頻譜（自相關結構）相同**、但與真實地震時間完全無關的「虛擬序列」，作為虛無假設（null hypothesis）分布的一種來源。
- **程式碼**：`scripts/stat_utils.py:81-101`
- **數學式**：

  $$
  X(f) = \mathcal{F}\{x_t-\bar{x}\}, \quad A(f)=|X(f)|
  $$
  $$
  X'(f) = A(f)\,e^{i\theta(f)}, \quad \theta(f)\sim \text{Uniform}(0,2\pi)\ \text{（獨立同分布，保留 Hermitian 對稱與 DC/Nyquist 相位為 0）}
  $$
  $$
  x_t' = \mathcal{F}^{-1}\{X'(f)\} + \bar{x}
  $$
- **優點**：屬非線性時間序列分析的標準「替代資料」方法，嚴格保留功率譜（等價於自相關函數），比單純打亂順序（shuffle）更貼近真實序列的統計性質。
- **缺點**：僅保留二階統計量（功率譜），若真實序列有非線性結構（如高階矩、局部突發特徵）不會被替代資料重現，可能低估某些非線性訊號的顯著性；相位完全隨機化也可能破壞真實訊號中「合理但非隨機」的局部相位耦合。

### D3. 循環區塊拔靴法（Circular Block Bootstrap Surrogate）
- **用途**：另一種虛無假設替代序列生成法，透過重組「連續區塊」而非單點打亂，保留序列的短期自相關（而非只保留頻譜）。
- **程式碼**：`scripts/stat_utils.py:40-78`
- **數學式**：

  區塊長度 $\ell$ 由自相關函數（ACF）首次降到 $1/e$ 以下的落後期（lag）決定（下限 3、上限 $n/5$）：

  $$
  \rho(k) = \frac{\sum_t (x_t-\bar x)(x_{t+k}-\bar x)}{n\cdot\operatorname{Var}(x)}, \qquad
  \ell = \min\{k : |\rho(k)| < e^{-1}\}
  $$

  替代序列由隨機起點的長度-$\ell$ 連續區塊（首尾循環相接）拼接而成，直到補滿原長度。
- **優點**：保留區塊內的局部（短期）自相關結構，較替代資料法更貼近「重排時間但維持局部型態」的直覺；區塊長度用資料驅動的方式估計，非任意選定。
- **缺點**：區塊長度估計為簡化的經驗法則（非如 Politis–White 的最適區塊長度理論），區塊邊界處會人為引入不連續；區塊間仍視為獨立可交換，若序列存在更長程的相依結構（如季節性）則無法保留。

### D4. 重抽樣式 p-value（觀測極值 vs. 虛無分布）
- **用途**：以 D2/D3 產生的 2000 組替代序列各自算出「該序列最負的 MAD z-score」，形成虛無分布，再看觀測值/固定閾值落在此分布的哪個百分位，得到經驗 p-value。
- **程式碼**：`scripts/surrogate_test.py:90-101`
- **數學式**：

  $$
  p = \frac{1}{N}\sum_{k=1}^{N} \mathbb{1}\bigl[z^{(k)}_{\min} \le z_{\text{obs}}\bigr]
  $$

  （$N=2000$，同時對「觀測極值」與「固定規則閾值 −4.1」各算一次）
- **優點**：不依賴任何母數分布假設，p-value 直接由資料本身的重抽樣分布給出，對非常態、自相關的地磁序列特別合適。
- **缺點**：p-value 精度受 $N$ 限制（$N=2000$ 時最小可解析 p-value 約 0.0005）；「取序列中最負一日」本身是一種事後挑選最極端值的做法（multiple comparison / look-elsewhere effect），此重抽樣法有部分緩解但未完全校正（跨群組、跨頻帶多重比較的校正在別處以其他方式處理，見 D5）。

### D5. 二項式檢定（Binomial Test，跨群組再現性）
- **用途**：檢驗「候選異常日落在震前視窗內」的命中率，是否顯著高於各群組自身的滑動窗基準率——用來評估整體規則跨 23 個獨立地震序列的再現性，而非單一事件的巧合。
- **程式碼**：`scripts/cross_group_analysis.py:92-201`，`scipy.stats.binomtest`
- **數學式**：設 $n$ 為受測群組數、$k$ 為命中群組數、$p_0$ 為平均基準命中率，單尾檢定：

  $$
  P(K\ge k \mid n, p_0) = \sum_{j=k}^{n}\binom{n}{j}p_0^j(1-p_0)^{n-j}
  $$
- **優點**：精確檢定（非常態近似），小樣本（群組數少，如本專案 23 群）下仍有效；直接回答「這不只是單一事件湊巧」的問題，是跨事件驗證的核心統計工具。
- **缺點**：把每個群組視為一次獨立的白努利試驗，忽略了群組間可能的地理/時間相依性（如同一斷層帶的連續地震）；基準率 $p_0$ 本身依賴滑動窗方法估計，若基準率估計有偏誤會直接傳導到檢定結果。

### D6. 標籤置換檢定（Label-Permutation Test）
- **用途**：在 coseismic 分析中，比較「與震動同步」事件組與「領先/持續於震動」事件組之間的堆疊統計量差異（$\Delta_{\text{peak}}$、$\Delta_{\text{tail}}$），檢驗此差異是否顯著大於隨機分組所致。
- **程式碼**：`scripts/coseismic_joint_analysis.py`（`N_PERM=2000`）
- **數學式**：

  $$
  p = \frac{1+\#\{k : |\Delta^{(k)}_{\text{perm}}| \ge |\Delta_{\text{obs}}|\}}{N+1}
  $$

  （隨機打亂事件的兩組標籤 $N=2000$ 次，重新計算組間差異）
- **優點**：置換檢定不需任何分布假設，天生適合處理小樣本、非常態的分組比較；分子/分母各加 1 的寫法（"add-one" 校正）避免 p-value 恰好為 0，是標準穩健做法。
- **缺點**：事件總數少（本專案僅 23 組），置換檢定的解析度（可達到的最小 p-value）受限；分組本身（noise_arm vs. signal_arm）依賴前一步驟的人工分類規則，並非統計上獨立產生。

### D7. 地震對照組重抽樣 p-value（Coseismic 隨機參考時刻）
- **用途**：在地震發生時刻附近偵測「階躍/尖峰」統計量後，另外抽取同一群組中大量與真實地震保持一定緩衝距離的隨機參考時刻，計算相同統計量在「無地震」情境下有多極端，藉此得到經驗 p-value。
- **程式碼**：`scripts/coseismic_step_analysis.py:292-338`（`_extremum_in_window`、`_draw_null_centers`），`N_NULL=2000`
- **數學式**：

  $$
  p = \frac{1+\#\{r : |\text{null}_r| \ge |\text{obs}|\}}{N_{\text{null}}+1}
  $$

  隨機參考時刻需與所有真實地震事件保持 `exclude_buffer_sec` 以上距離，避免虛無分布被真實異常污染。
- **優點**：虛無分布直接來自同一測站/同一時期的真實雜訊特性，不需假設雜訊分布形式；`_effective_half_sec` 機制確保搜尋窗不會跨越到鄰近的另一起真實地震，避免污染。
- **缺點**：搜尋窗越寬，「窗內最大值」本身就越容易偶然偏大（look-elsewhere effect）；本專案在秒級堆疊分析（見 E1）中特別改用「虛無序列自身峰值分布」而非逐點虛無帶來緩解此問題，但單事件層級的 p-value 仍需謹慎解讀多重比較的影響。

---

## E. 疊加時間分析與規則回測

### E1. 疊加時間分析（Superposed Epoch Analysis, SEA）+ Bootstrap 信賴區間 + Null Band
- **用途**：將多個獨立地震事件的異常序列，依「距地震發生日/秒的相對時間（lag）」對齊堆疊，檢驗是否存在跨事件一致的異常型態（而非單一事件的偶然現象）。日尺度版本用於震前 ULF 差分序列；秒尺度版本（`coseismic_stacking_analysis.py`）用於 coseismic 階躍/尖峰統計量。
- **程式碼**：`scripts/superposed_epoch_analysis.py:89-169`；`scripts/coseismic_stacking_analysis.py:157-407`
- **數學式**：設 $M$ 為 $n_{\text{events}}\times n_{\text{lags}}$ 矩陣，各列為單一事件對齊後的序列：

  $$
  \bar{S}(\text{lag}) = \frac{1}{n}\sum_{i=1}^n M_{i,\text{lag}}
  $$

  Bootstrap 90% CI（對事件重抽樣 $B$ 次，$B=2000$）：

  $$
  \bar{S}^{(b)}(\text{lag}) = \frac{1}{n}\sum_{i\in I_b} M_{i,\text{lag}}, \quad I_b \sim \text{均勻放回抽樣}
  $$
  $$
  \text{CI}_{90\%} = \bigl[P_5(\bar{S}^{(1..B)}),\ P_{95}(\bar{S}^{(1..B)})\bigr]
  $$

  Null band：以同群組但earthquake-unrelated 的隨機參考日期重複整個堆疊流程 $R$ 次（$R=1000$），取其 5/50/95 百分位作為「純巧合下堆疊結果應落在的範圍」。
- **優點**：直接檢驗「跨事件一致性」，是區分「單一事件的雜訊巧合」與「真正物理前兆」最有力的證據型態之一；Bootstrap CI 與 Null band 皆為非母數方法，適合小樣本、非常態資料；秒尺度版本额外用「虛無序列自身峰值分布」而非逐點百分位判斷顯著性，避免了 look-elsewhere 問題（見模組內文件字串說明，`coseismic_stacking_analysis.py:366-376`）。
- **缺點**：事件數仍偏少（23 群/49 事件），bootstrap 對總體變異的估計在小樣本下可能偏窄；不同事件的資料品質/測站覆蓋不一致，堆疊時以 `nanmean`/`nanmedian` 處理缺值，可能讓「有效樣本數」隨 lag 而變動，邊緣 lag 的統計力較弱。

### E2. 規則回測（Precision / Recall / False-Alarm Rate）與滑動窗基準率
- **用途**：把「z ≤ −4.1」這條固定規則當作實際的地震前兆警報規則，對照完整地震目錄回測其實務表現：抓到的天數中有多少真的在觸發窗、真正抓到的事件比例多少、非事件期間誤報率多高。
- **程式碼**：`scripts/backtest_rule.py:69-136`；`sliding_baseline_rate`（`scripts/cross_group_analysis.py`）
- **數學式**：

  $$
  \text{Precision} = \frac{\text{hit\_days}}{\text{flagged\_days}}, \qquad
  \text{Recall} = \frac{\text{events\_with\_hit}}{n_{\text{events}}}, \qquad
  \text{FAR} = \frac{\text{false\_alarm\_days}}{\text{non\_precursor\_days}}
  $$
- **優點**：直接以「作為警報系統的實際效用」評估規則，比單純的統計顯著性更貼近應用場景；同時報告 precision/recall/FAR 三者，避免只看單一指標造成誤導（例如極寬鬆的規則 recall 高但 precision/FAR 很差）。
- **缺點**：固定閾值（−4.1）本身是從單一案例（G10 2024-03-30）事後選定，用同一批資料回測有資料窺探（data snooping）之嫌，需要獨立樣本外驗證才能真正評估規則的泛化能力；小樣本下 precision/recall 的點估計變異度高，文件中應搭配信賴區間或至少樣本數一併報告。

---

## F. 地震對照驗證（Seismometer Cross-check）

### F1. STA/LTA 地震觸發偵測（Classic Short-Term/Long-Term Average Ratio）
- **用途**：地震學標準的震動起訖偵測法，用於獨立標記地震波實際到達測站的時間窗，以便和地磁異常的時間點做比對，判斷地磁訊號是否只是「儀器被震動干擾」而非真正的磁場變化。
- **程式碼**：`scripts/seismometer_comparison.py:358-380`，`obspy.signal.trigger.classic_sta_lta` / `trigger_onset`
- **數學式**：

  $$
  \text{STA}(t) = \frac{1}{n_s}\sum_{i=t-n_s+1}^{t} |x_i|, \qquad
  \text{LTA}(t) = \frac{1}{n_l}\sum_{i=t-n_l+1}^{t} |x_i| \quad (n_s \ll n_l)
  $$
  $$
  R(t) = \frac{\text{STA}(t)}{\text{LTA}(t)}
  $$

  當 $R(t)$ 超過觸發閾值即判定震動開始，低於解除閾值則判定結束。
- **優點**：地震學界數十年驗證過的標準方法，對突發性強震動極為敏感、計算成本低，適合即時/大量事件處理。
- **缺點**：需要事件前有足夠長的「安靜期」估計 LTA 基準，對本專案中部分**短時觸發式**加速度計紀錄（`triggered_short_trace`，記錄本身就是被觸發後才開始存的）不適用，此時退回 F2 的方法。

### F2. 包絡閾值震動窗偵測
- **用途**：作為 F1 的替代方案，用於沒有足夠事件前基準期的短觸發式紀錄，判斷震動的起訖時間窗。
- **程式碼**：`scripts/seismometer_comparison.py:383-411`，`obspy.signal.filter.envelope`（先做 1–20 Hz 帶通）
- **數學式**：

  $$
  \text{env}(t) = |\mathcal{H}(x(t))|, \qquad \text{threshold} = 0.1 \times \max_t \text{env}(t)
  $$

  震動窗 = envelope 超過此閾值的最早/最晚時間點。
- **優點**：不需要事件前基準期，適用於任何長度的觸發式紀錄；實作簡單、對強震動訊號穩健。
- **缺點**：閾值取「自身峰值的 10%」是相對而非絕對標準，不同事件間的震動窗定義口徑不完全一致；對訊噪比低或波形逐漸衰減不明顯的紀錄，起訖時間判定會較模糊。

---

## 附註：資料品管輔助步驟（非核心統計推論）

- **Kp/Dst 磁暴日旗標**（`scripts/fetch_space_weather.py`）：以 Kp ≥ 5.0 或 Dst ≤ −30 nT 判定磁暴日（含 1–2 天恢復期），排除在 B1/B2 的「乾淨日」基線之外；官方指數不可得時，退回內部代理指標（遠測站 H/F 全距的 MAD z-score > 3.0）。屬資料清理步驟，不直接產生前兆判定。
- **時區驗證（日變訊號診斷）**（`scripts/timezone_check.py`）：利用 Sq 日變曲線極值與高頻雜訊功率的日週期，反推確認 `.sec` 檔 TIME 欄位確實為 UTC。屬一次性資料驗證，非分析方法。

---

## 本專案「未採用」但同類文獻常見的方法，及取捨原因

| 方法 | 本專案未採用原因 |
|---|---|
| 傳統 t 檢定 / ANOVA | 假設常態、獨立同分布，與地磁序列的自相關、厚尾特性不符（見 D1 的 Shapiro–Wilk 檢定結果） |
| Pearson / Spearman 相關係數 | 本專案的「近-遠站關係」改用對離群值更穩健的 Theil–Sen 回歸（B2），而非相關係數 |
| Kolmogorov–Smirnov 檢定 | 常態性檢驗改採檢定力更高的 Shapiro–Wilk（D1）；分布比較需求由重抽樣法（D2-D4）取代 |
| PCA（主成分分析） | 本專案的多測站降維採用「近/遠測站群中位數」這種穩健、可解釋性更高的簡單聚合，而非依賴共變異矩陣（對離群值敏感）的 PCA |

整體取捨邏輯：地震前兆研究的樣本量天生受限於「M≥6 地震發生頻率」，難以用大樣本漸近理論（母數統計）獲得可靠的信賴區間，因此本專案系統性地改用「穩健統計量估計 + 電腦密集重抽樣式顯著性檢定」路線，這也是目前地球物理前兆研究文獻中，樣本量有限情境下較常見、較被接受的做法。
