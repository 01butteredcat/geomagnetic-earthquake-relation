# 互動式 pipeline 流程圖

把 `docs/statistical_methods.md` 的 18 個統計方法（A1–F2，含數學式）與 `scripts/` 的實際 pipeline 合成一張可點擊的流程圖：
點任一步驟看簡介、輸入輸出、數學式（KaTeX）、優缺點、程式碼位置；點方法徽章可在圖上標出所有用到該方法的步驟。

成品：`output/pipeline_flowchart.html`（單檔、離線，雙擊即可開；亮/暗色自動）。

## 怎麼改內容

**只編輯 `pipeline_flow.yaml`**，不要碰成品 HTML。

| 想做的事 | 改哪裡 |
|---|---|
| 改某方法的說明、公式、優缺點 | `methods.<ID>`（LaTeX 寫在 `\|-` 區塊內，反斜線不需跳脫） |
| 新增 pipeline 步驟 | `steps:` 加一項，填 `lane`、`kind`、`title`、`summary`、`scripts`、`after`（資料相依）、`methods`。**不需要給座標**，位置由 `after` 自動推導 |
| 新增統計方法 | `methods:` 加一項（id 需與 `docs/statistical_methods.md` 的 `### X#.` 標題一致），再由步驟的 `methods:` 引用 |
| 寫出假設檢定的虛無／對立假說 | `methods.<ID>.hypotheses: {h0, h1}`（兩者都要有；可用行內 `$...$`）。目前 D1–D7、E1 有，`docs/statistical_methods.md` 對應各節也要同步改 |
| 顯示程式中的參數值 | `consts: [{name: N_PERM, file: scripts/xxx.py}]`，數值在 build 時自動從程式碼讀，不會過期 |

## 怎麼檢查、怎麼建置

```bash
cd geomag_precursor
.venv/bin/pip install -r requirements.txt            # 需要 PyYAML
.venv/bin/python3 docs/flowchart/check_flowchart.py  # 對帳；有 ERROR 時結束碼為 1
.venv/bin/python3 docs/flowchart/build_flowchart.py  # 產生 output/pipeline_flowchart.html
```

`check_flowchart.py` 會驗證：結構（id、相依、無環）；每個 `code` 的檔案／行號／symbol 是否仍在（行號漂移只警告）；
`consts` 是否還找得到；方法 id 與 `statistical_methods.md` 標題是否一致；`scripts/` 有沒有新腳本沒進圖；
以及（有 node 時）所有公式的 KaTeX 語法。**改了 scripts/ 之後建議跑一次。**

## 檔案

- `pipeline_flow.yaml` — 唯一內容來源
- `template.html` — 版面與互動（不含內容；佔位符 `/*__KATEX_CSS__*/`、`/*__KATEX_JS__*/`、`/*__DATA__*/null`）
- `build_flowchart.py`、`check_flowchart.py` — 建置與對帳
- `vendor/katex/` — KaTeX 0.16.11（js、css、woff2），build 時內嵌，所以成品不需網路

建置輸出是決定性的（不含時間戳，只帶 YAML 的短 hash），YAML 沒變則 HTML 不變，方便 git diff。

## 與 `statistical_methods.md` 的關係

兩者目前是**各自維護**：md 是給人閱讀的完整文字版，YAML 是流程圖的資料來源。
`check_flowchart.py` 只保證方法 id 對得上，不比對內文；若改了某方法的公式，兩邊都要改。
