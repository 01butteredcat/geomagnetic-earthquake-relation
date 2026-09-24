#!/usr/bin/env python3
"""檢查 pipeline_flow.yaml 是否仍與專案對得上。

    .venv/bin/python3 docs/flowchart/check_flowchart.py

檢查項目（有 ERROR 時結束碼為 1，只有 WARN 時為 0）：
  1. 結構：id 唯一、after/methods 參照存在、無環、必要欄位齊全。
  2. 程式碼參照：檔案存在、行號範圍沒超出檔案、symbol 仍在檔內；
     symbol 的定義行離標註範圍超過 10 行時警告「行號可能漂移」。
  3. consts：常數仍可在該檔案頂層找到（build 會從這裡讀值）。
  4. 與 docs/statistical_methods.md 對帳：md 內的每個 `### X#.` 方法都要在 YAML 有對應 id，反之亦然。
  5. steps 的 scripts 都存在；scripts/ 內未出現在任何步驟的 .py 給警告（提醒補圖）。
  6. （需要 node）用 vendored KaTeX 試算每一條公式與行內 $...$，語法錯誤視為 ERROR。
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent          # geomag_precursor/
YAML_PATH = HERE / "pipeline_flow.yaml"
KATEX_JS = HERE / "vendor" / "katex" / "katex.min.js"
DRIFT_TOLERANCE = 10


def load(path: Path | None = None) -> dict:
    return yaml.safe_load((path or YAML_PATH).read_text(encoding="utf-8"))


# ---------------------------------------------------------------- 結構
def validate_structure(d: dict) -> tuple[list[str], list[str]]:
    errs: list[str] = []
    warns: list[str] = []
    for key in ("meta", "lanes", "kinds", "methods", "steps"):
        if key not in d:
            errs.append(f"缺少頂層欄位 {key}")
    if errs:
        return errs, warns

    lane_ids = [l["id"] for l in d["lanes"]]
    if len(set(lane_ids)) != len(lane_ids):
        errs.append("lanes 內有重複 id")
    for mid, m in d["methods"].items():
        for f in ("name", "purpose", "group"):
            if not m.get(f):
                errs.append(f"方法 {mid} 缺少 {f}")
        if not m.get("code"):
            warns.append(f"方法 {mid} 沒有 code 位置")
        hyp = m.get("hypotheses")
        if hyp is not None and not (isinstance(hyp, dict) and hyp.get("h0") and hyp.get("h1")):
            errs.append(f"方法 {mid} 的 hypotheses 需同時有 h0 與 h1")

    ids = [s["id"] for s in d["steps"]]
    for i in {x for x in ids if ids.count(x) > 1}:
        errs.append(f"步驟 id 重複：{i}")
    idset = set(ids)
    for s in d["steps"]:
        sid = s["id"]
        for f in ("lane", "kind", "title", "summary"):
            if not s.get(f):
                errs.append(f"步驟 {sid} 缺少 {f}")
        if s.get("lane") not in lane_ids:
            errs.append(f"步驟 {sid} 的 lane「{s.get('lane')}」不存在")
        if s.get("kind") not in d["kinds"]:
            errs.append(f"步驟 {sid} 的 kind「{s.get('kind')}」不存在")
        for a in s.get("after") or []:
            if a not in idset:
                errs.append(f"步驟 {sid} 的 after 參照不存在的步驟「{a}」")
            if a == sid:
                errs.append(f"步驟 {sid} 依賴自己")
        for m in s.get("methods") or []:
            if m not in d["methods"]:
                errs.append(f"步驟 {sid} 引用了未定義的方法「{m}」")

    # 環偵測（DFS 三色）
    color: dict[str, int] = {}
    by = {s["id"]: s for s in d["steps"]}

    def dfs(n: str, path: list[str]) -> None:
        color[n] = 1
        for p in by[n].get("after") or []:
            if p not in by:
                continue
            if color.get(p) == 1:
                errs.append("依賴有環：" + " -> ".join(path + [n, p]))
            elif color.get(p) is None:
                dfs(p, path + [n])
        color[n] = 2

    for n in by:
        if color.get(n) is None:
            dfs(n, [])

    used = {m for s in d["steps"] for m in s.get("methods") or []}
    for mid in d["methods"]:
        if mid not in used:
            warns.append(f"方法 {mid} 沒有被任何步驟引用（圖上看不到它）")
    return errs, warns


# ---------------------------------------------------------------- 程式碼參照
_file_cache: dict[Path, list[str]] = {}


def _lines(path: Path) -> list[str]:
    if path not in _file_cache:
        _file_cache[path] = path.read_text(encoding="utf-8").splitlines()
    return _file_cache[path]


def _def_lines(lines: list[str], sym: str) -> list[int]:
    pat = re.compile(rf"^\s*(?:def|class)\s+{re.escape(sym)}\b|^\s*{re.escape(sym)}\s*(?::[^=]+)?=(?!=)")
    return [i + 1 for i, ln in enumerate(lines) if pat.search(ln)]


def _word_lines(lines: list[str], sym: str) -> list[int]:
    pat = re.compile(rf"\b{re.escape(sym)}\b")
    return [i + 1 for i, ln in enumerate(lines) if pat.search(ln)]


def resolve_const(file: str, name: str) -> str | None:
    """從程式碼頂層讀出 `NAME = value`，去掉行尾註解；找不到回 None。"""
    path = PROJECT / file
    if not path.exists():
        return None
    for ln in _lines(path):
        m = re.match(rf"^{re.escape(name)}\s*(?::[^=]+)?=\s*(.+)$", ln)
        if m:
            val = m.group(1)
            val = re.split(r"\s+#", val)[0].strip()
            return val
    return None


def _iter_consts(d: dict):
    for mid, m in d["methods"].items():
        for c in m.get("consts") or []:
            yield f"方法 {mid}", c
    for s in d["steps"]:
        for c in s.get("consts") or []:
            yield f"步驟 {s['id']}", c


def validate_code_refs(d: dict) -> tuple[list[str], list[str]]:
    errs: list[str] = []
    warns: list[str] = []
    for mid, m in d["methods"].items():
        for ref in m.get("code") or []:
            where = f"方法 {mid} 的 {ref['file']}:{ref['lines']}"
            path = PROJECT / ref["file"]
            if not path.exists():
                errs.append(f"{where}：檔案不存在")
                continue
            lines = _lines(path)
            try:
                a, b = (int(x) for x in str(ref["lines"]).split("-"))
            except ValueError:
                errs.append(f"{where}：lines 格式應為 'a-b'")
                continue
            if b > len(lines) or a < 1 or a > b:
                errs.append(f"{where}：行號超出檔案（共 {len(lines)} 行）")
                continue
            sym = ref.get("symbol")
            if not sym:
                continue
            words = _word_lines(lines, sym)
            if not words:
                errs.append(f"{where}：symbol「{sym}」不在該檔內")
                continue
            defs = _def_lines(lines, sym) or words
            if not any(a - DRIFT_TOLERANCE <= n <= b + DRIFT_TOLERANCE for n in defs):
                warns.append(f"{where}：symbol「{sym}」實際在第 {defs[0]} 行，標註範圍可能已漂移")
    for who, c in _iter_consts(d):
        if resolve_const(c["file"], c["name"]) is None:
            errs.append(f"{who}：找不到常數 {c['name']}（{c['file']}）")
    return errs, warns


# ---------------------------------------------------------------- 文件 / scripts 對帳
def validate_md(d: dict) -> tuple[list[str], list[str]]:
    errs: list[str] = []
    md = PROJECT / d["meta"].get("source_doc", "docs/statistical_methods.md")
    if not md.exists():
        return [f"找不到 {md}"], []
    md_ids = set(re.findall(r"^###\s+([A-Z]\d+)\.", md.read_text(encoding="utf-8"), flags=re.M))
    yaml_ids = set(d["methods"])
    for i in sorted(md_ids - yaml_ids):
        errs.append(f"{md.name} 有方法 {i}，但 YAML 沒有")
    for i in sorted(yaml_ids - md_ids):
        errs.append(f"YAML 有方法 {i}，但 {md.name} 沒有對應標題")
    return errs, []


def validate_scripts(d: dict) -> tuple[list[str], list[str]]:
    errs: list[str] = []
    warns: list[str] = []
    referenced: set[str] = set()
    for s in d["steps"]:
        for sc in s.get("scripts") or []:
            referenced.add(sc)
            if not (PROJECT / "scripts" / sc).exists():
                errs.append(f"步驟 {s['id']} 的腳本 scripts/{sc} 不存在")
    helpers = set(d["meta"].get("helper_scripts") or [])
    for p in sorted((PROJECT / "scripts").glob("*.py")):
        if p.name not in referenced and p.name not in helpers:
            warns.append(f"scripts/{p.name} 沒有出現在任何步驟（新增腳本了嗎？若只是共用模組，請加進 meta.helper_scripts）")
    return errs, warns


# ---------------------------------------------------------------- 公式語法（需 node）
def _collect_tex(d: dict) -> list[tuple[str, str, bool]]:
    out: list[tuple[str, str, bool]] = []
    inline = re.compile(r"\$([^$]+)\$")

    def scan_text(where: str, text) -> None:
        for m in inline.finditer(str(text or "")):
            out.append((where, m.group(1), False))

    for mid, m in d["methods"].items():
        for i, f in enumerate(m.get("formula") or []):
            out.append((f"方法 {mid} formula[{i}]", f, True))
        for key in ("purpose", "where"):
            scan_text(f"方法 {mid} {key}", m.get(key))
        for key in ("h0", "h1"):
            scan_text(f"方法 {mid} hypotheses.{key}", (m.get("hypotheses") or {}).get(key))
        for key in ("notes", "pros", "cons"):
            for t in m.get(key) or []:
                scan_text(f"方法 {mid} {key}", t)
    for s in d["steps"]:
        for key in ("summary", "scope"):
            scan_text(f"步驟 {s['id']} {key}", s.get(key))
        for key in ("inputs", "outputs", "notes"):
            for t in s.get(key) or []:
                scan_text(f"步驟 {s['id']} {key}", t)
    return out


_NODE_SCRIPT = r"""
const katex = require(process.argv[1]);
const items = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const bad = [];
for (const [where, src, display] of items) {
  try { katex.renderToString(src, {displayMode: display, throwOnError: true, strict: 'ignore'}); }
  catch (e) { bad.push([where, String(e.message).split('\n')[0]]); }
}
process.stdout.write(JSON.stringify(bad));
"""


def validate_tex(d: dict) -> tuple[list[str], list[str]]:
    node = shutil.which("node")
    items = _collect_tex(d)
    if not node or not KATEX_JS.exists():
        return [], ["找不到 node 或 vendored KaTeX，略過公式語法檢查（頁面仍會顯示，但錯誤公式會以紅字呈現）"]
    res = subprocess.run([node, "-e", _NODE_SCRIPT, str(KATEX_JS)], input=json.dumps(items),
                         capture_output=True, text=True, timeout=60)
    if res.returncode != 0:
        return [], [f"node 執行失敗，略過公式檢查：{res.stderr.strip()[:200]}"]
    bad = json.loads(res.stdout or "[]")
    return [f"{w}：KaTeX 語法錯誤 — {msg}" for w, msg in bad], []


def run_all(d: dict) -> tuple[list[str], list[str]]:
    errs: list[str] = []
    warns: list[str] = []
    for fn in (validate_structure, validate_code_refs, validate_md, validate_scripts, validate_tex):
        e, w = fn(d)
        errs += e
        warns += w
    return errs, warns


def main() -> int:
    d = load()
    errs, warns = run_all(d)
    for w in warns:
        print(f"WARN  {w}")
    for e in errs:
        print(f"ERROR {e}")
    n_formula = sum(len(m.get("formula") or []) for m in d["methods"].values())
    print(f"\n{len(d['steps'])} 個步驟、{len(d['methods'])} 個方法、{n_formula} 條顯示公式 — "
          f"{len(errs)} 個錯誤、{len(warns)} 個警告")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
