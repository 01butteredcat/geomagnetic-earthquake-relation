#!/usr/bin/env python3
"""pipeline_flow.yaml + template.html + vendored KaTeX -> 單檔離線 HTML。

    .venv/bin/python3 docs/flowchart/build_flowchart.py [--out output/pipeline_flowchart.html]

只做結構檢查（有錯就拒絕建置）；完整的程式碼行號/文件對帳請跑 check_flowchart.py。
輸出是決定性的（不含時間戳，只含 YAML 的 hash），方便 git diff。
"""
from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_flowchart as chk  # noqa: E402

KATEX_DIR = HERE / "vendor" / "katex"
DEFAULT_OUT = chk.PROJECT / "output" / "pipeline_flowchart.html"


def inline_katex_css() -> str:
    css = (KATEX_DIR / "katex.min.css").read_text(encoding="utf-8")
    # 只保留 woff2 並以 base64 內嵌（現代瀏覽器都支援），移除 woff/ttf 備援來源
    css = re.sub(r',url\(fonts/[^)]*\.woff\) format\("woff"\)', "", css)
    css = re.sub(r',url\(fonts/[^)]*\.ttf\) format\("truetype"\)', "", css)

    def repl(m: re.Match) -> str:
        data = (KATEX_DIR / m.group(1)).read_bytes()
        return f"url(data:font/woff2;base64,{base64.b64encode(data).decode()}) format(\"woff2\")"

    return re.sub(r'url\((fonts/[^)]*\.woff2)\) format\("woff2"\)', repl, css)


def resolve_consts(d: dict) -> list[str]:
    missing: list[str] = []
    for who, c in chk._iter_consts(d):
        c["value"] = chk.resolve_const(c["file"], c["name"])
        if c["value"] is None:
            missing.append(f"{who}: {c['name']} ({c['file']})")
    return missing


def build(out: Path) -> int:
    raw = chk.YAML_PATH.read_bytes()
    d = chk.load()
    errs, warns = chk.validate_structure(d)
    for w in warns:
        print(f"WARN  {w}")
    if errs:
        for e in errs:
            print(f"ERROR {e}")
        print("結構檢查失敗，未建置。")
        return 1

    data = copy.deepcopy(d)
    for m in resolve_consts(data):
        print(f"WARN  找不到常數 {m}")
    data["meta"]["hash"] = hashlib.sha1(raw).hexdigest()[:8]
    data["meta"].pop("helper_scripts", None)

    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    html = (HERE / "template.html").read_text(encoding="utf-8")
    for marker, value in (
        ("/*__KATEX_CSS__*/", inline_katex_css()),
        ("/*__KATEX_JS__*/", (KATEX_DIR / "katex.min.js").read_text(encoding="utf-8")),
        ("/*__DATA__*/null", payload),
    ):
        if marker not in html:
            print(f"ERROR template.html 缺少佔位符 {marker}")
            return 1
        html = html.replace(marker, value)

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size / 1024:.0f} KB, {len(d['steps'])} steps, {len(d['methods'])} methods)")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    sys.exit(build(ap.parse_args().out))
