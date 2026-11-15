"""反向測試：驗證流程必須「抓得到失敗」。

過去 run_all_groups 的 verify 解析不管結果都標成 ok，全部 ok 的報告其實沒有意義。
這裡確認三層都會把失敗往外傳：
  1. 單一檢查：餵壞資料，pass 是 False
  2. verify_pipeline.main()：有任何 pass=False 就 exit 1；只有 None（無法判定）不算失敗
  3. run_pipeline.sh：verify 失敗時整支以非 0 結束、不往下建報告（run_all_groups.sh 靠這個結束碼判斷 ok / verify_failed）
"""
import json
import os
import stat
import subprocess
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

import verify_pipeline as vp

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


@pytest.fixture(autouse=True)
def _fresh_checks(monkeypatch):
    monkeypatch.setattr(vp, "checks", {})


def _cfg(tmp_path):
    interim = tmp_path / "interim"
    interim.mkdir()
    return types.SimpleNamespace(interim_dir=interim)


def _write_quiet_day_inputs(cfg, n_spikes):
    (cfg.interim_dir / "timezone_check.json").write_text(
        json.dumps({"station": "LNU", "quiet_days_used": ["20240101", "20240102"]}))
    pd.DataFrame({
        "station": ["LNU", "LNU"], "date": ["20240101", "20240102"],
        "n_spikes": [0, n_spikes], "max_abs_jump": [3.0, 250.0],
    }).to_csv(cfg.interim_dir / "daily_features.csv", index=False)


# --- 1. 單一檢查 -----------------------------------------------------

def test_quiet_day_smoothness_fails_on_spike(tmp_path):
    cfg = _cfg(tmp_path)
    _write_quiet_day_inputs(cfg, n_spikes=3)
    vp.check_quiet_day_smoothness(cfg)
    assert vp.checks["quiet_day_smoothness"]["pass"] is False


def test_quiet_day_smoothness_passes_on_clean(tmp_path):
    cfg = _cfg(tmp_path)
    _write_quiet_day_inputs(cfg, n_spikes=0)
    vp.check_quiet_day_smoothness(cfg)
    assert vp.checks["quiet_day_smoothness"]["pass"] is True


# --- 2. main() 的結束碼 ---------------------------------------------

CHECK_FUNCS = [
    "check_raw_vs_parsed_spotcheck", "check_known_outage_ground_truth", "check_timezone_conclusion",
    "check_storm_cancellation", "check_candidates_not_in_outage_windows", "check_quiet_day_smoothness",
    "check_no_spikes_adjacent_to_missing_data", "check_baseline_window_excludes_storms",
]


def _run_main(tmp_path, monkeypatch, results):
    """把每個檢查換成直接寫入指定的 pass 值，跑 main()，回傳 (結束碼, 報告)。"""
    cfg = _cfg(tmp_path)
    monkeypatch.setattr(vp, "load_group_config", lambda g: cfg)
    for name, value in zip(CHECK_FUNCS, results):
        monkeypatch.setattr(vp, name, lambda c, n=name, v=value: vp.checks.__setitem__(n, {"pass": v}))
    monkeypatch.setattr(sys, "argv", ["verify_pipeline.py", "--group", "GX"])
    code = 0
    try:
        vp.main()
    except SystemExit as e:
        code = e.code
    report = json.loads((cfg.interim_dir / "verification_report.json").read_text())
    return code, report


def test_main_exits_1_when_any_check_fails(tmp_path, monkeypatch, capsys):
    results = [True] * len(CHECK_FUNCS)
    results[4] = False
    code, report = _run_main(tmp_path, monkeypatch, results)
    assert code == 1
    assert report["summary"]["fail"] == 1


def test_main_exits_0_when_all_pass(tmp_path, monkeypatch, capsys):
    code, report = _run_main(tmp_path, monkeypatch, [True] * len(CHECK_FUNCS))
    assert code == 0
    assert report["summary"] == {"pass": len(CHECK_FUNCS), "fail": 0, "inconclusive": 0}


def test_inconclusive_is_reported_not_failed(tmp_path, monkeypatch, capsys):
    """None 不擋流程，但要出現在 summary 裡，不能被算成 pass。"""
    results = [True] * len(CHECK_FUNCS)
    results[0] = None
    code, report = _run_main(tmp_path, monkeypatch, results)
    assert code == 0
    assert report["summary"]["inconclusive"] == 1
    assert report["summary"]["pass"] == len(CHECK_FUNCS) - 1


# --- 3. run_pipeline.sh 傳遞結束碼 ----------------------------------

def _fake_python(tmp_path, verify_rc):
    """假的 PYTHON：其他步驟都成功，只有 verify_pipeline.py 回傳指定結束碼；記錄被呼叫的腳本。"""
    log = tmp_path / "calls.log"
    fake = tmp_path / "fakepython"
    fake.write_text(
        "#!/usr/bin/env bash\n"
        f'echo "$1" >> "{log}"\n'
        f'if [ "$1" = verify_pipeline.py ]; then exit {verify_rc}; fi\n'
        "exit 0\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    return fake, log


@pytest.mark.parametrize("verify_rc, expect_ok", [(1, False), (0, True)])
def test_run_pipeline_propagates_verify_result(tmp_path, verify_rc, expect_ok):
    fake, log = _fake_python(tmp_path, verify_rc)
    r = subprocess.run(
        ["bash", str(SCRIPTS / "run_pipeline.sh"), "--group", "G10", "--full-report"],
        env={**os.environ, "PYTHON": str(fake)}, capture_output=True, text=True, timeout=60)
    calls = log.read_text().split()
    assert (r.returncode == 0) is expect_ok, r.stderr
    assert "verify_pipeline.py" in calls
    # verify 失敗時，報告建置步驟一個都不能跑
    built = {"prepare_report_data.py", "build_artifact.py"} & set(calls)
    assert built == (set() if not expect_ok else {"prepare_report_data.py", "build_artifact.py"})
