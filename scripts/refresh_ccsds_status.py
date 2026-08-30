#!/usr/bin/env python3
"""Refresh an auditable AFF3CT-era CCSDS status table.

Only complete shard counters produced by the corrected executable are used for
the measured columns.  Seeded/mixed aggregate.csv files and pre-AFF3CT runs are
listed as excluded metadata, never folded into the final totals.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts" / "audit_single_ccsds_dir.py"
OUT = ROOT / "output" / "29_ccsds_status"

# Numeric FER reference points used for the Green Book curve comparison.  The
# Green Book figures themselves are plots; the exact tabulated points are from
# NASA DSN 810-005 Module 207 Tables A-3/A-4, the numerical source used by the
# project validation reports.
GREEN = {
    (1784, "1/2"): {0.6: 7.5000e-1, 0.8: 3.8931e-1, 1.0: 7.5529e-2,
                    1.2: 7.9605e-3, 1.4: 3.2503e-4, 1.6: 1.1620e-5,
                    1.8: 3.7500e-6},
    (1784, "1/3"): {-0.4: 9.9020e-1, -0.2: 9.0090e-1, 0.0: 6.8493e-1,
                    0.2: 2.9762e-1, 0.4: 4.7174e-2, 0.6: 4.4583e-3,
                    0.8: 9.2350e-5},
    (8920, "1/2"): {0.5: 1.0000e0, 0.6: 8.8496e-1, 0.7: 5.9524e-1,
                    0.8: 1.8939e-1, 0.9: 2.0309e-2, 1.0: 8.4691e-4,
                    1.1: 3.5650e-5, 1.2: 1.5510e-5, 1.3: 1.2280e-5,
                    1.4: 6.7700e-6},
    (8920, "1/3"): {0.0: 8.3333e-1, 0.1: 4.9505e-1, 0.2: 9.7752e-2,
                    0.3: 8.9847e-3, 0.4: 2.0755e-4, 0.5: 2.8730e-5,
                    0.6: 1.4360e-5},
}

POINTS = {
    (1784, "1/2"): [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8],
    (1784, "1/3"): [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
    (8920, "1/2"): [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3],
    (8920, "1/3"): [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
}

COUNTERS = ("Bit_Errors", "Total_Bits", "Frame_Errors", "Total_Frames")


def key_snr(x: float) -> str:
    return f"{x:.1f}"


def audit_local(path: Path) -> dict:
    if not path.exists():
        return {"error": f"missing local directory: {path}"}
    last = None
    for _ in range(3):
        try:
            proc = subprocess.run([sys.executable, str(AUDIT), str(path)],
                                  capture_output=True, text=True, check=True)
            return json.loads(proc.stdout)
        except Exception as exc:
            last = exc
            time.sleep(0.5)
    return {"error": f"local audit failed after retries: {last}"}


def parse_json_blob(text: str) -> dict:
    # SSH may print a login banner; retain the JSON object emitted by the audit.
    start, end = text.find("{"), text.rfind("}")
    if start < 0 or end <= start:
        raise ValueError(text[-500:])
    return json.loads(text[start:end + 1])


def audit_remote(host: str, key: str, user: str, helper: str, campaign: str) -> dict:
    cmd = ["ssh.exe", "-i", key, "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
           f"{user}@{host}", f"python3 {helper} {campaign}"]
    last = None
    for _ in range(2):
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=45)
            if proc.returncode:
                raise RuntimeError(f"ssh/audit exit {proc.returncode}: {proc.stderr[-500:]}")
            return parse_json_blob(proc.stdout)
        except Exception as exc:
            last = exc
            time.sleep(1)
    return {"error": f"remote audit failed after retries: {last}"}


def add_audit(store: dict, audit: dict, label: str, config: tuple[int, str]) -> None:
    store["audits"][label] = {
        "config": f"K={config[0]}, R={config[1]}",
        "files": audit.get("files"), "valid": audit.get("valid"),
        "bad": len(audit.get("bad", [])), "matches": audit.get("matches", {}),
        "error": audit.get("error"),
    }
    if audit.get("error"):
        store.setdefault("progress", {})[label] = {
            "config": f"K={config[0]}, R={config[1]}", "total_frames": 0,
            "by_snr": {}, "error": audit["error"]}
        return
    by_snr = audit.get("by_snr", {})
    store.setdefault("progress", {})[label] = {
        "config": f"K={config[0]}, R={config[1]}",
        "total_frames": sum(int(row.get("Total_Frames", 0)) for row in by_snr.values()),
        "by_snr": {str(snr): int(row.get("Total_Frames", 0)) for snr, row in by_snr.items()},
        "error": None}
    for snr, row in by_snr.items():
        dst = store["counters"].setdefault((config, snr), {name: 0 for name in COUNTERS})
        for name in COUNTERS:
            dst[name] += int(row.get(name, 0))
        store["sources"].setdefault((config, snr), []).append(label)


def fmt(x: object, digits: int = 6) -> str:
    if x is None or x == "":
        return "—"
    if isinstance(x, float):
        if x == 0:
            return "0"
        return f"{x:.{digits}g}"
    return str(x)


def build_snapshot() -> dict:
    snap = {"generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
            "counters": {}, "sources": {}, "audits": {}, "progress": {}, "excluded": [
                "output/20_ccsds_8920_r12_local: aggregate seeded/mixed; not used",
                "output/19_ccsds_8920_r13_trimmed: pre-AFF3CT/old campaign snapshot; not used",
                "output/24_ccsds_supplement_1784_r12_mini: old build-ccsds executable; not used",
                "output/25_ccsds_supplement_1784_r12_mini and output/26_ccsds_supplement_1784_r12_mini: old executable; not used",
            ]}

    # Clean local campaign currently running, plus the independent valid
    # supplement shards for K=8920 R=1/2.
    local_specs = [
        ("local_27_clean", ROOT / "output" / "27_clean_1784_r12_low_local", (1784, "1/2")),
        ("local_24_new_shards", ROOT / "output" / "24_ccsds_supplement_8920_r12_local", (8920, "1/2")),
        ("local_30_high", ROOT / "output" / "30_clean_1784_r12_high_local", (1784, "1/2")),
    ]
    for label, path, config in local_specs:
        add_audit(snap, audit_local(path), label, config)

    # Current final mini-host campaign.  It is intentionally audited from
    # complete shards rather than its live aggregate file.
    mini = audit_remote("192.168.17.2", r"C:\Users\Lenovo\.ssh\id_ed25519",
                        "zhangzw0170", "/home/zhangzw0170/Lab/chencode-hw/audit_single_ccsds_dir.py",
                        "/home/zhangzw0170/Lab/chencode-hw/output/28_clean_8920_r12_low_mini")
    add_audit(snap, mini, "mini_28_clean", (8920, "1/2"))
    mini_next = audit_remote("192.168.17.2", r"C:\Users\Lenovo\.ssh\id_ed25519",
                             "zhangzw0170", "/home/zhangzw0170/Lab/chencode-hw/audit_single_ccsds_dir.py",
                             "/home/zhangzw0170/Lab/chencode-hw/output/30_clean_1784_r12_high_mini")
    add_audit(snap, mini_next, "mini_30_high", (1784, "1/2"))

    # Current final cloud campaign.  Its aggregate is seeded, so only shards
    # are audited and summed here.
    cloud = audit_remote("8.138.144.13", r"E:\Main\Career\CHENLi\ccsds-test.pem",
                         "ecs-user", "/home/ecs-user/ccsds-cloud-20260830/audit_single_ccsds_dir.py",
                         "/home/ecs-user/ccsds-cloud-20260830/campaign_r12_full")
    add_audit(snap, cloud, "cloud_r12_new_shards", (8920, "1/2"))

    rows = []
    for config, snrs in POINTS.items():
        k, rate = config
        for snr in snrs:
            skey = key_snr(snr)
            c = snap["counters"].get((config, skey))
            ref = GREEN[config].get(snr)
            if c and c.get("Total_Frames", 0):
                bits, frames = c["Total_Bits"], c["Total_Frames"]
                be, fe = c["Bit_Errors"], c["Frame_Errors"]
                ber, fer = be / bits if bits else None, fe / frames
                log_ratio = math.log10(fer / ref) if ref and fer > 0 else None
                status = "AFF3CT后有效（持续追加）"
                if fe == 0:
                    status += "；0个错误帧，log误差不可定"
            else:
                bits = frames = be = fe = 0
                ber = fer = log_ratio = None
                status = "暂无AFF3CT后有效分片"
            rows.append({"K": k, "rate": rate, "snr_db": snr, "status": status,
                         "bit_errors": be, "total_bits": bits, "ber": ber,
                         "frame_errors": fe, "total_frames": frames, "fer": fer,
                         "green_fer": ref, "log10_fer_ratio": log_ratio,
                         "log10_measured_fer": math.log10(fer) if fer and fer > 0 else None,
                         "log10_green_fer": math.log10(ref) if ref and ref > 0 else None,
                         "abs_log10_fer_error": abs(log_ratio) if log_ratio is not None else None,
                         "source": "; ".join(snap["sources"].get((config, skey), [])) or "—"})
    snap["rows"] = rows
    # JSON has no tuple-keyed dictionaries; expose the internal audit maps in
    # a readable form while retaining the point table as the primary output.
    snap["counters"] = {
        f"K={cfg[0]},R={cfg[1]},SNR={snr}": value
        for (cfg, snr), value in snap["counters"].items()
    }
    snap["sources"] = {
        f"K={cfg[0]},R={cfg[1]},SNR={snr}": value
        for (cfg, snr), value in snap["sources"].items()
    }
    return snap


def write_outputs(snap: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ccsds_status.json").write_text(json.dumps(snap, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = ["# CCSDS Turbo 数据有效性与 Green Book 误差状态",
             "", f"最后刷新：`{snap['generated_at']}`（本机时间）", "",
             "判定规则：只累加 AFF3CT 修正版可执行文件产生的完整 `shards/*.csv`；" 
             "旧程序、预置/混合 `ccsds_aggregate.csv` 不进入实测列。",
             "Green Book 图 7-6/7-9 是曲线图；数值参考采用项目已核对的 NASA DSN 810-005 Module 207 Tables A-3/A-4。",
             "误差定义：`log10_fer_ratio = log10(实测 FER / Green FER)`；绝对值列为其绝对对数误差。", "",
             "有效程序来源：本机 `chencode_sim_corrected.exe`；迷你主机 `build-final/chencode_sim`；云主机 `build/chencode_sim`。",
             "## 当前有效性与实测 BER/FER", "",
             "| K | R | Eb/N0 (dB) | 有效性 | Bit errors | Total bits | BER | Frame errors | Total frames | FER | Green FER | log10(FER/Green) | 来源 |",
             "|---:|:---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|"]
    for r in snap["rows"]:
        lines.append("| {K} | {rate} | {snr_db:.1f} | {status} | {bit_errors} | {total_bits} | {ber} | {frame_errors} | {total_frames} | {fer} | {green_fer} | {log10_fer_ratio} | {source} |".format(
            **{**r, "ber": fmt(r["ber"]), "fer": fmt(r["fer"]), "green_fer": fmt(r["green_fer"]),
               "log10_fer_ratio": fmt(r["log10_fer_ratio"])}))
    lines += ["", "## 主机进度（完整分片）", "",
              "目标帧数只作为当前运行的上限参考；达到 300 个错误帧后可能提前停止，因此完成度不是任务成功判据。", "",
              "| 数据源 | 配置 | 已采样帧 | 目标帧上限 | 上限完成度 | 各 Eb/N0 帧数 |", "|:---|:---|---:|---:|---:|:---|"]
    target_by_label = {
        "local_27_clean": 5 * 3_000_000,
        "mini_28_clean": 6 * 3_000_000,
        "cloud_r12_new_shards": 3 * 3_000_000,
        "local_30_high": 3_000_000,
        "mini_30_high": 3_000_000,
    }
    for label, p in snap.get("progress", {}).items():
        total = int(p.get("total_frames", 0)); target = target_by_label.get(label)
        pct = f"{100 * total / target:.2f}%" if target else "—"
        by = ", ".join(f"{k}: {v}" for k, v in sorted(p.get("by_snr", {}).items(), key=lambda kv: float(kv[0]))) or "—"
        lines.append(f"| {label} | {p['config']} | {total} | {target or '—'} | {pct} | {by} |")
    lines += ["", "## Green 数值误差摘要", "",
              "| K | R | 有 Green 参考且有实测 FER 的点 | 平均绝对 log10 误差 | 最大绝对 log10 误差 |", "|---:|:---:|---:|---:|---:|"]
    for config in POINTS:
        rs = [r for r in snap["rows"] if (r["K"], r["rate"]) == config and r["abs_log10_fer_error"] is not None]
        vals = [r["abs_log10_fer_error"] for r in rs]
        lines.append(f"| {config[0]} | {config[1]} | {len(rs)} | {fmt(sum(vals)/len(vals) if vals else None)} | {fmt(max(vals) if vals else None)} |")
    lines += ["", "## 审计与排除项", "", "| 数据源 | 配置 | 分片文件 | 完整分片 | 坏/未完成 | aggregate 是否匹配（仅供诊断） | 错误 |", "|:---|:---|---:|---:|---:|:---|:---|"]
    for label, a in snap["audits"].items():
        match = ", ".join(f"{k}:{'是' if v else '否'}" for k, v in a.get("matches", {}).items()) or "—"
        lines.append(f"| {label} | {a['config']} | {a.get('files', '—')} | {a.get('valid', '—')} | {a.get('bad', '—')} | {match} | {a.get('error') or '—'} |")
    lines += ["", "排除的历史/混合目录："] + [f"- {x}" for x in snap["excluded"]]
    (OUT / "ccsds_status.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--watch", action="store_true", help="refresh repeatedly")
    ap.add_argument("--interval", type=int, default=900, help="seconds between refreshes")
    args = ap.parse_args()
    while True:
        snap = build_snapshot()
        write_outputs(snap)
        print(f"[{snap['generated_at']}] wrote {OUT / 'ccsds_status.md'}", flush=True)
        if not args.watch:
            return
        time.sleep(max(30, args.interval))


if __name__ == "__main__":
    main()
