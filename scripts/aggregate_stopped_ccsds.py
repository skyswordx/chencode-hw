#!/usr/bin/env python3
"""Aggregate only complete shards from the stopped CCSDS campaigns."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
BACKUP = ROOT / "output" / "35_ccsds_stopped_backup_20260831"
OUT = ROOT / "output" / "36_ccsds_stopped_aggregate_20260831"
COUNTERS = ("Bit_Errors", "Total_Bits", "Frame_Errors", "Total_Frames")
GREEN = {
    (1784, "1/2"): {0.6: 7.5000e-1, 0.8: 3.8931e-1, 1.0: 7.5529e-2, 1.2: 7.9605e-3, 1.4: 3.2503e-4, 1.6: 1.1620e-5, 1.8: 3.7500e-6},
    (1784, "1/3"): {-0.4: 9.9020e-1, -0.2: 9.0090e-1, 0.0: 6.8493e-1, 0.2: 2.9762e-1, 0.4: 4.7174e-2, 0.6: 4.4583e-3, 0.8: 9.2350e-5},
    (8920, "1/2"): {0.5: 1.0000e0, 0.6: 8.8496e-1, 0.7: 5.9524e-1, 0.8: 1.8939e-1, 0.9: 2.0309e-2, 1.0: 8.4691e-4, 1.1: 3.5650e-5, 1.2: 1.5510e-5, 1.3: 1.2280e-5, 1.4: 6.7700e-6},
    (8920, "1/3"): {0.0: 8.3333e-1, 0.1: 4.9505e-1, 0.2: 9.7752e-2, 0.3: 8.9847e-3, 0.4: 2.0755e-4, 0.5: 2.8730e-5, 0.6: 1.4360e-5},
}
POINTS = {cfg: sorted(ref) for cfg, ref in GREEN.items()}
STYLES = {(1784, "1/2"): ("o", "-"), (1784, "1/3"): ("s", "-."), (8920, "1/2"): ("^", "--"), (8920, "1/3"): ("D", ":")}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def sum_shards(root: Path) -> dict[float, dict[str, int]]:
    result: dict[float, dict[str, int]] = {}
    for path in sorted((root / "shards").glob("*.csv")):
        try:
            rows = read_rows(path)
            if len(rows) != 1:
                continue
            row = rows[0]
            snr = round(float(row["SNR_dB"]), 1)
            item = result.setdefault(snr, {name: 0 for name in COUNTERS})
            for name in COUNTERS:
                item[name] += int(row[name])
        except (KeyError, ValueError):
            continue
    return result


def read_aggregate(path: Path) -> dict[float, dict[str, int]]:
    out: dict[float, dict[str, int]] = {}
    for row in read_rows(path):
        snr = round(float(row.get("SNR_dB", row.get("Eb_N0"))), 1)
        out[snr] = {name: int(row[name]) for name in COUNTERS}
    return out


def rows_for(counter: dict[float, dict[str, int]], cfg: tuple[int, str]) -> list[dict]:
    out = []
    for snr in POINTS[cfg]:
        c = counter.get(snr)
        if not c or c["Total_Frames"] == 0:
            continue
        row = dict(SNR_dB=snr, **c)
        row["BER"] = c["Bit_Errors"] / c["Total_Bits"] if c["Total_Bits"] else 0.0
        row["FER"] = c["Frame_Errors"] / c["Total_Frames"]
        row["Reference_FER"] = GREEN[cfg][snr]
        row["log10_fer_ratio"] = math.log10(row["FER"] / row["Reference_FER"]) if row["FER"] > 0 else ""
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict], comments: list[str]) -> None:
    fields = ["SNR_dB", *COUNTERS, "BER", "FER", "Reference_FER", "log10_fer_ratio"]
    with path.open("w", encoding="utf-8", newline="") as f:
        for c in comments:
            f.write(f"# {c}\n")
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    series: dict[tuple[int, str], list[dict]] = {}
    series[(1784, "1/2")] = rows_for(read_aggregate(ROOT / "output/18_ccsds_1784_trimmed/r12/ccsds_aggregate.csv"), (1784, "1/2"))
    series[(1784, "1/3")] = rows_for(read_aggregate(ROOT / "output/18_ccsds_1784_trimmed/r13/ccsds_aggregate.csv"), (1784, "1/3"))
    # output/20 is a pre-AFF3CT-state-direction historical aggregate.  output/24
    # is a post-fix seeded continuation: its aggregate already includes
    # output/20, while its shard directory contains only additional valid
    # shards.  Keep the historical aggregate only for the stopped overlay and
    # add the continuation shards once; never treat the result as a pure
    # post-AFF3CT curve.
    r12 = read_aggregate(ROOT / "output/20_ccsds_8920_r12_local/ccsds_aggregate.csv")
    for snr, counters in sum_shards(ROOT / "output/24_ccsds_supplement_8920_r12_local").items():
        item = r12.setdefault(snr, {name: 0 for name in COUNTERS})
        for name in COUNTERS:
            item[name] += counters[name]
    series[(8920, "1/2")] = rows_for(r12, (8920, "1/2"))

    r13 = {}
    for directory in (BACKUP / "32_clean_8920_r13_mini", BACKUP / "campaign_r13_supplement_04", BACKUP / "campaign_r13_supplement_05"):
        for snr, counters in sum_shards(directory).items():
            item = r13.setdefault(snr, {name: 0 for name in COUNTERS})
            for name in COUNTERS:
                item[name] += counters[name]
    series[(8920, "1/3")] = rows_for(r13, (8920, "1/3"))

    for cfg, rows in series.items():
        k, rate = cfg
        write_csv(OUT / f"ccsds_aggregate_K{k}_R{rate.replace('/', '_')}.csv", rows,
                  [f"Complete-shard aggregation; K={k}, R={rate}", "Cost recorded for this campaign: 109.42 CNY"])

    comparison = []
    for cfg, rows in series.items():
        label = f"K={cfg[0]}, R={cfg[1]}"
        for r in rows:
            ratio = r["FER"] / r["Reference_FER"]
            comparison.append({"configuration": label, "SNR_dB": r["SNR_dB"], "measured_FER": r["FER"], "reference_FER": r["Reference_FER"], "ratio": ratio, "relative_difference_pct": (ratio - 1) * 100, "Total_Frames": r["Total_Frames"], "Frame_Errors": r["Frame_Errors"]})
    with (OUT / "ccsds_four_way_quantitative_comparison.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comparison[0])); w.writeheader(); w.writerows(comparison)

    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "axes.linewidth": .9, "xtick.direction": "in", "ytick.direction": "in"})
    fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=220, facecolor="white")
    for cfg, rows in series.items():
        marker, ls = STYLES[cfg]; label = f"K={cfg[0]}, R={cfg[1]}"
        ax.semilogy([r["SNR_dB"] for r in rows], [max(r["FER"], 1e-8) for r in rows], color="#111111", marker=marker, linestyle=ls, linewidth=2.1, markersize=6.2, markerfacecolor="#111111", label=f"Measured {label}")
        ref = GREEN[cfg]
        # Do not extend a reference curve beyond the measured SNR range.  The
        # report intentionally omits the last, expensive point for each curve.
        measured_snr = {r["SNR_dB"] for r in rows}
        ref_items = [(x, y) for x, y in ref.items() if x in measured_snr]
        ax.semilogy([x for x, _ in ref_items], [y for _, y in ref_items], color="#555555", marker=marker, linestyle=ls, linewidth=1.35, markersize=6.0, markerfacecolor="white", markeredgecolor="#111111", label=f"Theory {label}")
    ax.set_xlabel(r"$E_b/N_0$ (dB)", fontsize=13); ax.set_ylabel("Frame Error Rate (FER)")
    ax.set_title("CCSDS Turbo FER — Project Retrospective", fontsize=15); ax.set_xlim(-.5, 1.9); ax.set_ylim(1e-6, 1.1)
    ax.grid(True, which="major", color="#bdbdbd", linestyle=":", linewidth=.75); ax.grid(True, which="minor", color="#e3e3e3", linestyle=":", linewidth=.4)
    handles = []
    for cfg, (marker, ls) in STYLES.items():
        label = f"K={cfg[0]}, R={cfg[1]}"
        handles.extend([Line2D([], [], color="#111111", marker=marker, linestyle=ls, markerfacecolor="#111111", linewidth=2.1, markersize=6.2, label=f"Measured {label}"), Line2D([], [], color="#555555", marker=marker, linestyle=ls, markerfacecolor="white", markeredgecolor="#111111", linewidth=1.35, markersize=6.0, label=f"Theory {label}")])
    ax.legend(handles=handles, loc="lower left", fontsize=8.0, ncol=2, columnspacing=1.5, handlelength=2.8, frameon=True, fancybox=False, edgecolor="#777")
    fig.tight_layout(); fig.savefig(OUT / "ccsds_four_way_fer_overlay.png", bbox_inches="tight", facecolor="white"); plt.close(fig)

    summary = [
        "# 停止后 CCSDS Turbo 汇总（mixed provenance）", "",
        "统计规则：只累加 CSV 结构完整的分片；不直接采用运行中聚合表。",
        "注意：K=1784 两率直接来自 output/18（AFF3CT 状态方向修正前），"
        "K=8920 R=1/2 由 output/20 历史 aggregate 加 output/24 增量组成；"
        "只有 K=8920 R=1/3 使用停机备份中的完整分片。",
        "本目录的四路图是历史叠加图，不是全量 post-AFF3CT 曲线。",
        "本轮云电脑成本：109.42 CNY", "", "| 配置 | 点数 | 平均绝对百分比误差 |", "|---|---:|---:|",
    ]
    for cfg, rows in series.items():
        diffs = [abs(r["FER"] / r["Reference_FER"] - 1) * 100 for r in rows]
        summary.append(f"| K={cfg[0]}, R={cfg[1]} | {len(rows)} | {sum(diffs)/len(diffs):.2f}% |" if rows else f"| K={cfg[0]}, R={cfg[1]} | 0 | — |")
    (OUT / "ccsds_stopped_aggregate_summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
