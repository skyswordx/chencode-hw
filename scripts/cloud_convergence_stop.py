#!/usr/bin/env python3
"""Stop a CCSDS runner after a conservative BER/FER plateau.

The monitor never edits or deletes shards.  It sums only complete shard rows,
tracks the two-decimal scientific-notation BER and FER strings, and sends
SIGINT to the matching runner once every requested SNR has remained unchanged
for MIN_STABLE_FRAMES frames and is close to the Green FER reference.
"""
from __future__ import annotations

import csv
import json
import math
import os
import signal
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CAMPAIGN = ROOT / "campaign_r12_full"
SNRS = ("1.1", "1.2", "1.3")
GREEN = {"1.1": 3.5650e-5, "1.2": 1.5510e-5, "1.3": 1.2280e-5}
MIN_STABLE_FRAMES = 200_000
MIN_FRAMES = 1_000_000
MAX_ABS_LOG10_ERROR = 0.30
POLL_SECONDS = 60


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def snapshot() -> dict[str, dict[str, int | float | str | None]]:
    result = {snr: {"frames": 0, "bits": 0, "bit_errors": 0,
                    "frame_errors": 0} for snr in SNRS}
    for path in sorted((CAMPAIGN / "shards").glob("*.csv")):
        try:
            with path.open(encoding="utf-8-sig", newline="") as stream:
                rows = list(csv.DictReader(line for line in stream if not line.startswith("#")))
            if len(rows) != 1:
                continue
            row = rows[0]
            snr = f"{float(row['SNR_dB']):.1f}"
            if snr not in result:
                continue
            item = result[snr]
            item["bits"] += int(row["Total_Bits"])
            item["bit_errors"] += int(row["Bit_Errors"])
            item["frames"] += int(row["Total_Frames"])
            item["frame_errors"] += int(row["Frame_Errors"])
        except (OSError, ValueError, KeyError):
            # A shard being written is simply ignored until the next poll.
            continue
    for snr, item in result.items():
        bits, frames = int(item["bits"]), int(item["frames"])
        ber = int(item["bit_errors"]) / bits if bits else 0.0
        fer = int(item["frame_errors"]) / frames if frames else 0.0
        item["ber"] = ber; item["fer"] = fer
        item["ber_repr"] = f"{ber:.2e}" if bits else "—"
        item["fer_repr"] = f"{fer:.2e}" if frames else "—"
        item["abs_log10_error"] = (abs(math.log10(fer / GREEN[snr]))
                                    if fer > 0 and GREEN[snr] > 0 else None)
    return result


def runner_pid() -> int | None:
    for entry in Path("/proc").glob("[0-9]*"):
        try:
            pid = int(entry.name)
            cmdline = (entry / "cmdline").read_bytes().replace(b"\0", b" ").decode(errors="ignore")
            if "ccsds_parallel_runner.py" in cmdline and "campaign_r12_full" in cmdline and pid != os.getpid():
                return pid
        except (OSError, ValueError):
            continue
    return None


def main() -> None:
    state: dict[str, dict[str, object]] = {}
    log_path = ROOT / "convergence_stop.log"
    state_path = ROOT / "convergence_stop_state.json"
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"[{now()}] monitor started; stable={MIN_STABLE_FRAMES} frames, "
                  f"max_log10_error={MAX_ABS_LOG10_ERROR}\n")
        while True:
            data = snapshot()
            candidates = []
            for snr in SNRS:
                item = data[snr]
                current = (item["ber_repr"], item["fer_repr"])
                entry = state.setdefault(snr, {"last_repr": current, "last_change_frames": item["frames"]})
                if current != entry["last_repr"]:
                    entry["last_repr"] = current
                    entry["last_change_frames"] = item["frames"]
                stable = int(item["frames"]) - int(entry["last_change_frames"])
                ok = (int(item["frames"]) >= MIN_FRAMES and stable >= MIN_STABLE_FRAMES
                      and item["abs_log10_error"] is not None
                      and float(item["abs_log10_error"]) <= MAX_ABS_LOG10_ERROR)
                entry.update({"frames": item["frames"], "ber": item["ber_repr"],
                              "fer": item["fer_repr"], "stable_frames": stable,
                              "log10_error": item["abs_log10_error"], "eligible": ok})
                candidates.append(ok)
            state_path.write_text(json.dumps({"updated": now(), "points": state}, indent=2), encoding="utf-8")
            summary = "; ".join(f"{snr}:F={data[snr]['frames']},FER={data[snr]['fer_repr']},stable={state[snr]['stable_frames']}"
                                for snr in SNRS)
            log.write(f"[{now()}] {summary}\n"); log.flush()
            if all(candidates):
                pid = runner_pid()
                log.write(f"[{now()}] plateau condition met; runner_pid={pid}; sending SIGINT\n")
                log.flush()
                if pid:
                    os.kill(pid, signal.SIGINT)
                return
            if runner_pid() is None:
                log.write(f"[{now()}] runner not found; monitor exits without stopping anything\n")
                return
            time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    main()
