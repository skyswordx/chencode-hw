#!/usr/bin/env python3
"""Allocate the next R=1/2 points after the current campaigns finish.

This small coordinator is intentionally conservative: it waits for the
current runner's explicit ``complete:`` log line, then starts one new point on
the freed host.  It does not stop, kill, or alter an active simulation.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "ccsds_parallel_runner.py"
LOCAL_CURRENT = ROOT / "output" / "27_clean_1784_r12_low_local"
LOCAL_NEXT = ROOT / "output" / "30_clean_1784_r12_high_local"
MINI_CURRENT = "/home/zhangzw0170/Lab/chencode-hw/output/28_clean_8920_r12_low_mini"
MINI_NEXT = "/home/zhangzw0170/Lab/chencode-hw/output/30_clean_1784_r12_high_mini"
MINI_ROOT = "/home/zhangzw0170/Lab/chencode-hw"
MINI_KEY = r"C:\Users\Lenovo\.ssh\id_ed25519"
POLL = 60


def local_complete(path: Path) -> bool:
    log = path / "overnight.log"
    if not log.exists():
        return False
    try:
        return "complete:" in log.read_text(encoding="utf-8", errors="ignore").splitlines()[-1]
    except (OSError, IndexError):
        return False


def mini_complete() -> bool:
    command = ["ssh.exe", "-i", MINI_KEY, "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
               "zhangzw0170@192.168.17.2",
               f"tail -n 1 {MINI_CURRENT}/overnight.log 2>/dev/null || true"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        return "complete:" in result.stdout
    except (OSError, subprocess.SubprocessError):
        return False


def launch_local() -> None:
    LOCAL_NEXT.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(RUNNER), "--exe", str(ROOT / "chencode_sim_corrected.exe"),
               "--output-dir", str(LOCAL_NEXT), "--hours", "24", "--workers", "12",
               "--frames-per-shard", "1000", "--target-frames-per-snr", "3000000",
               "--target-frame-errors-per-snr", "300", "--minimum-frames-per-snr", "1000",
               "--k", "1784", "--rate", "1", "--snr", "1.6"]
    subprocess.Popen(command, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def launch_mini() -> None:
    command = (f"cd {MINI_ROOT} && mkdir -p {MINI_NEXT} && "
               f"nohup python3 ccsds_parallel_runner.py --exe build-final/chencode_sim "
               f"--output-dir {MINI_NEXT} --hours 24 --workers 16 --frames-per-shard 1000 "
               f"--target-frames-per-snr 3000000 --target-frame-errors-per-snr 300 "
               f"--minimum-frames-per-snr 1000 --k 1784 --rate 1 --snr 1.8 "
               f"> {MINI_NEXT}/runner.stdout.log 2> {MINI_NEXT}/runner.stderr.log < /dev/null &")
    commandline = ["ssh.exe", "-i", MINI_KEY, "-o", "BatchMode=yes", "-o", "ConnectTimeout=12",
                   "zhangzw0170@192.168.17.2", command]
    subprocess.run(commandline, check=True, timeout=45)


def main() -> None:
    marker_local = ROOT / "output" / "29_ccsds_status" / "continuation_local.started"
    marker_mini = ROOT / "output" / "29_ccsds_status" / "continuation_mini.started"
    log = ROOT / "output" / "29_ccsds_status" / "continuation.log"
    while True:
        if not marker_local.exists() and local_complete(LOCAL_CURRENT):
            launch_local(); marker_local.write_text(time.strftime("%Y-%m-%d %H:%M:%S"), encoding="utf-8")
            log.open("a", encoding="utf-8").write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] launched local K1784 R=1/2 SNR=1.6\n")
        if not marker_mini.exists() and mini_complete():
            launch_mini(); marker_mini.write_text(time.strftime("%Y-%m-%d %H:%M:%S"), encoding="utf-8")
            log.open("a", encoding="utf-8").write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] launched mini K1784 R=1/2 SNR=1.8\n")
        if marker_local.exists() and marker_mini.exists():
            return
        time.sleep(POLL)


if __name__ == "__main__":
    main()
