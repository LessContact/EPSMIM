#!/usr/bin/env python3
import subprocess
import re
import argparse

def run(cmd):
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc.stdout

def parse_perfctr(output):
    dp = re.search(r'\|\s*DP \[MFLOP/s\]\s*\|\s*([\d.]+)', output)
    oi = re.search(r'Operational intensity.*\|\s*([\d.]+)', output)
    return (dp.group(1) if dp else 'N/A', oi.group(1) if oi else 'N/A')

def parse_bench(output, metric):
    # metric: 'MFlops/s' or 'MByte/s'
    pat = rf'{re.escape(metric)}:\s*([\d.]+)'
    m = re.search(pat, output)
    return m.group(1) if m else 'N/A'

def main():
    # parse CLI args
    parser = argparse.ArgumentParser(description="Collect likwid metrics")
    parser.add_argument("binary", help="Path to the executable (e.g. ./v2.1)")
    parser.add_argument("--bench", action="store_true",
                        help="Enable likwid‑bench tests (cache/RAM)")
    args = parser.parse_args()

    # 1) likwid-perfctr
    perf_tests = [
        ('MEM_DP',   f"sudo likwid-perfctr -m -C E:N:16 -g MEM_DP {args.binary}"),
        ('FLOPS_DP', f"sudo likwid-perfctr -m -C E:N:16 -g FLOPS_DP {args.binary}"),
    ]
    print("== likwid-perfctr metrics ==")
    for name, cmd in perf_tests:
        out = run(cmd)
        dp, oi = parse_perfctr(out)
        print(f"{name:10} | DP [MFLOP/s]: {dp:>8} | OpInt [FLOP/Byte]: {oi:>6}")

    # 2) likwid-bench tests (optional)
    # ===========================================================================
    # NOTE: THESE ARE ARCHITECTURE SPECIFIC AND NEED TO BE ADJUSTED FOR YOUR CPU
    # ===========================================================================
    if args.bench:
        ws = {
            'L1':  'N:32kB:1',
            'L2':  'N:512kB:1',
            'L3':  'N:4MB:1',
            'RAM': 'N:512MB:1'
        }
        print("\n== likwid-bench cache/RAM bandwidth & flops ==")
        for level, param in ws.items():
            pf = run(f"likwid-bench -t peakflops_avx_fma -W {param}")
            ld = run(f"likwid-bench -t load_avx          -W {param}")
            mf = parse_bench(pf, 'MFlops/s')
            mb = parse_bench(ld, 'MByte/s')
            print(f"{level:4} | PeakFLOP/s: {mf:>8} | Load BW [MB/s]: {mb:>8}")
    else:
        print("\nSkipping likwid-bench tests (use --bench to enable)")

if __name__ == '__main__':
    main()