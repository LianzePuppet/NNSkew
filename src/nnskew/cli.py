import argparse
import os
import shutil
import subprocess as sp
import sys
import tempfile
from typing import Dict, Iterable, Tuple

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def which_or_die(name: str):
    path = shutil.which(name)
    if path is None:
        eprint(f"[ERROR] Required executable not found: {name}")
        sys.exit(127)
    return path

def fasta_sizes_from_fai(fa: str) -> Dict[str, int]:
    fai = fa + ".fai"
    sizes = {}
    if not os.path.exists(fai):
        eprint("[ERROR] Fasta index (.fai) not found. Run `samtools faidx <fasta>` first.")
        sys.exit(1)
    with open(fai) as fh:
        for line in fh:
            parts = line.split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes

def iter_bins(chrom: str, length: int, bin_size: int) -> Iterable[Tuple[int,int]]:
    start = 0
    while start < length:
        end = min(start + bin_size, length)
        yield start, end
        start = end

def fetch_seq(fasta: str, chrom: str, start: int, end: int) -> str:
    faidx = which_or_die("samtools")
    cmd = [faidx, "faidx", fasta, f"{chrom}:{start+1}-{end}"]
    try:
        out = sp.check_output(cmd, text=True)
    except sp.CalledProcessError:
        return ""
    return "".join([l.strip() for l in out.splitlines() if not l.startswith(">")]).upper()

def count_dinuc(seq: str, xy: str) -> Tuple[int, int]:
    x, y = xy[0], xy[1]
    xy_count = 0
    yx_count = 0
    for i in range(len(seq)-1):
        di = seq[i:i+2]
        if di == xy:
            xy_count += 1
        elif di == y + x:
            yx_count += 1
    return xy_count, yx_count

def write_bedgraph(fasta: str, sizes: Dict[str,int], bin_size: int, nn: str, out_bed: str):
    with open(out_bed, "w") as out:
        for chrom, ln in sizes.items():
            for s, e in iter_bins(chrom, ln, bin_size):
                seq = fetch_seq(fasta, chrom, s, e)
                if not seq:
                    continue
                xy, yx = count_dinuc(seq, nn)
                denom = xy + yx
                val = 0 if denom == 0 else (xy - yx)/denom
                out.write(f"{chrom}\t{s}\t{e}\t{val:.6f}\n")

def main():
    ap = argparse.ArgumentParser(description="Compute dinucleotide skew (NNskew) and output BigWig")
    ap.add_argument("-f", "--fasta", required=True, help="Reference fasta (indexed with samtools faidx)")
    ap.add_argument("-o", "--out", required=True, help="Output bigWig file")
    ap.add_argument("-b", "--bin", type=int, default=50, help="Bin size (default 50)")
    ap.add_argument("--nn", required=True, help="Dinucleotide to compute (e.g. GC, AT, AG)")
    args = ap.parse_args()

    sizes = fasta_sizes_from_fai(args.fasta)

    with tempfile.TemporaryDirectory() as td:
        bedgraph = os.path.join(td, "out.bedGraph")
        write_bedgraph(args.fasta, sizes, args.bin, args.nn.upper(), bedgraph)
        sizes_file = os.path.join(td, "chrom.sizes")
        with open(sizes_file, "w") as f:
            for k,v in sizes.items():
                f.write(f"{k}\t{v}\n")
        bedGraphToBigWig = which_or_die("bedGraphToBigWig")
        sp.check_call([bedGraphToBigWig, bedgraph, sizes_file, args.out])
        eprint(f"[OK] wrote {args.out}")
