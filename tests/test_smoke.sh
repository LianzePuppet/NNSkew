#!/usr/bin/env bash
set -euo pipefail

FA=mini.fa
cat > $FA <<'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT
EOF
samtools faidx $FA
nnskew -f $FA -o test.bw -b 4 --nn GC
ls -lh test.bw
