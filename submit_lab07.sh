#!/bin/sh
set -eu

make clean
make
./siggen

for file in sig-[0-3].bin; do
  echo "== $file =="
  ./band_scan bin "$file" 400000 32 10 | tail -n 1
  ./p_band_scan bin "$file" 400000 32 10 10 10 | tail -n 1
done

echo "After identifying the alien file number, run:"
echo "./scripts/seti-run --eval -a <alien-file-number>"
