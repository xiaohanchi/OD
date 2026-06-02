#!/bin/bash

mkdir empty_dir

for i in $(seq 1 42); do
  if [ -d "jobtype${i}" ]; then
    rsync -a --delete empty_dir/ jobtype${i}/
    rmdir jobtype${i}
    if (( i % 10 == 0 )); then
      echo "Progress: deleted up to jobtype${i}"
    fi
  fi
done

rmdir empty_dir