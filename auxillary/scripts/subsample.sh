#!/bin/bash

for f in *.fq.gz; do
  echo $f
  zcat -d $f | head -4000 | gzip -c -  > $PWD/sub/$f
done
