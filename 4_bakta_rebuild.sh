mamba activate BaktaEnv
nohup nice find . -maxdepth 1 -name '*.bakta.json' -print0   | parallel -0 --jobs 94 '
      base={= s/\.bakta\.json$// =};
      bakta_io --output "/home/def-jacquesp/jeanneau/9-Pangenomics/sample_annotation_total/${base}_complete" --prefix "$base" "{}"
' > nohup_alljson2complete.out &
