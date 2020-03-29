# Alternative primers for the ARTIC Network's nCov2019 multiplex PCR

![GIF](https://raw.githubusercontent.com/ItokawaK/Alt_nCov2019_primers/master/nCoV_coverage.gif)

Here we provide some resources regarding the alternative primers for the [ARTIC Network's mupltiplex PCR for SARS-CoV-2](https://github.com/artic-network/artic-ncov2019).

See below article for detail of the modifications:

[**A proposal of alternative primers for the ARTIC Network's multiplex PCR to improve coverage of SARS-CoV-2 genome sequencing**](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v3) (ver. 3)

 Kentaro Itokawa, Tsuyoshi Sekizuka, Masanori Hashino, Rina Tanaka, Makoto Kuroda

doi: https://doi.org/10.1101/2020.03.10.985150

https://www.biorxiv.org/content/10.1101/2020.03.10.985150v3

### Tools
-------
 Those tools are still beta and have no warranty to properly work.

- tools/plot_depth_v0.3.py

   Generates depth plots in sigle PDF file from multiple BAM files to briefly check coverages.

   Requires:

     - samtools in $PATH
     - python3
     - matplotlib
     - numpy
     - pandas

  ```
  Usage: plot_depth_v0.3.py [-h] [-i [BAMS [BAMS ...]]] [-p PRIMER] [-o OUT]

      -i [BAMS [BAMS ...]], --bams [BAMS [BAMS ...]]
                      Paths for input BAMs
      -p PRIMER, --primer PRIMER
                      primer_region in BED format
      -o OUT, --out OUT     Output PDF file name

  ```

- tools/trim_primers/trim_primer_parts.py

    Trim susceptive primer parts from reads. Input a bwa mem output directly via PIPE.

    ```
    Usage:

      bwa mem nCov_bwadb untrimmed_R1.fq untrimmed_R2.fq |
         trim_primer_parts.py [--gziped] primer.bed trimmed_R1 trimmed_R2
    ```
