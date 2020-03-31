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

- tools/plot_depth.py

   Generates depth plots in sigle PDF file from multiple BAM files to briefly check coverages.

   Requires:

     - samtools in $PATH
     - python3
     - matplotlib
     - numpy
     - pandas

  ```
  Usage: plot_depth.py [-h] [-i [BAMS [BAMS ...]]] [-p PRIMER] [-o OUT]

      -i [BAMS [BAMS ...]], --bams [BAMS [BAMS ...]]
                      Paths for input BAMs
      -p PRIMER, --primer PRIMER
                      primer_region in BED format
      -o OUT, --out OUT     Output PDF file name
      -r REF_FA, --ref_fa REF_FA
                      Reference fasta file [optional]
  ```
    If `-r` option is set, mismatches found on >80% reads will be highlighted. This, however, takes additional time.

    Output image

  ![plot_depth_out_image](https://user-images.githubusercontent.com/38896687/77901776-1ed5cb80-72bb-11ea-9b48-fa62a8bbc86a.png)


- tools/trim_primers/trim_primer_parts.py

    Trim primer parts of paired-end reads obtained from illumina machines.

    This program work as:
       1. Look alignments of mapped fragments from paired reads.
       1. Find fragment ends *contained* in a primer region.
       1. Trim sequence overlapping the primer region.

    Note that this program work conservatively, it throw away all soft masked parts and reads not mapped in proper pair (non 0x2 bit flag).

![trimming_image](https://user-images.githubusercontent.com/38896687/78016726-2a41f900-7386-11ea-8dfd-a3960ee3283f.PNG)

    Send the bwa mem output or name sorted SAM to this script *via* PIPE.
    ```
    Usage:

      bwa mem nCov_bwadb untrimmed_R1.fq untrimmed_R2.fq |
         trim_primer_parts.py [--gziped] primer.bed trimmed_R1 trimmed_R2
    ```
    

![trimming_image](https://user-images.githubusercontent.com/38896687/77902160-b89d7880-72bb-11ea-9ef6-9beaa33310bb.png)
