# Alternative primers for the ARTIC Network's nCov2019 multiplex PCR


Here we provide some resources regarding the alternative primers for the [ARTIC Network's mupltiplex PCR for SARS-CoV-2](https://github.com/artic-network/artic-ncov2019).

See below article for detail of the modifications:

Itokawa K, Sekizuka T, Hashino M, Tanaka R, Kuroda M (2020) Disentangling primer interactions improves SARS-CoV-2 genome sequencing by multiplex tiling PCR. PLoS ONE 15(9): [e0239403](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0239403).

### Primer version (updated on 2020/5/18)
- Primers/ver_niid-200325/: Including the 5 primer exchanges as described in the [preprint ver.3](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v3).
- Primers/ver_niid-200407/: Including the 5 primer exchanges as described in the [preprint ver.3](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v3) + an alternative for the 13_RIGHT.
- Primers/ver_N1/: Including the 12 primer exchanges as described in the [preprint ver.4](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v4.full.pdf) and [peer-reviewed paper in PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0239403)


## Tools
-------
### tools/plot_depth.py

   Generates depth plots in sigle PDF file from multiple BAM files to briefly check coverages.

   Requires:

     - samtools available in $PATH
     - python3
     - matplotlib
     - numpy
     - pandas

  ```
  usage: plot_depth.py [-h] [-i [BAMS [BAMS ...]]] [-o OUT] [-p PRIMER] [-l HIGHLIGHTS] [-r REF_FA] [-t THREADS] [-m MISMATCHES_THRESH] [-s]
                       [--min_readlen MIN_READLEN] [--out_consensus]

  Output depth plot in PDF. Ver: 0.10

  optional arguments:
    -h, --help            show this help message and exit
    -i [BAMS [BAMS ...]], --bams [BAMS [BAMS ...]]
                          Paths for input BAMs
    -o OUT, --out OUT     Output PDF file name
    -p PRIMER, --primer PRIMER
                          Primer regions in BED format [optional]
    -l HIGHLIGHTS, --highlights HIGHLIGHTS
                          Add highlights on selected amplicons. Give amplicon numbers delimited by comma (e.g. 18,76,...) Can only be used with the -p --primer    
                          option. [optional]
    -r REF_FA, --ref_fa REF_FA
                          Reference fasta file [optional]
    -t THREADS, --threads THREADS
                          Num tasks to process concurrently [optional]
    -m MISMATCHES_THRESH, --mismatches_thresh MISMATCHES_THRESH
                          Show mismatches higher than this ratio (default=0.8). Only effective with the -r option [optional]
    -s, --ignore_softclipped
                          Ignore softclipped reads (default=False). [optional]
    --min_readlen MIN_READLEN
                          Minumum length of read (default=0). [optional]
    --dump_consensus       Output consensus to STDOUT. Experimental.
  ```
    If `-r` option is set, mismatches found on >80% reads (parsed from *mpileup*'s output) will be highlighted. This, however, takes additional time. Yellow and red lines indicate mismatches out of and inside of a primer target region, respecitively.

    **Important Note 1**: This program is designed to use the Wuhan-Hu-1 assembly ([MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)) as a reference (-r, --reference). With other assemblies, this program would not work properly.  

    **Important Note 2**: The program now output consensus sequences to STDOUT with --dump_consensus flag (used with the -r, --reference option). This function was implemented so that I can briefly check lineage of viruses by throwing the fasta file to [pangolin](https://github.com/cov-lineages/pangolin). At the moment, the resulted sequencies should not be used for more productive purposes such as depositing to GISAID. For such purposes, there are more sophisticated software, [iVar](https://github.com/andersen-lab/ivar) for instance, to obtain consensuses.

    Output image

![plot_depth_out_image_2](https://user-images.githubusercontent.com/38896687/104012314-ec1a1580-51f2-11eb-97f1-a5c4b21a19c7.png)



### tools/trim_primers/trim_primer_parts.py

    Trim primer parts of paired-end reads obtained from illumina machines.

    This program works as:

 1. Looks alignments of mapped fragments from paired reads.
 1. Finds fragment ends *contained* in a primer region.
 1. Trims sequence overlapping the primer region.

    Note that this program work conservatively, it may throw away all soft masked ends if subjected to trimming.

![trimming_image](https://user-images.githubusercontent.com/38896687/78016726-2a41f900-7386-11ea-8dfd-a3960ee3283f.PNG)

 Send the bwa mem output or name sorted SAM to this script *via* PIPE.
 ```
   Usage:

      bwa mem nCov_bwadb untrimmed_R1.fq untrimmed_R2.fq | \
         trim_primer_parts.py [--gziped] primer.bed trimmed_R1 trimmed_R2
  ```

  Then remap the reads.
  ```
  bwa mem nCov_bwadb trimmed_R1.fq trimmed_R2.fq > ...

  ```

  If you do not want to write trimmed read files on the disk, use `trim_primer_parts2.py`.

  ```
    bwa mem nCov_bwadb untrimmed_R1.fq untrimmed_R2.fq | \
       trim_primer_parts2.py primer.bed trimmed_R1 trimmed_R2 | \
       bwa mem nCov_bwadb -p - | \
       samtools sort > out.bam
  ```

![trimming_image](https://user-images.githubusercontent.com/38896687/77902160-b89d7880-72bb-11ea-9ef6-9beaa33310bb.png)
