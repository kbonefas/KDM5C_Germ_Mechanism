# Analysis of bisulfite sequencing

1) Make bismark conda env and install all necessary packages in this environment

2) Fastqc on files (FASTQC_240722.py)

3) Trim the reads to remove adapter contamination (trimREADS_240722.py)

4) Make the bismark genome (only need to do once) - (bismark_genome_prep.py)

5) Run bismark to generate the bam alignment (bismark_align.py)

6) Deduplicate and merge the bam files (bismark_deduplicate.py)

7) Sort and index the bamfiles (sort_bams.py)

8) Read bam files into R and get differentially methylated CpGs (methylkit_mincov3.R)

9) Identify differentially methylated regions over a given input of regions using regionCounts (germline gene promoters or CGIs) (methyl_diff_germCGIpromo_exEpiLC_WTvsKO_min3.R)