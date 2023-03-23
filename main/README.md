# Selection of Intergenic Sites
- 'grna_off_target.py' identifies the gRNAs with more than X mismatches. We scan both the genomic strands for gRNAs with NGG PAM. gRNA sequences occurring once in the genome are then screened against the gRNA library for at least 6 mismatches (X = 6).
- 'grna_synbio' selects the intergenic gRNA and apply the cloning criteria for efficient plasmid assembly and experiments.
- 'site_selection' integrates information regarding gene density, transcriptomics, regulatory element disruption, gene essentiality and Rule Set 2 on-target gRNA activity to prioritize the selection of intergenic loci.
