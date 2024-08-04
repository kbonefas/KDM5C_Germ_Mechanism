#24.07.24 - Prepare the mm10 genome for bisulfite analysis. Only need to do once

import sys
import os



os.system("bismark_genome_preparation --parallel 4 --path_to_aligner ~/mambaforge/envs/bismark/bin/ --verbose /nfs/value/siwase2/GENOMES/MM10/bismark_genome/")