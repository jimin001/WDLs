version 1.0

import "bwa_mem2_ora.wdl" as ora

# Decompress-only test of the decompressOra task in bwa_mem2_ora.wdl.
# Uses the same task, so Terra call caching can reuse these outputs in alignShortReadsOra.
workflow decompressOraOnly {

    input {
        Array[File] ora_fastqs
    }

    scatter (ora_fastq in ora_fastqs) {
        call ora.decompressOra {
            input:
                ora_fastq = ora_fastq
        }
    }

    output {
        Array[File] fastq_gz = decompressOra.fastq_gz
        Array[File] summary = decompressOra.summary
    }

}
