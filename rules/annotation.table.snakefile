rule annotation_table:
    input:
        gtf = expand(REF_GENOME + "{file}", file = GENOME_FILENAMES["GTF"])
    output:
        RMD + "gene_annotation.txt"
    shell:
        "sed '/^[[:blank:]]*#/d;s/#.*//' {input.gtf} | awk '($3 == \"gene\")' | \
        awk -F';' '$1=$1' OFS='\\t' > {output}"
