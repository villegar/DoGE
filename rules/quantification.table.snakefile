rule quantification_table:
    input:
        counts = rules.feature_counts.output
    params:
        libs = "".join(el(["\\t"],LIBS)),
        cols = "1," + ",".join([str(i) for i in range(7, 7 + len(LIBS))])
    output:
        COUNTS + "counts.matrix"
    shell:
        "cat {input.counts} | grep -v '^#' | cut -f {params.cols} | \
        sed '1d' | sed '1i\Geneid{params.libs}' > {output}"
