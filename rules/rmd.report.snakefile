rule rmd_report:
    input:
        annotation = rules.annotation_table.output,
        counts = rules.quantification_table.output,
        experiment = "exp_design.csv"
    output:
        RMD + "doge_report_simple.html"
    shell:
        "cp html/* rmd/* " + RMD + " && " +
        "cp {input.counts} " + RMD + " && " +
        "cp {input.experiment} " + RMD + " && cd " + RMD + " && " +
        "Rscript -e \"rmarkdown::render(\'doge_report_simple.Rmd\'," +
        "output_dir=\'.\', clean = TRUE, quiet = TRUE)\" && " +
        "cd ../"
