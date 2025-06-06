repo = "~/repos/emseq"
data_dir=config["data_dir"]

# Explicitly select which references to build
index_targets = ["unmeth_lambda", "puc19"]

rule all:
    input:
        expand(f"{data_dir}/ref/{{name}}.fa.bwameth.c2t", name=index_targets)


include: f"{repo}/workflows/dev.smk"
