rule emseq_dedup:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}.bam",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        f"{log_dir}/{{library_id}}_emseq_dedup.log",
    output:
        bam = f"{emseq_bam_dir}/{{library_id}}_deduped.bam",
        index = f"{emseq_bam_dir}/{{library_id}}_deduped.bam.bai",
    shell:
        r"""
        rm -f {output.bam}.tmp.*.bam
        samtools view -h -f 0x2 {input.bam} \
        | samtools sort -n -@ 4 -O BAM -o /dev/stdout \
        | dupsifter \
            --add-mate-tags \
            --stats-output {log} \
            {input.fasta} - \
        | samtools sort -@ 8 -o {output.bam}
        samtools index -@ 8 {output.bam}
        """

rule emseq_pileup:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        bam = f"{emseq_bam_dir}/{{library_id}}_deduped.bam",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        f"{log_dir}/{{library_id}}_emseq_pileup.log",
    output:
        vcf = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf.gz",
        tsv = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf_meth_average.tsv",
    params:
        out_base = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf",
    shell:
        """
        biscuit pileup \
	-@ 8 \
	-o {params.out_base} \
        {input.fasta} {input.bam} \
        && bgzip -@ 8 {params.out_base}
        """
rule emseq_post_pileup:
    conda:
        "../config/biscuit-conda-env.yaml",
    input:
        vcf = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf.gz",
    log:
        f"{log_dir}/{{library_id}}_emseq_post_pileup.log",
    output:
        tbi = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.vcf.gz.tbi",
        bed = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_pileup.bed",
        bismark = f"{data_dir}/analysis/emseq/pileup/{{library_id}}_bismark_cov.bed",
    shell:
        """
        tabix -p vcf {input.vcf} \
        && biscuit vcf2bed \
	-t cg {input.vcf} > {output.bed} \
        && biscuit vcf2bed -c {input.vcf} > {output.bismark}
        """
