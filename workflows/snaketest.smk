rule emseq_fastp:
    input:
        r1 = f"{emseq_raw_fastq_dir}/{{libid}}_R1.fastq.gz",
    log:
        cmd = f"{log_dir}/{{libid}}-emseq-fastp.log",
        json = f"{log_dir}/{{libid}}-emseq-fastp.json",
        html = f"{log_dir}/{{libid}}-emseq-fastp.html",
    output:
        r1 = f"{emseq_trimmed_fastq_dir}/{{libid}}_R1.fastq.gz",
        r2 = f"{emseq_trimmed_fastq_dir}/{{libid}}_R2.fastq.gz",
        failed = f"{emseq_trimmed_fastq_dir}/{{libid}}-failed.fastq.gz",
    params:
        script = f"{emseq_script_dir}/fastp-emseq-wrapper.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input.r1} \
        {output.r1} \
        {output.r2} \
        {output.failed} \
        {log.cmd} \
        {log.json} \
        {log.html} \
        {params.threads}
        """
rule emseq_biscuit_align:
    input:
        r1 = f"{emseq_trimmed_fastq_dir}/{{libid}}_R1.fastq.gz",
        fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
    log:
        cmd = f"{log_dir}/{{libid}}_emseq_biscuit_align.log",
    output:
        bam = f"{emseq_unmerged_bam_dir}/{{libid}}_unmerged.bam",
    params:
        script = f"{emseq_script_dir}/emseq_biscuit_align_wrapper.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input.r1} \
        {input.fasta} \
        {output.bam} \
        {log.cmd} \
        {params.threads}
        """
  rule emseq_fastp:
      input:
          r1 = f"{emseq_raw_fastq_dir}/{{libid}}_R1.fastq.gz",
      log:
          cmd = f"{log_dir}/{{libid}}-emseq-fastp.log",
          json = f"{log_dir}/{{libid}}-emseq-fastp.json",
          html = f"{log_dir}/{{libid}}-emseq-fastp.html",
      output:
          r1 = f"{emseq_trimmed_fastq_dir}/{{libid}}_R1.fastq.gz",
          r2 = f"{emseq_trimmed_fastq_dir}/{{libid}}_R2.fastq.gz",
          failed = f"{emseq_trimmed_fastq_dir}/{{libid}}-failed.fastq.gz",
      params:
          script = f"{emseq_script_dir}/fastp-emseq-wrapper.sh",
          threads = threads,
      shell:
          """
          {params.script} \
          {input.r1} \
          {output.r1} \
          {output.r2} \
          {output.failed} \
          {log.cmd} \
          {log.json} \
          {log.html} \
          {params.threads}
          """
  rule emseq_biscuit_align:
      input:
          r1 = f"{emseq_trimmed_fastq_dir}/{{libid}}_R1.fastq.gz",
          fasta = f"{ref_dir}/biscuit/{emseq_ref_fasta}",
      log:
          cmd = f"{log_dir}/{{libid}}_emseq_biscuit_align.log",
      output:
          bam = f"{emseq_unmerged_bam_dir}/{{libid}}_unmerged.bam",
      params:
          script = f"{emseq_script_dir}/emseq_biscuit_align_wrapper.sh",
          threads = threads,
      shell:
          """
          {params.script} \
          {input.r1} \
          {input.fasta} \
          {output.bam} \
          {log.cmd} \
          {params.threads}
          """
