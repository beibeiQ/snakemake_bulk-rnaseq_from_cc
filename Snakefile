# ============================================================
# Bulk RNA-seq Snakemake Pipeline
# FastQC → fastp → FastQC → MultiQC → Salmon → MultiQC → DESeq2
# ============================================================

import os
import re
import pandas as pd
from pathlib import Path

# ── 加载配置 ──────────────────────────────────────────────
configfile: "config.yaml"

FASTQ_DIR     = config["fastq_dir"]
SALMON_INDEX  = config["genome"]["salmon_index"]
GENOME_FA     = config["genome"]["genome_fasta"]
TRANSCRIPTOME = config["genome"]["transcriptome_fasta"]
GTF           = config["genome"]["gtf"]

# ── 自动识别样本与 paired-end 文件 ────────────────────────
def detect_samples(fastq_dir):
    """
    自动识别 paired-end fastq 文件，支持以下命名规则：
      sample_R1.fastq.gz / sample_R2.fastq.gz
      sample_1.fastq.gz  / sample_2.fastq.gz
      sample_f.fastq.gz  / sample_r.fastq.gz  (或 _F/_R)
    返回 dict: {sample_name: {"R1": path, "R2": path}}
    """
    patterns = [
        (r"^(.+?)_R1(\.fastq\.gz|\.fq\.gz)$",  r"_R2"),
        (r"^(.+?)_1(\.fastq\.gz|\.fq\.gz)$",   r"_2"),
        (r"^(.+?)_f(\.fastq\.gz|\.fq\.gz)$",   r"_r"),
        (r"^(.+?)_F(\.fastq\.gz|\.fq\.gz)$",   r"_R"),
    ]
    samples = {}
    for fname in sorted(os.listdir(fastq_dir)):
        for pat, pair_suffix in patterns:
            m = re.match(pat, fname, re.IGNORECASE)
            if m:
                sample = m.group(1)
                ext    = m.group(2)
                r1     = os.path.join(fastq_dir, fname)
                r2     = os.path.join(fastq_dir, sample + pair_suffix + ext)
                if os.path.exists(r2):
                    samples[sample] = {"R1": r1, "R2": r2}
                break
    if not samples:
        raise ValueError(f"在 {fastq_dir} 中未找到任何 paired-end fastq 文件，请检查命名格式。")
    return samples

SAMPLES = detect_samples(FASTQ_DIR)
SAMPLE_NAMES = list(SAMPLES.keys())

# ── 辅助函数 ──────────────────────────────────────────────
def get_r1(wildcards):
    return SAMPLES[wildcards.sample]["R1"]

def get_r2(wildcards):
    return SAMPLES[wildcards.sample]["R2"]

# 获取所有 FASTQ 文件的完整路径列表
ALL_FASTQS = [path for sample_info in SAMPLES.values() for path in (sample_info["R1"], sample_info["R2"])]

# ── 读取 metadata ─────────────────────────────────────────
metadata = pd.read_csv(config["samples_file"], sep="\t", index_col=0)

# 获取所有唯一的 tissue 名称
TISSUES = metadata["tissue"].dropna().unique().tolist()

# 验证 metadata 中的样本与检测到的样本一致
missing = set(SAMPLE_NAMES) - set(metadata.index)
if missing:
    raise ValueError(f"以下样本在 samples.tsv 中缺失: {missing}")

# ── 辅助函数 ──────────────────────────────────────────────
def get_r1(wildcards):
    return SAMPLES[wildcards.sample]["R1"]

def get_r2(wildcards):
    return SAMPLES[wildcards.sample]["R2"]

# ── 最终目标 ──────────────────────────────────────────────
rule all:
    input:
        # 原始 QC
        expand("results/fastqc/raw/{sample}_{read}_fastqc.html",
               sample=SAMPLE_NAMES, read=["R1", "R2"]),
        # trim 后 QC
        expand("results/fastqc/trimmed/{sample}_{read}_fastqc.html",
               sample=SAMPLE_NAMES, read=["R1", "R2"]),
        # MultiQC 报告
        "results/multiqc/pre_trim/multiqc_report.html",
        "results/multiqc/post_trim/multiqc_report.html",
        "results/multiqc/salmon/multiqc_report.html",
        # Salmon 定量
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLE_NAMES),
        # DESeq2 结果
        "results/deseq2/pca_plot.pdf",
        expand("results/deseq2/{tissue}/{contrast}_DEG_results.csv",
               tissue=TISSUES,
               contrast=["_vs_".join(c) for c in config["deseq2"]["contrasts"]]),

# ════════════════════════════════════════════════════════════
# 0. 原始数据完整性校验 (MD5 或 gzip -t)
# ════════════════════════════════════════════════════════════
rule check_fastq_integrity:
    input:
        fastqs = ALL_FASTQS
    output:
        flag = "results/qc/integrity_passed.flag"
    log:
        "logs/check_fastq_integrity.log"
    threads: 8 # 可以根据本地 CPU 核心数进行调整，加速多文件的并行检测
    run:
        import subprocess
        import os
        from concurrent.futures import ThreadPoolExecutor

        corrupt_dir = "results/corrupted_files"
        os.makedirs(corrupt_dir, exist_ok=True)
        corrupted_list = []

        # 定义并行检测函数
        def check_gz(f):
            if f.endswith(".gz"):
                # 使用 gzip -t 快速检测压缩包尾部和 CRC
                res = subprocess.run(["gzip", "-t", f], capture_output=True)
                return f, res.returncode
            return f, 0 # 非 gz 文件默认跳过压缩校验

        with open(log[0], "w") as log_file:
            log_file.write("Starting FastQ integrity check...\n")
            
            # 使用多线程加速检测
            with ThreadPoolExecutor(max_workers=threads) as pool:
                results = pool.map(check_gz, input.fastqs)
                for f, retcode in results:
                    if retcode != 0:
                        corrupted_list.append(f)
                        bad_flag = os.path.join(corrupt_dir, os.path.basename(f) + ".corrupted")
                        with open(bad_flag, "w") as err:
                            err.write(f"File corrupted or truncated: {f}")
                        log_file.write(f"ERROR: {f} is corrupted.\n")

            if corrupted_list:
                err_msg = f"检测到损坏的 FASTQ 文件！已在 {corrupt_dir} 创建标记。流程终止。"
                log_file.write("\n" + err_msg + "\n")
                raise ValueError(err_msg)

        # 全部通过后，生成 flag 文件，解锁下游任务
        with open(output.flag, "w") as f:
            f.write("All files passed gzip integrity check.\n")




# ════════════════════════════════════════════════════════════
# 1. 构建 Salmon decoy-aware 索引
# ════════════════════════════════════════════════════════════
rule salmon_decoy_list:
    input:
        genome = GENOME_FA
    output:
        decoy = "refs/decoys.txt"
    shell:
        "grep '^>' {input.genome} | cut -d ' ' -f 1 | sed 's/>//' > {output.decoy}"

rule salmon_gentrome:
    input:
        transcriptome = TRANSCRIPTOME,
        genome        = GENOME_FA
    output:
        gentrome = "refs/gentrome.fa"
    shell:
        "cat {input.transcriptome} {input.genome} > {output.gentrome}"

rule salmon_index:
    input:
        gentrome = "refs/gentrome.fa",
        decoy    = "refs/decoys.txt"
    output:
        directory(SALMON_INDEX)
    threads: config["threads"]["salmon_index"]
    conda: "envs/salmon.yaml"
    log: "logs/salmon_index.log"
    shell:
        """
        salmon index \
            --transcripts {input.gentrome} \
            --decoys {input.decoy} \
            --index {output} \
            --threads {threads} \
            2> {log}
        """


# ════════════════════════════════════════════════════════════
# 2. 原始数据 FastQC
# ════════════════════════════════════════════════════════════
rule fastqc_raw:
    input:
        r1 = get_r1,
        r2 = get_r2,
        integrity = "results/qc/integrity_passed.flag"
    output:
        html_r1 = "results/fastqc/raw/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc/raw/{sample}_R2_fastqc.html",
        zip_r1  = "results/fastqc/raw/{sample}_R1_fastqc.zip",
        zip_r2  = "results/fastqc/raw/{sample}_R2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    conda: "envs/qc.yaml"
    log: "logs/fastqc/raw/{sample}.log"
    run:
        import subprocess, glob, re, os
        os.makedirs("results/fastqc/raw", exist_ok=True)
        with open(log[0], "w") as logf:
            subprocess.run(
                ["fastqc", "--threads", str(threads),
                 "--outdir", "results/fastqc/raw",
                 input.r1, input.r2],
                stderr=logf, check=True
            )
        # FastQC 输出文件名基于输入文件 stem，精确重命名为 _R1/_R2
        for src_path, dst_html, dst_zip in [
            (input.r1, output.html_r1, output.zip_r1),
            (input.r2, output.html_r2, output.zip_r2),
        ]:
            stem = Path(src_path).name
            for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
                stem = stem.replace(ext, "")
            src_html = f"results/fastqc/raw/{stem}_fastqc.html"
            src_zip  = f"results/fastqc/raw/{stem}_fastqc.zip"
            if src_html != dst_html:
                os.rename(src_html, dst_html)
            if src_zip != dst_zip:
                os.rename(src_zip, dst_zip)


# ════════════════════════════════════════════════════════════
# 3. MultiQC — 原始数据汇总
# ════════════════════════════════════════════════════════════
rule multiqc_raw:
    input:
        expand("results/fastqc/raw/{sample}_{read}_fastqc.zip",
               sample=SAMPLE_NAMES, read=["R1", "R2"])
    output:
        "results/multiqc/pre_trim/multiqc_report.html"
    conda: "envs/qc.yaml"
    log: "logs/multiqc_raw.log"
    shell:
        """
        multiqc results/fastqc/raw \
            --outdir results/multiqc/pre_trim \
            --filename multiqc_report.html \
            --force 2> {log}
        """


# ════════════════════════════════════════════════════════════
# 4. fastp 质控与过滤
# ════════════════════════════════════════════════════════════
rule fastp:
    input:
        r1 = get_r1,
        r2 = get_r2,
        integrity = "results/qc/integrity_passed.flag"
    output:
        r1      = "results/trimmed/{sample}_R1.fastq.gz",
        r2      = "results/trimmed/{sample}_R2.fastq.gz",
        html    = "results/fastp/{sample}_fastp.html",
        json    = "results/fastp/{sample}_fastp.json"
    threads: config["threads"]["fastp"]
    conda: "envs/qc.yaml"
    log: "logs/fastp/{sample}.log"
    params:
        min_len   = config["fastp"]["min_length"],
        qual      = config["fastp"]["qualified_quality_phred"],
        unqual_pct= config["fastp"]["unqualified_percent_limit"],
        extra     = config["fastp"]["extra"]
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.qual} \
            --unqualified_percent_limit {params.unqual_pct} \
            {params.extra} \
            2> {log}
        """


# ════════════════════════════════════════════════════════════
# 5. trim 后 FastQC
# ════════════════════════════════════════════════════════════
rule fastqc_trimmed:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz"
    output:
        html_r1 = "results/fastqc/trimmed/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc/trimmed/{sample}_R2_fastqc.html",
        zip_r1  = "results/fastqc/trimmed/{sample}_R1_fastqc.zip",
        zip_r2  = "results/fastqc/trimmed/{sample}_R2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    conda: "envs/qc.yaml"
    log: "logs/fastqc/trimmed/{sample}.log"
    shell:
        """
        fastqc --threads {threads} --outdir results/fastqc/trimmed \
            {input.r1} {input.r2} 2> {log}
        """


# ════════════════════════════════════════════════════════════
# 6. MultiQC — trim 后汇总（含 fastp 报告）
# ════════════════════════════════════════════════════════════
rule multiqc_trimmed:
    input:
        fastqc = expand("results/fastqc/trimmed/{sample}_{read}_fastqc.zip",
                        sample=SAMPLE_NAMES, read=["R1", "R2"]),
        fastp  = expand("results/fastp/{sample}_fastp.json", sample=SAMPLE_NAMES)
    output:
        "results/multiqc/post_trim/multiqc_report.html"
    conda: "envs/qc.yaml"
    log: "logs/multiqc_trimmed.log"
    shell:
        """
        multiqc results/fastqc/trimmed results/fastp \
            --outdir results/multiqc/post_trim \
            --filename multiqc_report.html \
            --force 2> {log}
        """


# ════════════════════════════════════════════════════════════
# 7. Salmon 定量（pseudo-alignment）
# ════════════════════════════════════════════════════════════
rule salmon_quant:
    input:
        r1    = "results/trimmed/{sample}_R1.fastq.gz",
        r2    = "results/trimmed/{sample}_R2.fastq.gz",
        index = SALMON_INDEX
    output:
        quant = "results/salmon/{sample}/quant.sf",
        lib   = "results/salmon/{sample}/lib_format_counts.json"
    threads: config["threads"]["salmon_quant"]
    conda: "envs/salmon.yaml"
    log: "logs/salmon/{sample}.log"
    params:
        lib_type = config["salmon"]["lib_type"],
        extra    = config["salmon"]["extra"],
        outdir   = "results/salmon/{sample}"
    shell:
        """
        salmon quant \
            --index {input.index} \
            --libType {params.lib_type} \
            -1 {input.r1} -2 {input.r2} \
            --output {params.outdir} \
            --threads {threads} \
            --numBootstraps 100 \
            {params.extra} \
            2> {log}
        """


# ════════════════════════════════════════════════════════════
# 8. MultiQC — Salmon 定量 QC 汇总
# ════════════════════════════════════════════════════════════
rule multiqc_salmon:
    input:
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLE_NAMES)
    output:
        "results/multiqc/salmon/multiqc_report.html"
    conda: "envs/qc.yaml"
    log: "logs/multiqc_salmon.log"
    shell:
        """
        multiqc results/salmon \
            --outdir results/multiqc/salmon \
            --filename multiqc_report.html \
            --force 2> {log}
        """


# ════════════════════════════════════════════════════════════
# 9. DESeq2 差异分析（tximport → PCA → DEG）
# ════════════════════════════════════════════════════════════
rule deseq2:
    input:
        quant_files = expand("results/salmon/{sample}/quant.sf", sample=SAMPLE_NAMES),
        metadata    = config["samples_file"],
        gtf         = GTF
    output:
        pca     = "results/deseq2/pca_plot.pdf",
        results = expand("results/deseq2/{tissue}/{contrast}_DEG_results.csv",
                         tissue=TISSUES, 
                         contrast=["_vs_".join(c) for c in config["deseq2"]["contrasts"]])
    conda: "envs/deseq2.yaml"
    log: "logs/deseq2.log"
    params:
        salmon_dir           = "results/salmon",
        outdir               = "results/deseq2",
        contrasts            = config["deseq2"]["contrasts"],
        padj_cutoff          = config["deseq2"]["padj_cutoff"],
        lfc_cutoff           = config["deseq2"]["log2fc_cutoff"],
        count_type           = config["deseq2"]["count_type"],
        use_tissue_covariate = config["deseq2"]["use_tissue_covariate"],
        split_by_tissue      = config["deseq2"]["split_by_tissue"]
    script:
        "scripts/deseq2_analysis.R"
