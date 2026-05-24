params.raw = "${workDir}/raw/*{1,2}.fastq.gz"
params.outdir    = "results_genomics"
params.reference_folder = "ref"
params.reference        = "hg38.fasta"
params.reference_index  = "hg38.fasta.fai"
params.reference_dict   = "hg38.dict"
params.bqsr_1000g       = "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.bqsr_mills       = "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.bqsr_dbsnp       = "dbsnp_144.hg38.vcf.gz"
params.vqsr_omni        = "1000G_omni2.5.hg38.vcf.gz"
params.vqsr_hapmap      = "hapmap_3.3.hg38.vcf.gz"
params.tmp_dir       = "${workDir}/tmp/"
params.threads          = 8

process download_reference {
    conda "samtools"

    publishDir params.outdir, mode: 'symlink'

    input:
        val folder
        val reference
        val bqsr_1000g
        val bqsr_mills
        val bqsr_dbsnp

    output:
        path reference, emit: reference
        path bqsr_1000g, emit: bqsr_1000g
        path bqsr_mills, emit: bqsr_mills
        path bqsr_dbsnp, emit: bqsr_dbsnp

    script:
    """
    mkdir ${folder}
    wget -O ${reference} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/GRCh38.primary_assembly.genome.fa
    wget -O ${bqsr_1000g} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    wget -O ${bqsr_mills} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    wget -O ${bqsr_dbsnp} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/gdc/dbsnp_144.hg38.vcf.gz
    """
}

process index_reference {
    conda "samtools"

    publishDir params.outdir, mode: 'symlink'

    input:
        path reference
        val reference_index
        val reference_dict

    output:
        path reference, emit: reference
        path reference_index, emit: reference_index
        path reference_dict, emit: reference_dict

    script:
    """
    samtools faidx -o ${reference_index} ${reference}
    samtools dict -o ${reference_dict} ${reference}
    """
}

process fastqc {
    conda "fastqc"

    publishDir params.outdir, mode: 'symlink'

    input:
        tuple val(sample_id), path(reads)

    output:
        path "*_fastqc.{zip,html}"

    script:
    """
    fastqc $reads
    """
}

process multiqc {
    conda "multiqc"

    publishDir params.outdir, mode: 'symlink'
 
    input:
        path fastqc_result

    output:
        file "multiqc_report.html"
        file "multiqc_data"

    script:
    """
    multiqc .
    """
}

process soapnuke_filter {
    container "filisiaa/soapnuke:2.1.8"

    publishDir params.outdir, mode: 'symlink'

    input:
        tuple val(sample_id), path(reads)

    output:
        val sample_id, emit: sample_id
        path "${sample_id}_1.filtered.fq.gz", emit: reads1
        path "${sample_id}_2.filtered.fq.gz", emit: reads2

    script:
    """
    SOAPnuke filter \
    -1 ${reads[0]} -2 ${reads[1]} \
    -C ${sample_id}_1.filtered.fq.gz \
    -D ${sample_id}_2.filtered.fq.gz
    """
}

process bwa_index {
    conda "bwa"

    publishDir params.outdir, mode: 'symlink'

    input:
        path ref_fasta

    output:
        path "${ref_fasta}.{,amb,ann,bwt,pac,sa}"

    script:
    """
    bwa index ${ref_fasta}
    """
}

process bwa_mem {
    conda "bwa samtools"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path reads_1
        path reads_2
        path ref_fasta
        path ref_fasta_index
        val threads

    output:
        val(sample_id), emit: sample_id
        path("${sample_id}.raw.bam"), emit: bam

    script:
    """
    id=\$(zcat < ${reads_1} | head -n 1 | cut -f 3-4 -d":" | sed 's/^@//')
    bwa mem -M -Y -R "@RG\\tID:\${id}\\tSM:${sample_id}\\tPL:COMPLETE\\tCN:BGI" -t ${threads} ${ref_fasta} ${reads_1} ${reads_2} | samtools view -Sb -o ${sample_id}.raw.bam -
    """
}

process samtools_sort {
    conda "samtools"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam
        val threads

    output:
        path "${sample_id}.sorted.bam"

    script:
    """
    samtools sort '$input_bam' -o '${sample_id}.sorted.bam'
    """
}

process gatk_mark_duplicates {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam

    output:
        path "${sample_id}.marked.bam", emit: bam
        path "${sample_id}.marked.bai", emit: bai
        path "${sample_id}.metrics.txt", emit: metrics

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def java_options = [
        avail_mem ? "-Xmx${avail_mem}G" : "",
        "-Djava.io.tmpdir='\${PWD}/tmp'",
        "-XX:+UseSerialGC",
    ]
    """
    gatk \
        --java-options "${java_options.join(' ')}" \
        MarkDuplicates \
        -I ${input_bam} \
        -O ${sample_id}.marked.bam \
        -M ${sample_id}.metrics.txt \
        --CREATE_INDEX true
    """
}

process gatk_index_known {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        path bqsr_1000g
        path bqsr_mills
        path bqsr_dbsnp
        path vqsr_omni
        path vqsr_hapmap

    output:
        path "${bqsr_1000g}.tbi", emit: bqsr_1000g
        path "${bqsr_mills}.tbi", emit: bqsr_mills
        path "${bqsr_dbsnp}.tbi", emit: bqsr_dbsnp
        path "${vqsr_omni}.tbi", emit: vqsr_omni
        path "${vqsr_hapmap}.tbi", emit: vqsr_hapmap

    script:
    """
    gatk IndexFeatureFile -I ${bqsr_1000g} -O ${bqsr_1000g}.tbi
    gatk IndexFeatureFile -I ${bqsr_mills} -O ${bqsr_mills}.tbi
    gatk IndexFeatureFile -I ${bqsr_dbsnp} -O ${bqsr_dbsnp}.tbi
    gatk IndexFeatureFile -I ${vqsr_omni} -O ${vqsr_omni}.tbi
    gatk IndexFeatureFile -I ${vqsr_hapmap} -O ${vqsr_hapmap}.tbi
    """
}

process gatk_base_recalibrator {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam
        path ref_fasta
        path ref_index
        path ref_dict
        path bqsr_1000g
        path bqsr_mills
        path bqsr_dbsnp
        path bqsr_1000g_index
        path bqsr_mills_index
        path bqsr_dbsnp_index
        path tmp_dir

    output:
        path "${sample_id}.recal_data.table", emit: table

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def java_options = [
        avail_mem ? "-Xmx${avail_mem}G" : "",
        "-Djava.io.tmpdir='\${PWD}/tmp'",
        "-XX:+UseSerialGC",
    ]
    """
    gatk \
        --java-options "${java_options.join(' ')}" \
        BaseRecalibrator \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${sample_id}.recal_data.table \
        --known-sites ${bqsr_1000g} \
        --known-sites ${bqsr_mills} \
        --known-sites ${bqsr_dbsnp} \
        --tmp-dir ${tmp_dir}
    """
}

process gatk_apply_bqsr {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam
        path recal_table
        path ref_fasta
        path ref_index
        path ref_dict
        path tmp_dir

    output:
        path "${sample_id}.bqsr.bam", emit: bam
        path "${sample_id}.bqsr.bai", emit: bai

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def java_options = [
        avail_mem ? "-Xmx${avail_mem}G" : "",
        "-Djava.io.tmpdir='\${PWD}/tmp'",
        "-XX:+UseSerialGC",
    ]
    """
    gatk \
        --java-options "${java_options.join(' ')}" \
        ApplyBQSR \
        -R ${ref_fasta} \
        -I ${input_bam} \
        --bqsr-recal-file ${recal_table} \
        -O ${sample_id}.bqsr.bam \
        --create-output-bam-index true \
        --tmp-dir ${tmp_dir}
    """
}

process gatk_haplotype_caller {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path tmp_dir

    output:
        path "${sample_id}.vcf"     , emit: vcf
        path "${sample_id}.vcf.idx" , emit: idx

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def java_options = [
        avail_mem ? "-Xmx${avail_mem}G" : "",
        "-Djava.io.tmpdir='\${PWD}/tmp'",
        "-XX:+UseSerialGC",
    ]
    """
    gatk \
        --java-options "${java_options.join(' ')}" \
        HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${sample_id}.vcf \
        --tmp-dir ${tmp_dir}
    """
}

process deepvariant {

    container "google/deepvariant:1.6.1"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index

    output:
        path "${sample_id}.deepvariant.vcf.gz", emit: vcf
        path "${sample_id}.deepvariant.vcf.gz.tbi", emit: tbi
        path "${sample_id}.deepvariant.gvcf.gz", emit: gvcf

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${ref_fasta} \
        --reads=${input_bam} \
        --output_vcf=${sample_id}.deepvariant.vcf.gz \
        --output_gvcf=${sample_id}.deepvariant.gvcf.gz \
        --num_shards=\$(nproc)
    """
}

process freebayes {

    container "staphb/freebayes:1.3.10"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index

    output:
        path "${sample_id}.freebayes.vcf", emit: vcf

    script:
    """
    freebayes \
        -f ${ref_fasta} \
        ${input_bam} \
        > ${sample_id}.freebayes.vcf
    """
}

process download_vqsr_resources {
    conda "samtools"

    publishDir params.outdir, mode: 'symlink'

    input:
        val vqsr_omni
        val vqsr_hapmap

    output:
        path vqsr_omni, emit: vqsr_omni
        path vqsr_hapmap, emit: vqsr_hapmap

    script:
    """
    wget -O ${vqsr_omni} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz
    wget -O ${vqsr_hapmap} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
    """
}

process gatk_variant_recalibrator {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_vcf
        path ref_fasta
        path ref_index
        path ref_dict
        path tmp_dir
        path vqsr_1000g
        path vqsr_omni
        path vqsr_dbsnp
        path vqsr_hapmap
        path vqsr_1000g_index
        path vqsr_omni_index
        path vqsr_dbsnp_index
        path vqsr_hapmap_index

    output:
        path "${sample_id}.recal", emit: table
        path "${sample_id}.tranches", emit: tranches

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def java_options = [
        avail_mem ? "-Xmx${avail_mem}G" : "",
        "-Djava.io.tmpdir='\${PWD}/tmp'",
        "-XX:+UseSerialGC",
    ]
    """
    gatk \
        --java-options "${java_options.join(' ')}" \
        VariantRecalibrator \
        -R ${ref_fasta} \
        -V ${input_vcf} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${vqsr_hapmap} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${vqsr_omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${vqsr_1000g} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${vqsr_dbsnp} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O ${sample_id}.recal \
        --tranches-file ${sample_id}.tranches \
        --tmp-dir ${tmp_dir}
    """
}

process gatk_index_recal {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        path recal

    output:
        path "${recal}.idx"

    script:
    """
    gatk IndexFeatureFile -I ${recal} -O ${recal}.idx
    """
}

process gatk_apply_vqsr {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path input_vcf
        path ref_fasta
        path ref_index
        path ref_dict
        path tranches
        path recal
        path recal_index
        path tmp_dir

    output:
        path "${sample_id}.vqsred.vcf.gz"

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0
    def java_options = [
        avail_mem ? "-Xmx${avail_mem}G" : "",
        "-Djava.io.tmpdir='\${PWD}/tmp'",
        "-XX:+UseSerialGC",
    ]
    """
    gatk \
        --java-options "${java_options.join(' ')}" \
        ApplyVQSR \
        -R ${ref_fasta} \
        -V ${input_vcf} \
        -O ${sample_id}.vqsred.vcf.gz \
        --ts_filter_level 99.0 \
        --tranches-file ${tranches} \
        --recal-file ${recal} \
        -mode SNP \
        --tmp-dir ${tmp_dir}
    """
}

process filter_freebayes_vcf {
    conda "bcftools py-bgzip tabix"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.freebayes.filtered.vcf"

    script:
    """
        bgzip ${vcf}
        tabix -p vcf ${vcf}.gz
        bcftools filter -i 'QUAL>5 & INFO/DP>2' ${vcf}.gz > ${sample_id}.freebayes.filtered.vcf
    """
}

process filter_deepvariant_vcf {
    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.deepvariant.filtered.vcf"

    script:
    """
        awk -F'\\t' '\$7=="PASS"' ${vcf} | grep "^chr" > ${sample_id}.deepvariant.filtered.vcf
    """
}

process filter_gatk_vcf {
    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.gatk.filtered.vcf"

    script:
    """
        awk -F'\\t' '\$7=="PASS"' ${vcf} | grep "^chr" > ${sample_id}.gatk.filtered.vcf
    """
}

process merge_vcf {
    conda "sklearn pandas numpy"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path haplotype_caller_vcf
        path deepvariant_vcf
        path freebayes_vcf

    output:
        path "${sample_id}.merged.vcf", emit: merged

    script:
    """
    python isolation_forest_predict.py \
        ${haplotype_caller_vcf} \
        ${deepvariant_vcf} \
        ${freebayes_vcf} \
        -o binary_vector.csv
    
    python construct_vcf.py \
        ${haplotype_caller_vcf} \
        ${deepvariant_vcf} \
        ${freebayes_vcf} \
        binary_vector.csv \
        ${sample_id}.merged.vcf
    """
}

process snpeff {
    conda "snpeff"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.snpeff.vcf"

    script:
    """
        snpeff  -v -stats variants_report.html hg38 -canon ${vcf} > ${sample_id}.snpeff.vcf
    """
}

process annotate_with_alphamissence {
    conda "snpsift tabix"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.snpeff.am.vcf"

    script:
    """
        wget https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz
        tabix -p vcf AlphaMissense_hg38.tsv.gz
        snpsift annotate -v AlphaMissense_hg38.tsv.gz ${vcf} > ${sample_id}.snpeff.am.vcf
    """
}

process annotate_with_clinvar {
    conda "snpsift"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.snpeff.am.clinvar.vcf"

    script:
    """
        wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
        wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
        snpsift annotate -v clinvar.vcf.gz ${vcf} > ${sample_id}.snpeff.am.clinvar.vcf
    """
}

process annotate_with_dbsnp {
    conda "snpsift"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}.snpeff.am.clinvar.dbsnp.vcf"

    script:
    """
        wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
        wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi
        snpsift annotate -v GCF_000001405.40.gz ${vcf} > ${sample_id}.snpeff.am.clinvar.dbsnp.vcf
    """
}

process make_gpt_call {
    conda "openai"

    publishDir params.outdir, mode: 'symlink'

    input:
        val sample_id
        path vcf

    output:
        path "${sample_id}_report.txt"

    script:
    """
        python gpt_call.py ${vcf} > ${sample_id}_report.txt
    """
}

workflow {

    reads_ch = channel.fromFilePairs(params.raw, checkIfExists: true)

    files = download_reference(
        params.reference_folder,
        params.reference,
        params.bqsr_1000g,
        params.bqsr_mills,
        params.bqsr_dbsnp,
    )

    reference = index_reference(
        files.reference,
        params.reference_index,
        params.reference_dict,
    )

    ref_file        = reference.reference
    ref_index       = reference.reference_index
    ref_dict_file   = reference.reference_dict
    bqsr_1000g      = files.bqsr_1000g
    bqsr_mills      = files.bqsr_mills
    bqsr_dbsnp      = files.bqsr_dbsnp
    tmp_dir         = file(params.tmp_dir)
    threads         = params.threads

    fastqc_result = fastqc(reads_ch)
    multiqc(fastqc_result)

    filtered_reads = soapnuke_filter(reads_ch)

    bwa_indices = bwa_index(ref_file)

    result = bwa_mem(
        filtered_reads.sample_id, 
        filtered_reads.reads_1, 
        filtered_reads.reads_2, 
        ref_file, 
        bwa_indices, 
        threads,
    )

    sorted = samtools_sort(result.sample_id, result.bam, threads)

    marked = gatk_mark_duplicates(result.sample_id, sorted)

    vqsr_files = download_vqsr_resources(
        params.vqsr_omni,
        params.vqsr_hapmap,
    )

    gatk_indices = gatk_index_known(bqsr_1000g, bqsr_mills, bqsr_dbsnp, vqsr_files.vqsr_omni, vqsr_files.vqsr_hapmap)

    recal_table = gatk_base_recalibrator(
        result.sample_id,
        marked.bam,
        ref_file,
        ref_index,
        ref_dict_file,
        bqsr_1000g,
        bqsr_mills,
        bqsr_dbsnp,
        gatk_indices.bqsr_1000g,
        gatk_indices.bqsr_mills,
        gatk_indices.bqsr_dbsnp,
        tmp_dir,
    )

    bqsr_bam = gatk_apply_bqsr(
        result.sample_id,
        marked.bam,
        recal_table,
        ref_file,
        ref_index,
        ref_dict_file,
        tmp_dir,
    )

    hc_result = gatk_haplotype_caller(
        result.sample_id,
        bqsr_bam.bam,
        bqsr_bam.bai,
        ref_file,
        ref_index,
        ref_dict_file,
        tmp_dir,
    )

    dv_result = deepvariant(
        result.sample_id,
        bqsr_bam.bam,
        bqsr_bam.bai,
        ref_file,
        ref_index,
    )

    fb_result = freebayes(
        result.sample_id,
        bqsr_bam.bam,
        bqsr_bam.bai,
        ref_file,
        ref_index,
    )

    data = gatk_variant_recalibrator(
        result.sample_id,
        hc_result.vcf,
        ref_file,
        ref_index,
        ref_dict_file,
        tmp_dir,
        bqsr_mills,
        vqsr_files.vqsr_omni,
        bqsr_dbsnp,
        vqsr_files.vqsr_hapmap,
        gatk_indices.bqsr_mills,
        gatk_indices.vqsr_omni,
        gatk_indices.bqsr_dbsnp,
        gatk_indices.vqsr_hapmap,
    )

    recal_index = gatk_index_recal(
        data.recal
    )

    vqsred = gatk_apply_vqsr(
        result.sample_id,
        hc_result.vcf,
        ref_file,
        ref_index,
        ref_dict_file,
        data.tranches,
        data.recal,
        recal_index,
        tmp_dir
    )

    hc_filtered = filter_gatk_vcf(
        result.sample_id,
        vqsred,
    )

    dv_filtered = filter_deepvariant_vcf(
        result.sample_id,
        dv_result.vcf,
    )

    fb_filtered = filter_freebayes_vcf(
        result.sample_id,
        fb_result.vcf,
    )

    merged_vcf = merge_vcf(
        result.sample_id,
        hc_filtered,
        dv_filtered,
        fb_filtered,
    )

    annotated_vcf = snpeff(
        result.sample_id,
        merged_vcf,
    )

    annotated_vcf = annotate_with_alphamissence(
        result.sample_id,
        annotated_vcf,
    )

    annotated_vcf = annotate_with_clinvar(
        result.sample_id,
        annotated_vcf,
    )

    annotated_vcf = annotate_with_dbsnp(
        result.sample_id,
        annotated_vcf,
    )

    make_gpt_call(result.sample_id, annotated_vcf)
}