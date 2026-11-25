nextflow.enable.dsl=2


process PREPARE {
    tag "${params.sample_name}"
    container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'

    input:
    path vcf
    path index
    path fasta
    path fai

    output:
    path 'proband.vcf.gz', emit: vcf
    path 'proband.vcf.gz.tbi', emit: index
    path 'inheritance.tsv.gz', emit: inheritance
    path 'original_trio_evidence.tsv.gz', emit: original_evidence
    path 'prepare.txt', emit: report

    script:
    def regionArg = params.region ? "-r '${params.region}'" : ''

    """
    set -euo pipefail

    bcftools query -l '${vcf}' > samples.txt

    for sample in '${params.proband}' '${params.mother}' '${params.father}'; do
        grep -Fxq "\$sample" samples.txt
    done

    bcftools view ${regionArg} -Oz \
        -o trio.input.vcf.gz '${vcf}'

    bcftools index -t trio.input.vcf.gz

    {
        printf 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\t${params.proband}\\t${params.mother}\\t${params.father}\\n'

        bcftools query \
            -s '${params.proband},${params.mother},${params.father}' \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT:%DP:%GQ:%AD:%PL:%PGT:%PID:%PS]\\n' \
            trio.input.vcf.gz

    } | gzip -c > original_trio_evidence.tsv.gz

    bcftools norm \
        -f '${fasta}' \
        --check-ref e \
        -m -any \
        --multi-overlaps . \
        --old-rec-tag ORIGINAL_RECORD \
        -Ou trio.input.vcf.gz |

    bcftools sort -Oz \
        -o trio.normalised.vcf.gz

    bcftools index -t trio.normalised.vcf.gz

    {
        printf 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\t${params.proband}\\t${params.mother}\\t${params.father}\\n'

        bcftools query \
            -s '${params.proband},${params.mother},${params.father}' \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT:%DP:%GQ:%AD]\\n' \
            trio.normalised.vcf.gz

    } | gzip -c > inheritance.tsv.gz

    bcftools view \
        -s '${params.proband}' \
        -Ou trio.normalised.vcf.gz |

    bcftools view \
        -i 'GT="alt" && ALT!="*"' \
        -v snps,indels,mnps \
        -Ou |

    bcftools annotate \
        -x INFO/MLEAC,INFO/MLEAF,INFO/AF,FORMAT/PGT,FORMAT/PID \
        -Oz -o proband.vcf.gz

    bcftools index -t proband.vcf.gz

    {
        bcftools --version
        echo 'Samples:'
        cat samples.txt
        echo 'Input records:'
        bcftools index -n trio.input.vcf.gz
        echo 'Normalised trio records:'
        bcftools index -n trio.normalised.vcf.gz
        echo 'Proband records:'
        bcftools index -n proband.vcf.gz
    } > prepare.txt

    test "\$(bcftools index -n proband.vcf.gz)" -gt 0
    """
}

process ANNOTATE {
    tag "${params.proband}"
    container 'ensemblorg/ensembl-vep:release_116.0'
    input:
    path vcf
    path idx
    path fasta
    path fai
    path cache
    path clinvar
    path clinvar_idx
    output:
    path 'annotated.vcf.gz', emit: vcf
    path 'annotated.vcf.gz.tbi', emit: index
    path 'vep.html', emit: stats
    path 'vep.log', emit: log
    script:
    def extra = params.vep_extra ?: ''
    """
    set -euo pipefail
    vep --input_file '${vcf}' --output_file annotated.vcf.gz --vcf --compress_output bgzip \\
      --offline --cache --dir_cache '${cache}' --cache_version '${params.vep_release}' \\
      --species homo_sapiens --assembly GRCh38 --fasta '${fasta}' --check_ref \\
      --everything --allele_number --flag_pick_allele --fork ${task.cpus} \\
      --custom 'file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNDN%CLNREVSTAT' \\
      --stats_file vep.html ${extra} > vep.log 2>&1
    tabix -p vcf annotated.vcf.gz
    """
}

process QC {
    tag "${params.proband}"
    container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_1'
    publishDir params.outdir, mode: 'copy'
    input:
    path original
    path annotated
    path annotated_idx
    path inheritance
    path original_evidence
    path prep_report
    path stats
    path vep_log
    output:
    path "${params.sample_name}.clinical.GRCh38.vcf.gz"
    path "${params.sample_name}.clinical.GRCh38.vcf.gz.tbi"
    path "${params.sample_name}.clinical.PASS.vcf.gz*"
    path "${params.sample_name}.inheritance.tsv.gz"
    path "${params.sample_name}.original_trio_evidence.tsv.gz"
    path 'annotation_qc.txt'
    path 'clinical_vep.html'
    path 'clinical_vep.log'
    path 'clinical_prepare.txt'
    script:
    """
    set -euo pipefail
    test "\$(bcftools query -l '${annotated}')" = '${params.proband}'
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ\\t%AD]\\n' '${original}' > before.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ\\t%AD]\\n' '${annotated}' > after.tsv
    cmp before.tsv after.tsv
    bcftools view -h '${annotated}' > header.txt
    grep -q '^##INFO=<ID=CSQ' header.txt
    for field in ALLELE_NUM MANE_SELECT CLIN_SIG ClinVar; do grep -q "\$field" header.txt; done
    total=\$(bcftools index -n '${annotated}')
    test "\$total" -gt 0
    missing=\$(bcftools view -H -i 'INFO/CSQ="."' '${annotated}' | wc -l)
    { echo "records=\$total"; echo "missing_CSQ=\$missing"; echo 'alleles_genotypes_filters=identical'; echo 'inheritance_columns=CHROM POS REF ALT FILTER ${params.proband}_GT DP GQ AD ${params.mother}_GT DP GQ AD ${params.father}_GT DP GQ AD'; } > annotation_qc.txt
    echo "Records without CSQ: \$missing"
    cp -L '${annotated}' '${params.sample_name}.clinical.GRCh38.vcf.gz'
    cp -L '${annotated_idx}' '${params.sample_name}.clinical.GRCh38.vcf.gz.tbi'
    bcftools view -f PASS -Oz -o '${params.sample_name}.clinical.PASS.vcf.gz' '${annotated}'
    bcftools index -t '${params.sample_name}.clinical.PASS.vcf.gz'
    cp -L '${inheritance}' '${params.sample_name}.inheritance.tsv.gz'
    cp -L '${original_evidence}' '${params.sample_name}.original_trio_evidence.tsv.gz'
    cp -L '${stats}' clinical_vep.html
    cp -L '${vep_log}' clinical_vep.log
    cp -L '${prep_report}' clinical_prepare.txt
    """
}

workflow {
    for (key in ['vcf','vcf_index','fasta','fai','cache','clinvar','clinvar_index','proband','mother','father','sample_name']) {
        if (!params[key]) error "Missing --${key}"
    }
    if ([params.proband,params.mother,params.father].toSet().size() != 3) error 'Sample IDs must be distinct'
    PREPARE(file(params.vcf, checkIfExists:true), file(params.vcf_index, checkIfExists:true), file(params.fasta, checkIfExists:true), file(params.fai, checkIfExists:true))
    ANNOTATE(PREPARE.out.vcf, PREPARE.out.index, file(params.fasta, checkIfExists:true), file(params.fai, checkIfExists:true), file(params.cache, checkIfExists:true), file(params.clinvar, checkIfExists:true), file(params.clinvar_index, checkIfExists:true))
    QC(PREPARE.out.vcf, ANNOTATE.out.vcf, ANNOTATE.out.index, PREPARE.out.inheritance, PREPARE.out.original_evidence, PREPARE.out.report, ANNOTATE.out.stats, ANNOTATE.out.log)
}
