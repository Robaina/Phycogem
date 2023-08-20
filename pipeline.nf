#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.genomes_dir = './input_mags'
params.media_file = './marine_media/media_db.tsv'
params.outdir = './results'

mags = Channel.fromPath("${params.genomes_dir}/*.fasta")

process carveme {
    publishDir "${params.outdir}/gems", mode: 'copy'

    input:
    path fasta 

    output:
    path("${fasta.baseName}.xml")

    script:
    """
    carve --solver gurobi \\
     -o ${fasta.baseName}.xml \\
     -u bacteria \\
     --init M9[marine] \\
     --gapfill M9[marine] \\
     --mediadb ${params.media_file} \\
     ${fasta}
    """
}

process memote {
    publishDir "${params.outdir}/memote_reports", mode: 'copy'

    input:
    path gem_file 

    output:
    path "${gem_file.baseName}.json"

    script:
    """
    memote run --filename ${gem_file.baseName}.json --ignore-git $gem_file
    """
}

process remove_duplicated_reactions

process merge_community {
    publishDir "${params.outdir}/merged_community", mode: 'copy'

    input:
    path xml_files from gems_ch

    output:
    path "merged.xml"

    script:
    """
    carve --solver gurobi \\
     --init M9[marine] \\
     --gapfill M9[marine] \\
     --mediadb ${params.media_file} \\
     --output merged.xml \\
     --fbc2 \\
     $xml_files
    """
}

workflow {
    mags.view().set { mags_ch }
    carveme( mags_ch ).view().set { gems_ch }
    merge_community( gems_ch )
    memote( gems_ch )
}
