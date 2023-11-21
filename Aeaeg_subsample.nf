#!/usr/bin/env nextflow

// Params from config files (system-dependent)
input=params.input
output=params.output
work=params.work
aux=params.aux

huge=params.huge
big=params.big
small=params.small

// Global Params
params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

params.rlen = null
if( !params.rlen ) error "Missing length (average read length) parameter"
println "rlen: $params.rlen"


////////////////////////////////////////////////
// ** - Pull in fq files
////////////////////////////////////////////////

Channel.fromFilePairs(input + "/${params.dir}/*_{1,2}.fq.gz", flat: true)
          .set { fqs }


////////////////////////////////////////////////
// ** - Subsample reads
////////////////////////////////////////////////

process sample_reads {

  publishDir "${output}/${params.dir}_sub", mode: 'copy'

  cpus small
  tag { id }

  input:
    tuple val(id), file(forward), file(reverse) from fqs

  output:
    tuple id, file("${id}_1.fq.gz"), file("${id}_2.fq.gz") into subsampled_fqs

  """
    seqtk sample -s 10 $forward 10000 > ${id}_1.fq.gz
    seqtk sample -s 10 $reverse 10000 > ${id}_2.fq.gz
  """

}