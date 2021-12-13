#!/usr/bin/env nextflow

/*
================================================================================
                                Workflow Help
================================================================================
*/

// Re-usable component for adding a helpful help message in our Nextflow script

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    
    nextflow run main.nf --alnformat bam --ref /path/to/reference.fa
    
    Mandatory arguments:
      --alnformat     [enum] Alignment file format ('bam', 'cram')
      --ref           [file] Reference FASTA

    Optional arguments:
      --wgs           [int] indicate the memory factor for WGS
      --help          [flag] Show help messages

    """.stripIndent()
}

/* Unused arguments
  --index         [file] Index in tab format (index_tab.txt)
  --binlist       [file] A txt file with paths to all bin files (one per row)
  --index         [file] Index in tab format (index_tab.txt)
  --bindir        [path] Path to the directory of all bin files
  --gender        [file] Gender file from part 2
  --batch         [int]  Batch size for references
  --samples       [file] Samples to process
  --rbindir
  --cordir
  --index
  --gender
  --cov
 */

// Show help message
if (params.help) exit 0, helpMessage()

batch = 1000
if (params.batch) {
  batch = params.batch
}

mem_factor = 1
if (params.wgs) {
  mem_factor = params.wgs
}

/*
================================================================================
                                Pre-process Setup
================================================================================
*/

// Pre-process: create design.csv file and create map

design_header = "name,bam,bai"
if (params.alnformat == "cram") {
  design_header = "name,cram,crai"
}

i = 0
iterate = true
design_body = ""
while(iterate) {
  name_param = "name_" + i
  aln_param = "aln_" + i
  idx_param = "idx_" + i

  if (params.containsKey(name_param)) {
    if (params.containsKey(aln_param) && params.containsKey(idx_param)) {
      design_body = design_body.concat(params[name_param] + "," + params[aln_param] + "," + params[idx_param] + "\n")
      
    } else {
      iterate = false
      print("Invalid parameter set: " + aln_param + " and/or " + idx_param + " not specified, but " + name_param + " was specified, exiting")
      exit 1, helpMessage()
    }

  } else {
    iterate = false
    if (i == 0) {
      print("No inputs to populate design file, exiting")
      exit 1, helpMessage()
    }
  }

  i++
}
File design_file = new File("design.csv")
design_file.write design_header + "\n"
design_file << design_body

Channel.fromPath('design.csv')
  .splitCsv(sep: ',', skip: 1)
  .map { name, file_path, index_path -> [ name, file(file_path), file(index_path) ] }
  .set { ch_files_sets }

// Pre-process: channels for common handled files

if (params.bed) ch_bed = Channel.value(file(params.bed))
if (params.ref) ch_ref = Channel.value(file(params.ref))

/*
================================================================================
                                Main Processes
================================================================================
*/

ch_bedgz = Channel.value(file("$baseDir/data/hg38.1kb.baits.bed.gz"))

if (!params.bed) {
  process step0 {
    input:
    file(bedgz) from ch_bedgz

    output:
    file("hg38.1kb.baits.bed") into ch_bed

    when:
    !params.bed

    script:
    if (params.test)
      """
      gzip -cd ${bedgz} | head -1000 > "hg38.1kb.baits.bed"
      """
    else
      """
      gzip -cd ${bedgz} > "hg38.1kb.baits.bed"
      """
  }
}

process make_subdirs {
  
  output:
  path "indexes" into ch_index_dir
  path "bin" into ch_bin_dir

  script:
  """
  mkdir indexes
  mkdir bin
  """
}

// Step1 create work directory
process step1 {

  input: 
  path ch_index_dir
  file(bed) from ch_bed

  output: 
  path "indexes/index_tab.txt" into ch_index_tab
  path "indexes/index.txt" into ch_index
  path "indexes/index.bed" into ch_index_bed

  script:
  """
  cnest.py step1 --project indexes --bed ${bed}
  """
}

process step2 {
  tag "id:${name}-file:${file_path}-index:${index_path}"

  input:
  set val(name), file(file_path), file(index_path) from ch_files_sets
  file("genome.fa") from ch_ref
  path "project/bin" from ch_bin_dir
  path "project/index.bed" from ch_index_bed

  output:
  path "project/bin/$name" into bin_paths 
  val(name) into ch_names_sets

  script:
  """
  mkdir -p project/tmp/
  cnest.py step2 --project project --sample ${name} --input ${file_path} --fasta genome.fa --fast
  """
}

process gender_qc {
  time '10h'

  input:
  val a from bin_paths.collect()
  path ch_index_tab
  path ch_bin_dir

  output:
  path "gender_qc.txt"
  path "gender_classification.txt" into ch_gender_classification
  path "mean_coverage.txt"

  script:
  """
  cnest.py step3 \
    --indextab $ch_index_tab \
    --bindir $ch_bin_dir \
    --qc gender_qc.txt \
    --gender gender_classification.txt \
    --cov mean_coverage.txt
  """
}


process logR_ratio {
  tag "${name}"
  echo true
  // memory { 1.GB * params.batch * mem_factor / 100 }
  // time { 40.m * params.batch * mem_factor / 100  }

  input:
  path ch_bin_dir
  path ch_index_tab
  path ch_gender_classification
  val(name) from ch_names_sets
  // path bin_dir from ch_bin
  // path index from ch_index
  // path gender from ch_gender
  // val sample_name from ch_sample_names

  output:
  path "project/cor/$name"
  path "project/logr/$name"
  path "project/rbin/$name"

  script:
  """
    mkdir -p project/cor/ project/logr/ project/rbin/
    cnest.py step4 \
      --bindir $ch_bin_dir \
      --indextab $ch_index_tab \
      --gender $ch_gender_classification \
      --sample $name \
      --batch $batch \
      --cordir project/cor/ \
      --logrdir project/logr/ \
      --rbindir project/rbin/
  """
}

/*
if (params.part == 4){
  process hmm_call {
    tag "${sample_name}"
    echo true
    publishDir "results/", mode: "move"
    memory { 5.GB * params.batch * mem_factor / 100 }
    time { 40.m * params.batch * mem_factor / 100  }

    input:
    path rbin_dir from ch_rbin
    path cor_dir from ch_cor
    path index from ch_index
    path gender_file from ch_gender
    path cov_file from ch_cov
    val sample_name from ch_sample_names

    output:
    path "${params.project}/cnv/${sample_name}/${sample_name}_mixed_calls.txt"
    path "${params.project}/cnv/${sample_name}/${sample_name}_mixed_states.txt"

    script:
    """
      echo "Processing sample $sample_name"
      mkdir -p ${params.project}/cnv/
      cnest.py step5 \
        --indextab $index \
        --rbindir $rbin_dir \
        --cordir $cor_dir \
        --cnvdir ${params.project}/cnv/ \
        --cov    $cov_file \
        --sample $sample_name \
        --gender $gender_file \
        --batch $params.batch
    """
  }
}
*/
