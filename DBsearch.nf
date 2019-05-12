/*
vim: syntax=groovy
-*- mode: groovy;-*-

==============================
a pipeline to search MS data in a protein database
==============================

*/

nf_required_version = '0.26.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


/* SET DEFAULT PARAMS */

params.outdir = 'result'

protdb = file(params.protdb)
params.MS1tolerance = '10ppm'
params.FragmentMethodID = 0 // 0 use the FragMethod as written in spectra,  3 for HCD, 1 for CID, 2 for ETD
params.inst = 1 // 0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR/Lumos, 2: TOF, 3: Q-Exactive
params.protocol = 0
params.qval = 0.01 // q-value cutoff

config = file(params.config) // instead of giving all parameters through command-line, use one config files for all search parameters


///////////////////

println "PrecursorMassTolerance used in MSGF+:${params.MS1tolerance}"
println "q-value cutoff:${params.qval}"
println "InstrumentID used in MSGF+: ${params.inst}"
println "FragmentMethodID used in MSGF+: ${params.FragmentMethodID}"
println "ProtocolID used in MSGF+: ${params.protocol}"


/* PIPELINE START */

// feed an mzmldef file (tab separated lines with absolute filepath and setname)
// If one MS experiment contains multiple MS spectra files, use the same setname
Channel
  .from(file("${params.mzmldef}").readLines())
  .map { it -> it.tokenize('\t') }
  .set { mzml_in }


mzml_in
  .tap { sets }
  .map { it -> [ it[1], (it[0]=~/.*\/(.+)\..*$/)[0][1], file(it[0])]} // extract set, file name and file object from input
  .tap{ mzml_msgf } // create channels for downstream process
  .count()
  .set{ amount_mzml } //count the number of MS spectra file per MS experiment


process DatabaseSearch {

  input:
  set val(setname), val(sample), file(x) from mzml_msgf

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  
  """
  msgf_plus -Xmx12G -thread 1 -d $protdb -s $x -o "${sample}.mzid" -mod $mods -tda 0 -t ${params.MS1tolerance} -ti -1,2 -m ${params.FragmentMethodID} -inst ${params.inst} -e 1 -protocol ${params.protocol} -ntt 2 -minLength 8 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1
  """
}

process ConvertmzidTotsv {

  input:
  set val(setname), file(x) from mzids
  
  output:
  set val(setname), file('out.mzid.tsv') into mzidtsvs 

  """
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i $x -o out.mzid.tsv
  """
}

mzidtsvs
  .groupTuple()
  .set {mzidtsvs_byset}

process MergeTSVbyset {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), file('tsv?') from mzidtsvs_byset

  output:
  set val(setname), file("${setname}_allpsms.sorted.txt") into merged_tsv

  script:
  """
  head -1 tsv1 |cut -f 1-14 > psmheader
  sort -g -s -t \$'\\t' -k 13,13 <(tail -q -n +2 tsv*) |cut -f 1-14 >psmmerge.sorted
  cat psmheader psmmerge.sorted > ${setname}_allpsms.sorted.txt
  """
}

merged_tsv
  .tap{ globalFDRtsv; subgroupFDRtsv }


process GlobalFDR {
 publishDir "${params.outdir}", mode: 'copy', overwrite: true

 input:
 set val(setname), file('sorted.psm') from globalFDRtsv
 
 output:
 set val(setname), file("${setname}_psms.FDR0.01.txt") into all_psm
 
 script:
 """
 GlobalFDR.py --input sorted.psm --output ${setname}_psms.FDR0.01.txt --decoy_prefix XXX_ --psm_qval 0.01
 """

}

process subGroupFDR {
 publishDir "${params.outdir}", mode: 'copy', overwrite: true

 input:
 set val(setname), file('sorted.psm') from subgroupFDRtsv

 output:
 set val(setname), file("${setname}_psms.txt") into novel_psm

 """
 subgroupFDR.py --input sorted.psm --group_tag ${group_tag} --output ${setname}_${group_tag}.txt --psm_qval ${params.qval}
 """
}


process ProteinInference {
  

}


Process ProteinFDR {



}




