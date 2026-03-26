include { ANNOTATE_MAGS } from './modules/local/prodigal.nf'
include { RUN_HMMSEARCH } from './modules/local/hmmer.nf'
//include { HMM_MOTIF_VALIDATE } from './modules/local/hmmotifvalidation.nf'


workflow {
    // 1. Create a channel of all fasta files
    // The .map creates the 'meta' object for naming
    ch_inputs = Channel
        .fromPath("${params.input_genome_folder}/*.fasta")
        .map { file -> [ [id: file.baseName], file ] }

    // 2. Point to the accessory script in your bin/ or repo folder
    ch_gff_script = file("${baseDir}/bin/Accessory_scripts/gff2fasta_mdf.pl")

    // 3. Run the annotation
    ANNOTATE_MAGS(ch_inputs, ch_gff_script)

    // 4. Prepare HMM database sources
    ch_hmms = Channel.fromPath("${params.kofam_db_path}/*.hmm")
    ch_faa = Channel.fromPath("results/annotations/*.faa")

    // 5. Run HMM search for every genome and each HMM file (this is heavy by design)
    ch_faa_hmm = ch_faa.combine(ch_hmms).map { faa_file, hmm_file ->
        def meta_id = faa_file.baseName   // e.g., Ga0485169_metabat2_ours.127_sub
        def hmm_name = hmm_file.baseName  // e.g., K00001
        [meta_id, faa_file, hmm_file, hmm_name]
    }
    
    RUN_HMMSEARCH(ch_faa_hmm)
}