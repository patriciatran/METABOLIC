include { ANNOTATE_MAGS } from './modules/local/prodigal.nf'
include { RUN_HMMSEARCH } from './modules/local/hmmer.nf'
//include { HMM_MOTIF_VALIDATE } from './modules/local/hmmotifvalidation.nf'


workflow {
    // 0. Load KO thresholds from ko_list (Kxxxxx -> threshold)
    // Prefer CLI --kofam_list_path when provided, then fallback to --ko_list_path/default
    def ko_list_input = params.kofam_list_path ?: params.ko_list_path
    if( !ko_list_input ) {
        error "Missing KO list path. Provide --ko_list_path or --kofam_list_path"
    }

    def ko_list_file = file(ko_list_input)
    if( ko_list_file.isDirectory() ) {
        ko_list_file = file("${ko_list_file}/ko_list")
    }
    if( !ko_list_file.exists() ) {
        error "KO list file not found: ${ko_list_file}"
    }

    def ko_threshold = [:]
    ko_list_file.eachLine { line ->
        if (!line || line.startsWith('#')) return

        def cols = line.split(/\t/)
        if (cols.size() < 2) return

        def ko_id = cols[0].trim()
        def raw = cols[1].trim()
        def threshold = (raw && raw != '-') ? raw : '50'

        // Support lookup by both "K00001" and "K00001.hmm"
        ko_threshold[ko_id] = threshold
        ko_threshold["${ko_id}.hmm"] = threshold
    }

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
    ch_hmms = Channel
        .fromPath("${params.kofam_db_path}/*.hmm")
        .map { hmm_file ->
            def hmm_name = hmm_file.baseName
            def threshold = ko_threshold[hmm_file.name] ?: ko_threshold[hmm_name] ?: '50'
            [hmm_file, hmm_name, threshold]
        }
    ch_faa = Channel.fromPath("${params.outdir}/annotations/*.faa")

    // 5. Run HMM search for every genome and each HMM file (this is heavy by design)
    ch_faa_hmm = ch_faa.combine(ch_hmms).map { row ->
        def faa_file = row[0]
        def hmm_file = row[1]
        def hmm_name = row[2]             // e.g., K00001
        def threshold = row[3]
        def meta_id = faa_file.baseName   // e.g., Ga0485169_metabat2_ours.127_sub
        [meta_id, faa_file, hmm_file, hmm_name, threshold]
    }
    
    RUN_HMMSEARCH(ch_faa_hmm)
}