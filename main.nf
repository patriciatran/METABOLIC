include { ANNOTATE_MAGS } from './modules/local/prodigal.nf'
include { RUN_HMMSEARCH } from './modules/local/hmmer.nf'
include { SETUP_DBCAN } from './modules/local/setup_dbcan.nf'
include { RUN_DBCAN } from './modules/local/run_dbcan.nf'
include { HMM_MOTIF_VALIDATE } from './modules/local/hmm_motif_validation.nf'
include { MAKE_RESULT_TABLES } from './modules/local/make_result_tables.nf'
include { MAKE_RESULT_TABLE2 } from './modules/local/make_result_table2.nf'
include { MAKE_RESULT_TABLE3 } from './modules/local/make_result_table3.nf'
include { MAKE_RESULT_TABLE4 } from './modules/local/make_result_table4.nf'
include { MAKE_KEGG_IDENTIFIER_RESULTS } from './modules/local/make_kegg_identifier_results.nf'
include { MAKE_HMM_SEQUENCE_COLLECTIONS } from './modules/local/make_hmm_sequence_collections.nf'
include { MAKE_RESULT_TABLE5 } from './modules/local/make_result_table5.nf'
include { SETUP_MEROPS } from './modules/local/setup_merops.nf'
include { RUN_MEROPS } from './modules/local/run_merops.nf'
include { MAKE_RESULT_TABLE6 } from './modules/local/make_result_table6.nf'
include { MAKE_EXCEL_WORKBOOK } from './modules/local/make_excel_workbook.nf'
include { DRAW_ELEMENT_CYCLES } from './modules/local/draw_element_cycles.nf'


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

    // 4b. Set up dbCAN database, then run dbCAN per genome FAA
    SETUP_DBCAN()
    RUN_DBCAN(ANNOTATE_MAGS.out.faa, SETUP_DBCAN.out.hmm_db)

    // 4c. Set up MEROPS database, then run MEROPS per genome FAA
    SETUP_MEROPS()
    RUN_MEROPS(ANNOTATE_MAGS.out.faa, SETUP_MEROPS.out.merops_db, SETUP_MEROPS.out.merops_lib)

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

    //6. Run the HMM motif validation step using hmmsearch outputs
    HMM_MOTIF_VALIDATE(RUN_HMMSEARCH.out.results)

    def ch_motifcheck_all = HMM_MOTIF_VALIDATE.out.collect()
    def ch_faa_all = ANNOTATE_MAGS.out.faa.map { meta, faa -> faa }.collect()
    def ch_gene_all = ANNOTATE_MAGS.out.genes.map { meta, gene -> gene }.collect()

    // 7. Aggregate all motifcheck TSVs + FAA files -> worksheet 1
    MAKE_RESULT_TABLES(
        ch_motifcheck_all,
        ch_faa_all,
        file(params.hmm_table_template)
    )

    // 8. Build worksheet 2 (Function presence)
    MAKE_RESULT_TABLE2(
        ch_motifcheck_all,
        ch_faa_all,
        file(params.hmm_table_template),
        file(params.hmm_table_template_2)
    )

    // 9. Build per-HMM protein/gene sequence collections
    MAKE_HMM_SEQUENCE_COLLECTIONS(
        ch_motifcheck_all,
        ch_faa_all,
        ch_gene_all
    )

    // 10. Build worksheet 3 (KEGG module presence)
    MAKE_RESULT_TABLE3(
        ch_motifcheck_all,
        ch_faa_all,
        file(params.hmm_table_template), 
        file(params.ko_module_table), 
        file(params.ko_module_step_db), 
        params.m_cutoff 
    )

    // 11. Build worksheet 4 (KEGG module-step presence)
    MAKE_RESULT_TABLE4(
        ch_motifcheck_all,
        ch_faa_all,
        file(params.hmm_table_template),
        file(params.ko_module_table),
        file(params.ko_module_step_db)
    )

    // 12. Build per-genome KEGG identifier results
    MAKE_KEGG_IDENTIFIER_RESULTS(
        ch_motifcheck_all,
        ch_faa_all,
        file(params.hmm_table_template)
    )

    // 13. Build worksheet 5 (dbCAN results)
    MAKE_RESULT_TABLE5(
        RUN_DBCAN.out.results.map { meta, domtbl, out, ps -> ps }.collect()
    )

    // 14. Build worksheet 6 (MEROPS results)
    MAKE_RESULT_TABLE6(
        RUN_MEROPS.out.parsed.map { meta, parsed -> parsed }.collect()
    )

    // 15. Build final Excel workbook from worksheets 1-6
    MAKE_EXCEL_WORKBOOK(
        MAKE_RESULT_TABLES.out,
        MAKE_RESULT_TABLE2.out,
        MAKE_RESULT_TABLE3.out,
        MAKE_RESULT_TABLE4.out,
        MAKE_RESULT_TABLE5.out.worksheet5,
        MAKE_RESULT_TABLE6.out.worksheet6,
        file("${baseDir}/bin/R/create_excel_spreadsheet.R")
    )

    // 16. Draw element cycling diagrams
    DRAW_ELEMENT_CYCLES(
        ch_motifcheck_all,
        ch_faa_all,
        file(params.r_pathways),
        file("${baseDir}/bin/R/draw_biogeochemical_cycles.R")
    )
}