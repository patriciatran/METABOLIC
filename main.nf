include { ANNOTATE_MAGS } from './modules/local/prodigal.nf'

workflow {
    // 1. Create a channel of all fasta files
    // The .map creates the 'meta' object for naming
    ch_inputs = Channel
        .fromPath("${params.input_genome_folder}/*.fasta")
        .map { file -> [ [id: file.baseName], file ] }

    // 2. Point to the accessory script in your bin/ or repo folder
    ch_gff_script = file("${baseDir}/bin/Accessory_scripts/gff2fasta_mdf.pl")

    // 3. Run the annotation
    // Nextflow will automatically parallelize this across all CPUs available
    ANNOTATE_MAGS(ch_inputs, ch_gff_script)

    // 4. Access the results for the next step (e.g., HMM searching)
    // ANNOTATE_MAGS.out.faa.view() 
}
