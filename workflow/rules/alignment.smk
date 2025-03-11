# Get aa sequences from the pdb structures
rule get_seqs_aa:
    input:
        pdb_files=f"{input_dir}",
        done="{output_dir}" + "/done.txt"
    output: "{output_dir}/alignment/processed_sequences.fa"
    script: "../scripts/structs2fasta.py"

# Instead of mafft we use foldmason for both the aa and 3di alignment
# rule aln_aa:
#     input: rules.get_seqs_aa.output
#     output: "{output_dir}/alignment/aln_aa/alignment.alg"
#     log: "{output_dir}/alignment/aln_aa/logs/mafft_alignment.log"
#     benchmark: "{output_dir}/alignment/aln_aa/benchmarks/mafft_alignment.txt"
#     threads: 8
#     shell:'''
#     mafft --localpair --maxiterate 1000 --thread {threads} {input} > {output} 2> {log}
#     '''
#     # shell:'''
#     # mafft --auto --thread {threads} {input} > {output} 2> {log}
#     # '''

rule get_seqs_3di:
    input:
        ids="{output_dir}/foldseek/foldseek_ids.tsv",
        fa="{output_dir}/foldseek/db/all_seqs_fsdb_ss.fa"
    output: "{output_dir}/alignment/get_seqs_3di/selected_3di_sequences.fa"
    shell: '''
    seqkit grep -f {input.ids} {input.fa} > {output}
    '''

rule mask_seqs_3di:
    input: "{output_dir}/alignment/aln_aa/alignment.alg",
    output: "{output_dir}/alignment/mask_seqs_3di/masked_sequences.alg"
    params:
        structdir={input_dir},
        min_lddt=config['min_lddt']
    script: "../scripts/mask_structures.py"

# Run foldmason to align all pdbs found in the input folder and create an aa and 3di alignment, as well as a html report
# could potentially refine the alignment with e.g. --refine-iters 1000
rule foldmason:
    input:
        db=f"{input_dir}",
        done="{output_dir}" + "/done.txt"
    output:
        aa_alignment = "{output_dir}/alignment/foldmason/alignment_aa.fa",
        structure_alignment = "{output_dir}/alignment/foldmason/alignment_3di.fa"
    log: "{output_dir}/alignment/foldmason/logs/foldmason_3di.log"
    benchmark: "{output_dir}/alignment/foldmason/benchmarks/foldmason_3di.txt"
    threads: 8
    shell: '''
    foldmason easy-msa {input.db} {output_dir}/alignment/foldmason/alignment /tmp --report-mode 1 --threads {threads}
    '''

# Trim the alignment using clipkit with the -m gappy option which is fairly conservative overall
rule trim_aln:
    input: "{output_dir}/alignment/foldmason/alignment_{alphabet}.fa"
    output: "{output_dir}/alignment/foldmason/alignment_{alphabet}_trimmed.fa"
    wildcard_constraints:
        alphabet="aa|3di"
    shell: '''
    # Step 1: First Clipkit run to generate the log file and trimmed alignment
    clipkit {input} -o {output}.tmp -m gappy -l 

    # Step 2: Read trimmed alignment from file and process it
    cat {output}.tmp | seqtk seq -A | awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
	
    rm {output}.tmp 
    '''

# concatenate the aa and 3di alignment
rule concat_aln:
    input:
        aa="{output_dir}/alignment/foldmason/alignment_aa_trimmed.fa",
        msa_3di="{output_dir}/alignment/foldmason/alignment_3di_trimmed.fa"
    output:
        "{output_dir}/alignment/concat_aln/concatenated_alignment.alg"
    shell:"""
    seqkit concat {input.aa} {input.msa_3di} > {output}
    """


