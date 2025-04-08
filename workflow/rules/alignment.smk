# Get aa sequences from the pdb structures
rule get_seqs_aa:
    input:
        pdb_files=f"{input_dir}",
        done="{output_dir}" + "/done.txt"
    output: "{output_dir}/alignment/processed_sequences.fa"
    script: "../scripts/structs2fasta.py"

# Run foldmason to align all pdbs found in the input folder and create an aa and 3di alignment, as well as a html report
# could potentially refine the alignment with e.g. --refine-iters 1000
rule foldmason:
    input:
        db=f"{input_dir}",
        done="{output_dir}" + "/done.txt"
    output:
        aa_alignment = "{output_dir}/alignment/foldmason/alignment_aa.fa",
        structure_alignment = "{output_dir}/alignment/foldmason/alignment_3di.fa"
    params: min_plddt = config['min_lddt']
    threads: 8
    shell: '''
    foldmason easy-msa {input.db} {output_dir}/alignment/foldmason/alignment /tmp --report-mode 1 \
    --threads {threads} --mask-bfactor-threshold {params.min_plddt}
    '''

# Trim the alignment using clipkit with the -m gappy option which is fairly conservative overall
rule trim_aln_clipkit:
    input: "{output_dir}/alignment/foldmason/alignment_{alphabet}.fa"
    output: "{output_dir}/alignment/foldmason/alignment_{alphabet}_trimmed.fa"
    wildcard_constraints:
        alphabet="aa|3di"
    shell: '''
    clipkit {input} -o {output}.tmp -m gappy -l 
    
    cat {output}.tmp | seqtk seq -A | awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
	
    rm {output}.tmp 
    '''
# a small workaround here so i can get the logfile

# Also trim with trimal which will be used to create the IMD distance matrix
rule trim_aln_trimal:
    input: "{output_dir}/alignment/foldmason/alignment_{alphabet}.fa"
    output: 
    	msa = "{output_dir}/alignment/trimal/alignment_{alphabet}_trimmed.fa",
    	report = "{output_dir}/alignment/trimal/trimal_{alphabet}_report.txt"
    wildcard_constraints:
        alphabet="aa|3di"
    shell: '''
    trimal -in {input} -out {output.msa} -gappyout > {output.report} 2>&1
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


