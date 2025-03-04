rule make_foldseekdb:
    input:
        structs = config['input_dir'],
        done = config["outdir"] + "/done.txt"
    output: "{output_dir}/foldseek/db/all_seqs_fsdb"
    threads: 8
    conda: "../envs/sp_homology.yaml"
    shell: '''
    foldseek createdb --threads {threads} {input.structs} {output}
    '''

rule make_foldseekdb_seq:
    input: rules.make_foldseekdb.output
    output: "{output_dir}/foldseek/db/all_seqs_fsdb_ss.fa"
    conda: "../envs/sp_homology.yaml"
    shell: '''
    foldseek lndb {input}_h {input}_ss_h
    foldseek convert2fasta {input}_ss {output}
    sed -i 's/-model_v4.cif//g' {output}
    '''

# rule filter_foldseek:
#     input: "{output_dir}/foldseek/db/all_seqs_fsdb"
#     output: "{output_dir}/foldseek/results/filtered_homologs.tsv"
#     params:
#         eval_both=config['eval_both'],
#         coverage=config['coverage'],
#         max_seqs=config['max_seqs']
#     shell: '''
#     awk '$11<{params.eval_both} && $17*100>{params.coverage} && $18*100>{params.coverage}' {input} > {output}
#     '''

# rule foldseek:
#     input:
#         q=config['input_dir'],
#         db=rules.make_foldseekdb.output
#     output:
#         "{output_dir}/foldseek/results/foldseek.tsv"
#     params: config['target_seqs']
#     log: "{output_dir}/foldseek/logs/foldseek.log"
#     benchmark: "{output_dir}/foldseek/benchmarks/foldseek.txt"
#     threads: 24
#     conda: "../envs/sp_homology.yaml"
#     shell: '''
#     foldseek easy-search {input.q} {input.db} {output} $TMPDIR/foldseek --threads {threads} --max-seqs {params} \
#     --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
#     sed -i 's/-model_v4.cif//g' {output}
#     '''

rule foldseek_allvall:
    input:
        q=config['input_dir'],
        t=config['input_dir']
    output:
        q_t="{output_dir}/foldseek/results/foldseek_allvall_qt.tsv",
        t_q="{output_dir}/foldseek/results/foldseek_allvall_tq.tsv"
    params: config['max_seqs_brh']
    log: "{output_dir}/foldseek/logs/foldseek_allvall.log"
    benchmark: "{output_dir}/foldseek/benchmarks/foldseek_allvall.txt"
    threads: 4
    conda: "../envs/sp_homology.yaml"
    shell: '''
    foldseek easy-search {input.q} {input.t} {output.q_t} $TMPDIR/foldseek_allvall_qt --threads {threads} --max-seqs {params} \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
    sed -i 's/-model_v4.cif//g' {output.q_t}

    foldseek easy-search {input.t} {input.q} {output.t_q} $TMPDIR/foldseek_allvall_tq --threads {threads} --max-seqs {params} \
    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
    sed -i 's/-model_v4.cif//g' {output.t_q}
    '''

rule foldseek_brh:
    input: "{output_dir}/foldseek/results/foldseek.tsv"
    output: "{output_dir}/foldseek/results/foldseek_brh.tsv"
    params: config['eval_brh']
    conda: "../envs/sp_python.yaml"
    script: "../scripts/get_BRH.py"


