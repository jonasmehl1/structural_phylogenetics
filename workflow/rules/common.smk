rule get_seqs_aa:
    input:
        pdb_files=config['input_dir'],
        done=config["outdir"] + "/done.txt"
    output: "{output_dir}/alignment/processed_sequences.fa"
    conda: "../envs/sp_python.yaml"
    script: "../scripts/structs2fasta.py"

# rule aln_aa:
#     input: rules.get_seqs_aa.output
#     output: "{output_dir}/alignment/aln_aa/alignment.alg"
#     log: "{output_dir}/alignment/aln_aa/logs/mafft_alignment.log"
#     benchmark: "{output_dir}/alignment/aln_aa/benchmarks/mafft_alignment.txt"
#     threads: 8
#     conda: "../envs/sp_tree.yaml"
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
    conda: "../envs/sp_utils.yaml"
    shell: '''
    seqkit grep -f {input.ids} {input.fa} > {output}
    '''

rule mask_seqs_3di:
    input: "{output_dir}/alignment/aln_aa/alignment.alg",
    output: "{output_dir}/alignment/mask_seqs_3di/masked_sequences.alg"
    params:
        structdir=config['input_dir'],
        min_lddt=config['min_lddt']
    conda: "../envs/sp_python.yaml"
    script: "../scripts/mask_structures.py"


rule foldmason:
    input:
        db=config["input_dir"],
        done=config["outdir"] + "/done.txt"
    output:
        aa_alignment = "{output_dir}/foldmason/alignment_aa.fa",
        structure_alignment = "{output_dir}/foldmason/alignment_3di.fa"
    log: "{output_dir}/alignment/aln_3di/logs/foldmason_3di.log"
    benchmark: "{output_dir}/alignment/aln_3di/benchmarks/foldmason_3di.txt"
    threads: 8
    conda: "../envs/sp_tree.yaml"
    shell: '''
    foldmason easy-msa {input.db} {output_dir}/foldmason/alignment /tmp --report-mode 1 --threads {threads}
    '''
    #--refine-iters 10

# Actually we dont want to trim our alignments
# rule trim_aln:
#     input: "{output_dir}/alignment/aln_{alphabet}/alignment.alg"
#     output: "{output_dir}/alignment/aln_{alphabet}/alignment.alg.clean"
#     wildcard_constraints:
#         alphabet="aa|3di"
#     conda: "../envs/sp_tree.yaml"
#     shell: '''
#     trimal -in {input} -out /dev/stdout -gappyout | seqtk seq -A | \
#     awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
#     '''

rule concat_aln:
    input:
        aa="{output_dir}/foldmason/alignment_aa.fa",
        msa_3di="{output_dir}/foldmason/alignment_3di.fa"
    output:
        "{output_dir}/alignment/concat_aln/concatenated_alignment.alg"
    shell:"""
    seqkit concat {input.aa} {input.msa_3di} > {output}
    """

rule iqtree:
    input: "{output_dir}/foldmason/alignment_{alphabet}.fa"
    output:
        tree="{output_dir}/iqtree/{alphabet}/{model}.nwk",
        iqtree_output="{output_dir}/iqtree/{alphabet}/{model}.iqtree"
    wildcard_constraints:
        alphabet="aa|3di",
        model="LG|3di|AF"
    params:
        threedi_submat=config['subst_matrix']['3di'],
        AF_submat=config['subst_matrix']['AF'],
        ufboot=config['UF_boot']
    log: "{output_dir}/iqtree/logs/{alphabet}_{model}.log"
    benchmark: "{output_dir}/iqtree/benchmarks/{alphabet}_{model}.txt"
    threads: 16
    conda: "../envs/sp_tree.yaml"
    shell: '''
    tree_prefix={output_dir}/iqtree/{wildcards.alphabet}/{wildcards.model}

    model={wildcards.model}
    if [ "$model" == "3di" ]; then
        model="3di -mdef {params.threedi_submat}"
    elif [ "$model" == "AF" ]; then
        model={params.AF_submat}
    fi

    iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads} --quiet \
    --mem 16G --cmin 4 --cmax 10 --mset $model

    mv $tree_prefix.treefile {output.tree}
    mv $tree_prefix.log {log}
    '''

rule get_part:
    input:
        LG="{output_dir}/iqtree/aa/LG.iqtree",
        ThreeDI="{output_dir}/iqtree/3di/{combs_3di}.iqtree"
    output: "{output_dir}/iqtree_partitioned/{combs_3di}_combined_partition.part"
    shell: """
    model_LG=$(grep "^Model of substitution" {input.LG} | cut -f2 -d ':')
    model_3di=$(grep "^Model of substitution" {input.ThreeDI} | cut -f2 -d ':')
    length_LG=$(grep "^Input data" {input.LG} | cut -f6 -d ' ')
    length_3di=$(grep "^Input data" {input.ThreeDI} | cut -f6 -d ' ')

    end=$(expr $length_LG + $length_3di)
    start_3di=$(expr $length_LG + 1)

    echo "$model_LG, aa_partition = 1-$length_LG" > {output}
    echo "$model_3di, 3di_partition = $start_3di-$end" >> {output}
    """

rule iqtree_partitioned:
    input:
        fa="{output_dir}/alignment/concat_aln/concatenated_alignment.alg",
        part="{output_dir}/iqtree_partitioned/{combs_3di}_combined_partition.part"
    output:
        tree="{output_dir}/iqtree_partitioned/LG_{combs_3di}.nwk",
        iqtree="{output_dir}/iqtree_partitioned/LG_{combs_3di}.iqtree"
    wildcard_constraints:
        combs_3di="3di|AF"
    params:
        threedi_submat=config['subst_matrix']['3di'],
        AF_submat=config['subst_matrix']['AF'],
        ufboot=config['UF_boot']
    log: "{output_dir}/iqtree_partitioned/logs/LG_{combs_3di}.log"
    benchmark: "{output_dir}/iqtree_partitioned/benchmarks/LG_{combs_3di}.txt"
    threads: 16
    conda: "../envs/sp_tree.yaml"
    shell: """
    tree_prefix={output_dir}/iqtree_partitioned/LG_{wildcards.combs_3di}

    iqtree2 -s {input.fa} -p {input.part} --prefix $tree_prefix -mdef {params.threedi_submat} -B {params.ufboot} -T {threads} --quiet 

    mv $tree_prefix.treefile {output.tree}
    mv $tree_prefix.log {log}
    rm -f $tree_prefix.model.gz $tree_prefix.splits.nex $tree_prefix.contree $tree_prefix.ckp.gz
    """

##### FOLDTREE #####

# rule foldmason_report:
#     input:
#         fa = "{output_dir}/alignment/aln_3di/alignment.alg",
#         db = "{output_dir}/foldseek/db/all_seqs_fsdb"
#     output:
#         fa=temp("{output_dir}/foldtree/foldmason_report/foldmason_alignment_3di.alg.tmp"),
#         html="{output_dir}/foldtree/foldmason_report/foldmason_alignment_3di.html"
#     conda: "../envs/sp_tree.yaml"
#     shell: '''
#         cat {input.fa} > {output.fa}
#         foldmason msa2lddtreport {input.db} {output.fa} {output.html}
#     '''

rule foldseek_allvall_tree:
    input:
        structs = config["input_dir"],
        done=config["outdir"] + "/done.txt"
    output: "{output_dir}/foldtree/foldseek_allvall/foldseek_allvall.txt"
    log: "{output_dir}/foldtree/foldseek_allvall/logs/foldseek.log"
    benchmark: "{output_dir}/foldtree/foldseek_allvall/benchmarks/foldseek.txt"
    conda: "../envs/sp_homology.yaml"
    shell:'''
    foldseek easy-search {input.structs} {input.structs} {output} $TMPDIR/foldseek_tmp \
    --format-output 'query,target,fident,lddt,alntmscore' --exhaustive-search -e inf \
    --alignment-type 2 > {log}
    '''

rule foldseek_distmat:
    input: "{output_dir}/foldtree/foldseek_allvall/foldseek_allvall.txt"
    output: "{output_dir}/foldtree/foldseek_distmat/foldseek_3di_fident.txt"
    conda: "../envs/sp_python.yaml"
    script: "../scripts/foldtree/foldseekres2distmat_simple.py"

rule foldtree:
    input: "{output_dir}/foldtree/foldseek_distmat/foldseek_3di_fident.txt"
    output: "{output_dir}/foldtree/foldtree_3di_FT.nwk"
    wildcard_constraints:
        alphabet="3di"
    benchmark: "{output_dir}/foldtree/benchmarks/foldtree.txt"
    conda: "../envs/sp_tree.yaml"
    shell: '''
        quicktree -i m {input} | paste -s -d '' > {output}
        '''

# rule foldtree_py:
#     input:
#         struct_db = config["input_dir"],
#         foldseek_results = rules.foldseek_allvall_tree.output
#     output:
#         distmat="{output_dir}/foldtree/foldtree_py/foldtree_3di_FTPY.txt",
#         tree="{output_dir}/foldtree/foldtree_py/foldtree_3di_FTPY.nwk"
#     params: outdir="{output_dir}"
#     log: "{output_dir}/foldtree/foldtree_py/logs/foldtree_py.log"
#     conda: "../envs/sp_python.yaml"
#     shell: '''
#         indir=$(dirname {input.struct_db})
#         foldseek_output=$(dirname {input.foldseek_results})
# 
#         python workflow/scripts/foldtree/foldtree.py -i {input.struct_db} -o {params.outdir}/foldtree/foldtree_py_output \
#         -t $TMPDIR --outtree {output.tree} -c {params.outdir}/foldtree/foldtree_py_core \
#         --corecut --correction --kernel fident > {log}
# 
#         # Cleanup
#         rm -rf {params.outdir}/foldtree/foldtree_py_core
#         rm -f {output.distmat}_fastme_stat.txt {output.distmat}.tmp
#         rm -f {params.outdir}/foldtree/foldtree_py_allvall.tsv
#         rm -f {params.outdir}/foldtree/foldtree_py_core_allvall.tsv
#     '''


#Run FASTME on the aa alignment

# rule convert_phylip:
#     input: "{output_dir}/alignment/aln_aa/alignment.alg"
#     output: temp("{output_dir}/fastME/alignment.alg.phy")
#     conda: "../envs/sp_python.yaml"
#     script: "../scripts/fasta2phylip.py"
#
# rule fastme:
#     input: rules.convert_phylip.output
#     output: "{output_dir}/fastME/FM.nwk"
#     params: config["distboot"]
#     log: "{output_dir}/fastME/logs/fastme.log"
#     benchmark: "{output_dir}/fastME/benchmarks/fastme.txt"
#     threads: 4
#     conda: "../envs/sp_tree.yaml"
#     shell: '''
#         fastme -q -p -T {threads} -b {params} -i {input} -o {output} > {log}
#         rm {input}_fastme*
#         '''


    # rule quicktree:
#     input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
#     output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_QT.nwk"
#     wildcard_constraints:
#         alphabet="aa"
#     params: config["distboot"]
#     log: outdir+"/log/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.log"
#     benchmark: outdir+"/benchmarks/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.txt"
#     conda: "../envs/sp_tree.yaml"
#     shell:'''
# esl-reformat stockholm {input} | quicktree -boot {params} -in a -out t /dev/stdin | paste -s -d '' > {output} 2> {log}
# '''

