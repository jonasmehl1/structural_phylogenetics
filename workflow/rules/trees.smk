# Run iqtree with the specified ML options, so one tree with an aa substitution model and one or more 3di options
rule iqtree:
    input: "{output_dir}/alignment/foldmason/alignment_{alphabet}_trimmed.fa"
    output:
        tree="{output_dir}/iqtree/{alphabet}/{model}.iqtree",
    wildcard_constraints:
        alphabet="aa|3di",
        model="ML|3di|AF"
    params:
        threedi_submat=config['subst_matrix']['3di'],
        AF_submat=config['subst_matrix']['AF'],
        models = config['ML_model'],
        ufboot=config['UF_boot']
    log: "{output_dir}/iqtree/logs/{alphabet}_{model}.log"
    benchmark: "{output_dir}/iqtree/benchmarks/{alphabet}_{model}.txt"
    threads: 8
    shell: '''
    tree_prefix={output_dir}/iqtree/{wildcards.alphabet}/{wildcards.model}

    model={wildcards.model}
    if [ "$model" == "3di" ]; then
        model="--mset {params.threedi_submat}"
    elif [ "$model" == "AF" ]; then
        model="--mset {params.AF_submat}"
    elif [ "$model" == "ML" ]; then
        model="-m {params.models}"
    fi

    iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads} --quiet \
    --mem 16G --cmin 4 --cmax 12 $model

    mv $tree_prefix.treefile {output.tree}
    mv $tree_prefix.log {log}
    '''

# Create a partition file of the ML tree with the 3di options
rule get_part:
    input:
        ML="{output_dir}/iqtree/aa/ML.iqtree",
        ThreeDI="{output_dir}/iqtree/3di/{combs_3di}.iqtree"
    output: "{output_dir}/iqtree_partitioned/{combs_3di}_combined_partition.part"
    shell: """
    model_ML=$(grep "^Model of substitution" {input.ML} | cut -f2 -d ':')
    model_3di=$(grep "^Model of substitution" {input.ThreeDI} | cut -f2 -d ':')
    length_ML=$(grep "^Input data" {input.ML} | cut -f6 -d ' ')
    length_3di=$(grep "^Input data" {input.ThreeDI} | cut -f6 -d ' ')

    end=$(expr $length_ML + $length_3di)
    start_3di=$(expr $length_ML + 1)

    echo "$model_ML, aa_partition = 1-$length_ML" > {output}
    echo "$model_3di, 3di_partition = $start_3di-$end\n" >> {output}
    """

# run iqtree with the concatenated alignment and the specified partition
rule iqtree_partitioned:
    input:
        fa="{output_dir}/alignment/concat_aln/concatenated_alignment.alg",
        part="{output_dir}/iqtree_partitioned/{combs_3di}_combined_partition.part"
    output:
        tree="{output_dir}/iqtree_partitioned/ML_{combs_3di}.nwk",
    wildcard_constraints:
        combs_3di="3di|AF"
    params:
        threedi_submat=config['subst_matrix']['3di'],
        AF_submat=config['subst_matrix']['AF'],
        ufboot=config['UF_boot'],
    log: "{output_dir}/iqtree_partitioned/logs/ML_{combs_3di}.log"
    benchmark: "{output_dir}/iqtree_partitioned/benchmarks/ML_{combs_3di}.txt"
    threads: 16
    shell: """
    tree_prefix={output_dir}/iqtree_partitioned/ML_{wildcards.combs_3di}

    iqtree2 -s {input.fa} -p {input.part} --prefix $tree_prefix -B {params.ufboot} -T {threads} --quiet 

    mv $tree_prefix.treefile {output.tree}
    mv $tree_prefix.log {log}
    rm -f $tree_prefix.model.gz $tree_prefix.splits.nex $tree_prefix.contree $tree_prefix.ckp.gz
    """

##### FOLDTREE #####

# rule foldseek_allvall_tree:
#     input:
#         structs = config["input_dir"],
#         done=config["outdir"] + "/done.txt"
#     output: "{output_dir}/foldtree/foldseek_allvall/foldseek_allvall.txt"
#     log: "{output_dir}/foldtree/foldseek_allvall/logs/foldseek.log"
#     benchmark: "{output_dir}/foldtree/foldseek_allvall/benchmarks/foldseek.txt"
#     shell:'''
#     foldseek easy-search {input.structs} {input.structs} {output} $TMPDIR/foldseek_tmp \
#     --format-output 'query,target,fident,lddt,alntmscore' --exhaustive-search -e inf \
#     --alignment-type 2 > {log}
#     '''
#
# rule foldseek_distmat:
#     input: "{output_dir}/foldtree/foldseek_allvall/foldseek_allvall.txt"
#     output: "{output_dir}/foldtree/foldseek_distmat/foldseek_3di_fident.txt"
#     script: "../scripts/foldtree/foldseekres2distmat_simple.py"
#
# rule foldtree:
#     input: "{output_dir}/foldtree/foldseek_distmat/foldseek_3di_fident.txt"
#     output: "{output_dir}/foldtree/foldtree_3di_FT.nwk"
#     wildcard_constraints:
#         alphabet="3di"
#     benchmark: "{output_dir}/foldtree/benchmarks/foldtree.txt"
#     shell: '''
#         quicktree -i m {input} | paste -s -d '' > {output}
#         '''

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

