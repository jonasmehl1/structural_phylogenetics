# Run iqtree with the specified ML options, so one tree with an aa substitution model and one or more 3di options
rule iqtree:
    input: "{output_dir}/alignment/foldmason/alignment_{alphabet}_trimmed.fa"
    output:
        tree="{output_dir}/iqtree/{alphabet}/{model}.treefile",
        iqtree="{output_dir}/iqtree/{alphabet}/{model}.iqtree"
    wildcard_constraints:
        alphabet="aa|3di",
        model="ML|3di|AF"
    params:
        threedi_submat=config['subst_matrix']['3di'],
        AF_submat=config['subst_matrix']['AF'],
        models = config['ML_model'],
        ufboot=config['UF_boot']
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
    '''

# Create a partition file of the ML tree with the 3di options
rule get_part:
    input:
        ML="{output_dir}/iqtree/aa/ML.iqtree",
        ThreeDI="{output_dir}/iqtree/3di/{combs_3di}.iqtree"
    output:
        "{output_dir}/iqtree_partitioned/{combs_3di}_combined_partition.part"
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
        tree="{output_dir}/iqtree_partitioned/ML_{combs_3di}.treefile",
        iqtree="{output_dir}/iqtree_partitioned/ML_{combs_3di}.iqtree"
    wildcard_constraints:
        combs_3di="3di|AF"
    params:
        threedi_submat=config['subst_matrix']['3di'],
        AF_submat=config['subst_matrix']['AF'],
        ufboot=config['UF_boot'],
    benchmark: "{output_dir}/iqtree_partitioned/benchmarks/ML_{combs_3di}.txt"
    threads: 24
    shell: """
    tree_prefix={output_dir}/iqtree_partitioned/ML_{wildcards.combs_3di}

    iqtree2 -s {input.fa} -p {input.part} --prefix $tree_prefix -B {params.ufboot} -T {threads} --quiet 
    """

##### FastME IMD #####


# Use T-Coffee to get IMD matrices from the fasta sequences
rule tcoffee_templates:
    input:   
        fasta = "{output_dir}/alignment/trimal/alignment_aa_trimmed.fa",
        pdb_dir = f"{input_dir}"
    output: 
        templates = "{output_dir}/fastME/tc_templates.txt"
    script:
        "../scripts/t_coffee_template.py"  # Python script to run

# Get the IMD matrices
rule IMD_matrix:
    input:
        templates = "{output_dir}/fastME/tc_templates.txt",
        alignment = "{output_dir}/alignment/trimal/alignment_aa_trimmed.fa"
    output: "{output_dir}/fastME/IMD_matrices.txt"
    threads: 8
    params:
        boot = config['fastME_boot']
    shell: """
    
    export THREED_TREE_MODE=10
    
    t_coffee -other_pg seq_reformat \
        -in {input.alignment} \
        -in2 {input.templates} \
        -action +replicates {params.boot} +phylo3d +print_replicates \
        -output dm > {output}
    """
    
# replace all -1 values with 100
rule correct_matrix:
    input:
        raw_matrix="{output_dir}/fastME/IMD_matrices.txt"
    output:
        fixed_matrix="{output_dir}/fastME/IMD_matrices_corrected.txt",
        log="{output_dir}/fastME/IMD_matrices_corrected.log"
    shell: """
    fixed_max=100

    num_replacements=$(awk -v max=$fixed_max '{{count+=gsub("-1", max)}} END {{print count}}' {input.raw_matrix})
    num_matrices=$(awk 'BEGIN{{count=0}} NF==0{{count++}} END{{print count}}' {input.raw_matrix})

    # Replace -1 with the fixed maximum distance (20) and save the new matrix
    awk -v max=$fixed_max '{{gsub("-1", max)}}1' {input.raw_matrix} > {output.fixed_matrix}

    echo "Max Distance Used: $fixed_max" > {output.log}
    echo "Total -1 Replacements: $num_replacements" >> {output.log}
    echo "Number of Matrices (Replicates): $num_matrices" >> {output.log}
    """

# Extract individual replicates from the big matrix file and compute the reference tree
rule extract_matrices:
    input:
        matrices_all ="{output_dir}/fastME/IMD_matrices_corrected.txt"
    output:
        main= "{output_dir}/fastME/ref_IMD.txt",
        replicates=directory("{output_dir}/fastME/replicates/")

    params: outdir = "{output_dir}"
    shell: """
    # Ensure the replicates directory exists.
    mkdir -p {output.replicates}
    
    # Split the big matrix into individual replicate files.
    awk -v RS="" '{{print > ("{output.replicates}/replicate_" NR ".txt")}}' {input.matrices_all}

    mv {output.replicates}/replicate_1.txt {params.outdir}/fastME/ref_IMD.txt
    """

# Get the reference tree
rule fastME_reference:
    input: "{output_dir}/fastME/ref_IMD.txt"
    output: "{output_dir}/fastME/ref_IMD.nwk"
    log: "{output_dir}/fastME/logs/main_IMD.log"
    shell: """
    fastme -i {input} -o {output} -g 1.0 -s -n -z 5 -I {log}
    """

# Run FastME on the bootstrapped distance matrices to create a separate tree from each
rule fastME_bootstrap_trees:
    input: "{output_dir}/fastME/replicates/replicate_{replicate}.txt"
    output: "{output_dir}/fastME/replicate_trees/replicate_{replicate}.nwk"
    log: "{output_dir}/fastME/logs/replicate_{replicate}.log"
    shell: """
    fastme -i {input} -o {output} -g 1.0 -s -n -z 5 -I {log}
    """

# Combine all trees and get branch bootstrap support
rule fastME_combine_trees:
    input:
        ref_tree="{output_dir}/fastME/ref_IMD.nwk",
        replicates=expand("{{output_dir}}/fastME/replicate_trees/replicate_{replicate}.nwk",replicate=range(2,102))
    output: "{output_dir}/fastME/FastME_IMD_final.nwk"
    params:
        script="workflow/scripts/IMD_bootstrap.R",
        rep_dir="{output_dir}/fastME/replicate_trees"
    shell: """  
    Rscript {params.script} --ref_tree {input.ref_tree} --rep_dir {params.rep_dir} --out {output} 
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

