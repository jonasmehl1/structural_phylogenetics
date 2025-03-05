# Download pdb structures and get sequence information from uniprot
rule make_structures_dir:
    output:
        directory(f"{output_dir}/structs")
    shell:
        "mkdir -p {output}"
    
rule dl_ids_sequences:
    input:
        ids= "data/input/" + config["identifiers"]
    output:
        "{output_dir}" + "/sequence_dataset.csv"
    log:
        "{output_dir}" + "/logs/dlsequences.log"
    params:
        custom_structs=config["custom_structs"],
        clean_folder=config["clean_folder"]
    script:
        "../scripts/dl_sequences.py"

rule dl_ids_structs:
    input:
        "{output_dir}" + "/sequence_dataset.csv"
    output:
        "{output_dir}" + "/finalset.csv"
    log:
        "{output_dir}" + "/logs/dlstructs.log"
    params:
        filtervar=config["eval_both"],
        filtervar_min=config["low_confidence"],
        filtervar_avg=config["coverage"],
        custom_structs=config["custom_structs"],
        clean_folder=config["clean_folder"]
    script:
        "../scripts/dl_structs.py"

rule plddt:
    input:
        "{output_dir}" + "/finalset.csv"
    output:
        "{output_dir}" + "/plddt.json",
        touch("{output_dir}" + "/done.txt")
    log:
        "{output_dir}" + "/logs/plddt.log"
    script:
        "../scripts/grabplddt.py"



