#!/bin/bash

# =====================================================================
# CNVkit Pipeline in Paraganglioma Samples
# =====================================================================

# Main directory definitions - keeping consistent with original script
in_dir="/mnt/d/CNVkit"
ref_dir="/mnt/d/CNVkit/ref" #hg38, hg38-access, gnomAD, dgv, rmsk 

#48 TUMOR SAMPLES
def_tumor_dir() {
    in_dir="/mnt/d/CNVkit"
    ref_dir="/mnt/d/CNVkit/ref"
    base_dir="/mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789"
    targets_dir="/mnt/d/CNVkit/tumor/tumor_targets"
    out_dir="/mnt/d/CNVkit/tumor/tumor_out"
    res_dir="/mnt/d/CNVkit/tumor/tumor_res"
}

def_tumor_dir

#8 MODEL SAMPLES
def_model_dir_pooledref() {
    in_dir="/mnt/d/CNVkit"
    ref_dir="/mnt/d/CNVkit/ref"
    base_dir="/mnt/d/CNVkit/model/WES_modelli"
    targets_dir="/mnt/d/CNVkit/model/model_targets"
    out_dir="/mnt/d/CNVkit/model/with_pooledref/model_out"
    res_dir="/mnt/d/CNVkit/model/with_pooledref/model_res"
}

def_model_dir_pooledref

def_model_dir_ref() {
    in_dir="/mnt/d/CNVkit"
    ref_dir="/mnt/d/CNVkit/ref"
    base_dir="/mnt/d/CNVkit/model/WES_modelli"
    targets_dir="/mnt/d/CNVkit/model/model_targets"
    out_dir="/mnt/d/CNVkit/model/with_ref/model_out"
    res_dir="/mnt/d/CNVkit/model/with_ref/model_res"
}

def_model_dir_ref

def_model_dir_pooledref

def_model_dir_tumor_ref() {
    in_dir="/mnt/d/CNVkit"
    ref_dir="/mnt/d/CNVkit/ref"
    base_dir="/mnt/d/CNVkit/model/WES_modelli"
    targets_dir="/mnt/d/CNVkit/model/model_targets"
    out_dir="/mnt/d/CNVkit/model/with_tumor_ref/model_out"
    res_dir="/mnt/d/CNVkit/model/with_tumor_ref/model_res"
}

def_model_dir_tumor_ref

def_model_dir_flatref() {
    in_dir="/mnt/d/CNVkit"
    ref_dir="/mnt/d/CNVkit/ref"
    base_dir="/mnt/d/CNVkit/model/WES_modelli"
    targets_dir="/mnt/d/CNVkit/model/model_targets"
    out_dir="/mnt/d/CNVkit/model/with_flatref/model_out"
    res_dir="/mnt/d/CNVkit/model/with_flatref/model_res"
}

def_model_dir_flatref


# Log files
error_log="failed_files.log"
orphan_vcf_report="orphan_vcfs.txt"
vcf_missing_log="$res_dir/missing_vcf.log"
parameter_log="$res_dir/parameter_optimization.log"

# Create output directories if they don't exist
mkdir -p "$out_dir" "$res_dir" "$ref_dir"

# =====================================================================
# FUNCTIONS
# =====================================================================

# Function to find files by pattern and sample name
# Usage: find_sample_files <directory> <sample_name> <file_extension> [output_variable]
find_sample_files() {
    local directory="$1"
    local sample_name="$2"
    local file_extension="$3"
    local output_var_name="$4"
    
    local result
    
    # Handle different search patterns based on file extension
    case "$file_extension" in
        "cnr")
            result=$(find "$directory" -type f -name "${sample_name}.cnr" -print -quit)
            ;;
        "cns")
            result=$(find "$directory" -type f -name "${sample_name}.cns" -print -quit)
            ;;
        "vcf")
            result=$(find "$directory" -type f \( -path "*${sample_name}.hard-filtered.vcf*" -o -path "*${sample_name}.sv.vcf*" \) -print -quit)
            ;;
        "bam")
            result=$(find "$directory" -type f -name "${sample_name}.bam" -print -quit)
            ;;
        *)
            result=$(find "$directory" -type f -name "${sample_name}.${file_extension}" -print -quit)
            ;;
    esac
    
    # If output variable name is provided, set that variable
    if [ -n "$output_var_name" ]; then
        eval "$output_var_name=\"$result\""
    else
        # Otherwise just echo the result
        echo "$result"
    fi
}

# Function to process a single BAM file
#OUTPUT CNR CNS GENEMETRICS.TXT 
process_bam() {
    local file="$1"
    local ref_dir="$2"
    local targets_dir="$3"
    local out_dir="$4"
    local mode="$5" #with_pooledref with_flatref with_ref
    
    sample_name=$(basename "$file" .bam)
    sample_out="${out_dir}/${sample_name}"
    
    # Skip if output directory already exists
    if [ ! -d "$sample_out" ]; then
        echo "Creating directory ${sample_out}"
        mkdir -p "$sample_out"
        
        if [[ "$mode" == "with_pooledref" ]]; then 
            #Run CNVkit batch with pooled reference
            echo "Launching batch command for ${file} with ${targets_dir}/blood_pooled_reference.cnn"
                cnvkit.py batch "$file" \
                -r "${targets_dir}/blood_pooled_reference.cnn" \
                --output-dir "$sample_out" \
                --method hybrid \
                -p $(nproc) \
                --scatter 

            echo "Success: wrote ${sample_out}/${sample_name}.cnn"

        elif [[ "$mode" == "with_flatref" ]]; then 
            #Run CNVkit batch with flat reference
            echo "Launching batch command for ${file} with flat reference"
                cnvkit.py batch "$file" -r "${targets_dir}/flat_reference.cnn" \
                --output-reference "${sample_out}/${sample_name}.cnn" \
                --output-dir "$sample_out" \
                -p $(nproc) \
                --drop-low-coverage \
                --scatter 

            echo "Success: wrote ${sample_out}/${sample_name}.cnn"

        elif [[ "$mode" == "with_ref" ]]; then 
            #Run CNVkit batch with normal
            #normal=${file/tum-001/blood}
            #normal=${file/2D-001/blood}
            normal=${file/2D-001/tum-001}
            echo "Launching batch command for ${file} using ${normal} as reference"
            cnvkit.py batch "$file" -n "$normal" \
            --fasta "${ref_dir}/hg38.fa" \
            --targets "${targets_dir}/targets.bed" \
            --antitargets "${targets_dir}/antitargets.bed" \
            --annotate "${ref_dir}/refFlat.txt" \
            --output-dir "$sample_out" \
            -p $(nproc) \
            --scatter

            echo "Success: wrote ${sample_out}/${sample_name}.cnn"

        else 
            echo "mode not recognized, please type one of the following: [ with_pooledref, with_flatref, with_ref ]"
        fi 
            echo "Completed analysis for ${sample_name}"
    fi
}

# Function to analyze sample results
analyze_sample_results() {
    local sample_name="$1"
    local sample_dir="$2"
    local sample_res="$3"
    local vcf_pattern="$4"
    local vcf_missing_log="$5"

    # Skip if results directory already exists
    if [ ! -d "$sample_res" ]; then
        echo "Creating directory ${sample_res}"
        mkdir -p "$sample_res"
        
        # Use our new function to find files
        cnr_file=$(find_sample_files "$sample_dir" "$sample_name" "cnr")
        cns_file=$(find_sample_files "$sample_dir" "$sample_name" "cns")
        
        if [[ -f "$cnr_file" ]]; then
            # Identify breakpoints if cns file exists
            if [[ -f "$cns_file" ]]; then
                echo "Launching break command for ${sample_name}"
                cnvkit.py breaks "$cnr_file" "$cns_file" --min-probes 5 > "$sample_res/${sample_name}_breaks.txt"
            fi
            
            #Find vcf path based on pattern 
            sample_prefix=$(basename "$sample_dir" | awk -F'-' '{print $1}')
            vcf_path=$(find "${base_dir}" -type f \( -path "*${sample_prefix}*${vcf_pattern}*" \) -print -quit)

            # Call CNVs using VCF (if available)
            if [ -f "$vcf_path" ]; then
                echo "Found ${vcf_path} for $sample_name, launching cnvkit.py call"
                cnvkit.py call "$cns_file" --vcf "$vcf_path" -o "$sample_res/${sample_name}_call.cns"
            else
                echo "VCF not found for $sample_name. Skipping call command."
                echo "$sample_name" >> "$vcf_missing_log"
                cnvkit.py call "$cns_file" -o "$sample_res/${sample_name}_call.cns"
            fi
            
            # Run genemetrics on cns file
            if [ -f "$sample_res/${sample_name}_call.cns" ]; then
                echo "Running genemetrics for ${sample_name}"
                cnvkit.py genemetrics "$cnr_file" -s "$sample_res/${sample_name}_call.cns" -o "$sample_res/${sample_name}_genemetrics.txt"
            fi
        else
            echo "Missing .cnr file for $sample_name."
        fi
    fi
}

# Function to generate multi-sample heatmap
generate_multi_sample_heatmap() {
    local out_dir="$1"
    local res_dir="$2"
    local cns_list="$3"
    
    echo "Generating multi-sample heatmap..."
    
    # Create directory for multi-sample output
    mkdir -p "${res_dir}/multi_sample"
    
    # Generate heatmap for all samples from batch output
    if [ -s "${cns_list}" ]; then
        echo "Generating heatmap for all samples from batch output..."
        cnvkit.py heatmap $(cat "${cns_list}") -d \
            -o "${res_dir}/multi_sample/all_samples_batch_heatmap.pdf"
        
        # Generate heatmaps for chromosomes relevant to paraganglioma
        # These chromosomes are verified to be compatible with hg38
        for chr in chr1 chr3 chr4 chr11 chr17 chr22; do
            echo "Generating heatmap for chromosome $chr..."
            cnvkit.py heatmap $(cat "${cns_list}") -d \
                --chromosome $chr \
                -o "${res_dir}/multi_sample/all_samples_${chr}_heatmap.pdf"
        done
    else 
        echo "cns_list file  not present in ${res_dir}/multi_sample/"
    fi
    
    echo "Multi-sample heatmap generation completed."
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================

# Check available cores
num_cores=$(nproc)
echo "Number of available cores: $num_cores"

# Phase 2: BAM processing
echo "=== PROCESSING BAM FILES ==="
#find "$base_dir" -type f -name "*tum-001.bam" ! -name "tmp*.bam" | while read file; do
#find "$base_dir" -type f -name "*2D-001.bam" ! -name "tmp*.bam" | while read file; do
find "$base_dir" -type f -name "*blood.bam" ! -name "tmp*.bam" | while read file; do
    process_bam "$file" "$ref_dir" "$targets_dir" "$out_dir" "with_flatref"
done

# Phase 3: Results analysis
echo "=== ANALYZING RESULTS ==="
mkdir -p "$res_dir"
> "$vcf_missing_log" # Reset log file

sample_list=()
for sample_dir in "$out_dir"/*; do
    sample_name=$(basename "$sample_dir")
    sample_res="$res_dir/$sample_name"
    sample_list+=("$sample_name")
    analyze_sample_results "$sample_name" "$sample_dir" "$sample_res" ".hard-filtered.vcf" "$vcf_missing_log" 
done

#MODEL
#"tum.hard-filtered.vcf"
#"2D.hard-filtered.vcf"
#48 TUMOR
#".hard-filtered.vcf"

#Calculate genemetrics on cnr 
for sample_dir in "$out_dir"/*; do  
    sample_name=$(basename "$sample_dir")
    sample_res="$res_dir/$sample_name"
    cnr_file=$(find_sample_files "$sample_dir" "$sample_name" "cnr")
    if [[ -f "$sample_res/${sample_name}_CNR_genemetrics.txt" ]]; then
        echo "$cnr_file already exists skipping"
        continue 
    else 
        cnvkit.py genemetrics "$cnr_file" -o "$sample_res/${sample_name}_CNR_genemetrics.txt"
    fi 
done

# Phase 4: Multi-sample heatmap generation
# Find all .cns files (batch output)
find "$out_dir" -name "*.cns" ! -name "*call*" ! -name "*bintest*" > "${res_dir}/multi_sample/cns_list.txt"

echo "=== GENERATING MULTI-SAMPLE HEATMAP ==="
generate_multi_sample_heatmap "$out_dir" "$res_dir" "${res_dir}/multi_sample/sorted_cns_list.txt"


#Single sample heatmap 
patients=$(ls "$out_dir" | cut -d'-' -f1 | sort -u)
for patient in $patients; do
    #mkdir -p "${res_dir}/single_sample/$patient"
    
    #find "$out_dir" -name "${patient}-*.cns" ! -name "*call*" ! -name "*bintest*" > "${res_dir}/single_sample/$patient/cns_list.txt"

    generate_multi_sample_heatmap "$out_dir" "${res_dir}/single_sample/$patient" "${res_dir}/single_sample/$patient/cns_list.txt"
done

#sort by patient 
awk -F'/' '{split($NF, a, "-"); print a[1] "\t" $0}' cns_list.txt | \
sort -k1,1 -k2,2 | \
cut -f2- > sorted_cns_list

