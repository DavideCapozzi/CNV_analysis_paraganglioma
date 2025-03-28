#!/bin/bash

# =====================================================================
# CNVkit prep file for download, process needed files for cnvkit_pipeline
# =====================================================================

# Main directory definitions - keeping consistent with original script
in_dir="/mnt/d/CNVkit"
ref_dir="/mnt/d/CNVkit/ref" #hg38, hg38-access, gnomAD, dgv, rmsk 

#48 TUMOR SAMPLES
base_dir="/mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789"
targets_dir="/mnt/d/CNVkit/tumor/tumor_targets"
out_dir="/mnt/d/CNVkit/tumor/tumor_out"
res_dir="/mnt/d/CNVkit/tumor/tumor_res"

#8 MODEL SAMPLES
base_dir="/mnt/d/CNVkit/model/WES_modelli"
targets_dir="/mnt/d/CNVkit/model/model_targets"
out_dir="/mnt/d/CNVkit/model/model_out"
res_dir="/mnt/d/CNVkit/model/model_res"

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
            result=$(find "$directory" -type f \( -path "*${sample_name}.hard-filtered.vcf*" \) -print -quit)
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

# Function to index BAM files
index_bam_files() {
    local base_dir="$1"
    local error_log="$2"
    
    echo "Indexing BAM files..."
    > "$error_log" # Reset error log
    
    find "$base_dir" -type f -name "*.bam" ! -name "tmp*.bam" -print0 | while IFS= read -r -d $'\0' bam_file; do
        # Check if index already exists 
        if [[ -f "${bam_file}.bai" ]]; then
            echo "Index already exists for: $bam_file"
            continue
        fi
        
        echo "Processing: $bam_file"
        
        if samtools quickcheck -q "$bam_file" 2>/dev/null; then
            # Index file
            if samtools index "$bam_file"; then
                echo "Index created successfully"
            else
                echo "ERROR: Indexing failed for $bam_file" | tee -a "$error_log"
            fi
        else
            echo "ERROR: Quickcheck failed for $bam_file" | tee -a "$error_log"
        fi
    done
    
    echo "BAM indexing completed."
}


# Function to check correspondence between VCF and BAM
check_vcf_bam_correspondence() {
    local vcf_source="$1"
    local bam_base="$2"
    local orphan_vcf_report="$3"
    
    echo "Checking correspondence between VCF and BAM..."
    > "$orphan_vcf_report"  # Empty report
    
    find "$vcf_source" -type f -name "*hard-filtered.vcf.gz" | while read -r vcf_file; do
        sample_name=$(basename "$vcf_file" | awk -F'.' '{print $1}')
        
        # Search for corresponding BAM using our new function
        bam_path=$(find_sample_files "$bam_base" "$sample_name" "bam")
        
        if [ -z "$bam_path" ]; then
            echo "VCF without BAM: $vcf_file" | tee -a "$orphan_vcf_report"
        fi
    done
    
    echo "Correspondence check completed."
}


# Function to index VCF files
index_vcf_files() {
    local vcf_root="$1"
    
    echo "Indexing VCF files..."
    
    find "$vcf_root" -type f \( -name "*.hard-filtered.vcf.gz" -o -name "*.sv.vcf.gz" \) | while read vcf_gz; do
        tbi_file="${vcf_gz}.tbi"
        
        if [ ! -f "$tbi_file" ]; then
            echo "Indexing ${vcf_gz}"
            tabix -p vcf "$vcf_gz"
        else
            echo "Index already exists for ${vcf_gz}"
        fi
    done
    
    echo "VCF indexing completed."
}

# Function to download and prepare reference files
prepare_reference_files() {
    local ref_dir="$1"
    
    echo "Preparing reference files..."
    
    # Download gap file
    if [ ! -f "${ref_dir}/ucsc-gaps-hg38.bed" ]; then
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
        gunzip gap.txt.gz
        awk '{print $2"\t"$3"\t"$4"\t"$8}' gap.txt > "${ref_dir}/ucsc-gaps-hg38.bed"
        rm gap.txt
    fi
    
    # Download gnomAD structural variants
    if [ ! -f "${ref_dir}/gnomad.v4.1.cnv.hg38.bed" ]; then
        wget http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt
        
        # Convert to BED format - fixed AWK script with proper quoting
        if [ -f "gnomad.v4.1.cnv.all.vcf" ]; then
            awk '
            BEGIN {OFS="\t"} 
            /^#/ {next} 
            {
                split($8, info, ";"); 
                end_val = $2;  # Default if END does not exist
                for (i in info) {
                    if (info[i] ~ /^END=/) {
                        split(info[i], end, "="); 
                        end_val = end[2]
                    }
                } 
                print $1, $2-1, end_val, $5
            }' "gnomad.v4.1.cnv.all.vcf" > "${ref_dir}/gnomad.v4.1.cnv.hg38.bed"
            
            # Sort BED file
            sort -k1,1 -k2,2n "${ref_dir}/gnomad.v4.1.cnv.hg38.bed" > "${ref_dir}/gnomad.v4.1.cnv.hg38.sorted.bed"
        fi
    fi
    
    # Download and prepare repeats file
    if [ ! -f "${ref_dir}/repeats.bed" ] && [ -f "rmsk.txt.gz" ]; then
        gunzip rmsk.txt.gz
        awk -v OFS='\t' '{print $6, $7, $8, $11}' rmsk.txt > "${ref_dir}/repeats.bed"
        rm rmsk.txt
    fi
    
    # Prepare DGV clean file
    if [ -f "${ref_dir}/dgv_hg38.bed" ] && [ ! -f "${ref_dir}/dgv_hg38.clean.bed" ]; then
        grep -v -E "CNV|OTHER" "${ref_dir}/dgv_hg38.bed" > "${ref_dir}/dgv_hg38.clean.bed"
    fi
    
    # Create centromeres/telomeres file
    if [ -f "${ref_dir}/ucsc-gaps-hg38.bed" ] && [ ! -f "${ref_dir}/centromeres_telomeres.bed" ]; then
        awk '$4 ~ /centromere|telomere/ {print $1"\t"$2"\t"$3"\t"$4}' "${ref_dir}/ucsc-gaps-hg38.bed" > "${ref_dir}/centromeres_telomeres.bed"
    fi
    
    # Create file for filtering problematic regions (SINE/LINE)
    if [ -f "${ref_dir}/repeats.bed" ] && [ ! -f "${ref_dir}/problematic_repeats.bed" ]; then
        grep "SINE\|LINE" "${ref_dir}/repeats.bed" > "${ref_dir}/problematic_repeats.bed"
    fi
    
    echo "Reference file preparation completed."
}

# Function to create BAM file list
create_bam_list() {
    local base_dir="$1"
    local targets_dir="$2"
    
    echo "Creating BAM file list..."
    find "$base_dir" -name "*.bam" | grep -v "tmp" | grep -v "tmp2" > "${targets_dir}/bam_list.txt"
    echo "BAM list created: $output_file with $(wc -l < "${targets_dir}/bam_list.txt") files"
}

# Function to prepare accessibility target files and flat reference 
prepare_targets_flatreference() {
    local ref_dir="$1"
    local targets_dir="$2"
    local panel_file="$3"

    
    echo "Preparing accessibility and target files..."
    
    # Create base accessibility file
    if [ ! -f "${ref_dir}/access-hg38.bed" ]; then
        cnvkit.py access "${ref_dir}/hg38.fa" -x "${ref_dir}/ucsc-gaps-hg38.bed" -o "${ref_dir}/access-hg38.bed"
    fi
    
    # Divide accessibility file into bins
    if [ ! -f "${targets_dir}/targets.bed" ] || [ ! -f "${targets_dir}/antitargets.bed" ]; then
        cnvkit.py autobin $(cat "${targets_dir}/bam_list.txt") \
          --targets "${targets_dir}/${panel_file}" \
          --access "${ref_dir}/access-hg38.bed" \
          --target-min-size 100 \
          --target-max-size 500 \
          --antitarget-max-size 5000 \
          --target-output-bed "${targets_dir}/targets.bed" \
          --antitarget-output-bed "${targets_dir}/antitargets.bed"
    fi
    
    echo "Accessibility and target file preparation completed." ##Echoes not always accurate
    echo "Creating flat reference"
    
    if [ ! -f "${ref_dir}/flat_reference.cnn" ]; then 
        cnvkit.py reference -f "${ref_dir}/hg38.fa" -t "${targets_dir}/targets.bed" -a "${targets_dir}/antitargets.bed" -o "${targets_dir}/flat_reference.cnn"
    fi 

    echo "Flat reference created" ##Echoes not always accurate
}

# Function to create pooled reference
# Not used (don't have normal samples)
create_pooled_reference() {
    local ref_dir="$1"
    
    echo "Creating pooled reference..."
    
    if [ ! -f "${ref_dir}/pooled_reference.cnn" ]; then
        cnvkit.py batch \
            --method hybrid \
            --targets "${ref_dir}/targets.bed" \
            --fasta "${ref_dir}/hg38.fa" \
            --normal $(cat bam_list.txt) \
            --output-reference "${ref_dir}/pooled_reference.cnn" \
            --output-dir "${ref_dir}"
    fi
    
    echo "Pooled reference created."
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================

# Check available cores
num_cores=$(nproc)
echo "Number of available cores: $num_cores"

# Copy parameter optimization log to results directory
mkdir -p "$(dirname "$parameter_log")"
if [ -f "/home/ubuntu/parameter_optimization.log" ]; then
    cp "/home/ubuntu/parameter_optimization.log" "$parameter_log"
    echo "Parameter optimization log copied to $parameter_log"
fi

# Phase 1: File preparation
index_bam_files "$base_dir" "$error_log"
check_vcf_bam_correspondence "$base_dir" "$base_dir" "$orphan_vcf_report"
index_vcf_files "$base_dir"
create_bam_list "$base_dir" "$targets_dir"
prepare_reference_files "$ref_dir"
prepare_targets_flatreference "$ref_dir" "$targets_dir" "hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated.bed" #"xgen-exome-hyb-panel-v2-probes-hg38.bed"
create_pooled_reference "$ref_dir"

# =====================================================================
# SigProfiler GenerateMatrix
# =====================================================================

# prepare files  
find "$base_dir" -type f -name "*.bam" ! -name "tmp*.bam" -print0 | while IFS= read -r -d $'\0' bam_file; do
    sample_name=$(basename "$bam_file" | awk -F'.' '{print $1}')       
    find_sample_files $base_dir $sample_name "vcf" >> "${targets_dir}/vcf_list.txt"
done

target_dir="matrix_generator"
mkdir -p "$target_dir"

while read vcf_path; do
    index_path="${vcf_path/.gz/.gz.tbi}"
    rsync -av --progress "$vcf_path" "$target_dir/"
    rsync -av --progress "$index_path" "$target_dir/"
done < "${targets_dir}/vcf_list.txt"


ref_dir="/mnt/d/CNVkit/ref"
genome_name="GRCh38"
temp_dir="${ref_dir}/${genome_name}_temp"


mkdir -p "${temp_dir}/${genome_name}" || exit 1 


cp -v "${ref_dir}/hg38.fa" "${temp_dir}/${genome_name}/${genome_name}.fa"
cp -v "${ref_dir}/hg38.fa.fai" "${temp_dir}/${genome_name}/${genome_name}.fa.fai"
cp -v "${ref_dir}/hg38.dict" "${temp_dir}/${genome_name}/${genome_name}.dict"

# tarball creation
tar -czvf "${ref_dir}/GRCh38.tar.gz" -C "${temp_dir}" "${genome_name}"

# install GRCh38
SigProfilerMatrixGenerator install "${genome_name}" --local_genome "/mnt/d/CNVkit/ref/"

SigProfilerMatrixGenerator matrix_generator \
    "mutational_catalog_output" \
    "GRCh38" \
    "matrix_generator/vcf_files/" \
    --exome \
    --bed_file "${ref_dir}/access-hg38.bed" \
    --tsb_stat 