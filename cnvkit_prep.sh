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
out_dir="/mnt/d/CNVkit/model/without_ref/model_out"
res_dir="/mnt/d/CNVkit/model/without_ref/model_res"

out_dir="/mnt/d/CNVkit/model/with_ref/model_out"
res_dir="/mnt/d/CNVkit/model/with_ref/model_res"

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




check_vcf_bam_correspondence() {
    local vcf_source="$1"
    local bam_base="$2"
    local orphan_vcf_report="$3"
    local wes_modelli_base="/mnt/d/CNVkit/model/WES_modelli"  # Base path aggiunto
    
    echo "Checking correspondence between VCF and BAM..."
    > "$orphan_vcf_report"  # Svuota il report
    
    find "$vcf_source" -type f -name "*hard-filtered.vcf.gz" | while read -r vcf_file; do
        sample_name=$(basename "$vcf_file" | awk -F'.' '{print $1}')
        
        # 1. Controllo corrispondenza BAM
        bam_path=$(find_sample_files "$bam_base" "$sample_name" "bam")
        
        if [ -z "$bam_path" ]; then
            echo "VCF without BAM: $vcf_file" | tee -a "$orphan_vcf_report"
        fi

        # 2. Copia VCF nella cartella corrispondente
        base_sample=$(echo "$sample_name" | awk -F'-' '{print $1"-"$2}')  # Estrae PTJXXX-XXX
        dest_dir=$(find "$wes_modelli_base" -maxdepth 1 -type d -name "${base_sample}*" -print -quit)
        
        if [ -n "$dest_dir" ]; then
            echo "Copying VCF: $(basename "$vcf_file")"
            echo "Destination directory: $dest_dir"
            cp -v "$vcf_file" "$dest_dir/"  # Aggiunto verbose per verifica
        else
            echo "Warning: No destination folder found for $base_sample"
        fi
    done
    
    echo "Correspondence check and VCF copy completed."
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
    
    if [ ! -f "${ref_dir}/refFlat.txt*"]; then 
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
        gunzip refFlat.txt.gz
    fi 
    
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

# Function to create BAM file list with custom name
create_bam_list() {
    local base_dir="$1"
    local targets_dir="$2"
    local name="${3:-*.bam}" #default value of name "*.bam"
    
    echo "Creating BAM file list..."
    find "$base_dir" -name "$name" | grep -v "tmp" | grep -v "tmp2" > "${targets_dir}/normal_bam_list.txt"
    echo "BAM list created: $output_file with $(wc -l < "${targets_dir}/normal_bam_list.txt") files"
}


prepare_targets_flatreference() {
    local ref_dir="$1"
    local targets_dir="$2"
    local panel_file="$3"

    echo "Preparing accessibility and target files..."
    
    # Modifica critica 1: Usa SOLO campioni normali per autobin
    normal_bams="${targets_dir}/normal_bam_list.txt"  # File specifico per normali
    
    # Genera file access escludendo regioni non sequenziabili
    if [ ! -f "${ref_dir}/access-hg38.bed" ]; then
        cnvkit.py access "${ref_dir}/hg38.fa" -x "${ref_dir}/ucsc-gaps-hg38.bed" -o "${ref_dir}/access-hg38.bed"
    fi
    
    # Modifica critica 2: Usa esplicitamente i normali per autobin
    if [ ! -f "${targets_dir}/targets.bed" ] || [ ! -f "${targets_dir}/antitargets.bed" ]; then
        cnvkit.py autobin $(cat "$normal_bams") \
          --targets "${targets_dir}/${panel_file}" \
          --access "${ref_dir}/access-hg38.bed" \
          --target-min-size 200 \
          --target-max-size 1000 \
          --antitarget-max-size 2000 \
          --method hybrid \
          --target-output-bed "${targets_dir}/targets.bed" \
          --antitarget-output-bed "${targets_dir}/antitargets.bed"
    fi
    
    # Modifica critica 3: Verifica manuale degli antitargets
    if [ $(wc -l < "${targets_dir}/antitargets.bed") -lt 1000 ]; then
        echo "ERROR: Antitargets malformati. Rigenerare con:"
        echo "cnvkit.py antitarget ${targets_dir}/targets.bed -g ${ref_dir}/access-hg38.bed -o ${targets_dir}/antitargets.bed"
        exit 1
    fi
    
    #echo "Creazione flat reference..."
    #cnvkit.py reference -f "${ref_dir}/hg38.fa" -t "${targets_dir}/targets.bed" -a "${targets_dir}/antitargets.bed" -o "${targets_dir}/flat_reference.cnn"
}

generate_coverage() {
    local ref_dir="$1"
    local targets_dir="$2"

    echo "Generating target/antitarget coverage for normal BAMs..."
    mkdir -p "${targets_dir}/coverage"

    while read -r bam; do
        sample_name=$(basename "$bam" .bam)
        target_cnn="${targets_dir}/coverage/${sample_name}.targetcoverage.cnn"
        antitarget_cnn="${targets_dir}/coverage/${sample_name}.antitargetcoverage.cnn"

        if [ -f "$target_cnn" ] && [ -f "$antitarget_cnn" ]; then
            echo "Coverage files exist for $sample_name. Skipping..."
            continue
        fi

        echo "Processing coverage for $sample_name..."

        if [ ! -f "$target_cnn" ]; then
            cnvkit.py coverage "$bam" "${targets_dir}/targets.bed" \
                -o "$target_cnn" 
        fi

        if [ ! -f "$antitarget_cnn" ]; then
            cnvkit.py coverage "$bam" "${targets_dir}/antitargets.bed" \
                -o "$antitarget_cnn" 
        fi

    done < "${targets_dir}/bam_list.txt"

    echo "Coverage file generation completed."
}

create_pooled_reference() {
    local ref_dir="$1"
    local targets_dir="$2"

    target_cnns=($(ls "${targets_dir}/coverage/"*.targetcoverage.cnn))
    antitarget_cnns=($(ls "${targets_dir}/coverage/"*.antitargetcoverage.cnn))

    # CORREZIONE CHIAVE: Rimuovi --method hybrid (non esiste in reference)
    cnvkit.py reference \
        "${target_cnns[@]}" "${antitarget_cnns[@]}" \
        -f "${ref_dir}/hg38.fa" \
        -o "${targets_dir}/blood_pooled_reference.cnn"

    # Verifica manuale della reference (sostituisce --view)
    echo "Validazione reference:"
    awk 'NR > 1 {print $5}' "${targets_dir}/blood_pooled_reference.cnn" | sort -n | uniq -c
}

#Sort alignment, mark and remove duplicates, index output files 
remove_PCR_duplicates() {
    local alignment="$1" 

    echo "Processing file: $alignment"
    base_name=${alignment%.bam}

    #echo "Sorting BAM file before marking duplicates..."
    #samtools sort "$alignment" -o "${base_name}_sorted.bam"

    echo "Running picard MarkDuplicates..."
    picard MarkDuplicates \
        -I "$alignment" \
        -O ${base_name}_marked_duplicates.bam \
        -M ${base_name}_marked_dup_metrics.txt \
        -VALIDATION_STRINGENCY LENIENT \
        --REMOVE_DUPLICATES true

    echo "Indexing BAM file..."
    samtools index "${base_name}_marked_duplicates.bam"
    
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
check_vcf_bam_correspondence "/mnt/d/CNVkit/model/PPGLs_model_WES-37315281/" "$base_dir" "$orph
an_vcf_report" 
index_vcf_files "$base_dir"

create_bam_list "$base_dir" "$targets_dir" 
prepare_reference_files "$ref_dir"

# Phase 1: removing duplicates 
echo "=== Removing picard duplicates ==="
find "$base_dir" -type f -name "*PGL214-363-blood.bam" ! -name "tmp*.bam" | while read file; do
    remove_PCR_duplicates "$file" 
done

create_bam_list "$base_dir" "$targets_dir" "*blood.bam"

prepare_targets_flatreference "$ref_dir" "$targets_dir" "hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated.bed" #"xgen-exome-hyb-panel-v2-probes-hg38.bed"
generate_coverage "$ref_dir" "$targets_dir"
create_pooled_reference "$ref_dir" "$targets_dir"

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