#!/bin/bash
# =============================================================================
# 1. 配置区域
# =============================================================================

# --- 路径设置 ---
# 根工作目录
WORK_DIR="/home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/psmc"

# 输入数据路径
GENOME_REF="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
# 注意末尾必须加斜杠 /，否则拼接文件名会出错
BAM_DIR="/home/data/4.picard_dedup/markdup/" 

# 输出数据路径 (会自动创建文件夹)
FQ_DIR="${WORK_DIR}/fqgz/"
PSMCFA_DIR="${WORK_DIR}/psmcfa/"
PSMC_OUT_DIR="${WORK_DIR}/psmc_out/"

# --- 软件路径 ---
# 建议检查 vcfutils.pl 是否在环境变量中，如果没有，需要写绝对路径
BCFTOOLS="bcftools"
VCFUTILS="vcfutils.pl" 
FQ2PSMCFA="fq2psmcfa"
PSMC_BIN="psmc"
SAMTOOLS="samtools"

# --- PSMC 参数 ---
# -p: 时间间隔参数
PSMC_PARAMS="-N25 -t15 -r5 -p 4+25*2+4+6"

# --- 样本列表 ---
species_list=(CP_1 CX_9 DA_1 DA_2 DA_4 DA_6 DA_7 DA_8 DRS_1 DRS_2 DRS_4 DRS_5 DRS_6 DRS_7 DST_1 DST_2 DST_4 DST_7 DST_8 
JD_1 JMX_17 JMX_18 JMX_19 JMX_1 LCJ_10 LCJ_11 LCJ_13 LCJ_6 LCJ_7 LCJ_9 LQ_15 LS_6 LX_5 QDH_12 RX_2 RX_4 RX_5 RX_7 RX_8 
RX_9 SFZ_2 SK_1 SL_4 SLZ_3 SXK_4 WYL_5 XC_8 XNF_3 YK_4 YX_2 YZY_7 ZJS_10)

# =============================================================================
# 2. 准备工作
# =============================================================================

# 创建输出目录
mkdir -p $FQ_DIR $PSMCFA_DIR $PSMC_OUT_DIR

# =============================================================================
# 3. 处理函数
# =============================================================================

process_species() {
    species=$1
    # 路径拼接
    bam_file="${BAM_DIR}${species}.markdup.bam"
    
    echo "[${species}] Starting process..."
    echo "Processing BAM: $bam_file"

    # --- 步骤 0: 检查 BAM 文件 ---
    if [ ! -f "$bam_file" ]; then
        echo "Error: BAM file not found: $bam_file"
        # 尝试检查是否是因为加上了 .bam 后缀导致的路径错误 (有些文件名本身不带后缀)
        return
    fi

    # --- 步骤 1: 智能计算深度 (优化版) ---
    echo "[${species}] Calculating average depth (using sampling)..."
    
    #  增加 head -n 10000000 
    # 仅读取 BAM 文件的前 20000 万个位点计算深度，极大提高速度，且估算足够准确。
    avg_depth=$($SAMTOOLS depth "$bam_file" | head -n 200000000 | awk '{sum+=$3} END { print sum/NR }')
    
    # 如果深度计算失败 (例如 bam 为空)，设置默认值防止报错
    if [ -z "$avg_depth" ]; then
        echo "Error: Failed to calculate depth. Is the BAM valid?"
        return
    fi

    # 设置阈值: d = 1/3 * avg, D = 2 * avg
    min_d=$(echo "$avg_depth * 0.33" | bc)
    max_D=$(echo "$avg_depth * 2" | bc)
    
    # 取整
    min_d=${min_d%.*}
    max_D=${max_D%.*}
    
    # 防止 min_d 小于 3 (PSMC 对极低深度很敏感，建议至少保留 3 或 5)
    if [ "$min_d" -lt 3 ]; then min_d=3; fi
    
    echo "[${species}] Avg Depth: $avg_depth | Filter: -d $min_d -D $max_D"

    # --- 步骤 2: Bam -> Fastq (vcf2fq) ---
    # 管道命令：mpileup -> call -> vcf2fq -> gzip
    # -C50 是 BWA 推荐参数，降低 mapping quality 过高但由 reads 错误导致的假阳性
    $BCFTOOLS mpileup -C50 -f "$GENOME_REF" "$bam_file" | \
    $BCFTOOLS call -c | \
    $VCFUTILS vcf2fq -d "$min_d" -D "$max_D" | \
    gzip > "${FQ_DIR}${species}.fq.gz"

    # 检查上一步是否成功生成了非空文件
    if [ ! -s "${FQ_DIR}${species}.fq.gz" ]; then
        echo "Error: fq.gz file is empty. Check pipeline for $species"
        return
    fi

    # --- 步骤 3: Fastq -> psmcfa ---
    echo "[${species}] Generating psmcfa..."
    $FQ2PSMCFA -q20 "${FQ_DIR}${species}.fq.gz" > "${PSMCFA_DIR}${species}.psmcfa"

    # --- 步骤 4: Run PSMC ---
    echo "[${species}] Running PSMC..."
    $PSMC_BIN $PSMC_PARAMS -o "${PSMC_OUT_DIR}${species}.psmc" "${PSMCFA_DIR}${species}.psmcfa"

    echo "[${species}] All Done."
}

# 导出变量和函数供 parallel 使用
export -f process_species
export WORK_DIR GENOME_REF BAM_DIR FQ_DIR PSMCFA_DIR PSMC_OUT_DIR
export BCFTOOLS VCFUTILS FQ2PSMCFA PSMC_BIN SAMTOOLS PSMC_PARAMS

# 并行运行 (根据服务器负载调整 -j)
# --keep-order 可以保持输出日志的顺序，方便查看
parallel -j 1 --eta --verbose process_species ::: "${species_list[@]}"
