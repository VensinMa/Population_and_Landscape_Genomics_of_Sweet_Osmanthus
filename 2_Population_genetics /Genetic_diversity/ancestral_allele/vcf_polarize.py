#!/usr/bin/env python3
# usage: python vcf_polarize.py <input_vcf> <estsfs_output> <positions_file>

import sys
import argparse
import os
import gzip

def create_parser():
    parser = argparse.ArgumentParser(
        description="Polarize VCF based on EstSFS output."
    )
    parser.add_argument('input_vcf', type=str, help="Input VCF file")
    parser.add_argument('estsfs_output', type=str, help="EstSFS output file (p_anc.txt)")
    parser.add_argument('positions_file', type=str, help="Positions file from vcf_to_estsfs.py")
    return parser

args = create_parser().parse_args()

# 文件名处理
input_filename = os.path.basename(args.input_vcf)
base_name = input_filename.replace(".vcf.gz", "").replace(".vcf", "")
output_vcf_name = f"{base_name}_polarized.vcf"
ancestral_out_name = f"{base_name}_ancestral.txt"

# VCF 打开方式
open_func = gzip.open if args.input_vcf.endswith(".gz") else open
mode = "rt" if args.input_vcf.endswith(".gz") else "r"

try:
    print("Loading data files...")
    site_data = {}
    
    # 1. 读取并清洗 EstSFS 输出文件
    # 规则：表头行以 '0' 开头，数据行以 '1', '2', '3'... 开头
    with open(args.estsfs_output, "r") as est_f:
        clean_est_lines = []
        for line in est_f:
            parts = line.strip().split()
            if len(parts) > 0:
                # 核心修改：检查第一列
                if parts[0] == "0":
                    continue  # 第一列是 0，说明是表头/元数据，跳过
                else:
                    clean_est_lines.append(line)  # 第一列不是 0 (如 1, 2...)，是数据，保留

    # 2. 读取 Positions 文件
    with open(args.positions_file, "r") as pos_f:
        pos_lines = [l for l in pos_f.readlines() if l.strip()] 

    # 3. 检查行数是否匹配
    if len(clean_est_lines) != len(pos_lines):
        print(f"Warning: Line counts do not match after cleaning!")
        print(f"  Positions file lines: {len(pos_lines)}")
        print(f"  Valid EstSFS lines:   {len(clean_est_lines)}")
        print("  (Script will proceed with overlapping parts.)")

    # 4. 构建映射字典
    limit = min(len(pos_lines), len(clean_est_lines))
    for i in range(limit):
        pos_line = pos_lines[i].strip()
        est_line = clean_est_lines[i].strip()
        
        parts = pos_line.split("\t")
        chrom, pos = parts[0], parts[1]
        key = f"{chrom}_{pos}"
        
        est_parts = est_line.split()
        # est-sfs 输出格式: [SiteIndex, State, P-major, ...]
        # 第3列 (索引2) 是概率值
        if len(est_parts) > 2:
            prob = est_parts[2] 
            site_data[key] = float(prob)

    print(f"Loaded {len(site_data)} sites with ancestral probabilities.")

    # 5. 处理 VCF
    with open_func(args.input_vcf, mode) as in_vcf, \
         open(output_vcf_name, "w") as out_vcf, \
         open(ancestral_out_name, "w") as ancestral_out:

        ancestral_out.write("CHROM\tPOS\tref_allele\tancestral_allele\n")

        total_sites_read = 0
        polarized_sites = 0
        trans_table = str.maketrans("01", "10")

        for line in in_vcf:
            if line.startswith("#"):
                out_vcf.write(line)
                continue
            
            total_sites_read += 1
            x = line.strip().split("\t")
            CHROM = x[0]
            POS = x[1]
            key = f"{CHROM}_{POS}"

            if key not in site_data:
                out_vcf.write(line) 
                continue

            ancestral_probability = site_data[key]
            ref_allele = x[3]
            alt_allele = x[4]
            
            # 统计 Major Allele
            c0 = 0
            c1 = 0
            for field in x[9:]:
                gt = field.split(":")[0]
                c0 += gt.count("0")
                c1 += gt.count("1")
            
            if c0 >= c1:
                major_allele = ref_allele
                minor_allele = alt_allele
            else:
                major_allele = alt_allele
                minor_allele = ref_allele

            # 判定 Ancestral Allele (概率 > 0.9)
            if ancestral_probability > 0.9:
                ancestral_allele = major_allele
            else:
                ancestral_allele = minor_allele

            # 极性化
            if ancestral_allele == ref_allele:
                final_fields = x[9:] 
            else:
                # 交换 Ref/Alt
                x[3] = ancestral_allele
                x[4] = ref_allele
                
                # 翻转基因型
                new_genotypes = []
                for sample_field in x[9:]:
                    field_parts = sample_field.split(":")
                    field_parts[0] = field_parts[0].translate(trans_table)
                    new_genotypes.append(":".join(field_parts))
                final_fields = new_genotypes

            out_vcf.write("\t".join(x[0:9]) + "\t" + "\t".join(final_fields) + "\n")
            ancestral_out.write(f"{CHROM}\t{POS}\t{ref_allele}\t{ancestral_allele}\n")
            polarized_sites += 1

    print("VCF file has been polarized based on ancestral states.")
    print(f"Total sites processed in VCF: {total_sites_read}")
    print(f"Successfully polarized sites: {polarized_sites}")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
