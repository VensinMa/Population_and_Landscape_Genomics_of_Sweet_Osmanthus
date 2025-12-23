#!/usr/bin/env python3
# usage: python vcf_polarize.py <input_vcf/vcf.gz> <estsfs_output> <positions_file>

import sys
import argparse
import os
import gzip

# 创建解析器
def create_parser():
    parser = argparse.ArgumentParser(
        description="This script polarizes a VCF file based on ancestral allele probabilities from EstSFS output."
    )
    parser.add_argument('input_vcf', type=str, help="Input VCF file to be polarized")
    parser.add_argument('estsfs_output', type=str, help="EstSFS output file containing ancestral allele probabilities")
    parser.add_argument('positions_file', type=str, help="File containing the positions of sites to be polarized from vcf_to_estsfs.py output")

    return parser

# 解析命令行参数
parser = create_parser()
args = parser.parse_args()

# 获取不带路径的文件名，确保输出在当前目录
input_filename = os.path.basename(args.input_vcf)
base_name = input_filename.replace(".vcf.gz", "").replace(".vcf", "")

output_vcf_name = f"{base_name}_polarized.vcf"
ancestral_out_name = f"{base_name}_ancestral.txt"

# 自动判断 VCF 打开方式
open_func = gzip.open if args.input_vcf.endswith(".gz") else open
mode = "rt" if args.input_vcf.endswith(".gz") else "r"

try:
    # 同时读取 positions 和 estsfs 文件，构建查找字典
    # 字典格式: {'CHROM_POS': probability_value}
    # 这样可以防止行数对不齐的问题，且查找速度更快
    site_data = {}
    
    with open(args.positions_file, "r") as pos_f, open(args.estsfs_output, "r") as est_f:
        pos_lines = pos_f.readlines()
        est_lines = est_f.readlines()
        
        # 确保两个文件行数一致（通常由 vcf_to_estsfs 生成时是一致的）
        if len(pos_lines) != len(est_lines):
            print("Warning: positions file and estsfs output have different line counts.")
        
        # 遍历较短的那个长度，避免越界
        min_len = min(len(pos_lines), len(est_lines))
        
        for i in range(min_len):
            pos_line = pos_lines[i].strip()
            if not pos_line: continue
            
            est_line = est_lines[i].strip()
            if not est_line: continue
            
            # 解析位置
            parts = pos_line.split("\t")
            chrom, pos = parts[0], parts[1]
            key = f"{chrom}_{pos}"
            
            # 解析概率 (假设 P-major-ancestral 在第3列，即索引2)
            est_parts = est_line.split()
            # 某些行可能以 metadata 开头，确保数据列足够
            if len(est_parts) > 2:
                prob = est_parts[2]
                site_data[key] = float(prob)

    # 打开主 VCF 和输出文件
    with open_func(args.input_vcf, mode) as in_vcf, \
         open(output_vcf_name, "w") as out_vcf, \
         open(ancestral_out_name, "w") as ancestral_out:

        ancestral_out.write("CHROM\tPOS\tref_allele\tancestral_allele\n")

        # 计数器
        total_sites_read = 0
        polarized_sites = 0
        
        # 定义基因型转换表: 0->1, 1->0. 保持分隔符(/或|)不变
        trans_table = str.maketrans("01", "10")

        for line in in_vcf:
            if line.startswith("#"):
                out_vcf.write(line) # 直接写入注释行
                continue
            
            total_sites_read += 1
            x = line.strip().split("\t")
            CHROM = x[0]
            POS = x[1] # 保持字符串格式即可
            key = f"{CHROM}_{POS}"

            # 1. 检查该位点是否在需要处理的列表中
            if key not in site_data:
                # 如果不在列表里，通常直接原样输出，或者可以选择跳过
                # 这里选择原样输出未极化的行，保持 VCF 完整性
                out_vcf.write(line)
                continue

            # 获取概率
            ancestral_probability = site_data[key]

            # 2. 获取参考和变异等位基因
            ref_allele = x[3]
            alt_allele = x[4]
            
            # 3. 统计内群基因型以确定 Major Allele
            # 注意：这里仅统计 0 和 1，忽略 . (缺失)
            c0 = 0
            c1 = 0
            for field in x[9:]:
                gt = field.split(":")[0] # 获取 GT 部分，忽略 DP, GQ 等
                c0 += gt.count("0")
                c1 += gt.count("1")
            
            # 判定 Major Allele
            if c0 >= c1:
                major_allele = ref_allele
                minor_allele = alt_allele
            else:
                major_allele = alt_allele
                minor_allele = ref_allele

            # 4. 判断祖先等位基因
            # 逻辑：如果 est-sfs 计算出 Major 为祖先的概率 > 0.9，则 Ancestral = Major
            if ancestral_probability > 0.9:
                ancestral_allele = major_allele
            else:
                ancestral_allele = minor_allele

            # 5. 极性化操作
            if ancestral_allele == ref_allele:
                # 祖先就是参考，不需要改变基因型
                final_fields = x[9:] 
            else:
                # 祖先是 Alt，需要交换 Ref/Alt 和基因型
                # 交换 Ref 和 Alt 列
                x[3] = ancestral_allele
                x[4] = ref_allele
                
                # 转换所有样本的基因型
                # 使用 translate 高效互换 0 和 1，同时保留 : 后面的深度等信息
                new_genotypes = []
                for sample_field in x[9:]:
                    # 分割 GT 和其他格式字段 (如 0/1:10:99)
                    field_parts = sample_field.split(":")
                    gt_part = field_parts[0]
                    
                    # 仅转换 GT 部分 (0->1, 1->0)
                    new_gt = gt_part.translate(trans_table)
                    
                    # 重新组合
                    field_parts[0] = new_gt
                    new_genotypes.append(":".join(field_parts))
                
                final_fields = new_genotypes

            # 写入 VCF
            # 前9列 + 处理后的基因型列
            out_vcf.write("\t".join(x[0:9]) + "\t" + "\t".join(final_fields) + "\n")
            
            # 写入记录文件
            ancestral_out.write(f"{CHROM}\t{POS}\t{ref_allele}\t{ancestral_allele}\n")
            
            polarized_sites += 1

    # 输出统计信息
    print("VCF file has been polarized based on ancestral states.")
    print(f"Total sites processed in VCF: {total_sites_read}")
    print(f"Successfully polarized sites: {polarized_sites}")

except IOError as e:
    print(f"Error: {e}")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    sys.exit(1)
