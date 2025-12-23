#!/usr/bin/env python3

"""
indv_GT_stats_v2.py
-------------------
高效统计 VCF 文件中每个样本的基因型 (hom0, hom1, het, missing, count_1)。
支持 .vcf 和 .vcf.gz，使用流式读取，内存占用极低

Usage:
    python indv_GT_stats_v2.py <input_vcf> [<Type>] [--output <output_file>]
"""

import sys
import gzip
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate genotype statistics from a VCF file efficiently.")
    parser.add_argument("input_vcf", help="Input VCF file (supports .vcf and .vcf.gz)")
    parser.add_argument("Type", nargs="?", default="", help="Optional filter string to match in INFO field")
    parser.add_argument("--output", help="Output file path. Defaults to 'indv_GT_stats_res.txt'")
    return parser.parse_args()

def open_file(filename):
    """自动判断是否为 gzip 文件并打开"""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def main():
    args = parse_args()
    vcf_file = args.input_vcf
    variant_type = args.Type
    output_file = args.output if args.output else "indv_GT_stats_res.txt"

    # 初始化统计字典: stats[sample_index] = { 'hom0': 0, ... }
    stats = []
    samples = []
    
    print(f"Processing {vcf_file}...")
    
    try:
        with open_file(vcf_file) as f:
            for line in f:
                # 1. 跳过注释行 (##)
                if line.startswith("##"):
                    continue
                
                # 2. 处理表头行 (#CHROM)
                if line.startswith("#CHROM"):
                    header = line.strip().split("\t")
                    samples = header[9:]
                    # 初始化每个样本的统计数据
                    for _ in samples:
                        stats.append({'hom0': 0, 'hom1': 0, 'het': 0, 'missing': 0, 'count_1': 0, 'total': 0})
                    continue

                # 3. 处理数据行
                if not line.strip(): continue # 跳过空行
                
                parts = line.strip().split("\t")
                
                # 如果指定了 Type 过滤，检查 INFO 列 (第8列，索引7)
                if variant_type and variant_type not in parts[7]:
                    continue
                
                # 遍历每个样本的基因型 (从第10列开始，索引9)
                # FORMAT 在索引 8，样本数据从索引 9 开始
                # 通常 GT 是 FORMAT 的第一个字段，例如 GT:AD:DP:GQ...
                
                for i, sample_data in enumerate(parts[9:]):
                    # 获取 GT 部分 (通常是第一个字段，以 : 分割)
                    # 考虑到某些 VCF 可能只写了 GT，或者不规范，这里做简单处理
                    gt_str = sample_data.split(":")[0]

                    # 统计逻辑
                    # 标准化分隔符：将 | (phased) 替换为 / (unphased) 以便统一处理
                    gt_norm = gt_str.replace("|", "/")
                    
                    if gt_norm == "0/0":
                        stats[i]['hom0'] += 1
                        stats[i]['total'] += 1
                    elif gt_norm == "1/1":
                        stats[i]['hom1'] += 1
                        stats[i]['count_1'] += 2
                        stats[i]['total'] += 1
                    elif gt_norm in ("0/1", "1/0"):
                        stats[i]['het'] += 1
                        stats[i]['count_1'] += 1
                        stats[i]['total'] += 1
                    elif gt_norm == "./." or gt_norm == ".":
                        stats[i]['missing'] += 1
                        stats[i]['total'] += 1
                    else:
                        # 处理多等位基因 (如 1/2, 0/2) 或其他情况
                        # 这里简单处理：如果包含 ., 算 missing；否则视情况而定
                        if "." in gt_norm:
                            stats[i]['missing'] += 1
                        else:
                            # 简单的多等位基因处理逻辑 (可视需求修改)
                            # 假设只统计是否含 1 (Alt1)
                            alleles = gt_norm.split("/")
                            if len(alleles) == 2:
                                c1 = alleles.count("1")
                                if c1 == 2: # 1/1
                                     stats[i]['hom1'] += 1
                                     stats[i]['count_1'] += 2
                                elif c1 == 1: # 0/1, 1/2 等
                                     stats[i]['het'] += 1
                                     stats[i]['count_1'] += 1
                                else: # 0/0, 2/2 等 (非1的纯合或杂合)
                                     # 非0/0非1/1非0/1的暂归类为 hom0 或忽略
                                     # 严谨起见，建议只关注 0 和 1
                                     if alleles[0] == alleles[1]:
                                         stats[i]['hom0'] += 1 # 如 2/2 这种
                                     else:
                                         stats[i]['het'] += 1
                                stats[i]['total'] += 1

    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # 写入结果
    print(f"Writing results to {output_file}...")
    with open(output_file, "w") as out:
        out.write("indv\thom0\thom1\thet\tmissing\ttotal\tcount_1\n")
        for i, sample in enumerate(samples):
            s = stats[i]
            out.write(f"{sample}\t{s['hom0']}\t{s['hom1']}\t{s['het']}\t{s['missing']}\t{s['total']}\t{s['count_1']}\n")

    print("Done.")

if __name__ == "__main__":
    main()
