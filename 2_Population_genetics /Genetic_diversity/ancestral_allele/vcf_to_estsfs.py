#!/usr/bin/env python3

import sys
import gzip
import os

# Usage: python vcf_to_estsfs.py input.vcf.gz outgroupname1 [outgroupname2] ...

if len(sys.argv) < 3:
    print("Usage: python vcf_to_estsfs.py input.vcf(.gz) outgroupname1 [outgroupname2] ...")
    sys.exit(1)

input_vcf = sys.argv[1]
outgroup_names = sys.argv[2:]

# 使用 os.path.basename 获取不带路径的文件名，使输出文件保存在当前目录
input_filename = os.path.basename(input_vcf)

# 根据输入文件是否为 .gz 决定输出文件名
base_name = input_filename.replace(".vcf.gz", "").replace(".vcf", "")
output_file_estsfs = f"{base_name}_estsfs_input.txt"
output_file_positions = f"{base_name}_estsfs.positions.txt"
log_file = f"{base_name}_processing.log"

total_sites = 0
kept_sites = 0
skipped_sites = []

# 自动判断打开方式
open_func = gzip.open if input_vcf.endswith(".gz") else open
mode = "rt" if input_vcf.endswith(".gz") else "r"

try:
    with open_func(input_vcf, mode) as vcf, \
         open(output_file_estsfs, "w") as estsfs_out, \
         open(output_file_positions, "w") as positions_out, \
         open(log_file, "w") as log_out:

        # 预定义需要的碱基，est-sfs 只需要 ACGT
        target_bases = ["A", "C", "G", "T"]
        
        outgroup_indices = []
        header_processed = False

        for line in vcf:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.strip().split()
                    for name in outgroup_names:
                        try:
                            outgroup_indices.append(samples.index(name))
                        except ValueError:
                            print(f"Error: Outgroup name '{name}' not found in VCF header.")
                            sys.exit(1)
                    header_processed = True
                continue

            if not header_processed:
                continue

            total_sites += 1
            fields = line.strip().split("\t")
            
            # --- STEP 1：基础信息检查 ---
            if len(fields) < 10: continue # 防止空行或错误行
            
            CHROM, POS, ref_allele, alt_allele = fields[0], fields[1], fields[3], fields[4]
            site_info = f"{CHROM}:{POS}"

            # 1. 仅处理二等位位点，且必须是单碱基替换 (SNP)
            if len(ref_allele) != 1 or len(alt_allele) != 1:
                skipped_sites.append(f"Skipping {site_info}: Non-biallelic or Indel ({ref_allele}/{alt_allele})")
                continue
            
            # 2. 检查 Ref/Alt 是否在 ACGT 范围内 (排除 N)
            if ref_allele not in target_bases or alt_allele not in target_bases:
                skipped_sites.append(f"Skipping {site_info}: Non-ACGT bases")
                continue

            # --- STEP 2：外类群处理逻辑 ---
            
            # 1. 提取所有外类群的  GT 字段 (只取冒号前部分，并统一分隔符)
            og_gt_raw = [fields[idx].split(":")[0] for idx in outgroup_indices]
            
            # 2. 检查缺失，缺失位点会被跳过
            if any(gt in ["./.", ".", ".|."] for gt in og_gt_raw):
                skipped_sites.append(f"Skipping {site_info}: Missing data in outgroup")
                continue
                
            # 3. 标准化分隔符，GT 字段可能同时包含 "|" 和 "/" 
            og_gt_normalized = [gt.replace("|", "/") for gt in og_gt_raw]

            # 4. 检查杂合 (est-sfs 通常要求外类群为固定纯合状态，即纯合性和一致性)
            # 外类群基因型必须纯合 (所有外类群基因型必须完全纯合)
            if any(gt.split('/')[0] != gt.split('/')[1] for gt in og_gt_normalized):
                skipped_sites.append(f"Skipping {site_info}: Heterozygous outgroup")
                continue

            # 外类群基因型必须一致 (所有外类群基因型必须完全一致)
            if len(set(og_gt_normalized)) > 1:
                skipped_sites.append(f"Skipping {site_info}: Inconsistent outgroups")
                continue
            
            # 初始化计数字典
            ingroup_counts = {b: 0 for b in target_bases}
            
            # 遍历所有样本列
            for i in range(9, len(fields)):
                if i in outgroup_indices:
                    continue
                
                # 安全获取 GT
                gt_str = fields[i].split(":")[0]
                
                # 显式统计 '0' 和 '1'，不依赖字符串拼接
                # 这样可以处理 '0/0', '0|1', '1|0', '1/1', '0' (haploid) 等情况
                c0 = gt_str.count('0')
                c1 = gt_str.count('1')
                
                # 将计数加到对应的碱基上
                ingroup_counts[ref_allele] += c0
                ingroup_counts[alt_allele] += c1

            # 生成内群字符串 "A,C,G,T"
            ingroup_str = ",".join(str(ingroup_counts[b]) for b in target_bases)

            # --- 外类群输出生成 ---
            # 因为前面检查了基因型一致且纯合，取第一个外类群的基因型即可代表所有外类群
            representative_gt = og_gt_normalized[0]
            
            outgroups_output = []
            # 重新生成每个外类群的 string (虽然它们是一样的，但需要保持格式对应)
            # 或者如果你需要分别输出每个外类群列：
            for gt in og_gt_normalized:
                og_counts = {b: 0 for b in target_bases}
                if gt == "0/0":
                    og_counts[ref_allele] = 1
                elif gt == "1/1":
                    og_counts[alt_allele] = 1
                # 杂合情况已被上面 continue 过滤，此处不需处理
                
                outgroups_output.append(",".join(str(og_counts[b]) for b in target_bases))

            # 写入
            estsfs_out.write(f"{ingroup_str}\t" + " ".join(outgroups_output) + "\n")
            positions_out.write(f"{CHROM}\t{POS}\n")
            kept_sites += 1

        # 写入日志 summary
        log_out.write(f"Total sites processed: {total_sites}\n")
        log_out.write(f"Successfully kept sites: {kept_sites}\n")
        log_out.write("-" * 20 + "\nSkipped details:\n")
        # 限制日志大小，防止几百万行被跳过导致日志过大
        if len(skipped_sites) > 1000:
             log_out.write("\n".join(skipped_sites[:1000]))
             log_out.write(f"\n... and {len(skipped_sites) - 1000} more lines.")
        else:
             log_out.write("\n".join(skipped_sites))

    # 运行结果
    print("VCF file has been converted to est-sfs input file.")
    print(f"Total sites processed: {total_sites}")
    print(f"Successfully kept sites: {kept_sites}")

except IOError as e:
    print(f"Error: {e}")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    sys.exit(1)
