#!/bin/bash
# fastp_parallel_processing.sh

# 设置目录
input_dir="/home/data/1.rawdata"
output_dir="/home/vensin/workspace/snpcalling_wild/2.cleaned_data"
report_dir="$output_dir/reports"

# 创建输出目录
mkdir -p "$output_dir" "$report_dir"

echo "Starting parallel processing at $(date)"
echo "Input: $input_dir"
echo "Output: $output_dir"
echo "Reports: $report_dir"

# 运行 parallel.py
python /home/vensin/software/script/parallel.py \
    -i "$input_dir" \
    -o "$output_dir" \
    -r "$report_dir" \
    -a "-3 -5 -l 50 -e 20 --detect_adapter_for_pe -w 2" \
    -1 "_1" \
    -2 "_2" \
    -p 15

echo "Processing completed at $(date)"
