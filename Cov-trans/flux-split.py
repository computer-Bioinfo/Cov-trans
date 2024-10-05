input_file = "/home/bioinfo/guoxy/cov-trans-data/simulator/flux_non/mers/fluxSimulator.fasta"  # 输入的FASTA文件路径
output_file1 = "/home/bioinfo/guoxy/cov-trans-data/simulator/flux_non/mers/sim9-mers-1.fasta"  # 第一个输出文件路径
output_file2 = "/home/bioinfo/guoxy/cov-trans-data/simulator/flux_non/mers/sim9-mers-2.fasta"  # 第二个输出文件路径



read1_lines = []  # 存储以"/1"结尾的序列行
read2_lines = []  # 存储以"/2"结尾的序列行

with open(input_file, 'r') as input_file:
    current_read = None
    for line in input_file:
        if line.startswith('>'):
            current_read = line.strip().endswith("/1")
        if current_read:
            read1_lines.append(line)
        else:
            read2_lines.append(line)

# 写入以"/1"结尾的序列到第一个输出文件
with open(output_file1, 'w') as output1:
    output1.writelines(read1_lines)

# 写入以"/2"结尾的序列到第二个输出文件
with open(output_file2, 'w') as output2:
    output2.writelines(read2_lines)
