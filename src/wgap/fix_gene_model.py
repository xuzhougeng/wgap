import sys

def process_gff(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # 跳过注释行
            if line.startswith('#'):
                continue

            # 写入原始行
            outfile.write(line)

            # 处理 CDS 行
            if '\tCDS\t' in line:
                fields = line.strip().split('\t')
                fields[2] = 'exon'  # 替换类型为 exon
                fields[8] = fields[8].replace('CDS', 'exon')  # 替换 ID 中的 CDS 为 exon
                new_line = '\t'.join(fields) + '\n'
                outfile.write(new_line)



# 主函数
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.gff output.gff")
        sys.exit(1)

    input_gff = sys.argv[1]
    output_gff = sys.argv[2]

    process_gff(input_gff, output_gff)
