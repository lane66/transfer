from Bio import SeqIO
import pandas as pd

# 创建字典存储基因和对应的 type
gene_type_dict = {}

# 读取三个 txt 文件并合并数据
for file in ['multi-component_structure.txt', 'single-component_structure.txt', 'two-component_structure.txt']:
    df = pd.read_csv(file, delimiter='\t')  # 假设是制表符分隔
    for index, row in df.iterrows():
        gene = row['gene']
        type_val = row['type']
        gene_type_dict[gene] = type_val

# 创建字典存储 type 对应的序列记录
type_records = {}

# 遍历 FASTA 文件
for record in SeqIO.parse('sarg.fasta', 'fasta'):
    # 提取序列名第四部分作为基因名
    gene = record.id.split(' ')[0]
    # 获取对应的 type
    if gene in gene_type_dict:
        type_val = gene_type_dict[gene]
    else:
        type_val = 'unknown'  # 基因名未找到对应的 type 时的处理方式
    # 添加序列记录到对应的 type 组
    if type_val in type_records:
        type_records[type_val].append(record)
    else:
        type_records[type_val] = [record]

# 将每个 type 组的序列写入对应的 FASTA 文件
for type_val, records in type_records.items():
    filename = f'{type_val}.fasta'
    SeqIO.write(records, filename, 'fasta')
    print(f'已创建文件: {filename}')