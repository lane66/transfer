import cyvcf2
import csv
import sys

def process_vcf(input_vcf, output_csv):
    # 打开 VCF 文件
    vcf = cyvcf2.VCF(input_vcf)

    # 打开输出文件
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        
        # 写入表头
        writer.writerow(['START', 'STOP', 'REF', 'ALT', 'FREQ'])
        
        # 遍历 VCF 文件中的每个变异
        for variant in vcf:
            # 提取信息
            start = variant.POS
            ref = variant.REF
            alts = variant.ALT
            ao_info = variant.INFO.get('AO', 0)  # 默认值为0
            
            # 如果 AO 是整数，将其转换为列表
            if isinstance(ao_info, int):
                ao_list = [ao_info]
            else:
                ao_list = list(ao_info)
            
            # 确保 ao_list 的长度与 alts 的长度一致
            if len(ao_list) < len(alts):
                ao_list = ao_list + [0] * (len(alts) - len(ao_list))
            
            dp = 2298  # 总测序读数
            
            # 遍历每个 ALT 和对应的 AO
            for alt, ao in zip(alts, ao_list):
                # 计算终止位点
                stop = start + len(ref) - 1
                
                # 计算频率
                freq = ao / dp if dp != 0 else 0
                
                # 写入一行数据
                writer.writerow([start, stop, ref, alt, freq])

    print(f"转换完成！输出文件为 {output_csv}")

if __name__ == "__main__":
    # 检查命令行参数数量
    if len(sys.argv) != 3:
        print("用法: python read_vcf.py <input_vcf> <output_csv>")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_csv = sys.argv[2]
    
    # 调用处理函数
    process_vcf(input_vcf, output_csv)