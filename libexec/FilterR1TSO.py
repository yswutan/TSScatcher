import pysam
import editdistance
import re
import sys

def filter_bam(input_bam, output_bam_all_R1):
    # 定义正链和负链的目标序列
    TS_oligo_last5_pos = "ATGGG"  # 正链目标序列
    TS_oligo_last5_neg = "CCCAT"  # 负链目标序列 (反向互补)

    # 打开 BAM 文件
    samfile = pysam.AlignmentFile(input_bam, "rb")
    outfile1 = pysam.AlignmentFile(output_bam_all_R1, "wb", template=samfile)

    for read in samfile.fetch():

        cigar = read.cigarstring  # 获取 CIGAR
        sequence = read.query_sequence  # 获取序列

        # 使用正则解析 CIGAR 字符串为数值和字符列表
        cigar_parts = re.findall(r'\d+|\D+', cigar)

        if read.is_reverse:
            # **负链处理**：检查 CIGAR 是否以 softclip 结尾，且长度在 6-20 之间
            if cigar.endswith('S'):
                softclip_len = int(cigar.split('S')[-2].split('M')[-1])
                if 6 <= softclip_len <= 20:
                    # 获取 softclip 序列并反转顺序
                    softclip_seq = sequence[-softclip_len:]  # 反转序列

                    # 检查匹配区域的长度
                    match_len = int(cigar_parts[-4]) if cigar_parts[-3] == 'M' else 0
                    if match_len > 5:
                        # 计算编辑距离（负链用 "CCCAT"）
                        for i in range(len(softclip_seq) - 4):  # len(softclip_seq) - 4 保证不越界
                            sub_seq = softclip_seq[i:i+5]
                            dist = editdistance.eval(sub_seq, TS_oligo_last5_neg)
                            if dist < 4:
                                outfile1.write(read)
                                break  # 只要找到一个符合条件的子串就保留该 read

        else:
            # **正链处理**：检查 CIGAR 是否以 softclip 开头，且长度在 6-20 之间
            if cigar_parts[1] == 'S':
                softclip_len = int(cigar_parts[0])
                if 6 < softclip_len < 20:
                    softclip_seq = sequence[:softclip_len]  # 获取 softclip 序列

                    # 检查匹配区域的长度
                    match_len = int(cigar_parts[2]) if cigar_parts[3] == 'M' else 0
                    if match_len > 5:
                        # 计算编辑距离（正链用 "ATGGG"）
                        for i in range(len(softclip_seq) - 5, -1, -1):
                            sub_seq = softclip_seq[i:i+5]
                            dist = editdistance.eval(sub_seq, TS_oligo_last5_pos)
                            if dist < 4:
                                outfile1.write(read)
                                break  # 只要找到一个符合条件的子串就保留该 read

    # 关闭 BAM 文件
    samfile.close()
    outfile1.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_bam> <output_bam_all_R1>")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam_all_R1 = sys.argv[2]

    filter_bam(input_bam, output_bam_all_R1)