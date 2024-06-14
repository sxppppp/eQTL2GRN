import pandas as pd
import subprocess
import os
import argparse
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Run GRN analysis.")
    parser.add_argument('-v', '--vcf', required=True, help='Vcf_file snp.vcf')
    parser.add_argument('-e', '--eQTL_file', required=True, help='eQTL_file eQTL.txt: col1:gene\tcol2:snp\tcol3:p_value')
    parser.add_argument('-b', '--bed_file', required=True, help='Bed_file Gene.bed')
    parser.add_argument('-o1', '--output_file1',
                        default='GRN.txt',
                        help='The output file, default: ./GRN.txt')
    parser.add_argument('-o2', '--output_file2',
                        default='Visualize.txt',
                        help='The output file, default: ./Visualize.txt')
    return parser.parse_args()


# Step 1: 提取eQTL中的SNP
def extract_snp(eQTL_file, vcf):
    snp_eQTL = pd.read_csv(eQTL_file, delimiter='\t', usecols=[1]).squeeze()
    snp_extract = set(snp_eQTL.str.replace(':', '\t'))

    with open(vcf, "r") as in_vcf, open("eQTL_snp.vcf", "w") as out_vcf:
        for line in in_vcf:
            if line.startswith("#"):
                out_vcf.write(line)
            else:
                break

        in_vcf.seek(0)
        for line in in_vcf:
            if not line.startswith("#"):
                chrom_pos = "\t".join(line.split("\t")[:2])
                if chrom_pos in snp_extract:
                    out_vcf.write(line)


# Step 2: 使用plink进行LD分析
def run_plink():
    subprocess.run(
        ["plink", "--noweb", "--vcf", "eQTL_snp.vcf", "--allow-extra-chr", "--make-bed", "--out", "eQTL_snp"])
    subprocess.run(["plink", "--noweb", "--bfile", "eQTL_snp", "--allow-extra-chr", "--make-bed", "--out", "eQTL_snp"])
    subprocess.run(
        ["plink", "--bfile", "eQTL_snp", "--allow-extra-chr", "--ld-window-kb", "1000", "--ld-window-r2", "0.1",
         "--ld-window", "99999", "--r2"])


# Step 3: 处理plink.ld文件生成snp-snp.txt
def process_ld_file():
    with open("plink.ld", "r") as ld_file, open("snp-snp.txt", "w") as snp_snp_file:
        for line in ld_file:
            parts = line.strip().split()
            snp1 = parts[2]
            snp2 = parts[5]
            snp_snp_file.write(f"{snp1}\t{snp2}\n")


# Step 4: 处理snp-snp.txt文件生成output_snp_classified.txt
def process_snp_file(snp_snp_file, output_snp_file):
    snp_clusters = {}
    with open(snp_snp_file, "r") as f:
        for line in f:
            snp1, snp2 = line.strip().split()
            if snp1 not in snp_clusters:
                snp_clusters[snp1] = set()
            if snp2 not in snp_clusters:
                snp_clusters[snp2] = set()
            snp_clusters[snp1].add(snp2)
            snp_clusters[snp2].add(snp1)

    visited = set()
    clusters = []
    for snp in snp_clusters:
        if snp not in visited:
            cluster = []
            stack = [snp]
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    cluster.append(node)
                    stack.extend(snp_clusters[node] - visited)
            clusters.append(cluster)

    with open(output_snp_file, "w") as f:
        for cluster in clusters:
            f.write(",".join(cluster) + ";\n")

    with open(output_snp_file, "r") as f:
        lines = f.readlines()
    with open(output_snp_file, "w") as f:
        f.writelines(lines[1:])


# Step 5: 创建SNP-基因字典，包含p_value
def create_snp_gene_dict(eQTL_file):
    snp_gene_dict = {}
    with open(eQTL_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            gene = parts[0]
            snp = parts[1]
            p_value = parts[2]
            if snp in snp_gene_dict:
                snp_gene_dict[snp].append((gene, p_value))
            else:
                snp_gene_dict[snp] = [(gene, p_value)]
    return snp_gene_dict


# Step 6: 创建基因位置字典
def create_gene_position_dict(bed_file):
    gene_position_dict = {}
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3]
            gene_position_dict[gene] = (chrom, start, end)
    return gene_position_dict


# Step 7: 判断SNP与基因的位置关系
def classify_snp_gene(snp, gene, gene_position_dict):
    if gene not in gene_position_dict:
        return 'trans'
    chrom, start, end = gene_position_dict[gene]
    snp_chrom, snp_pos = snp.split(':')
    snp_pos = int(snp_pos)
    if snp_chrom == chrom and start - 1000000 <= snp_pos <= end + 1000000:
        return 'cis'
    else:
        return 'trans'


# Step 8: 处理snp聚类文件并输出基因调控网络GRN.txt
def process_snp_cluster_file1(snp_cluster_file1, snp_gene_dict, gene_position_dict, output_file1):
    with open(snp_cluster_file1, 'r') as f1, open(output_file1, 'w') as out_f1:
        out_f1.write("Cluster\tGene_Count\tGene_SNP_Pairs\n")

        cluster_number1 = 1
        for line in f1:
            snps1 = line.strip().split(',')
            gene_snp_dict = {}
            for snp in snps1:
                if snp in snp_gene_dict:
                    for gene, _ in snp_gene_dict[snp]:
                        classification1 = classify_snp_gene(snp, gene, gene_position_dict)
                        if gene not in gene_snp_dict:
                            gene_snp_dict[gene] = {'cis': [], 'trans': []}
                        gene_snp_dict[gene][classification1].append(f"{snp}({classification1})")

            gene_snp_strs = []
            for gene, snp_dict in gene_snp_dict.items():
                snp_list = snp_dict['cis'] + snp_dict['trans']
                snp_str = ','.join(snp_list)
                gene_snp_strs.append(f"{gene}[{snp_str}]")

            gene_count = len(gene_snp_dict)
            gene_snp_strs = ','.join(gene_snp_strs)
            out_f1.write(f"Cluster{cluster_number1}\t{gene_count}\t{gene_snp_strs}\n")
            cluster_number1 += 1



# Step 9: 处理snp聚类文件并输出网络可视化Visualize.txt
def process_snp_cluster_file2(snp_cluster_file2, snp_gene_dict, gene_position_dict, output_file2):
    with open(snp_cluster_file2, 'r') as f2, open(output_file2, 'w') as out_f2:
        out_f2.write("Cluster\tGene\tSNP\tType\tP_value\n")
        cluster_number2 = 1
        for line in f2:
            snps = line.strip().split(',')
            gene_snp_list2 = []
            for snp in snps:
                if snp in snp_gene_dict:
                    for gene, p_value in snp_gene_dict[snp]:
                        snp_type = classify_snp_gene(snp, gene, gene_position_dict)
                        gene_snp_list2.append((gene, snp, snp_type, p_value))

            gene_snp_list2.sort(key=lambda x: (x[0], x[1], 0 if x[2] == 'cis' else 1))
            
            for gene, snp, snp_type, p_value in gene_snp_list2:
                out_f2.write(f"Cluster{cluster_number2}\t{gene}\t{snp}\t{snp_type}\t{p_value}\n")
            cluster_number2 += 1


def main():
    args = parse_args()
    temp_snp_classified_file = "output_snp_classified.txt"

    # 提取SNP
    extract_snp(args.eQTL_file, args.vcf)

    # 运行plink
    run_plink()

    # 处理plink.ld文件生成snp-snp.txt
    process_ld_file()

    # 处理snp-snp.txt文件生成temp_snp_classified.txt
    process_snp_file("snp-snp.txt", temp_snp_classified_file)

    # 创建SNP-基因和基因位置字典
    snp_gene_dict = create_snp_gene_dict(args.eQTL_file)
    gene_position_dict = create_gene_position_dict(args.bed_file)

    # 处理SNP聚类文件并输出基因调控网络文件GRN.txt
    process_snp_cluster_file1(temp_snp_classified_file, snp_gene_dict, gene_position_dict, args.output_file1)

    # 处理SNP聚类文件并输出网络可视化Visualize.txt
    process_snp_cluster_file2(temp_snp_classified_file, snp_gene_dict, gene_position_dict, args.output_file2)

    # 清理临时文件
    for temp_file in glob.glob("plink*") + glob.glob("eQTL_snp*") + ["snp-snp.txt", temp_snp_classified_file]:
        if os.path.exists(temp_file):
            os.remove(temp_file)

    print("分析完成，基因调控网络文件结果保存到", args.output_file1)
    print("分析完成，网络可视化文件结果保存到", args.output_file2)


if __name__ == '__main__':
    main()
