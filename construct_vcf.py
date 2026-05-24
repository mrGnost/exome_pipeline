import sys
import pandas as pd
import gzip

def parse_vcf_line(line):
    parts = line.strip().split('\t')
    if len(parts) >= 2:
        chrom = parts[0]
        pos = parts[1]
        return chrom, pos
    return None, None

def is_gzipped(filepath):
    return filepath.endswith('.gz')

def open_file(filepath):
    if is_gzipped(filepath):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def construct_vcf(vcf_files, binary_csv, output_vcf):
    binary_df = pd.read_csv(binary_csv)
    binary_df['chrom'] = binary_df['id'].apply(lambda x: x.split(':')[0])
    binary_df['pos'] = binary_df['id'].apply(lambda x: x.split(':')[1])
    selected_variants = binary_df[binary_df['binary'] == 1]
    variant_set = set(zip(selected_variants['chrom'], selected_variants['pos']))
    
    header_lines = []
    with open_file(vcf_files[0]) as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line)
            elif line.startswith('#CHROM'):
                header_lines.append(line)
                break
    
    variant_lines = []
    seen_variants = set()
    for vcf_file in vcf_files:
        with open_file(vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                chrom, pos = parse_vcf_line(line)
                if chrom is not None and (chrom, pos) in variant_set and (chrom, pos) not in seen_variants:
                    variant_lines.append((chrom, int(pos), line))
                    seen_variants.add((chrom, pos))
    
    variant_lines.sort(key=lambda x: (x[0], x[1]))
    
    with open(output_vcf, 'w') as out_f:
        for header_line in header_lines:
            out_f.write(header_line)
        
        for chrom, pos, line in variant_lines:
            out_f.write(line)

if __name__ == '__main__':
    vcf1 = sys.argv[1]
    vcf2 = sys.argv[2]
    vcf3 = sys.argv[3]
    binary_csv = sys.argv[4]
    output_vcf = sys.argv[5]
    
    construct_vcf([vcf1, vcf2, vcf3], binary_csv, output_vcf)
