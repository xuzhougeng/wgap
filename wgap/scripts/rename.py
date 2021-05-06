# Nomenclature 
def maker_rename(old_gff, new_gff, prefix, non_chr = 'scaffold'):
    f = open(old_gff, "r")
    fo = open(new_gff, "w")
    prev_chrom = 'flag'
    count = 1
    sub_count = 1
    non_chr_count = 0

    for line in f.readlines():
        gene_name = line.strip().split()[0]
        chrom = gene_name.split("-")[1]
        if chrom == prev_chrom:
            if 'mRNA' in gene_name:
                if non_chr in chrom:
                    print(gene_name + "\t{}{}U".format(prefix, non_chr_count) + str(count).zfill(4) + '0.{}'.format(sub_count))
                else:
                    chrom_id = chrom[3] 
                    print(gene_name + "\t{}{}G".format(prefix, chrom_id) + str(count).zfill(4) + '0.{}'.format(sub_count))
                count += 1
                sub_count += 1
            else:
                sub_count = 1
                if non_chr in chrom:
                    fo.write(gene_name + "\t{}{}U".format(prefix, non_chr_count) + str(count).zfill(4) + '0')
                else:
                    chrom_id = chrom[3] 
                    fo.write(gene_name + "\t{}{}G".format(prefix, chrom_id) + str(count).zfill(4) + '0')

            prev_chrom = chrom
        else :
            count = 1 
            if non_chr in chrom:
                non_chr_count += 1
                print(gene_name + "\t{}{}U".format(prefix, non_chr_count) + str(count).zfill(4) + '0')
            else:
                chrom_id = chrom[3] 
                print(gene_name + "\t{}{}G".format(prefix, chrom_id) + str(count).zfill(4) + '0')
            prev_chrom = chrom

    f.close()
    fo.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("usage: {} old_gff, new_gff, prefix, non_chr".format(sys.argv[0]), file = sys.stderr)
        sys.exit(1)
    old_gff = sys.argv[1]
    new_gff = sys.argv[2]
    prefix  = sys.argv[3]
    non_chr = sys.argv[4]
    maker_rename(old_gff, new_gff, prefix, non_chr)
    