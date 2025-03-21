import sys
import re

def load_bim_file(bim_file):
    """Loads a PLINK .bim file into a dictionary mapping RSID to (chrom, pos, ref, alt)."""
    bim_data = {}
    with open(bim_file, "r", encoding="utf-8") as bim:
        for line in bim:
            parts = re.split(r'\s+', line.strip())  # Handles spaces or tabs
            if len(parts) < 6:
                print(f"⚠️ Warning: Skipping malformed line in .bim file: {line.strip()}")
                continue
            chrom, rsid, _, pos, ref, alt = parts[:6]
            bim_data[rsid] = (chrom, pos, ref, alt)
    return bim_data

def convert_bgl_to_vcf(bgl_file, bim_file, output_vcf):
    """Converts a phased BEAGLE v3.1 file to a VCF file while preserving allele order in heterozygous cases."""
    bim_data = load_bim_file(bim_file)

    with open(bgl_file, "r", encoding="utf-8") as infile, open(output_vcf, "w", encoding="utf-8") as outfile:
        lines = infile.readlines()
        if len(lines) < 6:
            print("⚠️ Error: BEAGLE file seems too short or malformed.")
            return

        # Extract sample IDs from the second line (skip first two columns)
        sample_ids = re.split(r'\s+', lines[1].strip())[2::2]

        # Debugging: print extracted sample IDs
        print(f"🔍 Extracted {len(sample_ids)} sample IDs.")

        # Write VCF header
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids) + "\n")

        snp_count, skipped_snps = 0, 0

        for line in lines[5:]:  # Skip metadata lines (first five lines)
            parts = re.split(r'\s+', line.strip())
            if len(parts) < 3 or parts[0] != "M":
                continue  # Skip if not a marker line
            
            rsid = parts[1]
            if rsid not in bim_data:
                print(f"⚠️ Warning: SNP {rsid} not found in .bim file, skipping.")
                skipped_snps += 1
                continue
            
            chrom, pos, ref, alt = bim_data[rsid]
            genotypes = parts[2:]
            phased_genotypes = []

            for i in range(0, len(genotypes), 2):
                try:
                    allele1, allele2 = genotypes[i], genotypes[i + 1]

                    # 🛠 FIX: Handle cases where alleles in .bgl might not match .bim
                    if allele1 not in (ref, alt) or allele2 not in (ref, alt):
                        print(f"⚠️ Warning: Unexpected alleles {allele1}, {allele2} at {rsid}. Assigning missing genotype.")
                        phased_genotypes.append("./.")
                        continue

                    # ✅ Preserve allele order in heterozygous genotypes
                    if allele1 == ref and allele2 == ref:
                        phased_genotypes.append("0|0")
                    elif allele1 == alt and allele2 == alt:
                        phased_genotypes.append("1|1")
                    elif allele1 == ref and allele2 == alt:  
                        phased_genotypes.append("0|1")  # `ref alt` → `0|1`
                    elif allele1 == alt and allele2 == ref:  
                        phased_genotypes.append("1|0")  # `alt ref` → `1|0`
                    else:
                        print(f"⚠️ Warning: Unexpected alleles ({allele1}, {allele2}) at SNP {rsid}. Marking as missing.")
                        phased_genotypes.append("./.")

                except IndexError:
                    print(f"⚠️ Error: Genotype parsing issue at SNP {rsid} (line: {line.strip()})")
                    phased_genotypes.append("./.")

            outfile.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(phased_genotypes) + "\n")
            snp_count += 1

        print(f"✅ Conversion complete! {snp_count} SNPs processed, {skipped_snps} SNPs skipped.")

# Run script from command line
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python convert_bgl_to_vcf.py <bgl_file> <bim_file> <output_vcf>")
        sys.exit(1)

    bgl_file = sys.argv[1]
    bim_file = sys.argv[2]
    output_vcf = sys.argv[3]

    convert_bgl_to_vcf(bgl_file, bim_file, output_vcf)
