🧬 BEAGLE to VCF Converter (convert_bgl_to_vcf.py)

This Python script converts phased BEAGLE v3.1 genotype files (.bgl) to the standardized VCF (Variant Call Format), using variant annotations from a corresponding PLINK .bim file.

⚠️ Heterozygous genotypes preserve allele order (e.g., 0|1 vs. 1|0), which is important for phased datasets.
📂 Input Files

BEAGLE v3.1 file (*.bgl): Phased genotype data in BEAGLE format.
PLINK .bim file (*.bim): Contains SNP metadata (Chromosome, RSID, Position, Alleles)
🧾 Output

A .vcf file compatible with downstream tools and pipelines.
Phased genotypes represented using VCF-style genotype encoding (0|0, 0|1, etc.).
Any SNPs not found in the .bim file or with inconsistent allele encoding will be skipped or marked as missing (./.).

Usage 

python convert_bgl_to_vcf.py <input.bgl> <input.bim> <output.vcf>

Example

python convert_bgl_to_vcf.py example.bgl example.bim output.vcf

🔧 Requirements

Python 3.x
No external dependencies (uses only built-in libraries: sys, re)
