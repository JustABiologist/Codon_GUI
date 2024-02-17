from Bio import SeqIO
from collections import Counter
from Bio.Seq import Seq

def aggregate_codon_counts(gb_file, keyword):
    """Aggregate codon counts for all entries with a specific keyword in their 'note', across all feature types."""
    codon_counts = Counter()
    found = False  # Debugging flag
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if 'note' in feature.qualifiers and keyword in feature.qualifiers['note'][0]:
                found = True  # Set flag if keyword is found
                # Only process features with a sequence (e.g., 'CDS', 'gene', but open to others with sequences)
                if len(feature.location) % 3 == 0:  # Check if divisible by 3 to ensure it's likely coding
                    sequence = str(feature.extract(record.seq))
                    for i in range(0, len(sequence), 3):
                        codon = sequence[i:i+3].upper()
                        if len(codon) == 3:
                            codon_counts[codon] += 1
    if not found:
        print(f"No features with keyword '{keyword}' in 'note' were found.")  # Debugging message
    return codon_counts

def compute_codon_usage_from_counts(codon_counts):
    """Compute codon usage frequencies from aggregated codon counts."""
    total_codons = sum(codon_counts.values())
    codon_usage = {codon: count / total_codons for codon, count in codon_counts.items()}
    return codon_usage

def get_codon_counts_for_gene_by_locus(gb_file, locus_tag):
    """Retrieve and compute codon counts for a gene by locus tag."""
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == locus_tag:
                sequence = str(feature.extract(record.seq))
                codon_counts = Counter()
                for i in range(0, len(sequence), 3):
                    codon = sequence[i:i+3].upper()
                    if len(codon) == 3:
                        codon_counts[codon] += 1
                return codon_counts
    return Counter()

def compare_codon_usage(base_codon_usage, gene_codon_usage):
    """Prints the comparison of codon usage between the baseline and a gene, including amino acid translations, sorted by amino acid."""
    print("Codon\tAmino Acid\tBaseline Usage\tGene Usage")
    all_codons = set(base_codon_usage.keys()) | set(gene_codon_usage.keys())

    # Create a list of tuples (amino_acid, codon) for sorting
    codon_amino_acid_pairs = [(Seq(codon).translate(), codon) for codon in all_codons]
    # Sort by amino acid first, then by codon
    sorted_codons = sorted(codon_amino_acid_pairs)

    for amino_acid, codon in sorted_codons:
        base_usage = base_codon_usage.get(codon, 0) * 100
        gene_usage = gene_codon_usage.get(codon, 0) * 100
        print(f"{codon}\t{amino_acid}\t\t{base_usage:.4f}\t\t{gene_usage:.4f}")

def validate_codon_usage_percentages(codon_usage):
    """Validate that the sum of codon usage percentages is close to 100%."""
    total_percentage = sum(codon_usage.values())
    if 99.9 <= total_percentage <= 100.1:
        print("Validation passed: Codon usage percentages sum to approximately 100%.")
    else:
        print(f"Validation warning: Codon usage percentages sum to {total_percentage}%, which is outside the expected range.")


# Example usage
#gb_file = "/Users/floriangrun/Desktop/Bachelorarbeit_Koch/genome/GCA_001672515.gb"
#keyword = "Effector"  # Change this to your specific keyword
#locus_tag = "CH63R_02802"
#
# Step 1 & 2: Aggregate codon counts and compute codon usage frequencies for the baseline
#baseline_codon_counts = aggregate_codon_counts(gb_file, keyword)
#baseline_codon_usage = compute_codon_usage_from_counts(baseline_codon_counts)

# Step 3: Retrieve the gene of interest by locus tag and compute its codon usage
#gene_codon_counts = get_codon_counts_for_gene_by_locus(gb_file, locus_tag)
#gene_codon_usage = compute_codon_usage_from_counts(gene_codon_counts)

# Step 4: Compare the baseline to your gene of interest
#compare_codon_usage(baseline_codon_usage, gene_codon_usage)

#validate_codon_usage_percentages(baseline_codon_usage)
#validate_codon_usage_percentages(gene_codon_usage)