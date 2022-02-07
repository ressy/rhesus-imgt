from csv import DictWriter
from Bio import SeqIO
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

# http://www.imgt.org/IMGTindex/Fasta.php
# The IMGT FASTA header of nucleotide IMGT reference sequences contains 15 fields separated by '|':
FIELDS = [
    "Accession",   # 1. IMGT/LIGM-DB accession number(s)
    "Allele",      # 2. IMGT gene and allele name
    "Species",     # 3. species
    "FuncGroup",   # 4. IMGT allele functionality
    "Region",      # 5. exon(s), region name(s), or extracted label(s)
    "StartAndEnd", # 6. start and end positions in the IMGT/LIGM-DB accession number(s)
    "NumNT",       # 7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
    "CodonStart",  # 8. codon start, or 'NR' (not relevant) for non coding labels
    "NT5P",        # 9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
    "NT3P",        # 10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
    "NumNTCorr",   # 11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
    "NumAA",       # 12. number of amino acids (AA): this field indicates that the sequence is in amino acids
    "NumChars",    # 13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
    "Partial",     # 14. partial (if it is)
    "RevCmp"]      # 15. reverse complementary (if it is)

FIELDS_EXTRA = ["Gene", "Family", "Segment", "Locus", "SeqDesc", "Seq"]

URL = "http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"
SEGMENTS = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"]

def parse_fields(txt):
    return {key: item.strip() for key, item in zip(FIELDS, txt.split("|"))}

def parse_allele(txt):
    return {
        "Gene": re.sub(r"\*.*$", "", txt),
        "Family": re.sub("(IG[HKL][VDJ][0-9]+).*$", r"\1", txt),
        "Segment": txt[0:4],
        "Locus": txt[0:3]}

rule tabulate_segments_rhesus:
    input: "output/Macaca_mulatta.csv"

rule tabulate_segments:
    output: "output/{organism}.csv"
    input: expand("from-imgt/{{organism}}/{segment}.fasta", segment = SEGMENTS)
    run:
        with open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=FIELDS + FIELDS_EXTRA, lineterminator="\n")
            writer.writeheader()
            for input_fasta in input:
                with open(input_fasta) as f_in:
                    for record in SeqIO.parse(f_in, "fasta"):
                        row = parse_fields(record.description)
                        row.update(parse_allele(row["Allele"]))
                        row["SeqDesc"] = record.description
                        row["Seq"] = str(record.seq)
                        writer.writerow(row)

rule get_fasta:
    output: "from-imgt/{organism}/{segment}.fasta"
    input: HTTP.remote(expand("{url}/{{organism}}/IG/{{segment}}.fasta", url=URL), keep_local=True)
    run:
        shell("mv {input} {output}")
