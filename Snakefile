from csv import DictWriter, DictReader
from bs4 import BeautifulSoup
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

URL_VDJ = "http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"
URL_C = "https://www.imgt.org/genedb/GENElect" # see https://www.imgt.org/vquest/refseqh.html
SEGMENTS = ["IGHV", "IGHD", "IGHJ", "IGHC", "IGKV", "IGKJ", "IGLV", "IGLJ"]
ORGANISMS = ["Macaca_mulatta", "Homo_sapiens"]

wildcard_constraints:
    segment="IG[HKL][VDJC]",
    organism="[A-Za-z_]+"

def parse_fields(txt):
    return {key: item.strip() for key, item in zip(FIELDS, txt.split("|"))}

def parse_allele(txt):
    attrs = {
        "Gene": re.sub(r"\*.*$", "", txt),
        "Family": re.sub("(IG[HKL][VDJ][0-9]+).*$", r"\1", txt),
        "Segment": txt[0:4],
        "Locus": txt[0:3]}
    # handle constant region entries with some exceptions
    if txt[:4] in ["IGHA", "IGHD", "IGHE", "IGHG", "IGHM"]:
        # careful, "IGHD*01" means delta constant region while IGHD1-1*01 etc
        # mean D segment alleles
        if txt[3] != "D" or (txt[3] == "D" and "-" not in txt):
            attrs["Family"] = ""
            attrs["Segment"] = "IGHC"
    return attrs

TARGET_FASTAS = expand("output/{organism}.{segment}.fasta", organism=ORGANISMS, segment=SEGMENTS)
TARGET_CSVS = expand("output/{organism}.csv", organism=ORGANISMS)
TARGET_IMGT_FASTAS = expand("from-imgt/{organism}/{segment}.fasta", organism=ORGANISMS, segment=SEGMENTS)

rule all:
    input: TARGET_FASTAS

rule clean:
    run:
        for path in TARGET_FASTAS + TARGET_CSVS:
            shell(f"rm -f '{path}'")

rule realclean:
    run:
        for path in TARGET_FASTAS + TARGET_CSVS + TARGET_IMGT_FASTAS:
            shell(f"rm -f '{path}'")

rule simple_segment_fastas_defaults:
    input: TARGET_FASTAS

rule tabulate_segments_defaults:
    input: TARGET_CSVS

rule simple_segment_fastas:
    """Split detailed CSV into simple per-segment FASTA files.

    Uses short, IgBLAST-friendly seq IDs and tries to collapse exons for
    constant region sequences into a single CDS sequence per allele.
    """
    output: expand("output/{{organism}}.{segment}.fasta", segment=SEGMENTS)
    input: "output/{organism}.csv"
    params:
        partial=True
    run:
        handles = {key: open(path, "wt") for key, path in zip(SEGMENTS, output)}
        try:
            outs = {}
            with open(input[0]) as f_in:
                reader = DictReader(f_in)
                # For VDJ there's one row per allele, but each C is split up
                # into exons.
                for row in reader:
                    seg = row["Segment"]
                    if seg not in outs:
                        outs[seg] = {}
                    # makeblastdb refuses to parse names with character like "/"
                    # (such as in the human sequence "IGHV1/OR15-1*01") so we'll
                    # use "_" in place of any unusual characters
                    seqid = re.sub("[^-A-Za-z0-9*.]", "_", row["Allele"])
                    # remove "." or whatever else we encounter other than letters
                    seq = re.sub("[^A-Z]", "", row["Seq"].upper())
                    partial = row["Partial"]
                    if seqid not in outs[seg]:
                        outs[seg][seqid] = []
                    outs[seg][seqid].append((partial, seq))
            for seg, seqset in outs.items():
                handle = handles.get(seg, handles["IGHC"])
                for seqid, pairs in seqset.items():
                    partials = [p[0] for p in pairs]
                    # skip seqs with partial interior regions.  5' at the start
                    # is potentially OK and 3' at the end, but that's it.
                    p_5p = ["in 5'" in txt for txt in partials]
                    p_3p = ["in 3'" in txt for txt in partials]
                    if len(pairs) > 1 and (p_3p[0] or p_5p[-1] or any(left or right for left, right in zip(p_5p[1:-1], p_3p[1:-1]))):
                        # definitely not doing this
                        sys.stderr.write(f"Skipping writing sequence \"{seqid}\" with partial interior regions to FASTA\n")
                        continue
                    if any(left or right for left, right in zip(p_5p, p_3p)) and not params.partial:
                        # this case can go either way
                        sys.stderr.write(f"Skipping writing incomplete sequence {seqid} to FASTA\n")
                        continue
                    partial_map = {
                        "partial in 5'partial in 3'": "partial in 5' and in 3'",
                        "partial in 5' and in 3'": "partial in 5' and in 3'",
                        "partial in 5'": "partial in 5'",
                        "partial in 3'": "partial in 3'",
                        "": ""}
                    partial_key = "".join([p for p in partials if p])
                    desc = partial_map[partial_key]
                    # NOTE this implicitly assumes that exons are sorted
                    # correctly in the input CSV
                    seq = "".join([p[1] for p in pairs])
                    SeqIO.write(
                        SeqRecord(Seq(seq), id=seqid, description=desc),
                        handle,
                        "fasta-2line")
        finally:
            for val in handles.values():
                val.close()

rule tabulate_segments:
    output: "output/{organism}.csv"
    input: ancient("from-imgt/{organism}.fasta")
    run:
        with open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=FIELDS + FIELDS_EXTRA, lineterminator="\n")
            writer.writeheader()
            with open(input[0]) as f_in:
                for record in SeqIO.parse(f_in, "fasta"):
                    row = parse_fields(record.description)
                    row.update(parse_allele(row["Allele"]))
                    row["SeqDesc"] = record.description
                    row["Seq"] = str(record.seq)
                    writer.writerow(row)

# Snakemake has its HTTPRemoteProvider object for this instead of curl, but it
# has some side effects so I'm just sticking with curl.
# (Snakemake does a HEAD request on HTTP inputs to try to establish their mtime
# even when that input isn't needed, and imgt.org takes the same (kind of long)
# amount of time to respond to either HEAD or GET, so something like `snakemake
# -n` can take a surprisingly long while even when nothing needs doing.)

rule gather_fasta:
    output: "from-imgt/{organism}.fasta"
    input: expand("from-imgt/{{organism}}/{segment}.fasta", segment = SEGMENTS)
    shell: "cat {input} > {output}"

ruleorder: get_constant_fasta > get_fasta

rule get_fasta:
    output: "from-imgt/{organism}/{segment}.fasta"
    params:
        url=lambda w: f"{URL_VDJ}/{w.organism}/IG/{w.segment}.fasta"
    shell: "curl -L '{params.url}' > '{output}'"

# I know, there has to be a better way than scraping HTML, right?  I don't see
# one though.
rule get_constant_fasta:
    output: "from-imgt/{organism}/IGHC.fasta"
    input: "from-imgt/{organism}/IGHC.html"
    run:
        with open(input[0]) as f_in, open(output[0], "wt") as f_out:
            html = BeautifulSoup(f_in, "html.parser")
            fasta = html.find_all("pre")[1].get_text()
            # remove extra newlines on either end, if present
            if fasta.startswith("\n"):
                fasta = fasta[1:]
            if fasta.endswith("\n\n"):
                fasta = fasta[:-1]
            f_out.write(fasta)

rule get_constant_fasta_html:
    output: ensure("from-imgt/{organism}/IGHC.html", non_empty=True)
    params:
        url=lambda w: f"{URL_C}?query=7.2+IGHC&species={w.organism}"
    shell: "curl -L '{params.url}' > '{output}'"
