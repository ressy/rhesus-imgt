# Rhesus Macaque V(D)J Gathering from IMGT

This Snakefile downloads IMGT's per-segment FASTA files for [rhesus] and
[human] antibody genes and combines them into
[one CSV for Macaca mulatta](output/Macaca_mulatta.csv) and
[one CSV for Homo sapiens](output/Homo_sapiens.csv), parsing out details from
the FASTA description lines (as described in the [IMGT FASTA format]).  It also
creates per-segment FASTA files for each species using ungapped sequences and
modified sequence IDs for use with IgBLAST.

[IMGT FASTA format]: https://www.imgt.org/IMGTindex/Fasta.php
[rhesus]: http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Macaca_mulatta/IG
[human]: http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG
