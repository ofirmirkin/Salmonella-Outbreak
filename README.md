# Salmonella Outbreak

Python tool for finding mutations in a genome sequence using k-mer analysis and BLAST search.
BLAST hits are saved to XML files for further analysis.

## Usage
```python3 main.py <wild_filename> <resistant_filename> [<BLAST (0 or 1)>]```
- `<wild_filename>`: Input FASTA sequencing file containing the wild-type genome.
- `<resistant_filename>`: Input FASTA sequencing file containing the resistant genome.
- `[<BLAST (0 or 1)>]`: Optional flag to run BLAST (1 to enable, 0 to disable). If disabled, the program will used the pre-generated BLAST XML files.
