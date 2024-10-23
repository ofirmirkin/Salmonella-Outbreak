# Salmonella Outbreak

Python tool for finding mutations in a genome sequence using k-mer analysis and BLAST search.
BLAST hits are saved to XML files for further analysis.

## Usage
```python3 main.py <wild_filename> <resistant_filename> [<BLAST (0 or 1)>]```
- `<wild_filename>`: Input file containing the wild-type genome sequence.
- `<resistant_filename>`: Input file containing the resistant genome sequence.
- `[<BLAST (0 or 1)>]`: Optional flag to run BLAST (1 to enable, 0 to disable). If disabled, the program will used the pre-generated BLAST XML files.