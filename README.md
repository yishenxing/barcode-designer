# BarcodeDesigner

**BarcodeDesigner** is a Python script for designing robust and diverse barcode pools used in high-throughput sequencing (NGS). It generates barcode sequences that meet constraints on GC content, homopolymers, repeats, secondary structures, and edit distances, ensuring compatibility and reliability for multiplexed sequencing.

## Features

- Multiprocessing-enabled parallel barcode generation
- Configurable GC content constraints
- Homopolymer filtering
- Custom repeat unit filtering
- Secondary structure (hairpin) avoidance
- Substring overlap exclusion (to prevent barcode cross-talk)
- Minimum edit distance filtering using Levenshtein distance

---

## Requirements

- Python 3.6+
- [python-Levenshtein](https://pypi.org/project/python-Levenshtein/)

Install via pip:

```bash
pip install python-Levenshtein
pip install biopython
```

---

## Usage

```bash
python barcode_designer.py \
    --barcodeLen 12 \
    --polymerLen 4 \
    --minGCRate 0.3 \
    --maxGCRate 0.7 \
    --dictRepeats "GC:3,CG:3,AT:3,TA:3" \
    --innerPairLen 4 \
    --outerPairLen 6 \
    --minEditDistance 3 \
    --numCandidateBarcode 1000 \
    --numGeneratedBarcode 96 \
    --numProcess 8 \
    --outputSeqFile output/barcodes.fa \
    --outputDistanceFile output/distances.txt \
    --logFile output/run.log
```

---

### Arguments

| Argument                | Type  | Description                                                                    |
| ----------------------- | ----- | ------------------------------------------------------------------------------ |
| `--barcodeLen`          | int   | Length of each barcode                                                         |
| `--polymerLen`          | int   | Maximum allowed homopolymer run (e.g., AAAA)                                   |
| `--minGCRate`           | float | Minimum allowed GC content (e.g., 0.3)                                         |
| `--maxGCRate`           | float | Maximum allowed GC content (e.g., 0.7)                                         |
| `--dictRepeats`         | str   | Comma-separated repeat units and counts to exclude (e.g., `"GC:3,AT:3"`)       |
| `--innerPairLen`        | int   | Minimum length of internal reverse-complementary sequences to avoid (hairpins) |
| `--outerPairLen`        | int   | Length of subsequences that must not overlap between barcodes                  |
| `--minEditDistance`     | int   | Minimum Levenshtein edit distance between any pair of barcodes                 |
| `--numCandidateBarcode` | int   | Total number of candidate barcodes to initially generate                       |
| `--numGeneratedBarcode` | int   | Final number of filtered barcodes to output                                    |
| `--numProcess`          | int   | Number of parallel processes to use                                            |
| `--outputSeqFile`       | str   | Output FASTA file of final barcode sequences                                   |
| `--outputDistanceFile`  | str   | Output file with pairwise edit distance information                            |
| `--logFile`             | str   | Log file to record progress and filter statistics                              |

---

## Output

| File             | Description                                                          |
| ---------------- | -------------------------------------------------------------------- |
| `*.fa`           | FASTA file containing the final set of selected barcodes             |
| `*_candidate.fa` | FASTA file of all candidate barcodes before filtering                |
| `*.txt`          | Tab-separated file of pairwise edit distances between final barcodes |
| `*.log`          | Log file with filtering summary and rejection counts by error type   |

---

## Example

Generate a set of 96 barcodes, each 12 bp long, with:

- GC content between 40% and 60%
- No homopolymers longer than 2 bp
- No specific 3-mer repeats like "GC" or "TA"
- No self-complementary sequences of length >= 4 (to avoid hairpins)
- No 8-bp overlapping subsequences across barcodes
- A minimum edit distance of 9 between any two barcodes
- 500 candidate barcodes initially generated
- at least 96 final barcodes after filtering
- Using 16 parallel processes
```bash
python barcode_designer.py \
    --barcodeLen 20 \
    --polymerLen 3 \
    --minGCRate 0.4 \
    --maxGCRate 0.6 \
    --dictRepeats "GC:3,CG:3,AT:3,TA:3" \
    --innerPairLen 4 \
    --outerPairLen 8 \
    --minEditDistance 9 \
    --numCandidateBarcode 500 \
    --numGeneratedBarcode 96 \
    --numProcess 16 \
    --outputSeqFile barcodes.fa \
    --outputDistanceFile distances.txt \
    --logFile run.log
```

---

## License

This project is licensed under the MIT License.
