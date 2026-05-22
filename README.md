# Soybean Primer-BLAST

Local primer design and specificity checking tool for _Glycine max_ (soybean), replicating NCBI Primer-BLAST functionality.

## Features

- **Primer design** via Primer3 with customizable parameters (Tm, GC%, product size)
- **Genome-wide specificity checking** via local BLAST against the Wm82.a4.v1 reference
- **GFF3-aware gene targeting** — design primers against specific genes, CDS regions, or genomic intervals
- **Web interface** (Flask) replicating the NCBI Primer-BLAST interactive experience
- **CLI** for batch/scripted use with text, JSON, and TSV output formats
- **Primer quality scoring** — GC balance, 3′ stability, self-complementarity, penalty evaluation
- **Automatic fallback** — retries with relaxed constraints when no primers pass default filters

## Requirements

- Python 3.8+
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (`blastn`, `makeblastdb` on PATH)
- Primer3 (`pip install primer3-py`)

## Installation

```bash
git clone https://github.com/LeoYYXX/soybean-primer-blast.git
cd soybean-primer-blast
pip install -r requirements.txt
```

## Data Setup

Download the soybean reference genome and annotation:

| File | Size | Source |
|------|------|--------|
| `Gmax_508_v4.0.fa` | ~945 MB | [Phytozome / SoyBase](https://phytozome-next.jgi.doe.gov/info/Gmax_Wm82_a4_v1) |
| `Gmax_508_Wm82.a4.v1.gene.gff3` | ~119 MB | [Phytozome / SoyBase](https://phytozome-next.jgi.doe.gov/info/Gmax_Wm82_a4_v1) |

Place both files in the project root. The GFF3 index is auto-built on first run.

Configure BLAST binary paths in `primer_blast/constants.py` if they are not on your system PATH.

## Quick Start

### CLI

```bash
# Design primers for a gene by name
python -m primer_blast --gene Glyma.13G357600

# Design primers for a genomic region
python -m primer_blast --region Chr13:3000000-3005000

# JSON output with full specificity results
python -m primer_blast --gene Glyma.13G357600 --format json

# Custom parameters
python -m primer_blast --gene Glyma.13G357600 --product-size 100-300 --num-return 3
```

### Web Server

```bash
python -m primer_blast --serve
# Open http://localhost:5000
```

### Python API

```python
from primer_blast.fasta_index import load_fasta_index
from primer_blast.gff_index import load_or_build_index
from primer_blast.sequence import extract_gene_sequence
from primer_blast.primer_design import design_primers

fai = load_fasta_index("Gmax_508_v4.0.fa")
gene_index, transcript_index = load_or_build_index("Gmax_508_Wm82.a4.v1.gene.gff3")
transcript = transcript_index["Glyma.13G357600"]
template = extract_gene_sequence("Gmax_508_v4.0.fa", transcript, fai)
pairs = design_primers(template, product_size_min=100, product_size_max=400)
```

## Project Structure

```
primer_blast/
    __main__.py        Entry point (orchestrates CLI + server)
    cli.py             Argument parsing
    constants.py       Default Primer3 and BLAST parameters
    fasta_index.py     FASTA index (.fai) loader
    gff_index.py       GFF3 annotation index
    sequence.py        Template extraction (gene, CDS, region)
    primer_design.py   Primer3 wrapper
    primer_score.py    Quality scoring and penalty evaluation
    specificity.py     BLAST-based off-target detection
    blast_db.py        Local BLAST database management
    output.py          Text/JSON/TSV formatting
    server.py          Flask web server
    templates/
        index.html     Web UI
batch_test_100.py      Batch validation script
generate_final_report.py  Test report generator
```

## Default Parameters

Primer3 defaults tuned for soybean genome specificity:

| Parameter | Default |
|-----------|---------|
| Primer size | 20–27 bp (opt 22) |
| Tm | 57–63 °C (opt 60) |
| GC content | 35–65% |
| Product size | 100–400 bp |
| Max self-complementarity penalty | 9.0 |

## Validation

v1.0 tested on 1,000 soybean genes:
- **99%** success rate for primer design
- **100%** specificity pass rate (coverage-filtered BLAST)
- Average 4.8 primer pairs per gene

## License

MIT
