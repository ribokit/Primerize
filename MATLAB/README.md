# NA_Thermo

**NA_Thermo** is the Das lab in-house primer design tool. It is based on *MATLAB*.


## Installation

To install **NA_Thermo**, simply:

- Download the zip or tar file of the repository and unpack; or `git clone https://github.com/DasLab/Primerize.git`.

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/Primerize/MATLAB/Scripts/`.

## Usage

To design primers for your sequence, just follow these easy steps:

- Define your sequence. For example:

```
sequence = 'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC';
tag = 'P4P6';
```

This sequences includes a 20-nucleotide _T7 promoter_ sequence at the beginning, and then a construct (starting with `GG...`) encoding the _P4-P6 domain_ of the Tetrahymena ribozyme along with flanking sequences.

- Run with:

```
primers = design_primers(sequence, tag);
```

This will compute primers with minimal length, annealing temperatures above a cutoff (**60 &deg;C**, by default), through a recursive strategy. An additional score term helps avoid primers that share multiple 3' nucleotides with other parts of the sequence, as a heuristic to reduce mispriming. The algorithm has similarities (but was developed independently) of Thachuk & Condon (2007), BIBE, Proc. of 7th IEEE International Conf., p. 123-130.

If you want to use our original script (in use from 2008-2011), use `design_primers_OLD`. (This was slower and did not rigorously optimize length.)

- We usually copy/paste these to a Word or Excel document for easy look-up later. If you add primer labels, you can copy/paste these to the IDT website or wherever.

## Contact

Send questions and comments to [_Rhiju Das_](mailto:rhiju@stanford.edu).

