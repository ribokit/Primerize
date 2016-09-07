# NA_Thermo

**NA_Thermo** is the Das lab in-house primer design tool. It is based on *MATLAB*.


## Installation

To install legacy MATLAB scripts for **Primerize** (previously called **NA_Thermo**), simply:

- Download the zip or tar file of the repository and unpack; or `git clone https://github.com/DasLab/Primerize.git`.

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/Primerize/MATLAB/Scripts/`.

## Usage

### design_primers
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

This pipeline has been moved into *Python* for general use and for deployment in the **Primerize** server.

### mutate_primers
Suppose you have an assembly of primers to make an RNA, and then want to make a mutant. Many of the primers that you already have are already in hand, so you just need to get a few primers that define the mutation location(s). Tou can use `mutate_primers`

- Check out use of `mutate_primers` in the `Examples` directory:
   `P4P6_afewmutants_order_script.m`, `P4P6_de209_R1_Jaeger_mutant_script.m`, and `FN_KTtest_muts_script.m`

- Example to install a GGAA/R(1) tetraloop receptor in the P4-P6 RNA, replacing the GAAA/11-nt tetraloop receptor:
```
primers= ...
    {'TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCG','GGCCATCTCAAAGTTTCCCCTGAGACTTGGTACTGAACGGCTGTTGACCCCTTTCCCG','GGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGAC','GTTGACTTAGGACTTGGCTGCGTGTTAGGACCATGTCCGTCAGCTTATTACCATAC','CAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTT','GTTGTTGTTGTTGTTTCTTTGGTTTGGTTTTGAACTGCATCCATATCAACAGAAG'};

mutate_primers( primers, 'R1align.txt' )
```

- In the above, the primers were a set of our 'favorite' primers to make P4-P6 DNA template with a T7 promoter. Here, `R1align.txt` is a sequence alignment file, with dashes representing gaps:

```
TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACACGCAGCCAAGT-CCTAAGTCAACAGATCTTCTGTTGATATGG--ATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC
TTCTAATACGACTCACTATAGGCCAAAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGGAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACACGCAGCCAAGT-CCTGTGTCAACAGATCTTCTGTTGAATCTGG-ATGCAGTTCAAAACCAAACCAAAGAAACAACAACAACAAC
```

- In this case 5 of the 6 primers need to be replaced, since mutations are present in both the tetraloop & receptor (much of the sequence is affected!).

## Contact

Send questions and comments to [_Rhiju Das_](mailto:rhiju@stanford.edu).

