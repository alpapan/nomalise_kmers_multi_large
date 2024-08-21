## Normalise reads using a kmer function

Takes a few Tb of paired end Illumina data and then removes redundancy. 
Pretty useful for RNA-Seq and a drop-in replacement for Trinity's normalisation (use `Trinity --no_normalize_reads`).

Example output:

```bash
Processing rate: 21,497 (+10.85%) sequences per second, processed 245,869,236 pairs, printed: 33,999,830 (+0.06%), skipped: 211,869,406 (+0.10%), Total unique kmers across all sequences: 304,792,232 (+0.04%)
File's statistics: Processed 2,458,692,374, Printed 349,654,749, Skipped 2,109,037,625
```

### Compile:

```bash
gcc -o normalise_kmers_multi_large normalise_kmers_multi_large.c -lpthreads
```

I will write this at some point, in the meantime here is the helpful help:

### Usage

```bash


                Mandatory:
                * --forward|-f file1 [file2+]   List of forward (read1) sequence files
                * --reverse|-r file1 [file2+]   List of reverse (read2) sequence files

                Optional:
                [--ksize|-k (integer 5-32; def. 25)]    Number of what size of K to use (must be between 5 and 32)
                [--depth|-d (integer; def. 100)]        Number determining when a kmer is tagged as high coverage (defaults to 100),
                                                        must be above 2xCPU count as each CPU calculates depth independently
                [--coverage|-g (float 0-1; def. 0.9)]   Proportion (0-1) of sequence that must be covered by high coverage kmers before tagging as redundant
                [--canonical|-c]                        Flag to ask the program to merge kmers from forward and reverse complement forms (e.g. for DNA-Seq or unstranded RNA-Seq)
                [--filetype|-t (fq|fa; def. fq)]        Whether the input files are fastq or fasta
                [--outformat|-o (fq|fq; def. fq)]       Whether you want the output files as fastq or fasta (e.g. for Trinity)
                [--memory_start|-m (integer; def. 1)]   Number in Gb of the total memory the program will initially allocate across all threads. The program may request more memory when needed but very small values will cause it to slow down
                [--cpu|-p (int; def 1)]                 Number of CPUs that will process the input files, each file is processed sequentially after distributing to the CPUs
                [--verbose|-e]                          Entertain the user
                [--debug|-b]                            Annoy the developer
                [--version|-v]                          Print version and exit

```

The --coverage controls how much redundancy you want to remove. It defaults to 0.9 and I've used --coverage 0.96. One day I'll survey the effect...

I would recommend it is used on at least 200 Gb worth of input FASTQ data, otherwise you might have to use --cpu 1 to remove enough redundancy.

### Performance and Warnings

- It is fast (but could be faster), it processes 1.4 Tb of FASTQ in 7h using 10 CPUs, scaling pretty linearly.
- Memory is dictated by the data. Each unique kmer takes 16 bytes of memory per CPU so that 350 million kmers would need ~ 5Gb of RAM per CPU.
- There is no shared kmer table, each thread is independent (hence the speed).
- Light on the I/O: Input files are memory mapped.
- Output files can be made Trinity friendly with --outformat fa, no more wasted I/O!

```bash
--- Final Report ---
Processed Records: 2,987,923,777
Printed Records: 352,574,553
Skipped Records: 2,635,349,224
Total runtime: 24569.00 seconds
Overall processing rate: 121,614 sequences per second
136874.68user 1179.63system 6:49:28elapsed 561%CPU (0avgtext+0avgdata 142760856maxresident)k
2713964272inputs+339317512outputs (170major+100520914minor)pagefaults 0swaps

Input size: 1.4T
Output size: 162G
```

NB to developers: Do not use normalise_kmers_multi. 
It is an older version that uses a shared kmer table.
It also needs some bugfixes to be backported to it.
It is far less performant due to the need to sync kmer tables and the gain is minimal for the large datasets you'd need normalisation for.
