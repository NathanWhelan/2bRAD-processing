# `paris_optimize.py`: choosing Stacks de novo parameters following Paris et al. (2017)

When you assemble RAD-seq data without a reference genome, Stacks has to guess which reads come from the same place in the genome. Three parameters (`m`, `M`, `n`) control that guess, and the right values differ between species and datasets. Paris et al. (2017) showed how to choose them:

1. Assemble the data several times with different values.
2. Measure the same summary numbers for each assembly.
3. Choose the value where those numbers say alleles have been joined, without going further than needed.

This script does steps 2 and 3, using **only Paris et al.'s metrics and reasoning**. It reads Stacks 2 output folders and writes three files:
- `paris_report.txt`: the recommendation and the main evidence;
- `paris_plots.pdf`: one page of plots for each series of runs;
- `paris_summary.xlsx`: the supporting tables (per run, per tested value, per individual, SNPs per locus, fixed SNPs, high-heterozygosity loci, coverage per sample).

With `--image png` or `--image jpeg`, it also saves each plot page as an image.

It works for any Stacks 2 de novo dataset, single-end or paired-end. Throughout this guide, a real dataset (Alabama pearlshell, *Margaritifera marrianae*, 44 samples) is the example. Section 6 walks through its results. Section 4.10 lists everything the script does that Paris et al. did not.

---

## Contents

1. [Background you need](#1-background-you-need)
2. [Running the script](#2-running-the-script)
3. [How Paris et al. chose parameters](#3-how-paris-et-al-chose-parameters)
4. [What the script measures and how it decides](#4-what-the-script-measures-and-how-it-decides)
5. [Reading the outputs](#5-reading-the-outputs)
6. [Worked example: Alabama pearlshell](#6-worked-example-alabama-pearlshell)
7. [Troubleshooting](#7-troubleshooting)
8. [Glossary](#8-glossary)
9. [References](#9-references)

---

## 1. Background you need

### 1.1 What de novo assembly is doing

RAD-seq sequences short stretches of DNA next to restriction-enzyme cut sites. Every individual is cut at (nearly) the same sites, so you get reads from the same thousands of places in everyone's genome. Each such place is a **locus**.

With a reference genome, you would map reads to it to find out which locus each read belongs to. Without one ("de novo"), Stacks groups reads by how similar their sequences are. `denovo_map.pl` runs the Stacks programs in this order:

| Step | Program | What it does | Parameter |
|---|---|---|---|
| 1 | `ustacks` | Within one individual, piles up identical reads into **stacks**. Each stack is a candidate allele. | `m`: minimum reads needed to make a stack |
| 2 | `ustacks` | Within one individual, merges stacks that differ at ≤ `M` positions into one locus, treating them as alleles of the same locus. | `M` |
| 3 | `cstacks` | Across individuals, merges loci that differ at ≤ `n` positions into one **catalog** locus. | `n` |
| 4 | `sstacks`, `tsv2bam`, `gstacks` | Matches each individual back to the catalog, assembles paired-end reads (if any), then calls SNPs and genotypes. | — |
| 5 | `populations` | Filters loci (for example, by how many samples have them) and computes statistics. | `-R` |

Steps 1 and 2 happen **within each individual**, before any other individual is looked at. That is why Paris et al. could watch `M` at work within individuals (Section 3), and why those numbers don't change with `n`.

### 1.2 The two ways merging goes wrong

Paris et al. put it this way: "set it too low and alleles from the same locus will not collapse, set it too high and paralogous or repetitive loci will incorrectly merge together."

**Under-merging (`M` or `n` too low).** The two alleles of one locus differ at more positions than `M` allows, so Stacks keeps them as two separate "loci". A heterozygous individual then looks homozygous at each piece. You lose real SNPs, and the genotypes at that locus are wrong.

**Over-merging (`M` or `n` too high).** Two *different* loci with similar sequence get merged into one. These are usually **paralogs** (duplicated sequences) or repeats. The fixed differences between the two copies look like SNPs, and nearly every individual looks heterozygous at them, because everyone carries both copies.

```
Under-merging (M too low)              Over-merging (M too high)
allele A:  ACGTACGTAC                  copy 1 (one chromosome):     ACGTACGTAC
allele B:  ACTTACGAAC  (2 differences) copy 2 (another chromosome): ACTTACGAAC
M = 1  ->  two "loci", no SNPs         M = 2  ->  one "locus" with 2 fake SNPs
```

The two pairs of sequences above are identical. **From sequence alone, Stacks cannot tell whether two similar sequences are alleles or paralogs.** That is why no single value of `M` is "correct", and why Paris et al. looked at several metrics together.

### 1.3 The parameter `m`

`m` sets how many identical reads it takes to believe a sequence is real rather than a sequencing error. In all three of Paris et al.'s datasets, the number of loci dropped sharply from m2 to m3 and then stabilised, and they chose `m = 3` (the Stacks default) for all three. They suggest testing m3 to m7 (their Table 2), do not suggest m above 5, and advise disabling secondary reads above m6. A common design, and the one in the example, is to fix `m = 3` and vary `M` and `n`. The script also analyses a series in which `m` varies, if you run one.

### 1.4 r80 loci, and `-R` versus `-r`

Paris et al. judge each assembly by its **r80 loci**: loci found in at least 80% of individuals. Their reasoning is that a locus present in nearly everyone is unlikely to come from paralogs, repeats or sequencing error, so r80 loci are a good stand-in for the true genome. You are **not** required to use only r80 loci in later analyses. They are just the yardstick for comparing assemblies.

In Stacks 2, two `populations` options sound alike but do different things:

| Option | Meaning | Use for this script? |
|---|---|---|
| `-R 0.8` (`--min-samples-overall`) | Keep a locus if it is in ≥ 80% of **all** samples | **Yes.** This is Paris's r80 |
| `-r 0.8` (`--min-samples-per-pop`) | Keep a locus if it is in ≥ 80% of samples **within a population** | No |

With several populations, `-r` can keep a locus that is present in just one small population. `denovo_map.pl -r 0.8` runs `populations` with `-r`, so the `populations` output inside that `denovo_map.pl` folder is **not** r80. There are two ways to get r80 output (Section 2.1):
- run `populations -R 0.8` separately on each assembly;
- give `denovo_map.pl` the option `-X "populations:-R 0.8"`, **without** `-r`. Its `populations` output in the assembly folder is then r80, and you can pass that folder directly.

The script refuses output that was not made with `-R`. If `-r` was applied as well, it warns that both filters shaped the loci.

If you use another `-R` value, the script still runs and labels the loci accordingly (for example, "r70" for `-R 0.7`), but Paris's rule is defined for 0.8.

### 1.5 Paired-end loci in Stacks 2

With paired-end data, `gstacks` assembles the paired-end reads of each locus into a **paired-end contig** and joins it to the forward read with a run of `N`s:

```
forward read (assembled by ustacks/cstacks under M and n)  NNNNNNNNNN  paired-end contig (assembled by gstacks)
|<----------------- e.g. 140 bp ------------------>|                  |<---------- variable --------->|
```

**Only the forward read is assembled under `M` and `n`.** The paired-end contig is built afterwards from whichever reads belong to that locus. Its SNPs still count towards the locus's SNP total, but its length and quality depend on things `M` and `n` don't control. Paris et al. used Stacks 1.42, whose de novo pipeline built loci from the forward reads alone, so their SNP counts covered only the read assembled under `M` and `n`. Where it matters (SNP-dense loci, Section 4.7), the script counts the forward read separately. It finds the forward read from the `N` spacer in `catalog.fa.gz`. Loci whose two parts overlap have no spacer, and for them it uses the usual forward-read length.

---

## 2. Running the script

### 2.1 Make the inputs

You need one de novo assembly per parameter set, and then `populations -R 0.8` on each. For example, to fix `m = 3` and scan `M = n` from 1 to 6:

```bash
# one assembly per parameter set
for i in 1 2 3 4 5 6; do
  mkdir M${i}
  denovo_map.pl --samples ./reads --popmap ./popmap.txt -m 3 -M $i -n $i \
      -T 16 -o ./M${i} --paired
done

# Paris's r80 filter: loci in at least 80% of ALL samples
for i in 1 2 3 4 5 6; do
  mkdir M${i}_R80
  populations -P ./M${i} -O ./M${i}_R80 -M ./popmap.txt -R 0.8 -t 16
done
```

Setting `M = n` and scanning `M` is a shortcut from the Stacks 2 protocol (Rochette & Catchen 2017). Paris et al. varied `M` with `n` held fixed. Their **n rule** then compares `n = M − 1`, `M` and `M + 1` at the chosen `M` (Section 4.5). The report lists the runs that step still needs. For example, if `M = 2` is chosen:

```bash
for n in 1 3; do
  mkdir M2n${n}
  denovo_map.pl --samples ./reads --popmap ./popmap.txt -m 3 -M 2 -n ${n} \
      -T 16 -o ./M2n${n} --paired
  mkdir M2n${n}_R80
  populations -P ./M2n${n} -O ./M2n${n}_R80 -M ./popmap.txt -R 0.8 -t 16
done
```

Then run the script on all the r80 folders together. Leave out `--paired` for single-end data. Use the same `populations` options for every assembly, apart from the folder names.

**Alternative: let `denovo_map.pl` run the r80 filter.** Instead of the separate `populations` step, add `-X "populations:-R 0.8"` to `denovo_map.pl` and leave out `-r`:

```bash
denovo_map.pl --samples ./reads --popmap ./popmap.txt -m 3 -M $i -n $i \
    -T 16 -o ./M${i} --paired -X "populations:-R 0.8"
```

The r80 output is then in `M${i}/` itself, and you pass `M1 M2 …` to the script instead of `M1_R80 M2_R80 …`. Its `populations` options (`-R 0.8 -M ./popmap.txt`) match those of the separate step above, so the two kinds of folder can be analysed together. Don't pass both an assembly folder and a separate r80 folder made from the same assembly. The script warns that they have the same parameters and uses only one.

In the example dataset, the assemblies are `1_mmar/` … `6_mmar/` (made with `-r 0.8`, so their own `populations` output is not used), and the r80 outputs are `1_mmar_R80/` … `6_mmar_R80/`.

### 2.2 Run it

```bash
python paris_optimize.py M1_R80 M2_R80 M3_R80 M4_R80 M5_R80 M6_R80 -o paris_results
```

- Give it the `populations -R 0.8` **output** folders, in any order and as many as you like. If `denovo_map.pl` ran the r80 filter itself (Section 2.1), these are the assembly folders.
- It finds each assembly folder from the `-P` path recorded in `populations.log`. It looks relative to where `populations` was run, then next to the output folder, then in the current folder.
- It reads `m`, `M` and `n` from `denovo_map.log`, not from folder names, so names don't matter.
- The screen shows only the recommendation: the `m M n` pick, one line of reasoning for each parameter, any runs Paris's n rule still needs, and any warnings. The full report is always saved as `paris_report.txt`. Add `--verbose` to print it to the screen as well.
- Six runs take about 40 seconds with 8 CPUs. Most of that is reading each sample's ustacks SNP file for the fixed-SNP count (Section 4.9), which is spread over up to 8 worker processes (`-t`).
- **Python:** 3.6 or later. Parsing needs only the standard library. The spreadsheet needs `openpyxl` and the plots need `matplotlib` (`conda install openpyxl matplotlib`; JPEG images also need Pillow, which comes with matplotlib). Without them, the script still writes the report, and says what it skipped.

For the example dataset:

```bash
python paris_optimize.py {1..6}_mmar_R80 -o paris_results
```

### 2.3 Which files it reads

| File | Used for | If missing |
|---|---|---|
| `populations.log`, `populations.log.distribs` | Filters, r80 loci, polymorphic r80 loci, SNPs per locus | **Required** |
| `populations.haplotypes.tsv` | The samples `populations` genotyped, so the per-sample numbers use the same samples; a check of the r80 locus count | Samples are taken from `denovo_map.log` |
| `populations.sumstats.tsv` | SNP positions; heterozygosity per site, pooled over populations | No forward-read SNP counts or heterozygous-site check |
| `denovo_map.log` | `m`, `M`, `n`; each sample's ustacks stacks, loci and coverage ("Coverage after assembling stacks") | Give `m`, `M`, `n` in a runs table (Section 2.4); no loci per individual, allele collapse or coverage |
| `SAMPLE.alleles.tsv.gz` (in the assembly folder) | Heterozygous loci per individual (ustacks), which decides `M` together with the r80 rule | No heterozygous loci per individual, so `M` falls back to the r80 rule alone (Section 4.4) |
| `SAMPLE.snps.tsv.gz`, `SAMPLE.matches.tsv.gz` (in the assembly folder) | Each individual's heterozygous SNPs (ustacks) and how sstacks lined its loci up with the catalog | No heterozygous SNPs per individual or fixed-SNP count |
| `catalog.snps.tsv.gz` | The catalog's SNPs (cstacks) | No fixed-SNP count |
| `catalog.fa.gz` | Where each locus's forward read ends | SNP-dense loci are counted over whole loci only |

All of these are written by a normal `denovo_map.pl` and `populations` run. **Every population-level number comes from the `populations -R 0.8` folders you pass.** From the assembly folder, the script reads only files written before `populations` runs: the ustacks and cstacks lines of `denovo_map.log` (it stops reading the log where the populations step begins), the per-sample ustacks and sstacks files (`SAMPLE.alleles.tsv.gz`, `SAMPLE.snps.tsv.gz`, `SAMPLE.matches.tsv.gz`) and the catalog (`catalog.snps.tsv.gz`, `catalog.fa.gz`). These are the same whether `denovo_map.pl` was given `-r` or `-R`, because those options only affect its final `populations` step. The `populations` files that `denovo_map.pl` writes into an assembly folder are read only if you pass that folder as a `populations -R` folder, and the script checks `populations.log` to make sure they are `-R` output. In the example they were made with `-r`, so they are never read. A missing optional file switches off only what needs it, and the report's Warnings say what was lost. Without the `SAMPLE.alleles.tsv.gz` files there is no within-individual level-off point, so the recommendation for `M` falls back to the r80 rule alone (or to Fig. 4, if polymorphic r80 loci are still rising at the highest value), and says so.

### 2.4 A runs table, for unusual layouts

If you ran the Stacks programs by hand, so there is no `denovo_map.log`, or the script can't find the assembly folder, describe the runs in a tab-separated table:

```
populations_dir	denovo_dir	m	M	n
M1_R80	asm/M1	3	1	1
M2_R80	asm/M2	3	2	2
```

and run `python paris_optimize.py --runs runs.tsv`. Only `populations_dir` is required. Paths are relative to the table's location. Lines starting with `#` are ignored. Values in the table override those in `denovo_map.log`, with a warning if the two disagree. Without `denovo_map.log` there are no loci per individual or allele-collapse numbers, but heterozygous loci per individual still come from the `SAMPLE.alleles.tsv.gz` files if they are in `denovo_dir`.

### 2.5 Options

| Option | Default | What it changes |
|---|---|---|
| `-o DIR` | `paris_results` | Where the outputs go |
| `--runs FILE` | — | Read runs from a table (Section 2.4) |
| `-v`, `--verbose` | off | Print the full report to the screen. By default the screen shows only the recommendation, any warnings and where the files went. `paris_report.txt` is always the full report |
| `--image png` or `--image jpeg` | off | Also save each plot page as `paris_plot_N.png` or `paris_plot_N.jpg` |
| `--plateau` | `0.01` | Paris Fig. 4 reading: a step that adds fewer new polymorphic r80 loci than this share (1%) adds "very few". Used only when polymorphic r80 loci are still rising at the highest value **and** heterozygous loci per individual are missing |
| `--dist-max-snps` | `6` | SNPs-per-locus plot: loci with 1 … 6 SNPs; loci with more are pooled as "7+" |
| `--tail-snps` | `10` | A locus with at least this many SNPs counts as SNP-dense |
| `--fwd-tail-snps` | scaled | Paired-end data: the SNP-dense threshold for the forward read. The default scales `--tail-snps` to the same SNP density: 10 × forward length ÷ median locus length |
| `--het-site` | `0.7` | Paris's paralog sign: a site where at least this share (70%) of all samples are heterozygous ("heterozygous in a large majority of individuals") |
| `--max-snp-col` | `12` | Last SNP-count column in the spreadsheet; loci with more SNPs are pooled into the next column |
| `-t`, `--threads` | CPUs available, up to 8 | Worker processes for reading the per-sample files of the fixed-SNP count. Changes only the speed |

None of these options changes the recommendation in the example. The recommendation for `M` depends on the r80 counts and on the two-standard-error test for heterozygous loci per individual (Section 4.8). It does not use `--tail-snps` or `--het-site` at all. The default thresholds of `--plateau`, `--tail-snps`, `--fwd-tail-snps` and `--het-site` are this script's, not Paris et al.'s (Section 4.10).

---

## 3. How Paris et al. chose parameters

Paris et al. assembled three datasets (brown trout, king penguin, red earthworm) many times. They "plotted ... visually assessed and optimal sets of parameters were chosen, specific to each dataset". They didn't pick on one number. Each choice was the value where several of their metrics agreed:

| | Trout | Penguin | Earthworm |
|---|---|---|---|
| Most polymorphic r80 loci (Fig. 3b) | M5 | M2 | still rising at M11 |
| New polymorphic r80 loci per step (Fig. 4) | peak at M3/M4 | peak at M1/M2 | "very few" added above M3 |
| Polymorphic loci within individuals | "began to flatten out at M5" | "relatively uniform" at M1 (our reading: see below) | "began to plateau at approximately M3" |
| Allele collapse levels out (Table S6a) | between M3 and M4 | between M1 and M2 | between M2 and M3 |
| SNPs as `M` rises (Fig. 3c) | steep rise | much smaller rise | steep rise |
| **Chosen M** | **M5** | **M2** | **M3** |
| **Chosen n** (`n = M − 1`, `M`, `M + 1`) | n4 | n2 | n3 |

Their reasoning, in their words:
- **Trout:** "the average number of polymorphic loci identified within individuals began to flatten out at M5, and M5 also provided the highest amount of polymorphism across 80% of the population."
- **Penguin:** "At M1, the majority of polymorphism and SNPs were already captured in PGN, and the amount of polymorphic loci was relatively uniform; the highest polymorphism detected across 80% of the population was at M2." They also note that "the PGN dataset was relatively monomorphic". The sentence does not say whether "the amount of polymorphic loci" is within individuals. This guide reads it that way, as the penguin counterpart of the trout and earthworm statements; that is an inference, not something Paris et al. state. Either way, the penguin choice is the r80 maximum (M2).
- **SNPs (Fig. 3c):** trout and earthworm showed a steep rise in SNPs with `M`, while "a much smaller increase was detected in PGN". "The steeper increase in SNPs in TRT reflects the higher value of M required to detect the maximum amount of polymorphism across 80% of the population."
- **Earthworm:** "although polymorphic loci within individuals began to plateau at approximately M3, polymorphism at 80% continued to increase with incremental values of M (even up to unrealistically high levels of M11 ...), suggesting that the dataset contained a small number of loci with a high-density of SNPs". "Fig. 4 showed that very few new polymorphic loci were added at M values greater than 3 ... M3 also provides a balance between obtaining true polymorphism and introducing sequencing error."
- **Two maxims:** "(i) the 80% polymorphic (r80) loci rule, and (ii) n = M plus or minus one iteration".
- **n and fixed SNPs:** they counted the catalog SNPs at which no individual is heterozygous ("fixed SNPs"). "All datasets continued to show more fixed SNPs with increasing values of n (Fig. S3) ... This alone, however, does not provide enough information to choose a value for n. We found that the highest polymorphism of r80 loci resulted from n = M, n = M − 1 or n = M + 1 (Table S9) and suggest this as the best method for obtaining the optimal value for n." At the `n` they chose, 15% (trout), 13% (penguin) and 21% (earthworm) of catalog SNPs were fixed.
- **Paralogs:** they "will appear as heterozygous in a large majority of individuals; they can be filtered out by the populations module after the main pipeline runs ... (--max_obs_het option)". For their three datasets, Paris et al. filtered them after assembly. Their Table 2 also marks `M` "↓ if repetitive or polyploid genome", that is, test lower values of `M` for such genomes.
- **SNP-dense loci:** higher `M` adds "a small number of loci that contained many SNPs" at each increment. "With high coverage, these are likely to be true loci ... However, with low coverage these are more likely to be erroneously assembled loci."

Two patterns run through all three choices:
1. The chosen `M` was **never below** the point where polymorphism within individuals levelled off (for penguin, on the reading above).
2. Where the r80 count had a maximum at or above that point (trout, penguin), they took the maximum. Where it kept rising (earthworm), they took the point where polymorphism within individuals levelled off and Fig. 4 showed few new loci.

None of their datasets had the most polymorphic r80 loci at M1, below where within-individual polymorphism levels off. The example dataset does (Section 6).

---

## 4. What the script measures and how it decides

The script groups runs into **series** in which only one parameter changes:
- `M = n` (the design in Section 2.1);
- `M` with `n` fixed (three or more values);
- `n` with `M` fixed;
- `m` with `M` and `n` fixed.

A run can belong to more than one series. Runs with different sample counts or `-R` values are never compared. In this section, *P* is the number of polymorphic r80 loci (r80 loci with at least one SNP), and a "step" is going from one tested value to the next (for example M2n2 → M3n3).

### 4.1 r80 rule

- **What:** the value with the most polymorphic r80 loci (*P*). For an `n` series, only `n = M − 1`, `M` and `M + 1` are allowed (Section 4.5).
- **Paris:** their first maxim, and the deciding metric for trout and penguin.
- **In the report:** "r80 max". Panel (b).

### 4.2 Fig. 4: new polymorphic r80 loci per step

- **What:** the change in *P* at each step. Paris et al. read Fig. 4 for the point after which "very few new polymorphic loci were added".
- **How:** the first value after which no step adds `--plateau` (1%) or more new polymorphic r80 loci. A step that loses loci adds none.
- **Used for the pick** only when *P* is still rising at the highest value tested **and** the within-individual metrics are missing (Section 4.4). Otherwise it is reported alongside.
- **In the report:** "Fig. 4". Panel (d).

### 4.3 Within individuals

These come from each sample's **ustacks** output, before the catalog is built, so they show what `M` does inside each individual and don't change with `n`.

| Metric | Source | Paris's reading |
|---|---|---|
| **Loci per individual** | "Final number of stacks built" per sample in `denovo_map.log` | Falls as `M` merges alleles. At M0, "a high number of loci was assembled ... simply because putative alleles were not correctly merged" |
| **Allele collapse:** share of each individual's stacks merged into loci, 1 − (loci after `M` ÷ stacks) | "Assembled X stacks into Y" per sample in `denovo_map.log` | "the proportion of putative alleles that collapsed levelled out between M3 and M4" (trout; their Table S6a) |
| **Heterozygous loci per individual** | loci with two or more haplotypes in `SAMPLE.alleles.tsv.gz` | "polymorphic loci within individuals ... began to flatten out" (trout), "began to plateau" (earthworm). **Decides `M`** |
| **Heterozygous SNPs per individual** | SNPs of type E in `SAMPLE.snps.tsv.gz` | Count (i) of Paris et al.'s fixed-SNP comparison (Section 4.9). Reported only |

Each is averaged over individuals. A step is tested **paired by individual**: the script takes each individual's change between two values, then the mean change and its standard error. A step **rises** if its mean change is more than two standard errors above zero. Pairing removes the large differences between individuals (coverage, heterozygosity), so a consistent change shows up even when it is small.

**Where polymorphism within individuals first levels off** ("plateau" in the report) is the first value after which the next step of heterozygous loci per individual does **not** rise beyond two standard errors. This is the script's reading of Paris et al.'s "began to flatten out" (trout) and "began to plateau" (earthworm), which describe where the levelling starts. Later steps can rise again, as they do in the example (Section 6.2).

Paris et al. don't say whether their within-individual counts were limited to loci found in 40/60/80% of individuals. The script takes the literal reading: all of each individual's ustacks loci.

Allele collapse rises at every step in practice, and with dozens of individuals even tiny steps are significant. So the script reports how much of the collapse across the tested range has happened by each value, but doesn't take a pick from it.

### 4.4 The decision rule for `M`

For an `M` or `M = n` series:

1. If *P* is still rising at the highest value tested (the earthworm case), the pick is **where heterozygous loci per individual first level off**. If they are missing, the Fig. 4 reading is used instead.
2. Otherwise, the pick is the **higher of the r80 maximum and where heterozygous loci per individual first level off**.
3. If heterozygous loci per individual rise at every step, the tested range never reaches the point at or above which Paris et al. chose. There is no pick; test higher values.
4. If heterozygous loci per individual are missing and *P* is not still rising, only the r80 rule applies.

This reproduces all three of Paris et al.'s choices from their own descriptions:

| Dataset | r80 maximum | Within individuals first levels off | Rule | Paris chose |
|---|---|---|---|---|
| Trout | M5 | M5 | higher of the two: M5 | M5 |
| Penguin | M2 | M1 ("relatively uniform"; our reading, Section 3) | higher of the two: M2 | M2 |
| Earthworm | still rising at M11 | about M3 | still rising, so the level-off point: M3 | M3 |
| Pearlshell (Section 6) | M1 | M2 | higher of the two: M2 | — |

The rule is a written-down version of what Paris et al. did by eye. It is not their wording. For their three datasets, it gives what they chose. In each of them the r80 maximum was at or above the level-off point, or still rising. When the r80 maximum lies **below** the level-off point, as in the example, rule 2 applies their judgement to a case they did not meet. Taken alone, their r80 rule would give the r80 maximum. The report says this whenever it happens. Section 4.10 lists this and the other parts that are not from Paris.

**A level-off at the last tested step rests on that one step.** With no real change, a step still comes out "rising" about 2–3% of the time, so if the pick is the highest value tested, test one value beyond it.

### 4.5 `n`: Paris's n rule

For an `n` series at fixed `m` and `M`, the pick is the value with the most polymorphic r80 loci among `n = M − 1`, `M` and `M + 1`. The report notes any of the three that weren't tested, and any value outside that range with more loci.

The **Recommendation** section at the end of the report puts `m`, `M` and `n` together. When `M` comes from an `M = n` series and there is no `n` series at that `M`, it reports `n = M` "for now" and lists the runs still needed for the n rule.

### 4.6 `m`

For an `m` series, the pick is the r80 maximum, as for Paris et al., who chose m3 because it gave "the highest amount of polymorphism" at r40, r60 and r80. With no `m` series, the report says so, and notes Paris et al.'s advice (Section 1.3).

### 4.7 Over-merging signs Paris names (cautions, not picks)

Paris et al. name two signs of `M` being too high. The script reports both, from the pick upwards, as cautions:

- **SNP-dense loci:** the share of polymorphic r80 loci with ≥ `--tail-snps` (10) SNPs. With paired-end data, only SNPs in the forward read count (Section 1.5), with the threshold scaled to the same density: ≥ 5 in a 140 bp forward read, for a median 282 bp locus. These loci usually rise with `M`. Paris et al. found them "likely to be true loci" at high coverage and "more likely to be erroneously assembled" at low coverage, so the report gives your coverage alongside (Section 5).
- **Loci "heterozygous in a large majority of individuals":** the share of polymorphic r80 loci with a site where at least `--het-site` (70%) of all samples are heterozygous. A real single-copy site with two alleles is expected to have at most 50% heterozygotes (2*pq* at *p* = *q* = 0.5). For their datasets, Paris et al. filtered such loci **after** assembly (`populations --max_obs_het`). Their Table 2 also suggests testing lower `M` for a repetitive or polyploid genome. In Stacks 2, `--max-obs-het` removes the suspicious **SNP**, not the whole locus. To drop whole loci, list the loci that have such a site and pass them to `populations --blacklist`.

Paris et al. give no numeric thresholds for either sign. The 10 SNPs, the forward-read scaling and the 70% are this script's (Section 4.10).

The plots (panel e) and the spreadsheet also show the SNPs-per-locus distribution, where the "small number of loci that contained many SNPs" that Paris et al. saw added with each increment of `M` appears. No pick is drawn from it.

### 4.8 Is a difference real?

Paris et al. compared their plots by eye. The tests below are this script's (Section 4.10).

- **Within-individual steps** are tested paired by individual (Section 4.3): mean change ± 2 SE, and * when it is beyond 2 SE. The report also gives how many individuals went up.
- **Shares of loci** (SNP-dense, heterozygous sites): for a share *p* from *N* loci, the standard error is √(*p*(1 − *p*)/*N*). Two values count as different only if they are more than two standard errors apart. The runs share most loci, so this is cautious: "within two standard errors" means *not shown to differ*, not *shown to be the same*.
- **r80 counts** have no error term. Paris et al. compared them directly, and so does the script.

### 4.9 Fixed SNPs in the catalog (reported, not a pick)

To see what `n` does, Paris et al. compared the SNPs in each individual with those in the catalog (their `count_fixed_catalog_snps.py`, written for Stacks 1). It counts:

1. the heterozygous SNPs in each individual;
2. the catalog SNPs at which at least one individual is heterozygous;
3. all catalog SNPs;
4. **fixed SNPs** = 3 − 2: catalog SNPs at which no individual is heterozygous.

Paris et al. read a fixed SNP as a place where individuals differ but each is homozygous. The difference only appears when cstacks joins their loci into one catalog locus, which `n` allows. The script does the same count for Stacks 2. It reads the catalog's SNPs (`catalog.snps.tsv.gz`) and each sample's heterozygous SNPs (`SAMPLE.snps.tsv.gz`). It then places each sample's SNPs on the catalog through sstacks's alignment of that sample's locus to its catalog locus (the CIGAR string in `SAMPLE.matches.tsv.gz`). Like Paris's script, it uses only the forward reads, the part cstacks assembles. The report gives count 1 as the mean per individual ("het SNPs/ind"), next to the mean that lands on a catalog SNP ("placed"). The rest are in sample loci that sstacks did not match back to the catalog.

As Paris et al. did, the script **reports** fixed SNPs and doesn't choose from them: they found fixed SNPs rose with `n` in all three datasets without levelling off, and that this "alone ... does not provide enough information to choose a value for n".

**Reading the count in Stacks 2.** sstacks does not match every sample locus back to the catalog. Its lines in `denovo_map.log` count, for each sample, the loci it could not match or excluded. If a catalog SNP's second allele is only in such loci, no individual can be found heterozygous there, so it counts as fixed. To show how much of the count that is, the report splits the fixed SNPs into three groups. The split is this script's, not Paris's; it is only for reading the count, and the count itself is Paris's.

| Group | Meaning |
|---|---|
| **no-match** | In catalog loci that cstacks built but sstacks matched no sample locus back to. These never reach `gstacks` or `populations` |
| **between** | Two or more bases are seen at the SNP among the haplotypes of the individuals that sstacks matched to the locus, although none of them has a heterozygous SNP call there: a difference between individuals |
| **other** | At most one base is seen among those haplotypes. The second allele is in sample loci that sstacks did not match back, or is an N |

"between" is a subset of Paris's count, so it can't be set against their 13–21%.

With an `M = n` design, `M` and `n` change together, so changes in fixed SNPs along that series are not due to `n` alone. Paris et al. plotted fixed SNPs against `n` with `M` fixed (Fig. S3); the n-rule runs (Section 4.5) give that series.

### 4.10 What is not from Paris et al.

Everything the script chooses from is one of Paris et al.'s metrics. The way it reads them, and some of what it reports, is its own:

- **The two-standard-error test, paired by individual,** for steps within individuals (Sections 4.3, 4.8). Paris et al. judged their plots by eye.
- **"Began to flatten out" as a rule:** the first value after which heterozygous loci per individual don't rise beyond 2 SE. These are counted over all of each individual's ustacks loci, because Paris et al. don't say whether theirs were limited to r40/r60/r80 loci.
- **The higher of the r80 maximum and the first level-off point** (Section 4.4, rule 2), when the r80 maximum lies below the level-off point. None of Paris et al.'s datasets did this. The rule extends what they did in all three (never choosing `M` below that point) to a case they did not meet. Taken alone, their r80 rule gives the r80 maximum (M1n1 in the example).
- **Penguin's "relatively uniform"** read as within individuals (Section 3).
- **The binomial standard error on shares of loci** (Section 4.8).
- **The thresholds:** 70% for "heterozygous in a large majority of individuals" (`--het-site`); 10 SNPs for "loci that contained many SNPs" (`--tail-snps`), scaled to the forward read (≥ 5 in the example; `--fwd-tail-snps`); and 1% for Fig. 4's "very few new polymorphic loci" (`--plateau`).
- **Counting SNP-dense loci in the forward read** of paired-end Stacks 2 loci (Section 1.5). Stacks 1.42, which Paris et al. used, built loci from forward reads only.
- **The fixed-SNP count in Stacks 2:** placing each sample's SNPs on the catalog through sstacks's CIGAR alignment, and the no-match / between / other split (Section 4.9).
- **Allele collapse as the share of the tested range's change** reached by each value (Section 4.3).
- **The `M = n` design** of the example, which comes from the Stacks 2 protocol (Rochette & Catchen 2017). Paris et al. varied `M` with `n` held at 0 in their primary runs.

---

## 5. Reading the outputs

Everything goes to the `-o` folder: `paris_report.txt`, `paris_plots.pdf` and `paris_summary.xlsx`, plus `paris_plot_N.png` or `.jpg` with `--image`.

### `paris_report.txt`: the recommendation and the main evidence (printed to the screen with `--verbose`)

Without `--verbose`, the screen shows a short version: the pick (for example `m3 M2 n2 (n provisional)`), one line for each of `m`, `M` and `n` giving the reason, the runs Paris's n rule still needs, any warnings, and the path to this file.

1. **Runs table:** one line per run, with ustacks coverage (mean and lowest sample).
2. **Data line:** paired-end or single-end, forward-read length, and the SNP-dense thresholds used.
3. **For each series:**
   - **Paris et al.'s metrics on r80 loci:** r80 loci, polymorphic r80 loci, the change at each step (Fig. 4), SNPs. A `Paris` column marks the r80 maximum, the Fig. 4 point, where polymorphism within individuals first levels off ("plateau") and the **PICK**.
   - **Within individuals:** loci per individual, allele collapse and heterozygous loci per individual, with paired steps ± 2 SE.
   - **Over-merging signs Paris names:** SNP-dense loci (whole locus and, for paired-end data, forward read) and loci with a ≥ 70% heterozygous site.
   - **The verdict:** "Pick following Paris et al.", each metric's reading, and how the pick was made.
   - **Notes:** the evidence, with Paris's own comparisons.
4. **Recommendation following Paris et al. (2017):** `m`, `M` and `n` together, each with the reason for it, and any runs Paris's n rule still needs.
5. **Fixed SNPs in the catalog (Paris Fig. S3):** Paris et al.'s four counts per run, with count (i) as heterozygous SNPs per individual and the part placed on the catalog. Also their own percentages for comparison, and where the fixed SNPs are (Section 4.9). Reported, not used to choose.
6. **Coverage (Paris Fig. 1):** ustacks's "Coverage after assembling stacks" per sample: after merging alleles, before secondary reads are added (Paris Fig. 1's "merged coverage"). Also samples below 10× (genotyping error rises below 10×, Fountain et al. 2016, cited by Paris et al.).
7. **Warnings:** files that disagree, missing files, runs that aren't part of any series.

### `paris_plots.pdf`: one page per series

| Panel | Shows | Paris |
|---|---|---|
| (a) | r80 loci | Fig. 3a |
| (b) | Polymorphic r80 loci, with the r80 maximum and the Paris pick | Fig. 3b |
| (c) | SNPs in r80 loci | Fig. 3c |
| (d) | New polymorphic r80 loci per step; each bar sits between the two values it compares | Fig. 4 |
| (e) | SNPs per polymorphic r80 locus (log scale), one line per value | loci with many SNPs |
| (f) | Loci per individual (one box per value, one point per individual) | within individuals |
| (g) | Allele collapse: share of each individual's stacks merged into loci | Table S6a |
| (h) | Heterozygous loci per individual, paired steps ± 2 SE along the top, and where it first levels off | within individuals; decides `M` |
| (i) | Heterozygous SNPs per individual (ustacks), with paired steps the same way | fixed-SNP count (i) |
| (j) | Coverage per sample (ustacks, after merging alleles), with 10× and 25× lines | Fig. 1 |
| (k) | SNP-dense share (forward read for paired-end data), ± 2 SE, with the Paris pick | over-merging sign |
| (l) | Share with a ≥ 70% heterozygous site, ± 2 SE, with the Paris pick | over-merging sign |
| (m) | Mean heterozygous SNPs per individual placed on a catalog SNP | fixed-SNP count (i), placed |
| (n) | All catalog SNPs, and those at which at least one individual is heterozygous | fixed-SNP counts (ii), (iii) |
| (o) | Fixed SNPs, with their share of catalog SNPs, split into the three groups of Section 4.9 | Fig. S3 |

Rows: r80 loci on top, within individuals and coverage in the middle, the over-merging signs and the catalog comparison at the bottom. The verdict and as many notes as fit are printed down the right side; the report has them all.

### `paris_summary.xlsx`: supporting tables

| Sheet | Contents |
|---|---|
| `Runs` | One row per run, with every number the script computes |
| `Series` | One row per value in each series: the change from the previous value, the within-individual means and paired steps (± 2 SE, how many individuals went up), the over-merging shares, catalog and fixed SNPs, and which value is the r80 maximum, Fig. 4 point, plateau and pick |
| `Individuals` | One row per run and sample: ustacks stacks, loci, share of stacks merged, heterozygous loci and heterozygous SNPs (all, and those placed on the catalog), and ustacks coverage |
| `SNPs_per_locus` | Counts of r80 loci with 0, 1, 2, … SNPs and their proportions, with a line chart |
| `Fixed_SNPs` | One row per run: Paris et al.'s four counts (count (i) as heterozygous SNPs per individual, all and placed on the catalog), the fixed share, and the three groups of Section 4.9 |
| `High_het_loci` | For each run, the IDs of polymorphic r80 loci with a site where ≥ 70% of all samples are heterozygous: a ready-made `populations --blacklist` (Section 4.7) |
| `Coverage` | ustacks coverage per sample (rows) and run (columns) |

---

## 6. Worked example: Alabama pearlshell

The design was m3 with `M = n` from 1 to 6 on paired-end data (forward read 140 bp), filtered with `populations -R 0.8`, 44 samples in 7 populations. The full report is `paris_results/paris_report.txt`.

### 6.1 Paris et al.'s metrics on r80 loci

| | r80 loci | Polymorphic r80 loci (*P*) | Change in *P* (Fig. 4) | SNPs |
|---|---|---|---|---|
| M1n1 | 18,299 | **12,314** | — | 25,767 |
| M2n2 | 16,048 | 11,100 | −1,214 (−9.9%) | 23,884 |
| M3n3 | 15,251 | 10,603 | −497 (−4.5%) | 23,898 |
| M4n4 | 14,921 | 10,408 | −195 (−1.8%) | 23,838 |
| M5n5 | 14,713 | 10,263 | −145 (−1.4%) | 24,002 |
| M6n6 | 14,576 | 10,193 | −70 (−0.7%) | 24,138 |

*P* falls at every step, so the r80 maximum is M1n1 and no step adds new polymorphic r80 loci. None of Paris et al.'s datasets did this. SNPs fall from M1n1 to M2n2 and are flat after that: even less of a rise than in Paris et al.'s "relatively monomorphic" penguin data, where they chose M2.

### 6.2 Within individuals

Means over 44 individuals; steps are paired by individual, ± 2 SE, * = beyond 2 SE.

| | Loci per individual | Stacks merged into loci | Heterozygous loci per individual | Step | Individuals up |
|---|---|---|---|---|---|
| M1n1 | 40,420 | 13.3% | 2,360.6 | | |
| M2n2 | 35,407 | 23.8% | 2,426.8 | +66 ± 33* | 37 of 44 |
| M3n3 | 33,534 | 27.7% | 2,434.3 | +8 ± 20 | 22 of 44 |
| M4n4 | 32,779 | 29.2% | 2,494.2 | +60 ± 18* | 37 of 44 |
| M5n5 | 32,376 | 29.9% | 2,587.5 | +93 ± 12* | 44 of 44 |
| M6n6 | 32,121 | 30.2% | 2,701.7 | +114 ± 20* | 44 of 44 |

- **M1n1 is too low by Paris et al.'s own measure.** Going to M2n2, heterozygous loci per individual rise in 37 of 44 individuals (+66 ± 33), and another 10 percentage points of stacks merge into loci.
- **Polymorphism within individuals first levels off at M2n2.** The step from M2n2 to M3n3 is +8 ± 20, with 22 of 44 individuals up.
- **Allele collapse:** 62% of the change across the tested range has happened by M2n2, 85% by M3n3 and 94% by M4n4. Paris et al. saw collapse level out between M1 and M2 (penguin), M2 and M3 (earthworm), and M3 and M4 (trout).
- **Why *P* falls.** From M1n1 to M2n2, each individual loses about 5,000 loci but gains only 66 heterozygous loci. Most stacks merged at that step don't form a heterozygous locus, so raising `M` removes loci rather than turning monomorphic loci polymorphic. *P* can therefore fall while polymorphism within individuals still rises. By Paris et al.'s measure, alleles are still being joined at the r80 maximum (M1n1).
- **Above M2n2, heterozygous loci per individual rise again,** faster at each step (+60, +93, +114; 37, 44 and 44 of 44 individuals up). The level-off at M2n2 is therefore the first one, not a lasting plateau. Paris et al.'s wording, "began to flatten out" (trout) and "began to plateau" (earthworm), refers to where the levelling starts. They don't say whether polymorphism within individuals later rose again in their data. They did find that higher `M` adds "a small number of loci that contained many SNPs" at each increment. In the earthworm data, polymorphism across 80% of the population kept rising, which they put down to "a small number of loci with a high-density of SNPs". They chose M3 there, where polymorphism within individuals began to plateau.

### 6.3 The pick

The r80 maximum (M1n1) lies below where polymorphism within individuals first levels off (M2n2). In all three of Paris et al.'s datasets they chose at or above that point, so rule 2 of Section 4.4 picks **M2n2**. None of their datasets had the r80 maximum below the level-off point, so this applies their judgement to a case they did not meet. Taken alone, their r80 rule gives M1n1. Until the n-rule runs are done (Section 6.6), `n = 2` is provisional.

### 6.4 Over-merging signs Paris names

| | ≥ 5 SNPs in the forward read | ≥ 70% heterozygous site |
|---|---|---|
| M1n1 | 0.19% | 1.68% |
| M2n2 | 0.35% | 1.50% |
| M3n3 | 0.80% | 1.53% |
| M4n4 | 1.20% | 1.51% |
| M5n5 | 1.91% | 1.68% |
| M6n6 | 2.25% | 1.70% |

From M2n2 to M6n6, the share of loci with ≥ 5 SNPs in the forward read rises from 0.35% to 2.25%, more than two standard errors. The share with a ≥ 70% heterozygous site stays at 1.5–1.7%, within two standard errors. Mean ustacks coverage after merging alleles is 41–46× per run. That is between Paris et al.'s trout (28×) and their best-covered dataset, penguin (49×); earthworm had 21×. By their reasoning, SNP-dense loci are more likely to be real at high coverage. Paris et al. filtered likely paralogs after assembly. The `High_het_loci` sheet lists, for each run, the loci with a ≥ 70% heterozygous site, ready to pass to `populations --blacklist` (Section 4.7). The 70% is this script's threshold.

### 6.5 Fixed SNPs in the catalog

Paris et al.'s counts:

| | Heterozygous SNPs per individual (i) | placed on the catalog | Catalog SNPs (iii) | Heterozygous in an individual (ii) | Fixed (iv) |
|---|---|---|---|---|---|
| M1n1 | 2,722 | 2,028 | 64,824 | 30,295 | 34,529 (53.3%) |
| M2n2 | 3,465 | 2,441 | 74,995 | 32,466 | 42,529 (56.7%) |
| M3n3 | 4,040 | 2,815 | 74,676 | 35,198 | 39,478 (52.9%) |
| M4n4 | 4,672 | 3,234 | 78,784 | 38,504 | 40,280 (51.1%) |
| M5n5 | 5,395 | 3,654 | 86,095 | 41,912 | 44,183 (51.3%) |
| M6n6 | 6,230 | 4,142 | 93,350 | 44,863 | 48,487 (51.9%) |

Where the fixed SNPs are (Section 4.9; not from Paris), as shares of all catalog SNPs:

| | no-match | (catalog loci) | between | other |
|---|---|---|---|---|
| M1n1 | 3,388 (5.2%) | 200 | 19,157 (29.6%) | 11,984 (18.5%) |
| M2n2 | 17,047 (22.7%) | 421 | 16,267 (21.7%) | 9,215 (12.3%) |
| M3n3 | 15,248 (20.4%) | 365 | 16,410 (22.0%) | 7,820 (10.5%) |
| M4n4 | 14,723 (18.7%) | 407 | 17,655 (22.4%) | 7,902 (10.0%) |
| M5n5 | 15,887 (18.5%) | 450 | 18,771 (21.8%) | 9,525 (11.1%) |
| M6n6 | 18,056 (19.3%) | 491 | 20,177 (21.6%) | 10,254 (11.0%) |

- **Paris's count is about half of all catalog SNPs at every value,** well above the 13–21% Paris et al. report at the `n` they chose.
- **Much of the difference is in loci the data never use.** At M2n2, 17,047 fixed SNPs sit in 421 catalog loci that sstacks matched no sample back to (about 40 SNPs per locus on average), so they never reach `gstacks` or `populations`. For 9,215 more ("other"), the matched individuals' haplotypes show at most one base.
- **"between"** (two or more bases among the matched individuals, but no heterozygous call) is about 22% of catalog SNPs from M2n2 on, steady as `M = n` rises.
- **Heterozygous SNPs per individual** rise from 2,722 at M1n1 to 6,230 at M6n6. Between two-thirds and three-quarters of them land on a catalog SNP; the rest are in sample loci that sstacks did not match back.
- Because `M` and `n` change together here, these numbers can't show what `n` alone does. The n-rule runs (next section) give an `n` series at M2, which the script will report the way Paris et al. plotted it (Fig. S3). As in Paris et al., the count does not enter the choice of `n`.

### 6.6 Still to do: the n rule

The `M = n` series fixes `n = M`. Paris et al.'s n rule compares `n = 1, 2, 3` at `M = 2`, so **run m3 M2 n1 and m3 M2 n3**, each with the r80 filter. `denovo_map_n_rule.sh` in the project folder does both (`sbatch denovo_map_n_rule.sh`). It gives `denovo_map.pl` the option `-X "populations:-R 0.8"` and no `-r` (Section 2.1), so the r80 output is in `2_mmar_n1/` and `2_mmar_n3/` themselves. Then rerun the script on all eight r80 folders:

```bash
python paris_optimize.py {1..6}_mmar_R80 2_mmar_n1 2_mmar_n3 -o paris_results
```

It will add an `n` series at M2 and report the final `n`.

A **sensitivity check** is not part of Paris et al.'s method, but is worth doing downstream. The two obvious comparisons are M1n1, which the r80 rule alone gives, and M3n3, the next value, where heterozygous loci per individual are still flat. Both assemblies already exist. If your biological conclusions are the same, say so in your methods.

### 6.7 Coverage

ustacks coverage after merging alleles ("Coverage after assembling stacks") averages 41–46× per run. Three samples are below 10× in all six runs:
- A063: 7.0–7.4×;
- ALP_MBR_01: 7.9–8.0×;
- ALP_MBR_02: 8.6–8.8×.

AP6 is just under 10× at M1n1 (9.97×) and 10.2–10.5× in the other runs. The genotypes of these samples are less reliable (Fountain et al. 2016, cited by Paris et al.). MBR1 consists only of ALP_MBR_01 and ALP_MBR_02, so any conclusion about MBR1 rests on two low-coverage samples.

### 6.8 Example methods paragraph

Fill in the bracketed parts once the n-rule runs and any downstream checks are done:

> We optimised de novo assembly parameters following Paris et al. (2017). Holding m = 3, we assembled loci with M = n from 1 to 6, varying M and n together as in the Stacks 2 protocol (Rochette & Catchen 2017). We filtered each assembly to loci present in at least 80% of all individuals (populations -R 0.8; "r80 loci"). For each assembly we recorded the number of r80 loci, polymorphic r80 loci and SNPs. Within each individual (ustacks output), we recorded the number of loci, the proportion of stacks merged into loci and the number of heterozygous loci. Polymorphic r80 loci were highest at M = 1 and declined with every increase in M. However, heterozygous loci per individual rose from M = 1 to M = 2, and first levelled off at M = 2 (paired comparisons across 44 individuals; a step counted as a rise if its mean exceeded two standard errors). This indicated that alleles were not yet fully merged at M = 1. In all three of their datasets, Paris et al. chose M at or above the point where polymorphism within individuals levelled off. Following that judgement, we chose M = 2. We then compared n = 1, 2 and 3 at M = 2 and chose n = [n], which gave the most polymorphic r80 loci (Paris et al.'s n = M ± 1 rule). [Paralog filter, if used: loci with a site at which at least 70% of individuals were heterozygous were removed.] [Downstream results were (unchanged / changed as follows) with M = n = 1 and M = n = 3.]

---

## 7. Troubleshooting

| Message | Meaning and fix |
|---|---|
| `populations was not run with -R` | You passed a folder filtered with `-r` (or not filtered), such as a `denovo_map.pl` folder made with `-r 0.8`. Run `populations -R 0.8` on it and pass that output folder. For new assemblies, you can instead give `denovo_map.pl` `-X "populations:-R 0.8"` without `-r` (Section 2.1). |
| `a per-population -r … filter was applied as well as -R …` | `populations` used both filters, for example because `denovo_map.pl` was given `-r 0.8` and `-X "populations:-R 0.8"`. Those loci are not Paris's r80. Rerun `populations` with `-R 0.8` alone. |
| `could not find the assembly for this run` | The `-P` path in `populations.log` doesn't lead to an assembly folder from any of the places the script looked (they are listed in the message). Give it in a runs table (Section 2.4). |
| `no denovo_map.log, so m, M and n are unknown` | The assembly wasn't made with `denovo_map.pl`. Give `m`, `M`, `n` in a runs table. |
| `no per-sample ustacks files (SAMPLE.alleles.tsv.gz)` | The per-sample ustacks output was deleted or moved. Heterozygous loci per individual are skipped, so there is no within-individual level-off point: `M` comes from the r80 rule alone (or Fig. 4, if *P* is still rising), and the report says so. |
| `no populations.haplotypes.tsv` | The sample list is taken from `denovo_map.log` instead. Check that it matches the samples `populations` used. |
| `no 'Coverage after assembling stacks' lines` | An older or unusual ustacks log. Coverage falls back to ustacks's "Final coverage", which also counts secondary reads, so it runs a little higher than Paris's merged coverage. |
| `No pick following Paris et al. in the tested range` | Heterozygous loci per individual (or *P*, if it is still rising) haven't levelled off by the highest value. Test higher values (Paris et al. suggest M1 to M8). |
| `Runs …: no populations.sumstats.tsv …` (or another file) | The metrics named in the message can't be computed for those runs. The rest still work. |
| `… is not part of any series` | That run has no partner differing in just one parameter (an `M` series with `n` fixed needs three values), so it isn't used in any decision. |
| `runs differ in sample count / -R value / populations options` | You are comparing runs that were filtered differently. Re-run `populations` the same way for all of them. |
| `populations was run with -R 0.7 …` | Not Paris's 0.8. The script runs, but labels and interpretation follow your value. |
| `openpyxl not available` / `matplotlib not available` | Install them (`conda install openpyxl matplotlib`) or use a Python environment that has them. Without matplotlib there is no PDF and no `--image` output. |
| `… was given more than once; using it once` | The same folder appeared twice. It is analysed once. |
| `no catalog.snps.tsv.gz, so no fixed-SNP count` (or `no per-sample SAMPLE.snps.tsv.gz and SAMPLE.matches.tsv.gz`) | The catalog or the per-sample ustacks/sstacks files were deleted or moved. Only the fixed-SNP section is lost; the recommendation doesn't use it. If only some samples' files are missing, the count uses the others and says so. |
| `could not read the sstacks or ustacks file of SAMPLE` | That sample's file is damaged. It is left out of the fixed-SNP count. |
| `No series found` | You need at least two runs that differ in exactly one parameter. |

---

## 8. Glossary

| Term | Meaning |
|---|---|
| **Locus** | One place in the genome; in RAD-seq, the DNA next to one cut site. |
| **Allele** | One version of a locus's sequence. A diploid individual carries two. |
| **Stack** | A pile of identical reads within one individual; a candidate allele. |
| **Allele collapse** | Stacks merging into loci within an individual under `M`; Paris et al.'s "proportion of putative alleles that collapsed". |
| **Catalog** | Stacks' combined list of loci across all individuals. |
| **Fixed SNP** | A catalog SNP at which no individual is heterozygous (Paris et al.'s Fig. S3). Paris et al. read it as a place where individuals differ but each is homozygous. |
| **Forward read / paired-end contig** | The two parts of a paired-end Stacks 2 locus. Only the forward read is assembled under `M` and `n`. |
| **Paralog** | A duplicated copy of a sequence elsewhere in the genome. It looks like the original, but it is a different locus. |
| **Under-merging** | Alleles of one locus kept apart as separate loci (`M` or `n` too low). |
| **Over-merging** | Different loci (often paralogs) fused into one (`M` or `n` too high). |
| **r80 locus** | A locus present in at least 80% of individuals (`populations -R 0.8`). |
| **Polymorphic locus** | A locus with at least one SNP. |
| **Plateau (within individuals)** | The first value after which heterozygous loci per individual don't rise beyond two standard errors: where they "began to flatten out", in Paris et al.'s words. Later steps can rise again. |
| **Paired step** | The change from one value to the next computed within each individual, then averaged; its standard error comes from the spread of those changes. |
| **Coverage** | Average number of reads per locus in a sample. Here it is ustacks's "Coverage after assembling stacks": after merging alleles, before secondary reads are added (Paris Fig. 1's "merged coverage"). |
| **Standard error (SE)** | How much a number would vary from sampling noise alone. Values more than about 2 SE apart are unlikely to differ by chance. |

---

## 9. References

- Paris, J.R., Stevens, J.R. & Catchen, J.M. (2017) Lost in parameter space: a road map for STACKS. *Methods in Ecology and Evolution*, 8, 1360–1373. doi:10.1111/2041-210X.12775
- Rochette, N.C. & Catchen, J.M. (2017) Deriving genotypes from RAD-seq short-read data using Stacks. *Nature Protocols*, 12, 2640–2659.
- Fountain, E.D., Pauli, J.N., Reid, B.N., Palsbøll, P.J. & Peery, M.Z. (2016) Finding the right coverage: the impact of coverage and sequence quality on SNP genotyping error rates. *Molecular Ecology Resources*, 16, 966–978.
