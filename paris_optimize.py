#!/usr/bin/env python3
"""
paris_optimize.py -- choose Stacks 2 de novo parameters (m, M, n) following Paris et al. (2017),
"Lost in parameter space: a road map for STACKS", Methods in Ecology and Evolution 8:1360-1373.
Only Paris et al.'s metrics and reasoning are used. README_paris_optimize.md is the full guide.

Usage
  paris_optimize.py POPDIR [POPDIR ...] [-o OUTDIR] [--image png|jpeg] [-v|--verbose]
  paris_optimize.py --runs runs.tsv [-o OUTDIR] [-v|--verbose]

Each POPDIR is the output of `populations -R 0.8` (Paris's "r80" filter) run on one de novo
assembly. The assembly directory is found from the -P path in populations.log. A --runs table can
name it instead, and can give m, M and n when there is no denovo_map.log (pipelines run by hand).
POPDIR can be the assembly directory itself when denovo_map.pl was run with
-X "populations:-R 0.8" and without -r, because populations then writes its -R output there.

Files read (Stacks 2; developed on 2.68). Missing optional files switch off what needs them.
  POPDIR/populations.log, populations.log.distribs   r80 loci, polymorphic r80 loci, SNPs per locus
  POPDIR/populations.haplotypes.tsv                  the samples populations used; r80 locus IDs
  POPDIR/populations.sumstats.tsv                    SNP positions; heterozygosity per site
  DENOVO/denovo_map.log                              m, M, n; per-sample ustacks stacks, loci, coverage
                                                     (ustacks/cstacks/sstacks lines only)
  DENOVO/SAMPLE.alleles.tsv.gz                       heterozygous loci per individual (ustacks)
  DENOVO/SAMPLE.snps.tsv.gz, SAMPLE.matches.tsv.gz   each individual's SNPs, placed on the catalog
  DENOVO/catalog.snps.tsv.gz                         catalog SNPs (for fixed SNPs, Paris Fig. S3)
  DENOVO/catalog.fa.gz                               forward read vs paired-end contig of each locus
From DENOVO, only files written before populations runs are read. Its populations output is read only
when DENOVO is itself the POPDIR, and a POPDIR must be -R output (populations.log is checked).

Runs are grouped into series in which one parameter varies (M with n = M, M, n or m) and the rest,
the sample count and the -R value are shared. Paris et al.'s metrics:
  - r80 rule: the value with the most polymorphic r80 loci (for n, only n = M - 1, M, M + 1).
  - Fig. 4: new polymorphic r80 loci added at each step.
  - Within individuals (each sample's ustacks output): loci per individual, allele collapse (share
    of stacks merged into loci), and heterozygous loci per individual (Paris's "polymorphic loci
    within individuals").
  - Over-merging signs: SNP-dense loci, and loci heterozygous in a large majority of individuals.
  - Coverage (Fig. 1).
  - Fixed SNPs: catalog SNPs at which no individual is heterozygous (Fig. S3). Reported, not used to
    choose, as in Paris et al.
For M (and M = n), the pick is the higher of the r80 maximum and the first value after which
heterozygous loci per individual stop rising beyond two standard errors (paired by individual), the
point where Paris et al.'s polymorphic loci within individuals "began to flatten out". If polymorphic
r80 loci are still rising at the highest value, the pick is that level-off point instead. This
reproduces Paris et al.'s choices for their trout (M5), penguin (M2) and earthworm (M3) data; where
the r80 maximum lies below the level-off point, it applies their judgement to a case they did not
meet. For n and m the r80 rule decides. The report then combines m, M and n and lists the n = M - 1
and M + 1 runs Paris's n rule still needs. README_paris_optimize.md lists what is not from Paris.

Outputs, in -o (default paris_results/):
  paris_report.txt     the recommendation and the main evidence
  paris_plots.pdf      one page of plots per series (needs matplotlib)
  paris_summary.xlsx   supporting tables: per run, per series value, per individual, SNPs per
                       locus, fixed SNPs, high-heterozygosity loci, coverage per sample (needs openpyxl)
  paris_plot_N.png or paris_plot_N.jpg, one per page, with --image png or --image jpeg
The screen shows only the recommendation (one line per parameter, runs still to do) and any warnings;
--verbose prints the full paris_report.txt instead.
Parsing needs only the Python 3.6+ standard library.
"""

import argparse
import gzip
import multiprocessing
import math
import os
import re
import shlex
import statistics
import sys
import textwrap
import zlib
from collections import Counter, defaultdict

LOW_SAMPLE_COV = 10.0   # Fountain et al. 2016, cited by Paris et al.
LOW_MEAN_COV = 25.0     # Paris et al. suggest aiming above 25x
NRUN = re.compile(r"N{5,}")             # spacer between forward read and paired-end contig
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")

WARNINGS = []
MISSING = {}    # (what is missing, what is lost) -> directories


def warn(msg):
    if msg not in WARNINGS:
        WARNINGS.append(msg)


def show_path(p):
    """A path as the user would type it: relative below the working directory, else absolute."""
    rel = os.path.relpath(p)
    return os.path.abspath(p) if rel.startswith("..") else rel


def missing(directory, what, lost):
    """Record a missing optional file; reported once per file in flush_missing()."""
    MISSING.setdefault((what, lost), [])
    d = show_path(directory)
    if d not in MISSING[(what, lost)]:
        MISSING[(what, lost)].append(d)


def flush_missing(n_runs):
    for (what, lost), dirs in MISSING.items():
        if len(dirs) == n_runs and n_runs > 1:
            warn("No run has a {}, so {}.".format(what, lost))
        else:
            warn("{} {}: no {}, so {} for {}.".format(
                "Run" if len(dirs) == 1 else "Runs", ", ".join(dirs), what, lost,
                "this run" if len(dirs) == 1 else "these runs"))


def say(text):
    """Write to stdout even under an ASCII locale (Python 3.6 does not force UTF-8)."""
    try:
        sys.stdout.write(text)
    except UnicodeEncodeError:
        sys.stdout.flush()
        sys.stdout.buffer.write(text.encode("utf-8"))
        sys.stdout.buffer.flush()


class ParseError(Exception):
    pass


# --------------------------------------------------------------------------- small helpers

def rlabel(R):
    """Paris's name for the loci kept by populations -R: r80 for -R 0.8."""
    return "r{:g}".format(round(R * 100, 1))


def relabel(text, R):
    lab = rlabel(R)
    return text if lab == "r80" else text.replace("r80", lab)


def pct(x, digits=1):
    return "{:+.{d}%}".format(x, d=digits) if x == x else "n/a"


def share(k, n):
    return k / n if n else float("nan")


def shares_differ(p1, n1, p2, n2, z=2.0):
    """True if two proportions differ by more than z binomial standard errors.

    The runs are treated as independent samples. They share most loci, so this is cautious for a
    claim that two values differ; a real difference can still sit inside it.
    """
    if not n1 or not n2 or p1 != p1 or p2 != p2:
        return False
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    return abs(p1 - p2) > z * se


def se_share(p, n):
    return math.sqrt(p * (1 - p) / n) if n and p == p else 0.0


def join_list(items, last="and"):
    items = list(items)
    if len(items) < 2:
        return "".join(items)
    return "{} {} {}".format(", ".join(items[:-1]), last, items[-1])


def isare(items):
    return "is" if len(items) == 1 else "are"


# --------------------------------------------------------------------------- parsing

def read_text(path):
    if not os.path.isfile(path):
        raise ParseError("missing file: {}".format(path))
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def read_distribs(path):
    """Return {section_name: [row, ...]} for a Stacks .distribs file (comments dropped)."""
    sections = {}
    cur = None
    for line in read_text(path).splitlines():
        m = re.match(r"^BEGIN (\S+)\s*$", line)
        if m:
            cur = m.group(1)
            sections[cur] = []
            continue
        if re.match(r"^END (\S+)\s*$", line):
            cur = None
            continue
        if cur is not None and line.strip() and not line.startswith("#"):
            sections[cur].append(line.split("\t"))
    return sections


def as_table(rows):
    if not rows:
        return []
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:]]


def get_opt(tokens, names):
    for i, t in enumerate(tokens):
        if t in names and i + 1 < len(tokens):
            return tokens[i + 1]
        for nm in names:
            if nm.startswith("--") and t.startswith(nm + "="):
                return t.split("=", 1)[1]
    return None


def strip_opts(tokens, names):
    """Drop options (and their values) that legitimately differ between runs."""
    out = []
    skip = False
    for t in tokens:
        if skip:
            skip = False
            continue
        if t in names:
            skip = True
            continue
        if any(nm.startswith("--") and t.startswith(nm + "=") for nm in names):
            continue
        out.append(t)
    return " ".join(out)


def search(pattern, text, what, path, cast=float):
    m = re.search(pattern, text, re.M)
    if not m:
        raise ParseError("could not find {} in {}".format(what, path))
    return cast(m.group(1))


def header_cols(line):
    return [c.strip() for c in line.lstrip("#").rstrip("\n").split("\t")]


def col_of(header, name, default):
    if header and name in header:
        return header.index(name)
    return default


def parse_populations(popdir):
    log_path = os.path.join(popdir, "populations.log")
    text = read_text(log_path)
    lines = text.splitlines()
    if len(lines) < 2:
        raise ParseError("{} is too short to be a populations log".format(log_path))
    vm = re.match(r"populations v(\d+)\.(\d+)", lines[0])
    if vm and int(vm.group(1)) < 2:
        raise ParseError("{}: written by Stacks {}.{}; this script reads Stacks 2 output".format(
            log_path, vm.group(1), vm.group(2)))
    if not vm:
        warn("{}: could not read the Stacks version from the first line".format(log_path))
    tokens = shlex.split(lines[1])
    info = {
        "version": "{}.{}".format(vm.group(1), vm.group(2)) if vm else "",
        "command": lines[1].strip(),
        "in_path": get_opt(tokens, ["-P", "--in-path", "--in_path"]),
        "out_path": get_opt(tokens, ["-O", "--out-path", "--out_path"]),
        "other_flags": strip_opts(tokens[1:], ["-P", "--in-path", "--in_path", "-O", "--out-path", "--out_path",
                                               "-t", "--threads"]),
    }
    info["R"] = search(r"Percent samples overall:\s*([\d.]+)", text, "'Percent samples overall'", log_path)
    info["r"] = search(r"Percent samples limit per population:\s*([\d.]+)", text,
                       "'Percent samples limit per population'", log_path)
    if info["R"] <= 0:
        raise ParseError(
            "{}: populations was not run with -R (Percent samples overall = 0). Run "
            "`populations -R 0.8` on this assembly and pass that output directory.".format(popdir))
    if info["r"] > 0:
        warn("{}: a per-population -r {:g} filter was applied as well as -R {:g}; both filters "
             "shape these loci.".format(popdir, info["r"], info["R"]))
    info["n_samples"] = search(r"Working on (\d+) samples", text, "'Working on N samples'", log_path, int)
    m = re.search(r"Kept (\d+) loci, composed of (\d+) sites; (\d+) of those sites were filtered, "
                  r"(\d+) variant sites remained", text)
    if not m:
        raise ParseError("could not find the 'Kept ... loci' line in {} (did populations finish?)".format(log_path))
    info["kept_loci"] = int(m.group(1))
    info["kept_variant_sites"] = int(m.group(4))

    d = read_distribs(os.path.join(popdir, "populations.log.distribs"))
    rows = d.get("snps_per_loc_postfilters")
    if not rows:
        raise ParseError("no snps_per_loc_postfilters section in {}/populations.log.distribs".format(popdir))
    info["snp_dist"] = {int(r["n_snps"]): int(r["n_loci"]) for r in as_table(rows)}
    return info


def is_assembly_dir(d):
    return any(os.path.isfile(os.path.join(d, f)) for f in ("denovo_map.log", "gstacks.log", "catalog.fa.gz"))


def resolve_denovo_dir(popdir, in_path, out_path):
    """Find the assembly that populations read (-P), trying the directory populations ran from first.

    populations.log records -P and -O as typed. If -O was relative, the directory populations ran
    from is POPDIR with the -O path taken off its end, and -P is relative to that.
    """
    pop_abs = os.path.abspath(popdir)
    cands = []
    if in_path:
        if os.path.isabs(in_path):
            cands.append(in_path)
        else:
            if out_path and not os.path.isabs(out_path):
                parts = [p for p in os.path.normpath(out_path).split(os.sep) if p not in ("", ".")]
                here = pop_abs.split(os.sep)
                if parts and ".." not in parts and here[-len(parts):] == parts:
                    run_dir = os.sep.join(here[:-len(parts)]) or os.sep
                    cands.append(os.path.join(run_dir, in_path))
            cands.append(os.path.join(os.path.dirname(pop_abs), in_path))
            cands.append(os.path.abspath(in_path))
    cands.append(pop_abs)
    seen = []
    for c in cands:
        c = os.path.normpath(c)
        if c not in seen:
            seen.append(c)
    for c in seen:
        if is_assembly_dir(c):
            return c
    raise ParseError("{}: could not find the assembly for this run (populations -P was '{}'; looked in: {}). "
                     "Name it in a --runs table (column denovo_dir).".format(popdir, in_path, ", ".join(seen)))


def parse_denovo(denovo_dir, given):
    path = os.path.join(denovo_dir, "denovo_map.log")
    found = {"m": None, "M": None, "n": None}
    ustacks_cov, final_cov = {}, {}
    ustacks = {}
    if os.path.isfile(path):
        text = read_text(path)
        lines = text.splitlines()
        cmd_tokens = shlex.split(lines[1]) if len(lines) > 1 else []
        # Only the assembly stages (ustacks to gstacks) are parsed. The populations step that denovo_map.pl runs at
        # the end (often with -r, not Paris's -R) is cut off here and never read; r80 numbers come from POPDIR.
        cut_at = next((k for k, ln in enumerate(lines) if k > 1 and ln.strip() == "populations"), len(lines))
        asm_lines = lines[:cut_at]
        asm_text = "\n".join(asm_lines)

        def realised(pattern, label):
            allv = [int(v) for v in re.findall(pattern, asm_text, re.M)]
            if len(set(allv)) > 1:
                warn("{}: {} differs between samples ({}); using the most common value".format(
                    path, label, sorted(set(allv))))
            return max(set(allv), key=allv.count) if allv else None

        found["m"] = realised(r"Min depth of coverage to create a stack \(m\):\s*(\d+)", "ustacks m")
        found["M"] = realised(r"Max distance allowed between stacks \(M\):\s*(\d+)", "ustacks M")
        found["n"] = realised(r"Number of mismatches allowed between stacks:\s*(\d+)", "cstacks n")
        if found["n"] is None:
            mm = re.search(r"cstacks .*?-n\s+(\d+)", asm_text)
            found["n"] = int(mm.group(1)) if mm else None
        for key in ("m", "M", "n"):
            cv = get_opt(cmd_tokens, ["-" + key])
            if cv is not None and found[key] is not None and int(cv) != found[key]:
                warn("{}: denovo_map.pl was given -{} {} but the log shows {} was used".format(
                    path, key, cv, found[key]))
            if found[key] is None and cv is not None:
                found[key] = int(cv)
        if "denovo_map.pl is done" not in text and "completed at" not in "\n".join(lines[-3:]):
            warn("{}: denovo_map.pl does not appear to have finished".format(path))
        # Per-sample ustacks lines. sstacks also prints "Sample i of N" headers but none of these lines, so the
        # first of each after a header is the ustacks one.
        sample = None
        for line in asm_lines:
            mm = re.match(r"^Sample \d+ of \d+ '(.+)'", line)
            if mm:
                sample = mm.group(1)
                continue
            if sample is None:
                continue
            mm = re.match(r"^\s*Assembled (\d+) stacks into (\d+); blacklisted", line)
            if mm:
                st = ustacks.setdefault(sample, {})
                st.setdefault("stacks", int(mm.group(1)))
                st.setdefault("loci_M", int(mm.group(2)))
                continue
            mm = re.match(r"^Final number of stacks built: (\d+)", line)
            if mm:
                ustacks.setdefault(sample, {}).setdefault("final", int(mm.group(1)))
                continue
            # Paris Fig. 1's 'merged coverage' is 'the coverage after merging alleles': ustacks's coverage after
            # assembling stacks under M, before secondary reads are added.
            mm = re.match(r"^Coverage after assembling stacks: mean=([\d.]+)", line)
            if mm:
                ustacks_cov.setdefault(sample, float(mm.group(1)))
                continue
            mm = re.match(r"^Final coverage: mean=([\d.]+)", line)
            if mm:
                final_cov.setdefault(sample, float(mm.group(1)))
        if not ustacks_cov and final_cov:
            ustacks_cov = final_cov
            warn("{}: no 'Coverage after assembling stacks' lines, so coverage is ustacks's 'Final coverage', which "
                 "also counts secondary reads".format(path))

    params = {}
    for key in ("m", "M", "n"):
        g = given.get(key)
        if g is not None:
            if found[key] is not None and found[key] != g:
                warn("{}: the --runs table gives {} = {} but denovo_map.log shows {}; using {}".format(
                    denovo_dir, key, g, found[key], g))
            params[key] = g
        elif found[key] is not None:
            params[key] = found[key]
    unknown = [k for k in ("m", "M", "n") if k not in params]
    if unknown:
        if os.path.isfile(path):
            raise ParseError("{}: could not determine {}".format(path, ", ".join(unknown)))
        raise ParseError("{}: no denovo_map.log, so {} {} unknown. Give them in a --runs table (columns m, M, n)."
                         .format(denovo_dir, join_list(unknown), isare(unknown)))

    return {"m": params["m"], "M": params["M"], "n": params["n"], "ustacks_cov": ustacks_cov, "ustacks": ustacks}


def read_haplotypes(popdir):
    """From populations.haplotypes.tsv: the IDs of every locus populations kept, and the samples it genotyped.
    (None, None) if the file is missing."""
    path = os.path.join(popdir, "populations.haplotypes.tsv")
    if not os.path.isfile(path):
        return None, None
    ids = set()
    samples = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                cols = [c.strip() for c in line.lstrip("#").rstrip("\n").split("\t")]
                if samples is None and len(cols) > 2 and "Locus" in cols[0]:
                    samples = cols[2:]
                continue
            if line.strip():
                ids.add(line.split("\t", 1)[0].strip())
    return ids, samples or []


def read_ind_het(denovo_dir, samples):
    """Heterozygous loci per sample: loci with two or more haplotypes in SAMPLE.alleles.tsv.gz (ustacks)."""
    out = {}
    for smp in samples:
        path = os.path.join(denovo_dir, smp + ".alleles.tsv.gz")
        if not os.path.isfile(path):
            continue
        cnt = Counter()
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cnt[line.split("\t", 1)[0]] += 1
        out[smp] = sum(1 for v in cnt.values() if v >= 2)
    return out


# Paris et al.'s within-individual metrics, from each sample's ustacks output (before the catalog is built).
IND_METRICS = [
    ("ind_loci", lambda d: d["final"]),
    ("ind_collapse", lambda d: 1 - d["loci_M"] / d["stacks"] if d["stacks"] and d["loci_M"] is not None else None),
    ("ind_het", lambda d: d["het"]),                     # polymorphic loci within individuals: decides M
    ("ind_snps", lambda d: d["het_snps"]),               # count (i) of Paris's fixed-SNP comparison; reported only
]


def fasta_records(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lid, chunks = None, []
        for line in fh:
            if line.startswith(">"):
                if lid is not None:
                    yield lid, "".join(chunks)
                lid, chunks = line[1:].split(None, 1)[0], []
            else:
                chunks.append(line.strip())
        if lid is not None:
            yield lid, "".join(chunks)


def read_catalog(denovo_dir, keep=None):
    """Forward-read part of loci in catalog.fa.gz (only IDs in `keep`, if given). None if missing.

    Stacks 2 joins a locus's paired-end contig to its forward read with a run of Ns. Loci whose contigs
    overlap have no Ns; for them the forward read is taken as the dataset's usual forward length. If
    (almost) no locus has the spacer, the data are treated as single-end and the whole locus counts.
    """
    path = os.path.join(denovo_dir, "catalog.fa.gz")
    if not os.path.isfile(path):
        return None
    seqs = {}
    starts = Counter()
    total = 0
    for lid, s in fasta_records(path):
        total += 1
        m = NRUN.search(s)
        if m:
            starts[m.start()] += 1
        if keep is None or lid in keep:
            seqs[lid] = s
    paired = total > 0 and sum(starts.values()) >= 0.05 * total
    fwd_len = starts.most_common(1)[0][0] if paired else None
    out = {"paired": paired, "fwd_len": fwd_len, "n_loci": total, "fwd": {}, "bound": {}, "length": {}}
    for lid, s in seqs.items():
        m = NRUN.search(s)
        if m:
            b, ln = m.start(), len(s) - (m.end() - m.start())
        else:
            b, ln = (min(len(s), fwd_len) if paired else len(s)), len(s)
        out["fwd"][lid] = s[:b].upper()
        out["bound"][lid] = b
        out["length"][lid] = ln
    return out


# --------------------------------------------------------------------------- fixed SNPs (Paris's count_fixed_catalog_snps.py)

def cigar_to_catalog(cigar, col):
    """Catalog column of a column in a sample locus, from sstacks's alignment of that locus to its catalog locus
    (D = catalog bases the sample lacks, I = sample bases the catalog lacks). None if the column has no catalog
    counterpart. An empty CIGAR means the two line up base for base."""
    if not cigar:
        return col
    sp = cp = 0
    for n, op in CIGAR_OP.findall(cigar):
        n = int(n)
        if op in "M=X":
            if col < sp + n:
                return cp + col - sp
            sp += n
            cp += n
        elif op in "DN":
            cp += n
        elif op in "IS":
            if col < sp + n:
                return None
            sp += n
    return None


def read_matches(path):
    """{sample locus: (catalog locus, CIGAR, [haplotypes])} from sstacks's SAMPLE.matches.tsv.gz. Haplotypes are
    written over the catalog locus's SNPs. A sample locus matched to more than one catalog locus is left out, as
    sstacks itself excludes such loci."""
    out, bad = {}, set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            try:
                c, s = int(f[0]), int(f[1])
            except ValueError:
                continue
            cg = f[4] if len(f) > 4 else ""
            if s in out and out[s][0] != c:
                bad.add(s)
            out.setdefault(s, (c, cg, []))[2].append(f[2] if len(f) > 2 else "")
    for s in bad:
        del out[s]
    return out


def sample_het_snps(path):
    """(locus, column) of every heterozygous SNP (type E) in a ustacks SAMPLE.snps.tsv.gz.

    These files list every column of every locus, so they are scanned as bytes for the E type column."""
    out = []

    def scan(buf):
        pos = buf.find(b"\tE\t")
        while pos >= 0:
            start = buf.rfind(b"\n", 0, pos) + 1
            f = buf[start:pos].split(b"\t")
            if len(f) == 2 and f[0][:1] != b"#":
                try:
                    out.append((int(f[0]), int(f[1])))
                except ValueError:
                    pass
            pos = buf.find(b"\tE\t", pos + 3)

    rest = b""
    with gzip.open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 24)
            if not chunk:
                break
            buf = rest + chunk
            cut = buf.rfind(b"\n") + 1
            scan(buf[:cut])
            rest = buf[cut:]
    if rest:
        scan(rest)
    return out


_FX = {}


def _fx_init(cat, cols):
    _FX["cat"], _FX["cols"] = cat, cols


def _fx_sample(job):
    """One sample's part of the fixed-SNP count: its heterozygous SNPs, those placed on a catalog SNP, the bases it
    carries at each catalog SNP of the loci it matched, and those catalog loci. Runs in a worker process."""
    denovo_dir, smp = job
    mp = os.path.join(denovo_dir, smp + ".matches.tsv.gz")
    sp = os.path.join(denovo_dir, smp + ".snps.tsv.gz")
    if not (os.path.isfile(mp) and os.path.isfile(sp)):
        return smp, None
    try:
        matches, snps = read_matches(mp), sample_het_snps(sp)
    except (OSError, EOFError, zlib.error, ValueError) as e:
        return smp, str(e) or e.__class__.__name__
    cat, cols = _FX["cat"], _FX["cols"]
    placed = set()
    for loc, col in snps:
        mt = matches.get(loc)
        if mt is None:
            continue
        cc = cigar_to_catalog(mt[1], col)
        if cc is not None and (mt[0], cc) in cat:
            placed.add((mt[0], cc))
    carried, matched = set(), set()
    for c, _cg, haps in matches.values():
        matched.add(c)
        cs = cols.get(c)
        if not cs:
            continue
        for h in haps:
            if len(h) == len(cs):
                for col, base in zip(cs, h):
                    if base in "ACGT":
                        carried.add((c, col, base))
    return smp, (len(snps), placed, carried, matched)


def read_fixed_snps(denovo_dir, samples, threads):
    """Paris et al.'s comparison of each individual's SNPs with the catalog (their count_fixed_catalog_snps.py),
    for Stacks 2: (i) heterozygous SNPs per individual, (ii) catalog SNPs at which at least one individual is
    heterozygous, (iii) catalog SNPs, and (iv) fixed SNPs = (iii) - (ii), catalog SNPs at which no individual is
    heterozygous. Each sample's ustacks SNPs are placed on the catalog through sstacks's alignment. None if
    catalog.snps.tsv.gz is missing.

    The fixed SNPs are also split, to read them (not used to choose; not from Paris): 'nomatch' = in a catalog
    locus that sstacks matched no sample locus back to; 'between' = the haplotypes of the individuals sstacks
    matched to the locus show two or more bases there; 'unplaced' = the rest, where those haplotypes show at most
    one base (the other allele is in sample loci sstacks did not match back, or is an N)."""
    cat_path = os.path.join(denovo_dir, "catalog.snps.tsv.gz")
    if not os.path.isfile(cat_path):
        return None
    cat = set()
    cols = {}
    with gzip.open(cat_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 3)
            if len(f) > 2 and f[2] == "E":
                try:
                    k = (int(f[0]), int(f[1]))
                except ValueError:
                    continue
                cat.add(k)
                cols.setdefault(k[0], []).append(k[1])
    for c in cols:
        cols[c].sort()

    work = [(denovo_dir, smp) for smp in samples]
    results = None
    if threads > 1 and len(work) > 1:
        try:
            with multiprocessing.Pool(min(threads, len(work)), _fx_init, (cat, cols)) as pool:
                results = pool.map(_fx_sample, work, chunksize=1)
        except (OSError, ImportError, NotImplementedError):     # no process pool here: read the files one by one
            results = None
    if results is None:
        _fx_init(cat, cols)
        results = [_fx_sample(w) for w in work]
    found, carried, matched_cat = set(), set(), set()
    ind_all, ind_cat, skipped = {}, {}, []
    for smp, res in results:
        if res is None or isinstance(res, str):
            skipped.append(smp)
            if isinstance(res, str):
                warn("{}: could not read the sstacks or ustacks file of {} ({}); it is left out of the fixed-SNP "
                     "count".format(show_path(denovo_dir), smp, res))
            continue
        n_all, placed, carr, mcat = res
        ind_all[smp] = n_all
        ind_cat[smp] = len(placed)
        found |= placed
        carried |= carr
        matched_cat |= mcat
    if not ind_cat:
        return None
    fixed = cat - found
    n_bases = Counter((c, col) for c, col, _b in carried)
    between = sum(1 for k in fixed if n_bases.get(k, 0) >= 2)
    nomatch = sum(1 for k in fixed if k[0] not in matched_cat)
    return {"cat_snps": len(cat), "variable": len(found), "fixed": len(fixed),
            "fixed_frac": share(len(fixed), len(cat)), "between": between, "nomatch": nomatch,
            "unplaced": len(fixed) - between - nomatch,
            "nomatch_loci": sum(1 for c in cols if c not in matched_cat),
            "ind_all": ind_all, "ind_cat": ind_cat, "skipped": skipped}


def parse_sumstats(popdir, het_thr):
    """From populations.sumstats.tsv: the SNP columns of every polymorphic locus, and the polymorphic loci with a
    site where at least het_thr of all genotyped samples (pooled over populations) are heterozygous, Paris et al.'s
    'heterozygous in a large majority of individuals'."""
    out = {"het_sites": None, "het_loci": None, "het_set": None, "poly_set": None, "snp_cols": None,
           "het_poly_loci": None}
    path = os.path.join(popdir, "populations.sumstats.tsv")
    if not os.path.isfile(path):
        missing(popdir, "populations.sumstats.tsv", "no forward-read SNP counts or heterozygous-site check")
        return out
    pooled_n, pooled_h = defaultdict(float), defaultdict(float)
    cols = defaultdict(set)
    iL, iC, iN, iH = 0, 3, 7, 9
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("# Locus ID"):
                    hdr = header_cols(line)
                    iL, iC = col_of(hdr, "Locus ID", 0), col_of(hdr, "Col", 3)
                    iN, iH = col_of(hdr, "N", 7), col_of(hdr, "Obs Het", 9)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(iL, iC, iN, iH):
                continue
            try:
                loc, col = f[iL], int(f[iC])
                n, ho = float(f[iN]), float(f[iH])
            except ValueError:
                continue
            key = (loc, col)
            pooled_n[key] += n
            pooled_h[key] += ho * n
            cols[loc].add(col)
    sites = [k for k in pooled_n if pooled_n[k] > 0 and pooled_h[k] / pooled_n[k] >= het_thr]
    out["het_sites"] = len(sites)
    out["het_set"] = set(k[0] for k in sites)
    out["het_loci"] = len(out["het_set"])
    out["poly_set"] = set(cols)
    out["het_poly_loci"] = len(cols)
    out["snp_cols"] = {loc: sorted(cs) for loc, cs in cols.items()}
    return out


def load_run(spec, args):
    popdir = spec["popdir"].rstrip("/") or "/"
    pop = parse_populations(popdir)
    if spec.get("denovo_dir"):
        denovo_dir = os.path.normpath(os.path.abspath(spec["denovo_dir"]))
        if not os.path.isdir(denovo_dir):
            raise ParseError("{}: denovo_dir {} does not exist".format(popdir, spec["denovo_dir"]))
    else:
        denovo_dir = resolve_denovo_dir(popdir, pop["in_path"], pop["out_path"])
    dn = parse_denovo(denovo_dir, spec)
    run = dict(pop)
    run.update(dn)
    run["popdir"] = os.path.normpath(popdir)
    run["denovo_abs"] = denovo_dir
    run["denovo_dir"] = show_path(denovo_dir)
    run["label"] = "m{} M{} n{}".format(run["m"], run["M"], run["n"])

    dist = run["snp_dist"]
    run["loci"] = sum(dist.values())
    run["poly"] = sum(v for k, v in dist.items() if k > 0)
    run["snps"] = sum(k * v for k, v in dist.items())
    run["tail"] = sum(v for k, v in dist.items() if k >= args.tail_snps)
    run["tail_frac"] = share(run["tail"], run["poly"])
    run["pct_poly"] = share(run["poly"], run["loci"])
    run["snps_per_poly"] = share(run["snps"], run["poly"])

    # Coverage as in Paris Fig. 1: each sample's ustacks coverage (after merging alleles).
    covs = list(run["ustacks_cov"].values())
    run["cov_mean"] = statistics.mean(covs) if covs else float("nan")
    run["cov_median"] = statistics.median(covs) if covs else float("nan")
    run["cov_min"] = min(covs) if covs else float("nan")
    run["low_cov"] = sorted((s, c) for s, c in run["ustacks_cov"].items() if c < LOW_SAMPLE_COV)

    if run["loci"] != run["kept_loci"]:
        warn("{}: SNPs-per-locus distribution sums to {:,} loci but populations.log says {:,} were kept".format(
            popdir, run["loci"], run["kept_loci"]))
    if run["snps"] != run["kept_variant_sites"]:
        warn("{}: SNPs-per-locus distribution gives {:,} SNPs but populations.log says {:,} variant sites "
             "remained".format(popdir, run["snps"], run["kept_variant_sites"]))

    run.update(parse_sumstats(popdir, args.het_site))
    run["het_frac"] = share(run["het_loci"], run["poly"]) if run["het_loci"] is not None else None
    if run["het_poly_loci"] is not None and run["het_poly_loci"] != run["poly"]:
        warn("{}: populations.sumstats.tsv has variant sites in {:,} loci but the distribs file has {:,} "
             "polymorphic loci".format(popdir, run["het_poly_loci"], run["poly"]))

    run["r80_ids"], r80_samples = read_haplotypes(popdir)
    if run["r80_ids"] is None:
        missing(popdir, "populations.haplotypes.tsv", "the samples are taken from denovo_map.log")
    elif len(run["r80_ids"]) != run["loci"]:
        warn("{}: populations.haplotypes.tsv lists {:,} loci but the distribs file has {:,}".format(
            popdir, len(run["r80_ids"]), run["loci"]))

    # Paris et al.'s within-individual metrics, over the samples populations used (else every ustacks sample).
    samples = sorted(r80_samples) if r80_samples else sorted(run["ustacks"])
    if not run["ustacks"]:
        missing(denovo_dir, "per-sample ustacks lines in denovo_map.log", "no loci per individual, allele "
                "collapse or coverage")
    hets = read_ind_het(denovo_dir, samples)
    if not hets:
        missing(denovo_dir, "per-sample ustacks files (SAMPLE.alleles.tsv.gz)", "no heterozygous loci per "
                "individual")
    elif len(hets) < len(samples):
        warn("{}: {} of {} samples have no SAMPLE.alleles.tsv.gz; heterozygous loci per individual use the "
             "other {}".format(show_path(denovo_dir), len(samples) - len(hets), len(samples), len(hets)))
    run["fixed"] = read_fixed_snps(denovo_dir, samples, args.threads)
    if run["fixed"] is None:
        if not os.path.isfile(os.path.join(denovo_dir, "catalog.snps.tsv.gz")):
            missing(denovo_dir, "catalog.snps.tsv.gz", "no fixed-SNP count")
        else:
            missing(denovo_dir, "per-sample SAMPLE.snps.tsv.gz and SAMPLE.matches.tsv.gz", "no fixed-SNP count")
    elif run["fixed"]["skipped"]:
        warn("{}: {} of {} samples have no SAMPLE.snps.tsv.gz or SAMPLE.matches.tsv.gz; the fixed-SNP count uses "
             "the other {}".format(show_path(denovo_dir), len(run["fixed"]["skipped"]), len(samples),
                                   len(samples) - len(run["fixed"]["skipped"])))
    fx = run["fixed"] or {"ind_all": {}, "ind_cat": {}}
    run["ind"] = {}
    for smp in samples:
        u = run["ustacks"].get(smp, {})
        run["ind"][smp] = {"stacks": u.get("stacks"), "loci_M": u.get("loci_M"), "final": u.get("final"),
                           "het": hets.get(smp), "het_snps": fx["ind_all"].get(smp),
                           "het_snps_cat": fx["ind_cat"].get(smp)}
    for key, fn in IND_METRICS:
        vals = dict((smp, fn(d)) for smp, d in run["ind"].items())
        run[key] = dict((smp, v) for smp, v in vals.items() if v is not None)
        run[key + "_mean"] = statistics.mean(run[key].values()) if run[key] else None
    cat_snps = [d["het_snps_cat"] for d in run["ind"].values() if d["het_snps_cat"] is not None]
    run["ind_snps_cat_mean"] = statistics.mean(cat_snps) if cat_snps else None

    run["catalog"] = read_catalog(denovo_dir, run["poly_set"] if run["poly_set"] is not None else set())
    if run["catalog"] is None:
        missing(denovo_dir, "catalog.fa.gz", "SNP-dense loci are counted over whole loci only")
    run["fwd_snps"] = run["pe_snps"] = None
    if run["catalog"] is not None and run["snp_cols"] is not None:
        bound = run["catalog"]["bound"]
        fs, ps, miss = {}, {}, 0
        for loc, cs in run["snp_cols"].items():
            b = bound.get(loc)
            if b is None:
                miss += 1
                continue
            nf = sum(1 for c in cs if c < b)
            fs[loc], ps[loc] = nf, len(cs) - nf
        if miss:
            warn("{}: {:,} polymorphic loci are missing from catalog.fa.gz".format(popdir, miss))
        run["fwd_snps"], run["pe_snps"] = fs, ps
    return run


def read_runs_table(path):
    """Rows of a --runs table: populations_dir plus optional denovo_dir, m, M, n."""
    base = os.path.dirname(os.path.abspath(path))
    lines = [ln for ln in read_text(path).splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise ParseError("{} is empty".format(path))
    hdr = [h.strip() for h in lines[0].split("\t")]
    if "populations_dir" not in hdr:
        raise ParseError("{}: the first line must be a tab-separated header with a populations_dir column "
                         "(optional: denovo_dir, m, M, n)".format(path))
    specs = []
    for ln in lines[1:]:
        f = dict(zip(hdr, [x.strip() for x in ln.split("\t")]))
        spec = {}
        for key in ("populations_dir", "denovo_dir"):
            v = f.get(key, "")
            if v:
                spec["popdir" if key == "populations_dir" else key] = v if os.path.isabs(v) else os.path.join(base, v)
        if "popdir" not in spec:
            raise ParseError("{}: a row has no populations_dir".format(path))
        for key in ("m", "M", "n"):
            v = f.get(key, "")
            if v:
                try:
                    spec[key] = int(v)
                except ValueError:
                    raise ParseError("{}: {} = '{}' is not a whole number".format(path, key, v))
        specs.append(spec)
    return specs


def tail_settings(runs, args):
    """Thresholds for the SNP-dense tail: whole locus, forward read and paired-end part."""
    cats = [r["catalog"] for r in runs if r["catalog"] is not None]
    paired = bool(cats) and any(c["paired"] for c in cats)
    k = args.tail_snps
    info = {"paired": paired, "k": k, "k_fwd": k, "k_pe": None, "fwd_len": None, "med_len": None}
    if not paired:
        return info
    fl = Counter(c["fwd_len"] for c in cats if c["paired"]).most_common(1)[0][0]
    lens = []
    for r in runs:
        if r["catalog"] is not None and r["poly_set"] is not None:
            ln = r["catalog"]["length"]
            lens.extend(ln[l] for l in r["poly_set"] if l in ln)
    med = statistics.median(lens) if lens else None
    info["fwd_len"], info["med_len"] = fl, med
    if args.fwd_tail_snps:
        info["k_fwd"] = args.fwd_tail_snps
    elif med:
        info["k_fwd"] = max(2, int(round(k * fl / med)))
    if med and med > fl:
        info["k_pe"] = max(2, int(round(k * (med - fl) / med)))
    return info


def apply_tails(runs, ti):
    for r in runs:
        r["fwd_tail"] = r["pe_tail"] = None
        r["fwd_tail_frac"] = r["pe_tail_frac"] = None
        if r["fwd_snps"] is not None:
            r["fwd_tail"] = sum(1 for v in r["fwd_snps"].values() if v >= ti["k_fwd"])
            r["fwd_tail_frac"] = share(r["fwd_tail"], r["poly"])
            if ti["k_pe"]:
                r["pe_tail"] = sum(1 for v in r["pe_snps"].values() if v >= ti["k_pe"])
                r["pe_tail_frac"] = share(r["pe_tail"], r["poly"])
        if ti["paired"]:
            r["rule_tail_frac"] = r["fwd_tail_frac"]
        else:
            r["rule_tail_frac"] = r["tail_frac"]


# --------------------------------------------------------------------------- series

VAR_NAME = {"M=n": "M = n", "M": "M", "n": "n", "m": "m"}
KIND_ORDER = ["M=n", "M", "n", "m"]


def value_label(kind, run):
    if kind == "M=n":
        return "M{}n{}".format(run["M"], run["n"])
    return "{}{}".format(kind, run[kind])


def find_series(runs):
    groups = defaultdict(list)
    for r in runs:
        groups[(r["R"], r["n_samples"])].append(r)
    out = []
    for (R, ns), rs in sorted(groups.items()):
        specs = [
            ("M=n", lambda r: r["M"] == r["n"], lambda r: (r["m"],), "M", lambda k: "m{} fixed".format(*k)),
            ("M", lambda r: True, lambda r: (r["m"], r["n"]), "M", lambda k: "m{} n{} fixed".format(*k)),
            ("n", lambda r: True, lambda r: (r["m"], r["M"]), "n", lambda k: "m{} M{} fixed".format(*k)),
            ("m", lambda r: True, lambda r: (r["M"], r["n"]), "m", lambda k: "M{} n{} fixed".format(*k)),
        ]
        for kind, keep, key, var, fixed_txt in specs:
            buckets = defaultdict(list)
            for r in rs:
                if keep(r):
                    buckets[key(r)].append(r)
            for k, lst in sorted(buckets.items()):
                by_x = {}
                for r in lst:
                    x = r[var]
                    if x in by_x:
                        warn("{} and {} have the same parameters ({}); using {}".format(
                            by_x[x]["popdir"], r["popdir"], r["label"], by_x[x]["popdir"]))
                        continue
                    by_x[x] = r
                # A series needs two values; an M series with n fixed needs three, so that the n = M - 1 and M + 1
                # runs made for Paris's n rule do not form two-point M series of their own.
                if len(by_x) < (3 if kind == "M" else 2):
                    continue
                srt = [by_x[x] for x in sorted(by_x)]
                out.append({
                    "kind": kind, "var": var, "fixed": dict(zip(["m", "n"] if kind == "M" else
                                                                ["m", "M"] if kind == "n" else
                                                                ["M", "n"] if kind == "m" else ["m"], k)),
                    "R": R, "n_samples": ns, "runs": srt,
                    "xs": [r[var] for r in srt], "labels": [value_label(kind, r) for r in srt],
                    "title": "{} series ({})".format(VAR_NAME[kind], fixed_txt(k)),
                    "filter": "populations -R {:g}: loci in ≥{:.0%} of all {} samples".format(R, R, ns),
                })
    out.sort(key=lambda s: (s["R"], s["n_samples"], KIND_ORDER.index(s["kind"])))
    return out


def snp_classes(dist, top, include0=False):
    """Counts per SNP class: 1..top (or 0..top) plus one pooled class for > top."""
    ks = list(range(0 if include0 else 1, top + 1))
    counts = [dist.get(k, 0) for k in ks] + [sum(v for k, v in dist.items() if k > top)]
    return [str(k) for k in ks] + ["{}+".format(top + 1)], counts


def shares_of(counts):
    total = sum(counts)
    return [c / total if total else float("nan") for c in counts]


# --------------------------------------------------------------------------- analysis

M_KINDS = ("M=n", "M")


def paired_step(a, b):
    """Mean change b - a over the samples present in both runs, with its standard error (paired by individual)."""
    common = [smp for smp in a if smp in b]
    if len(common) < 3:
        return None
    d = [b[smp] - a[smp] for smp in common]
    m = statistics.mean(d)
    se = statistics.stdev(d) / math.sqrt(len(d))
    return {"d": m, "se": se, "n": len(d), "up": sum(1 for x in d if x > 0), "rise": m > 2 * se,
            "fall": m < -2 * se}


def first_flat(steps):
    """Index of the first value whose next step does not rise by more than 2 SE (None if every step does)."""
    for k, st in enumerate(steps):
        if not st["rise"]:
            return k
    return None


def pm(st, digits, scale=1.0):
    """A paired step as '+67 ± 32*' (* = beyond two standard errors)."""
    return "{:+.{d}f} ± {:.{d}f}{}".format(st["d"] * scale, 2 * st["se"] * scale,
                                          "*" if st["rise"] or st["fall"] else "", d=digits)


def sig_change(v0, n0, v1, n1):
    if shares_differ(v0, n0, v1, n1):
        return "more than two standard errors {}".format("higher" if v1 > v0 else "lower")
    return "within two standard errors"


def by_step(lab, items):
    """'a at X→Y, b at Y→Z' for a list of per-step texts."""
    return join_list("{} at {}→{}".format(t, lab[j], lab[j + 1]) for j, t in items)


def analyse(s, args, ti):
    runs, xs, lab = s["runs"], s["xs"], s["labels"]
    R, var, kind = s["R"], s["var"], s["kind"]
    P = [r["poly"] for r in runs]
    polys = P
    s["class_labels"] = snp_classes({}, args.dist_max_snps)[0]
    s["shares"] = [shares_of(snp_classes(r["snp_dist"], args.dist_max_snps)[1]) for r in runs]
    steps = []
    for i in range(1, len(runs)):
        d = P[i] - P[i - 1]
        steps.append({"from": i - 1, "to": i, "d": d, "rel": d / P[i - 1] if P[i - 1] else float("nan")})
    s["P"], s["steps"] = P, steps
    last = len(xs) - 1
    notes = []

    # Within-individual metrics: means over individuals, and steps paired by individual.
    ind = {}
    for key, _ in IND_METRICS:
        vals = [r.get(key) or {} for r in runs]
        st = [paired_step(vals[i - 1], vals[i]) for i in range(1, len(runs))] if all(vals) else None
        if st is not None and any(x is None for x in st):
            st = None
        ind[key] = {"ok": st is not None, "means": [r.get(key + "_mean") for r in runs], "steps": st}
    s["ind"] = ind

    # 1. r80 rule: the most polymorphic r80 loci (for an n series, only n = M - 1, M, M + 1).
    cand = list(range(len(xs)))
    if kind == "n":
        M = s["fixed"]["M"]
        allowed = [i for i, x in enumerate(xs) if abs(x - M) <= 1]
        untested = [v for v in (M - 1, M, M + 1) if v >= 0 and v not in xs]
        if untested:
            notes.append("Paris et al.'s n rule compares n = M - 1, M and M + 1. Not tested here: {}.".format(
                ", ".join("n{}".format(v) for v in untested)))
        if allowed:
            cand = allowed
            free = max(range(len(xs)), key=lambda i: (P[i], -xs[i]))
            if free not in allowed:
                notes.append("{} has more polymorphic r80 loci but lies outside n = M ± 1, so it is not "
                             "picked.".format(lab[free]))
    best = max(cand, key=lambda i: (P[i], -xs[i]))
    ties = [lab[i] for i in cand if P[i] == P[best] and i != best]
    if ties:
        notes.append("{} ties with {} on polymorphic r80 loci; the lower value is taken.".format(
            lab[best], ", ".join(ties)))
    falling = bool(steps) and all(st["d"] < 0 for st in steps)
    rising = bool(steps) and all(st["d"] > 0 for st in steps)
    rising_top = bool(steps) and best == last and steps[-1]["d"] > 0

    # 2. Paris Fig. 4: the first value after which no step adds --plateau or more new polymorphic r80 loci.
    fig4 = None
    for k in range(len(xs) - 1):
        if all(steps[j]["rel"] < args.plateau for j in range(k, len(steps))):
            fig4 = k
            break

    # 3. Within individuals. Paris et al. chose where polymorphic loci within individuals 'began to flatten out'
    # (trout) or 'began to plateau' (earthworm). Here: heterozygous loci per individual (ustacks), and the first
    # value after which the next step does not rise beyond 2 SE (paired by individual). Only for M: ustacks runs
    # before the catalog, so its numbers do not change with n.
    wpl = {}
    if ind["ind_het"]["ok"] and kind in M_KINDS:
        wpl["ind_het"] = first_flat(ind["ind_het"]["steps"])
    if not wpl:
        W, w_status = None, "missing"
    elif wpl["ind_het"] is None:
        W, w_status = None, "rising"
    else:
        W, w_status = wpl["ind_het"], "ok"

    # 4. The pick. For M: the higher of the r80 maximum and where polymorphic loci within individuals first level
    # off; if polymorphic r80 loci are still rising at the top of the range, the within-individual level-off point
    # (Fig. 4 only if the within-individual metric is missing).
    if kind in M_KINDS:
        if rising_top:
            if w_status == "ok":
                pick, basis = W, "earthworm"
            elif w_status == "rising":
                pick, basis = None, "rising"
            elif fig4 is not None:
                pick, basis = fig4, "fig4"
            else:
                pick, basis = None, "rising_nofig4"
        elif w_status == "ok":
            pick, basis = max(best, W), ("below" if best < W else "max")
        elif w_status == "rising":
            pick, basis = None, "w_rising"
        else:
            pick, basis = best, "r80_only"
    else:
        pick, basis = best, kind + "_rule"
    ref = pick if pick is not None else best

    # ---- notes
    if kind in M_KINDS:
        if falling:
            notes.append("Polymorphic r80 loci fall at every step, from {:,} at {} to {:,} at {}, so no step adds new "
                         "polymorphic r80 loci (Paris Fig. 4). None of Paris et al.'s datasets did this: polymorphic "
                         "r80 loci were highest at M2 (penguin) and M5 (trout), and still rising at M11 (earthworm)."
                         .format(P[0], lab[0], P[-1], lab[-1]))
        elif rising_top:
            notes.append("Polymorphic r80 loci are still rising at the highest value tested{}. Paris et al. saw this "
                         "in their earthworm data, up to M11, and put it down to a small number of loci with many "
                         "SNPs.".format(" (they rise at every step)" if rising else ""))
        if not falling:
            if fig4 is not None:
                notes.append("Paris Fig. 4: from {} on, no step adds {:.1%} or more new polymorphic r80 loci.".format(
                    lab[fig4], args.plateau))
            else:
                notes.append("Paris Fig. 4: the last step still adds {:.1%} or more new polymorphic r80 loci, so the "
                             "number of new loci has not levelled off in the tested range.".format(args.plateau))

        het = ind["ind_het"]
        if het["ok"]:
            st = het["steps"]
            k = wpl.get("ind_het")
            txt = ("Polymorphic loci within individuals: Paris et al. chose M5 for trout, where these 'began to "
                   "flatten out', and M3 for earthworm, where they 'began to plateau'. Here they are heterozygous loci "
                   "per individual, counted in each sample's ustacks output (before the catalog). They change by {} "
                   "(± 2 SE, paired over {} individuals).".format(
                       by_step(lab, [(j, "{} ({} of {} up)".format(pm(x, 0), x["up"], x["n"]))
                                     for j, x in enumerate(st)]), st[0]["n"]))
            txt += (" They first level off (the next step does not rise beyond two standard errors) after {}."
                    .format(lab[k]) if k is not None else " They rise beyond two standard errors at every step.")
            notes.append(txt)

        col = ind["ind_collapse"]
        if col["ok"]:
            cm = col["means"]
            txt = ("Allele collapse within individuals (Paris Table S6a, 'the proportion of putative alleles that "
                   "collapsed'): the share of each individual's stacks that ustacks merged into loci is {}.".format(
                       join_list("{:.1%} at {}".format(v, lab[i]) for i, v in enumerate(cm))))
            tot = cm[-1] - cm[0]
            if tot > 0 and len(xs) > 2:
                fr = ["{:.0%} by {}".format((cm[i] - cm[0]) / tot, lab[i]) for i in range(1, last)]
                txt += " Of the change across the tested range, {} has happened by {}{}.".format(
                    fr[0].split(" by ")[0], fr[0].split(" by ")[1],
                    ", " + join_list(fr[1:]) if len(fr) > 1 else "")
            txt += (" Paris et al. saw this level out between M3 and M4 (trout), M1 and M2 (penguin) and M2 and M3 "
                    "(earthworm), and chose M5, M2 and M3.")
            notes.append(txt)

        if basis == "below" and ind["ind_loci"]["ok"] and het["ok"]:
            dl = ind["ind_loci"]["steps"][best]["d"]
            dh = het["steps"][best]["d"]
            if dl < 0 < dh:
                notes.append("From {} to {}, each individual loses {:,.0f} loci on average but gains only {:.0f} "
                             "heterozygous loci, so most of the stacks merged at that step do not form a heterozygous "
                             "locus. Raising M is removing loci rather than turning monomorphic loci polymorphic, "
                             "which is how polymorphic r80 loci can fall while polymorphism within individuals still "
                             "rises.".format(lab[best], lab[best + 1], -dl, dh))

        if pick is not None and het["ok"]:
            later = [j for j in range(pick, len(steps)) if het["steps"][j]["rise"]]
            if later:
                ucov = [r["cov_mean"] for r in runs if r["cov_mean"] == r["cov_mean"]]
                notes.append("Above {}, heterozygous loci per individual rise again ({}), so the level-off is the "
                             "first one, not a lasting plateau. Paris et al.'s wording, 'began to flatten out' (trout) "
                             "and 'began to plateau' (earthworm), refers to where the levelling starts. They found that "
                             "higher M adds 'a small number of loci that contained many SNPs' at each increment: "
                             "'With high coverage, these are likely to be true loci', and with low coverage 'more "
                             "likely to be erroneously assembled loci'.{}".format(
                                 lab[pick], by_step(lab, [(j, pm(het["steps"][j], 0)) for j in later]),
                                 " Mean ustacks coverage (after merging alleles) here is {}–{}; Paris et al.'s datasets averaged "
                                 "21x (earthworm), 28x (trout) and 49x (penguin) after merging alleles at m3."
                                 .format(fmt_cov(min(ucov)), fmt_cov(max(ucov))) if ucov else ""))

        S = [r["snps"] for r in runs]
        notes.append("SNPs in r80 loci (Paris Fig. 3c) go from {:,} at {} to {:,} at {} ({}). Paris et al. saw a "
                     "steep rise in SNPs with M in trout and earthworm and 'a much smaller increase' in penguin. In "
                     "trout it 'reflects the higher value of M required'; in earthworm they put the continued rise in "
                     "r80 polymorphism down to 'a small number of loci with a high-density of SNPs'; in penguin 'the "
                     "majority of polymorphism and SNPs were already captured' at M1.".format(
                         S[0], lab[0], S[-1], lab[-1], pct((S[-1] - S[0]) / S[0]) if S[0] else "n/a"))

    # Over-merging signs Paris et al. name: SNP-dense loci and loci heterozygous in most individuals.
    paired = ti["paired"] and all(r["fwd_tail_frac"] is not None for r in runs)
    tvals = [r["rule_tail_frac"] if paired else r["tail_frac"] for r in runs]
    tail_what = ("≥{} SNPs in the forward read".format(ti["k_fwd"]) if paired else "≥{} SNPs".format(args.tail_snps))
    hvals = [r["het_frac"] for r in runs]
    het_ok = all(v is not None for v in hvals)
    if ref < last:
        parts = []
        if all(v is not None and v == v for v in tvals):
            parts.append("loci with {} go from {:.2%} at {} to {:.2%} at {} ({})".format(
                tail_what, tvals[ref], lab[ref], tvals[-1], lab[-1],
                sig_change(tvals[ref], polys[ref], tvals[-1], polys[-1])))
        if het_ok:
            parts.append("loci with a site where at least {:.0%} of all samples are heterozygous go from {:.2%} at {} "
                         "to {:.2%} at {} ({})".format(args.het_site, hvals[ref], lab[ref], hvals[-1], lab[-1],
                                                       sig_change(hvals[ref], polys[ref], hvals[-1], polys[-1])))
        if parts:
            txt = ("Over-merging signs Paris et al. name, as shares of polymorphic r80 loci from {} upwards: {}."
                   .format("the pick" if pick is not None else "the r80 maximum", "; ".join(parts)))
            if paired:
                txt += (" With paired-end data only the forward read is assembled under M and n, so only its SNPs are "
                        "counted; ≥{} there is the same density as ≥{} in a {:.0f} bp locus.".format(
                            ti["k_fwd"], args.tail_snps, ti["med_len"]))
            txt += (" Paris et al. note that SNP-dense loci are likely real at high coverage and likely misassembled "
                    "at low coverage. Over-merged paralogs 'will appear as heterozygous in a large majority of "
                    "individuals'; Paris et al. filter them after assembly with populations --max_obs_het, and their "
                    "Table 2 suggests a lower M for a repetitive or polyploid genome. In Stacks 2, --max-obs-het "
                    "removes SNPs, not loci; to drop whole loci, blacklist them. The thresholds used here are this "
                    "script's, not Paris et al.'s.")
            notes.append(txt)

    # Edge of the tested range.
    if pick is not None and len(xs) > 1 and kind != "n":
        if pick == last:
            notes.append("The pick ({}) is the highest value tested; test higher values to confirm it.".format(
                lab[pick]))
        elif pick == 0:
            if var == "M" and xs[0] == 1:
                notes.append("The pick ({}) is the lowest value tested. The only value below it is M0, which Paris "
                             "et al. reject because alleles cannot merge into loci.".format(lab[pick]))
            elif xs[0] > 0:
                notes.append("The pick ({}) is the lowest value tested; test {}{} to confirm it.".format(
                    lab[pick], "M{0}n".format(xs[0] - 1) if kind == "M=n" else var, xs[0] - 1))
    if kind in M_KINDS and xs[-1] < 8:
        notes.append("Paris et al. suggest testing M1 to M8 (Table 2); this series stops at {}.".format(lab[-1]))
    fxs = [r["fixed"] for r in runs]
    if kind == "n" and all(fxs):
        notes.append("Fixed SNPs in the catalog (Paris Fig. S3): {}. Paris et al. found that fixed SNPs rose with n in "
                     "all three datasets without levelling off, and that 'This alone, however, does not provide "
                     "enough information to choose a value for n', so they are reported, not used for the pick."
                     .format(join_list("{:,} ({:.1%}) at {}".format(f["fixed"], f["fixed_frac"], lab[i])
                                       for i, f in enumerate(fxs))))
    if kind == "m":
        notes.append("Paris et al. suggest testing m3 to m7 (Table 2), do not suggest m above 5, advise disabling "
                     "secondary reads above m6, and chose m3 for all three of their datasets.")

    # ---- verdict
    def wtxt(key):
        if key not in wpl:
            return "not available"
        k = wpl[key]
        return lab[k] if k is not None else "does not level off in the tested range"

    head = ("Pick following Paris et al.: {}.".format(runs[pick]["label"]) if pick is not None
            else "No pick following Paris et al. in the tested range.")
    verdict = [head]
    if kind in M_KINDS:
        verdict.append("  r80 rule (most polymorphic r80 loci): {} ({:,})".format(lab[best], P[best]))
        verdict.append("  Fig. 4 (later steps add <{:.1%} new polymorphic r80 loci): {}".format(
            args.plateau, lab[fig4] if fig4 is not None else "not reached"))
        verdict.append("  Polymorphic loci within individuals (heterozygous loci per individual) first level off: "
                       "{}".format(wtxt("ind_het")))
        how = {
            "max": "The r80 maximum ({}) is at or above where polymorphic loci within individuals first level off "
                   "({}), so the r80 rule decides, as in Paris et al.'s trout (M5) and penguin (M2) data."
                   .format(lab[best], lab[W] if W is not None else ""),
            "below": "The r80 maximum ({}) lies below where polymorphic loci within individuals first level off "
                     "({}). None of Paris et al.'s datasets did this; in all three they chose M at or above where "
                     "polymorphic loci within individuals levelled off (trout M5, earthworm M3, penguin M2), so the "
                     "pick is {}. This applies their judgement to a case they did not meet; it is not a rule they "
                     "state. Taken alone, their r80 rule gives {}."
                     .format(lab[best], lab[W] if W is not None else "", lab[pick] if pick is not None else "",
                             lab[best]),
            "earthworm": "Polymorphic r80 loci are still rising at the highest value, as in Paris et al.'s earthworm "
                         "data, where they took M3: where polymorphic loci within individuals began to plateau and "
                         "Fig. 4 showed few new polymorphic r80 loci. The same reading here gives {}."
                         .format(lab[pick] if pick is not None else ""),
            "fig4": "Polymorphic r80 loci are still rising at the highest value and heterozygous loci per individual "
                    "are not available, so the pick follows Fig. 4: {}.".format(lab[pick] if pick is not None else ""),
            "rising": "Polymorphic r80 loci and polymorphic loci within individuals are both still rising at the "
                      "highest value tested. Test higher values (Paris et al. suggest M1 to M8).",
            "rising_nofig4": "Polymorphic r80 loci are still rising at the highest value tested, heterozygous loci per "
                             "individual are not available, and Fig. 4 has not levelled off. Test higher values.",
            "w_rising": "Polymorphic loci within individuals rise beyond two standard errors at every step, so the "
                        "tested range does not reach the point at or above which Paris et al. chose. Test higher "
                        "values. The r80 maximum is {}.".format(lab[best]),
            "r80_only": "Heterozygous loci per individual are not available (see Warnings), so only the r80 rule "
                        "applies.",
        }[basis]
        verdict.append("  How: " + how)
    elif kind == "n":
        verdict.append("  n rule (most polymorphic r80 loci among n = M - 1, M, M + 1): {} ({:,})".format(
            lab[best], P[best]))
    else:
        verdict.append("  r80 rule (most polymorphic r80 loci): {} ({:,})".format(lab[best], P[best]))

    s.update({"best": best, "fig4": fig4, "wpl": wpl, "W": W, "w_status": w_status, "pick": pick, "basis": basis,
              "tail_vals": tvals, "tail_what": tail_what, "tail_paired": paired,
              "verdict": [relabel(v, R) for v in verdict], "notes": [relabel(n, R) for n in notes]})
    return s


def recommend(series, runs):
    """One Paris recommendation per group of comparable runs (same -R and sample count): m, then M, then n.
    Each entry is ((R, n_samples), report lines, label, brief lines); the brief lines are the short
    screen version."""
    out = []
    groups = defaultdict(list)
    for s in series:
        groups[(s["R"], s["n_samples"])].append(s)
    for (R, ns), ss in sorted(groups.items()):
        grp = [r for r in runs if r["R"] == R and r["n_samples"] == ns]
        lines = []
        brief = []
        mser = [s for s in ss if s["kind"] == "m" and s["pick"] is not None]
        m_vals = sorted(set(r["m"] for r in grp))
        if mser:
            picks = sorted(set(s["runs"][s["pick"]]["m"] for s in mser))
            m = picks[0]
            m_why = "{} (r80 rule over {}){}".format(
                m, join_list("m{}".format(x) for x in sorted(set(x for s in mser for x in s["xs"]))),
                "" if len(picks) == 1 else "; the m series disagree ({}), so the lowest is taken".format(
                    ", ".join("m{}".format(p) for p in picks)))
            lines.append("m: " + m_why)
            brief.append("m: " + m_why)
        else:
            m = 3 if 3 in m_vals else m_vals[0]
            m_why = "the only value tested" if len(m_vals) == 1 else "no m series"
            lines.append("m: {} ({}). Paris et al. suggest testing m3 to m7 and chose m3 for all three of their "
                         "datasets.".format(m, m_why))
            brief.append("m: {} ({})".format(m, m_why))
        Mser = [s for s in ss if s["kind"] in M_KINDS and s["fixed"]["m"] == m]
        Mser.sort(key=lambda s: (s["kind"] != "M=n", -len(s["xs"])))
        if not Mser:
            lines.append("M: no M or M = n series at m{}, so no recommendation for M.".format(m))
            brief.append("M: no recommendation (no M or M = n series at m{})".format(m))
            out.append(((R, ns), lines, None, brief))
            continue
        sM = Mser[0]
        if sM["pick"] is None:
            lines.append("M: no recommendation in the tested range ({}); see its verdict.".format(sM["title"]))
            brief.append("M: no recommendation in the tested range ({}); see the full report.".format(sM["title"]))
            out.append(((R, ns), lines, None, brief))
            continue
        pr = sM["runs"][sM["pick"]]
        M = pr["M"]
        lab = sM["labels"]
        best_l, w_l = lab[sM["best"]], (lab[sM["W"]] if sM["W"] is not None else "")
        why = {
            "below": "the r80 maximum ({0}) lies below where polymorphic loci within individuals (heterozygous loci "
                     "per individual) first level off ({1}). Paris et al. never chose M below that point (trout M5, "
                     "earthworm M3, penguin M2), so the higher of the two is taken. None of their datasets had the "
                     "r80 maximum below it, so this applies their judgement to a case they did not meet; taken "
                     "alone, their r80 rule gives {0}.".format(best_l, w_l),
            "max": "it has the most polymorphic r80 loci (Paris et al.'s r80 rule) and is at or above where "
                   "polymorphic loci within individuals first level off ({}).".format(w_l),
            "earthworm": "polymorphic r80 loci are still rising at the highest value tested, so, as for Paris et al.'s "
                         "earthworm data, the pick is where polymorphic loci within individuals first level off.",
            "fig4": "polymorphic r80 loci are still rising and heterozygous loci per individual are missing, so Paris "
                    "Fig. 4 decides.",
            "r80_only": "it has the most polymorphic r80 loci; heterozygous loci per individual are missing.",
        }.get(sM["basis"], "see the series verdict.")
        lines.append("M: {}, from the {}, because {}".format(M, sM["title"], why))
        short = {
            "below": "The r80 maximum ({0}) is below where heterozygous loci per individual first level off ({1}), "
                     "so the higher is taken; Paris's r80 rule alone gives {0}.".format(best_l, w_l),
            "max": "It has the most polymorphic r80 loci and is at or above where heterozygous loci per individual "
                   "first level off ({}).".format(w_l),
            "earthworm": "Polymorphic r80 loci are still rising at the highest value, so, as for Paris's earthworm "
                         "data, the pick is where heterozygous loci per individual first level off.",
            "fig4": "Polymorphic r80 loci are still rising and heterozygous loci per individual are missing, so "
                    "Paris Fig. 4 decides.",
            "r80_only": "It has the most polymorphic r80 loci (heterozygous loci per individual are missing).",
        }.get(sM["basis"], "See the series verdict in the full report.")
        series_name = "M = n series" if sM["kind"] == "M=n" else "M series"
        brief.append("M: {}, from the {}. {}".format(M, series_name, short))
        if sM["kind"] == "M=n":
            lines.append("That series changes M and n together. Paris et al. varied M with n fixed (n0 in their "
                         "primary runs), so the r80 counts here reflect n as well as M.")
        nser = [s for s in ss if s["kind"] == "n" and s["fixed"] == {"m": m, "M": M} and s["pick"] is not None]
        have = sorted(set(r["n"] for r in grp if r["m"] == m and r["M"] == M))
        want = [v for v in (M - 1, M, M + 1) if v >= 0]
        todo = [v for v in want if v not in have]
        if nser:
            n = nser[0]["runs"][nser[0]["pick"]]["n"]
            lines.append("n: {} (Paris et al.'s n rule over {})".format(
                n, ", ".join("n{}".format(v) for v in have if v in want)))
            n_brief = "n: {} (Paris's n rule over {})".format(
                n, join_list("n{}".format(v) for v in have if v in want))
        else:
            n = pr["n"]
            lines.append("n: {} for now (from the {}), until Paris et al.'s n rule is run.".format(
                n, "M = n series" if sM["kind"] == "M=n" else "M series"))
            n_brief = "n: {} for now.".format(n)
        if todo:
            lines.append("Paris et al.'s n rule compares n = M - 1, M and M + 1 at the chosen M and keeps the one "
                         "with the most polymorphic r80 loci. Still to run: {}, each followed by populations -R {:g} "
                         "with the same options as the other runs (or with -X \"populations:-R {:g}\" given to "
                         "denovo_map.pl, without -r).".format(
                             join_list("m{} M{} n{}".format(m, M, v) for v in todo), R, R))
            n_brief += " Still to run for Paris's n rule: {} (then populations -R {:g}).".format(
                join_list("m{} M{} n{}".format(m, M, v) for v in todo), R)
        brief.append(n_brief)
        out.append(((R, ns), lines, "m{} M{} n{}{}".format(m, M, n, " (n provisional)" if todo else ""), brief))
    return out


# --------------------------------------------------------------------------- report

def fmt_cov(x):
    return "{:.1f}x".format(x) if x == x else "n/a"


def count_share(k, p):
    if k is None or p is None or p != p:
        return "n/a"
    return "{:,} ({:.2%})".format(k, p)


def wrap_lines(lines, b, bullet="  - "):
    for nt in lines:
        b(textwrap.fill(nt, width=78, initial_indent=bullet, subsequent_indent=" " * len(bullet)))


def build_report(runs, series, recs, args, ti, footer=None):
    out = []
    w = out.append
    w("Stacks parameter selection following Paris et al. (2017)")
    w("=" * 57)
    w("Runs read: {}".format(len(runs)))
    w("")
    dw = max([len("populations dir")] + [len(show_path(r["popdir"])) for r in runs])
    w("{:<{dw}} {:<12} {:>8} {:>8} {:>8} {:>9} {:>8}".format(
        "populations dir", "params", "r80 loci", "poly", "SNPs", "mean cov", "min cov", dw=dw))
    for r in sorted(runs, key=lambda r: (r["m"], r["M"], r["n"])):
        w(relabel("{:<{dw}} {:<12} {:>8,} {:>8,} {:>8,} {:>9} {:>8}".format(
            show_path(r["popdir"]), r["label"], r["loci"], r["poly"], r["snps"], fmt_cov(r["cov_mean"]),
            fmt_cov(r["cov_min"]), dw=dw), r["R"]))
    w(textwrap.fill("poly = polymorphic loci. cov = each sample's ustacks coverage after merging alleles (Paris "
                    "Fig. 1): the mean and the lowest over samples.", width=78))
    inside = sorted(show_path(r["popdir"]) for r in runs
                    if os.path.realpath(r["popdir"]) == os.path.realpath(r["denovo_abs"]))
    about = ("Population-level numbers come only from these populations -R folders. From each assembly folder, "
             "the other files read are ones written before populations runs: the ustacks and cstacks lines of "
             "denovo_map.log, the per-sample ustacks and sstacks files (SAMPLE.alleles, .snps and .matches.tsv.gz) "
             "and the catalog (catalog.snps.tsv.gz, catalog.fa.gz).")
    if inside:
        about += (" {} {} the assembly folder itself: denovo_map.pl ran populations there with -R (its -X "
                  "option), so the populations output in {} is that run's -R output.".format(
                      join_list(inside), isare(inside), "it" if len(inside) == 1 else "each"))
    if not inside:
        about += " The populations output that denovo_map.pl writes into the assembly folder is never used."
    elif len(inside) < len(runs):
        about += (" For the other runs, the populations output that denovo_map.pl writes into the assembly folder "
                  "is never used.")
    w(textwrap.fill(about, width=78))
    w("")
    if ti["paired"]:
        w(textwrap.fill("Paired-end data: forward read {} bp; median polymorphic locus {:.0f} bp. Only the forward "
                        "read is assembled under M and n, so SNP-dense loci are also counted there: ≥{} SNPs in the "
                        "forward read is the same density as ≥{} in the whole locus.".format(
                            ti["fwd_len"], ti["med_len"] or float("nan"), ti["k_fwd"], args.tail_snps), width=78))
    elif any(r["catalog"] is not None for r in runs):
        w("Single-end data (no paired-end contigs found in catalog.fa.gz).")
    else:
        w("No catalog.fa.gz was read, so SNPs are counted over whole loci.")
    w("")
    if not series:
        w("No series found: a series needs at least two runs that differ in one parameter (M with n = M, M, n "
          "or m) and share the rest.")
    for s in series:
        R = s["R"]
        blk = []
        b = blk.append
        runs_s, lab = s["runs"], s["labels"]
        b("-" * 78)
        b("{}  |  {}".format(s["title"], s["filter"]))
        b("-" * 78)
        b("Paris et al. (2017) metrics on r80 loci (Figs. 3 and 4).")
        b("{:<8} {:>9} {:>11} {:>16} {:>8}  {}".format("value", "r80 loci", "r80 poly", "new poly (step)", "SNPs",
                                                     "Paris"))
        for i, r in enumerate(runs_s):
            st = s["steps"][i - 1] if i > 0 else None
            tags = []
            if s["best"] == i:
                tags.append("r80 max")
            if s["kind"] in M_KINDS and s["fig4"] == i:
                tags.append("Fig. 4")
            if s["W"] == i:
                tags.append("plateau")
            if s["pick"] == i:
                tags.append("PICK")
            b("{:<8} {:>9,} {:>11,} {:>16} {:>8,}  {}".format(
                lab[i], r["loci"], r["poly"], "{:+,} ({})".format(st["d"], pct(st["rel"])) if st else "",
                r["snps"], ", ".join(tags)))
        legend = "r80 max = most polymorphic r80 loci"
        if s["kind"] in M_KINDS:
            legend += ("; Fig. 4 = later steps add <{:.1%} new ones; plateau = polymorphic loci within individuals first "
                       "level off".format(args.plateau))
        b(textwrap.fill(legend + "; PICK = Paris recommendation.", width=78))
        b("")
        ind = s["ind"]
        if any(ind[k]["ok"] for k in ind):
            n_ind = next(ind[k]["steps"][0]["n"] for k in ind if ind[k]["ok"])
            b(textwrap.fill("Within individuals (Paris Results and Table S6): means over {} individuals. loci = loci per "
                            "individual; merged = share of each individual's stacks merged into loci (allele "
                            "collapse); het = heterozygous loci per individual (ustacks), which decides M. Steps are "
                            "paired by individual, ± 2 SE; * = beyond 2 SE.".format(n_ind), width=78))
            b("{:<8} {:>8} {:>7} {:>9} {:>12}".format("value", "loci", "merged", "het", "step"))

            def cell(key, i, fmt):
                v = ind[key]["means"][i]
                return fmt.format(v) if v is not None else "n/a"

            for i in range(len(runs_s)):
                hs = ind["ind_het"]["steps"][i - 1] if i > 0 and ind["ind_het"]["ok"] else None
                b("{:<8} {:>8} {:>7} {:>9} {:>12}".format(
                    lab[i], cell("ind_loci", i, "{:,.0f}"), cell("ind_collapse", i, "{:.1%}"),
                    cell("ind_het", i, "{:,.1f}"), pm(hs, 0) if hs else ""))
            b("")
        b(textwrap.fill("Over-merging signs Paris et al. name, as shares of polymorphic r80 loci: SNP-dense loci, and "
                        "loci with a site where at least {:.0%} of all samples are heterozygous ('heterozygous in a "
                        "large majority of individuals').".format(args.het_site), width=78))
        if ti["paired"]:
            b("{:<8} {:>18} {:>19} {:>16}".format(
                "value", "≥{} SNPs, locus".format(args.tail_snps), "≥{} forward read".format(ti["k_fwd"]),
                "≥{:.0%} het".format(args.het_site)))
            for i, r in enumerate(runs_s):
                b("{:<8} {:>18} {:>19} {:>16}".format(
                    lab[i], count_share(r["tail"], r["tail_frac"]), count_share(r["fwd_tail"], r["fwd_tail_frac"]),
                    count_share(r["het_loci"], r["het_frac"])))
        else:
            b("{:<8} {:>18} {:>16}".format("value", "≥{} SNPs".format(args.tail_snps),
                                           "≥{:.0%} het".format(args.het_site)))
            for i, r in enumerate(runs_s):
                b("{:<8} {:>18} {:>16}".format(lab[i], count_share(r["tail"], r["tail_frac"]),
                                               count_share(r["het_loci"], r["het_frac"])))
        b("")
        for line in s["verdict"]:
            ind_ = len(line) - len(line.lstrip())
            b(textwrap.fill(line.strip(), width=78, initial_indent=" " * ind_, subsequent_indent=" " * (ind_ + 4)))
        if s["notes"]:
            b("")
            b("Notes:")
            wrap_lines(s["notes"], b)
        b("")
        out.extend(relabel(x, R) for x in blk)

    w("-" * 78)
    w("Recommendation following Paris et al. (2017)")
    w("-" * 78)
    if not recs:
        w("  none: no series could be analysed")
    for (R, ns), lines, label, _brief in recs:
        if len(recs) > 1:
            w("Runs with -R {:g} and {} samples:".format(R, ns))
        if label:
            w("  {}".format(label))
        wrap_lines(lines, w, "  - ")
    w("")
    w("-" * 78)
    w("Fixed SNPs in the catalog (Paris Fig. S3; reported, not used to choose)")
    w("-" * 78)
    fxr = [r for r in sorted(runs, key=lambda r: (r["m"], r["M"], r["n"])) if r["fixed"]]
    if fxr:
        w(textwrap.fill(
            "Paris et al. compared the SNPs in each individual with those in the catalog (their "
            "count_fixed_catalog_snps.py, written for Stacks 1). A catalog SNP is fixed when no individual is "
            "heterozygous at it. Counted here on the catalog cstacks builds from the forward reads "
            "(catalog.snps.tsv.gz), with each sample's heterozygous ustacks SNPs (SAMPLE.snps.tsv.gz) placed on it "
            "through sstacks's alignment (SAMPLE.matches.tsv.gz). het in an ind = catalog SNPs at which at least "
            "one individual is heterozygous; het SNPs/ind = Paris et al.'s count (i), the mean heterozygous SNPs "
            "per individual (ustacks); placed = the mean that land on a catalog SNP (the rest are in loci sstacks "
            "did not match back to the catalog).", width=78))
        w("{:<9} {:>12} {:>14} {:>7} {:>7} {:>13} {:>7}".format("params", "catalog SNPs", "het in an ind", "fixed",
                                                              "% fixed", "het SNPs/ind", "placed"))
        for r in fxr:
            f = r["fixed"]
            w("{:<9} {:>12,} {:>14,} {:>7,} {:>7} {:>13} {:>7}".format(
                r["label"], f["cat_snps"], f["variable"], f["fixed"], "{:.1%}".format(f["fixed_frac"]),
                "{:,.0f}".format(r["ind_snps_mean"]) if r["ind_snps_mean"] is not None else "n/a",
                "{:,.0f}".format(r["ind_snps_cat_mean"]) if r["ind_snps_cat_mean"] is not None else "n/a"))
        w(textwrap.fill(
            "Paris et al. found that fixed SNPs rose with n in all three datasets without levelling off: 'This "
            "alone, however, does not provide enough information to choose a value for n.' The rise was gentlest in "
            "penguin, which has little divergence between colonies, and steepest in earthworm, which has highly "
            "divergent cryptic lineages. At the n they chose, 15% (trout, n4), 13% (penguin, n2) and 21% "
            "(earthworm, n3) of catalog SNPs were fixed.", width=78))
        w("")
        w(textwrap.fill(
            "Where the fixed SNPs are (a reading aid, not from Paris; not used to choose). between = the haplotypes "
            "of the individuals sstacks matched to the locus show two or more bases there, although no "
            "individual's SNP calls are heterozygous: a difference between individuals. no-match = in catalog loci "
            "that cstacks built but sstacks matched no sample locus back to, so they never reach gstacks or "
            "populations (loci = how many such catalog loci). other = those haplotypes show at most one base: the "
            "second allele is in sample loci that sstacks did not match back (sstacks's lines in denovo_map.log "
            "count these per sample), or is an N. Shares are of all catalog SNPs.", width=78))
        w("{:<12} {:>16} {:>16} {:>8} {:>16}".format("params", "between", "no-match", "loci", "other"))
        for r in fxr:
            f = r["fixed"]
            cs = f["cat_snps"]
            w("{:<12} {:>16} {:>16} {:>8,} {:>16}".format(
                r["label"], count_share(f["between"], share(f["between"], cs)),
                count_share(f["nomatch"], share(f["nomatch"], cs)), f["nomatch_loci"],
                count_share(f["unplaced"], share(f["unplaced"], cs))))
        nser = [s for s in series if s["kind"] == "n" and all(r["fixed"] for r in s["runs"])]
        for s in nser:
            w(textwrap.fill("Along n ({}): {}.".format(s["title"].split("(", 1)[1].rstrip(")"), join_list(
                "{:,} fixed ({:.1%}) at {}".format(r["fixed"]["fixed"], r["fixed"]["fixed_frac"], s["labels"][i])
                for i, r in enumerate(s["runs"]))), width=78))
        if not nser:
            w(textwrap.fill("There is no n series (n varied with m and M fixed) yet, which is how Paris et al. "
                            "plotted fixed SNPs (Fig. S3). {}".format(
                                "In an M = n series, M and n change together, so its changes in fixed SNPs are not "
                                "due to n alone." if any(s["kind"] == "M=n" for s in series) else ""), width=78))
    else:
        w("Not counted: catalog.snps.tsv.gz, SAMPLE.snps.tsv.gz or SAMPLE.matches.tsv.gz is missing (see "
          "Warnings).")
    w("")
    w("-" * 78)
    w("Coverage (Paris Fig. 1)")
    w("-" * 78)
    means = [r["cov_mean"] for r in runs if r["cov_mean"] == r["cov_mean"]]
    if means:
        w(textwrap.fill("ustacks coverage per sample after merging alleles ('Coverage after assembling stacks', before "
                        "secondary reads are added; Paris Fig. 1's merged coverage), averaged over "
                        "samples: {}–{} over the runs read. At m3, Paris et al.'s datasets averaged 19x–42x before "
                        "merging alleles and 21x–49x after, and they suggest aiming above {:.0f}x.".format(
                            fmt_cov(min(means)), fmt_cov(max(means)), LOW_MEAN_COV), width=78))
        if min(means) < LOW_MEAN_COV:
            w("  Below the {:.0f}x Paris et al. suggest.".format(LOW_MEAN_COV))
        low = sorted(set(smp for r in runs for smp, _ in r["low_cov"]))
        if low:
            w(textwrap.fill("Samples below {:.0f}x in at least one run (genotyping error rises below 10x, Fountain et "
                            "al. 2016, cited by Paris et al.), with their range over all runs:".format(LOW_SAMPLE_COV),
                            width=78))
            for smp in low:
                cv = [r["ustacks_cov"][smp] for r in runs if smp in r["ustacks_cov"]]
                w("  {}: {}–{}, below {:.0f}x in {} of {} runs".format(
                    smp, fmt_cov(min(cv)), fmt_cov(max(cv)), LOW_SAMPLE_COV,
                    sum(1 for c in cv if c < LOW_SAMPLE_COV), len(cv)))
        else:
            w("No sample is below {:.0f}x in any run.".format(LOW_SAMPLE_COV))
    else:
        w("No per-sample ustacks coverage was found in denovo_map.log.")
    w("")
    w("-" * 78)
    w("Warnings")
    w("-" * 78)
    if WARNINGS:
        for x in WARNINGS:
            w(textwrap.fill(x, width=78, initial_indent="  - ", subsequent_indent="    "))
    else:
        w("  none")
    if footer:
        w("")
        w(textwrap.fill(footer, width=78))
    return "\n".join(out) + "\n"


def build_brief(recs, out_dir):
    """The short screen output: the recommendation, one line per parameter, and any warnings."""
    out = []
    w = out.append
    if not recs:
        w("No recommendation: no series could be analysed.")
    for (R, ns), _lines, label, brief in recs:
        if len(recs) > 1:
            if out:
                w("")
            w("Runs with -R {:g} and {} samples:".format(R, ns))
        w("Paris et al. (2017) recommendation: {}".format(label or "none"))
        for line in brief:
            w(textwrap.fill(relabel(line, R), width=78, initial_indent="  ", subsequent_indent="     "))
    if WARNINGS:
        w("")
        w("Warnings:")
        for x in WARNINGS:
            w(textwrap.fill(x, width=78, initial_indent="  - ", subsequent_indent="    "))
    w("")
    w("Full report: {} (--verbose prints it)".format(os.path.join(out_dir, "paris_report.txt")))
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- spreadsheet tables

def rnd(v, d):
    if v is None or (isinstance(v, float) and v != v):
        return ""
    return round(v, d)


def summary_rows(runs, args, ti):
    hdr = ["run", "populations_dir", "denovo_dir", "m", "M", "n", "R", "samples", "r80_loci",
           "r80_polymorphic_loci", "pct_polymorphic", "snps", "snps_per_polymorphic_locus",
           "individuals", "loci_per_individual", "stacks_merged_share", "het_loci_per_individual",
           "het_snps_per_individual",
           "loci_ge{}_snps".format(args.tail_snps), "pct_poly_ge{}_snps".format(args.tail_snps)]
    if ti["paired"]:
        hdr += ["loci_ge{}_fwd_snps".format(ti["k_fwd"]), "pct_poly_ge{}_fwd_snps".format(ti["k_fwd"])]
        if ti["k_pe"]:
            hdr += ["loci_ge{}_pe_snps".format(ti["k_pe"]), "pct_poly_ge{}_pe_snps".format(ti["k_pe"])]
    hdr += ["high_het_loci", "pct_poly_high_het", "high_het_sites", "catalog_loci", "ustacks_cov_mean",
            "ustacks_cov_median", "ustacks_cov_min", "samples_lt{:.0f}x".format(LOW_SAMPLE_COV)]
    rows = []
    for r in sorted(runs, key=lambda r: (r["m"], r["M"], r["n"])):
        cat = r["catalog"]
        row = [r["label"], r["popdir"], r["denovo_dir"], r["m"], r["M"], r["n"], r["R"], r["n_samples"],
               r["loci"], r["poly"], rnd(r["pct_poly"], 4), r["snps"], rnd(r["snps_per_poly"], 3),
               len(r["ind_loci"]) or len(r["ind_het"]) or "", rnd(r["ind_loci_mean"], 1),
               rnd(r["ind_collapse_mean"], 5), rnd(r["ind_het_mean"], 2), rnd(r["ind_snps_mean"], 1),
               r["tail"], rnd(r["tail_frac"], 5)]
        if ti["paired"]:
            row += [rnd(r["fwd_tail"], 0), rnd(r["fwd_tail_frac"], 5)]
            if ti["k_pe"]:
                row += [rnd(r["pe_tail"], 0), rnd(r["pe_tail_frac"], 5)]
        row += [rnd(r["het_loci"], 0), rnd(r["het_frac"], 5), rnd(r["het_sites"], 0), cat["n_loci"] if cat else "",
                rnd(r["cov_mean"], 2), rnd(r["cov_median"], 2), rnd(r["cov_min"], 2), len(r["low_cov"])]
        rows.append(row)
    return hdr, rows


def series_rows(series):
    hdr = ["series", "filter", "value", "run", "r80_loci", "r80_polymorphic_loci", "snps", "new_polymorphic_loci",
           "pct_change", "loci_per_individual", "stacks_merged_share", "het_loci_per_individual", "het_loci_step",
           "het_loci_step_2se", "het_loci_step_individuals_up", "pct_poly_snp_dense",
           "pct_poly_high_het", "catalog_snps", "fixed_catalog_snps", "fixed_share",
           "het_snps_per_individual_on_catalog", "r80_max", "fig4_point", "within_ind_plateau", "paris_pick"]
    rows = []
    for s in series:
        ind = s["ind"]
        for i, r in enumerate(s["runs"]):
            st = s["steps"][i - 1] if i > 0 else None
            hs = ind["ind_het"]["steps"][i - 1] if i > 0 and ind["ind_het"]["ok"] else None
            yes = lambda k: "yes" if s[k] == i else ""
            rows.append([s["title"], s["filter"], s["labels"][i], r["label"], r["loci"], r["poly"], r["snps"],
                         st["d"] if st else "", rnd(st["rel"], 5) if st else "",
                         rnd(ind["ind_loci"]["means"][i], 1), rnd(ind["ind_collapse"]["means"][i], 5),
                         rnd(ind["ind_het"]["means"][i], 2), rnd(hs["d"], 2) if hs else "",
                         rnd(2 * hs["se"], 2) if hs else "", "{} of {}".format(hs["up"], hs["n"]) if hs else "",
                         rnd(s["tail_vals"][i], 5), rnd(r["het_frac"], 5),
                         r["fixed"]["cat_snps"] if r["fixed"] else "", r["fixed"]["fixed"] if r["fixed"] else "",
                         rnd(r["fixed"]["fixed_frac"], 5) if r["fixed"] else "",
                         rnd(r["ind_snps_cat_mean"], 1), yes("best"),
                         yes("fig4") if s["kind"] in M_KINDS else "", yes("W"), yes("pick")])
    return hdr, rows


def individual_rows(runs):
    """One row per run and sample: the per-individual numbers behind Paris et al.'s within-individual metrics."""
    hdr = ["run", "m", "M", "n", "sample", "ustacks_stacks", "ustacks_loci_after_M", "ustacks_final_loci",
           "stacks_merged_share", "het_loci",
           "het_snps", "het_snps_on_catalog", "ustacks_cov"]
    rows = []
    for r in sorted(runs, key=lambda r: (r["m"], r["M"], r["n"])):
        for smp in sorted(r["ind"]):
            d = r["ind"][smp]
            rows.append([r["label"], r["m"], r["M"], r["n"], smp, rnd(d["stacks"], 0), rnd(d["loci_M"], 0),
                         rnd(d["final"], 0), rnd(r["ind_collapse"].get(smp), 5), rnd(d["het"], 0),
                         rnd(d["het_snps"], 0), rnd(d["het_snps_cat"], 0), rnd(r["ustacks_cov"].get(smp), 2)])
    return hdr, rows


def fixed_rows(runs):
    """Paris et al.'s comparison of each individual's SNPs with the catalog, one row per run."""
    hdr = ["run", "m", "M", "n", "catalog_snps", "catalog_snps_het_in_an_individual", "fixed_snps", "fixed_share",
           "fixed_between_individuals", "fixed_in_loci_no_sample_matched", "loci_no_sample_matched",
           "fixed_other_at_most_one_base_seen", "individuals", "het_snps_per_individual",
           "het_snps_per_individual_on_catalog"]
    rows = []
    for r in sorted(runs, key=lambda r: (r["m"], r["M"], r["n"])):
        f = r["fixed"]
        if not f:
            continue
        rows.append([r["label"], r["m"], r["M"], r["n"], f["cat_snps"], f["variable"], f["fixed"],
                     rnd(f["fixed_frac"], 5), f["between"], f["nomatch"], f["nomatch_loci"], f["unplaced"],
                     len(f["ind_cat"]), rnd(statistics.mean(f["ind_all"].values()), 1),
                     rnd(statistics.mean(f["ind_cat"].values()), 1)])
    return hdr, rows


def high_het_rows(runs):
    """Polymorphic r80 loci with a site heterozygous in at least --het-site of all samples, per run."""
    def key(x):
        return (0, int(x), "") if x.isdigit() else (1, 0, x)
    rows = []
    for r in sorted(runs, key=lambda r: (r["m"], r["M"], r["n"])):
        for loc in sorted(r["het_set"] or (), key=key):
            rows.append([r["label"], int(loc) if loc.isdigit() else loc])
    return ["run", "locus_id"], rows


def relabel_table(hdr, rows, runs):
    Rs = set(r["R"] for r in runs)
    if len(Rs) == 1:
        R = Rs.pop()
        return [relabel(h, R) for h in hdr], rows
    return [h.replace("r80", "rR") for h in hdr], rows


def write_xlsx(path, runs, series, args, ti):
    """Supporting tables. The recommendation and the main evidence are in paris_report.txt."""
    from openpyxl import Workbook
    from openpyxl.chart import LineChart, Reference
    from openpyxl.styles import Font
    from openpyxl.utils import get_column_letter

    bold = Font(bold=True)
    wb = Workbook()
    srt = sorted(runs, key=lambda r: (r["m"], r["M"], r["n"]))
    R0 = srt[0]["R"]

    def put_table(ws, hdr, rows, note=None):
        for j, h in enumerate(hdr, 1):
            ws.cell(row=1, column=j, value=h).font = bold
        for i, row in enumerate(rows, 2):
            for j, v in enumerate(row, 1):
                ws.cell(row=i, column=j, value=v)
        for j, h in enumerate(hdr, 1):
            width = max([len(str(h))] + [len(str(r[j - 1])) for r in rows]) + 2
            ws.column_dimensions[get_column_letter(j)].width = min(width, 45)
        ws.freeze_panes = ws.cell(row=2, column=2)
        if note:
            ws.cell(row=len(rows) + 3, column=1, value=relabel(note, R0))

    ws = wb.active
    ws.title = "Runs"
    put_table(ws, *relabel_table(*summary_rows(runs, args, ti), runs),
              note="One row per run. Within-individual columns are means over individuals. pct_poly_* are shares of "
                   "polymorphic r80 loci. high_het = loci with a site where at least {:.0%} of all samples are "
                   "heterozygous. ustacks_cov = each sample's ustacks coverage after merging alleles (Paris Fig. 1)."
                   .format(args.het_site))

    ws = wb.create_sheet("Series")
    put_table(ws, *relabel_table(*series_rows(series), runs),
              note="One row per value in each series. Steps are from the previous value, paired by individual; 2se = "
                   "two standard errors. r80_max, fig4_point, within_ind_plateau and paris_pick mark the values the "
                   "report names.")

    ws = wb.create_sheet("Individuals")
    put_table(ws, *relabel_table(*individual_rows(runs), runs),
              note="Paris et al.'s within-individual metrics, per sample. ustacks numbers come from denovo_map.log and "
                   "SAMPLE.alleles.tsv.gz (het_loci = loci with two or more haplotypes); stacks_merged_share = 1 - "
                   "loci after M / stacks. het_snps = heterozygous "
                   "SNPs in SAMPLE.snps.tsv.gz; het_snps_on_catalog = those placed on a catalog SNP through "
                   "SAMPLE.matches.tsv.gz (the rest are in loci sstacks did not match to the catalog).")

    ws = wb.create_sheet("SNPs_per_locus")
    K = args.max_snp_col
    bins = [str(k) for k in range(K + 1)] + ["{}+".format(K + 1)]
    hdr = ["run"] + bins + ["Total", "Polymorphic", ""] + ["prop " + b for b in bins]
    rows = []
    for r in srt:
        counts = [r["snp_dist"].get(k, 0) for k in range(K + 1)]
        counts.append(sum(v for k, v in r["snp_dist"].items() if k > K))
        props = [round(c / r["loci"], 5) if r["loci"] else "" for c in counts]
        rows.append([r["label"]] + counts + [r["loci"], r["poly"], ""] + props)
    put_table(ws, hdr, rows, note="Counts are r80 loci by number of SNPs (snps_per_loc_postfilters). Proportions are "
                                  "of all r80 loci in the run.")
    total_col = len(bins) + 2
    ch = LineChart()
    ch.title = relabel("r80 loci by run", R0)
    ch.y_axis.title = "Loci"
    ch.x_axis.title = "Run"
    ch.x_axis.delete = False
    ch.y_axis.delete = False
    ch.add_data(Reference(ws, min_col=total_col, max_col=total_col + 1, min_row=1, max_row=len(rows) + 1),
                titles_from_data=True)
    ch.set_categories(Reference(ws, min_col=1, min_row=2, max_row=len(rows) + 1))
    ch.height, ch.width = 9, 18
    ws.add_chart(ch, "B{}".format(len(rows) + 5))

    ws = wb.create_sheet("Fixed_SNPs")
    put_table(ws, *fixed_rows(runs),
              note="Paris et al.'s comparison of each individual's SNPs with the catalog (count_fixed_catalog_snps.py; "
                   "Fig. S3). catalog_snps = SNPs in catalog.snps.tsv.gz (forward reads, from cstacks); "
                   "het_in_an_individual = catalog SNPs at which at least one individual is heterozygous; fixed = "
                   "catalog SNPs at which no individual is heterozygous. Paris et al. found fixed SNPs rose with n "
                   "without levelling off, and did not use them alone to choose n. The split of the fixed SNPs is a "
                   "reading aid, not from Paris: between individuals = the haplotypes of individuals sstacks matched "
                   "to the locus show two or more bases; in loci no sample matched = catalog loci sstacks matched no "
                   "sample locus back to; other = those haplotypes show at most one base (the second allele is in "
                   "sample loci sstacks did not match back, or is an N). het_snps_per_individual is Paris et al.'s "
                   "count (i); het_snps_per_individual_on_catalog is the part placed on a catalog SNP.")

    ws = wb.create_sheet("High_het_loci")
    put_table(ws, *high_het_rows(runs),
              note="Polymorphic r80 loci with a site where at least {:.0%} of all samples are heterozygous ('heterozygous "
                   "in a large majority of individuals', Paris et al.'s sign of merged paralogs). Paris et al. remove "
                   "such loci after assembly; in Stacks 2, --max-obs-het removes SNPs only, so pass these IDs to "
                   "populations --blacklist to remove whole loci.".format(args.het_site))

    ws = wb.create_sheet("Coverage")
    samples = sorted(set().union(*[set(r["ustacks_cov"]) for r in srt]))
    hdr = ["sample"] + [r["label"] for r in srt]
    rows = [[smp] + [r["ustacks_cov"].get(smp, "") for r in srt] for smp in samples]
    put_table(ws, hdr, rows, note="ustacks coverage per sample after merging alleles ('Coverage after assembling stacks' in "
                                  "denovo_map.log; Paris Fig. 1). Genotyping error rises below {:.0f}x (Fountain et al. "
                                  "2016, cited by Paris et al.).".format(LOW_SAMPLE_COV))
    wb.save(path)


# --------------------------------------------------------------------------- plots

# Reference palette from the dataviz skill (light mode; the PDF is for print and screen on white).
SURFACE = "#ffffff"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e6e5e1"
BLUE = "#2a78d6"
BLUE_LIGHT = "#cde2fb"
BLUE_DARK = "#104281"
ACCENT = "#eb6834"
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
        "#0d366b"]


def ramp(n):
    if n == 1:
        return [BLUE]
    return [RAMP[round(i * (len(RAMP) - 1) / (n - 1))] for i in range(n)]


def count_fmt(v, _pos=None):
    return "{:,.0f}".format(v)


def write_plots(path, series, args, image=None, image_dir=None):
    """One PDF page per series; with image = 'png' or 'jpeg', also one image file per page."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator

    plt.rcParams.update({
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.titlecolor": TEXT, "axes.labelcolor": TEXT2, "axes.edgecolor": GRID, "axes.linewidth": 0.8,
        "xtick.color": TEXT2, "ytick.color": TEXT2, "xtick.major.size": 0, "ytick.major.size": 0,
        "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.8,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "text.color": TEXT,
    })

    def line(ax, xs, ys, color=BLUE, label=None):
        ax.plot(xs, ys, color=color, lw=2, solid_joinstyle="round", solid_capstyle="round", zorder=2, label=label)
        ax.plot(xs, ys, "o", color=color, ms=7, mec=SURFACE, mew=1.5, zorder=3)

    def bars_err(ax, xs, ys, errs, color=BLUE):
        ax.errorbar(xs, ys, yerr=errs, fmt="none", ecolor=color, elinewidth=1, capsize=4, zorder=2)
        line(ax, xs, ys, color)

    def mark(ax, x, y, text):
        ax.plot([x], [y], "o", color=ACCENT, ms=11, mec=SURFACE, mew=2, zorder=4)
        ax.annotate(text, (x, y), xytext=(0, 14), textcoords="offset points", ha="center", va="bottom", fontsize=8,
                    color=TEXT, zorder=5)

    def xaxis(ax, s):
        ax.set_xticks(s["xs"])
        ax.set_xticklabels(s["labels"])
        ax.set_xlim(min(s["xs"]) - 0.5, max(s["xs"]) + 0.5)

    def headroom(ax, frac=0.18, low=0.05):
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - (hi - lo) * low, hi + (hi - lo) * frac)

    def midpoints(s):
        return [(s["xs"][st["from"]] + s["xs"][st["to"]]) / 2 for st in s["steps"]]

    def vpick(ax, x, text, top=True):
        """Vertical marker at a value, labelled at the top or bottom of the panel."""
        ax.axvline(x, color=ACCENT, lw=1.5, zorder=1)
        ax.annotate(text, (x, 1 if top else 0), xycoords=("data", "axes fraction"), xytext=(0, -2 if top else 2),
                    textcoords="offset points", ha="center", va="top" if top else "bottom", fontsize=8, color=TEXT,
                    zorder=5, bbox={"boxstyle": "square,pad=0.2", "fc": SURFACE, "ec": "none"})

    def missing_panel(ax, text):
        ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=TEXT2)

    def step_labels(ax, s, steps, digits, scale=1.0):
        """Paired steps (± 2 SE) written between the two values they compare, along the top of the panel."""
        for st, x in zip(steps, midpoints(s)):
            ax.annotate("{:+.{d}f}\n±{:.{d}f}{}".format(st["d"] * scale, 2 * st["se"] * scale,
                                                       "*" if st["rise"] or st["fall"] else "", d=digits),
                        (x, 1), xycoords=("data", "axes fraction"), xytext=(0, -2), textcoords="offset points",
                        ha="center", va="top", fontsize=7, color=TEXT2)

    def room(ax, top=0.35, bottom=0.18):
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - (hi - lo) * bottom, hi + (hi - lo) * top)

    def nonneg_ticks(ax):
        """Shares cannot be negative: keep the room below for labels, but drop ticks under zero."""
        lo, hi = ax.get_ylim()
        ax.set_yticks([t for t in ax.get_yticks() if -1e-12 <= t <= hi])
        ax.set_ylim(lo, hi)

    def share_panel(ax, s, vals, ns):
        v = [x * 100 for x in vals]
        bars_err(ax, s["xs"], v, [200 * se_share(x, n) for x, n in zip(vals, ns)])
        headroom(ax, 0.25, 0.25)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "{:.2f}%".format(v)))
        nonneg_ticks(ax)
        if s["pick"] is not None:
            vpick(ax, s["xs"][s["pick"]], "Paris pick", top=False)

    def boxes(ax, s, data):
        step = min(b2 - a2 for a2, b2 in zip(s["xs"], s["xs"][1:]))
        ax.boxplot(data, positions=s["xs"], widths=0.45 * step, patch_artist=True, showfliers=True,
                   medianprops={"color": TEXT, "lw": 1.5},
                   boxprops={"facecolor": BLUE_LIGHT, "edgecolor": BLUE, "lw": 1},
                   whiskerprops={"color": BLUE, "lw": 1}, capprops={"color": BLUE, "lw": 1},
                   flierprops={"marker": "o", "ms": 4, "mfc": BLUE, "mec": SURFACE, "mew": 0.8})

    with PdfPages(path) as pdf:
        for si, s in enumerate(series):
            xs, lab, runs, R = s["xs"], s["labels"], s["runs"], s["R"]
            rl = lambda t: relabel(t, R)
            step = min(b2 - a2 for a2, b2 in zip(xs, xs[1:]))
            mx = midpoints(s)
            ind = s["ind"]
            fig = plt.figure(figsize=(28, 15))
            gs = fig.add_gridspec(3, 6, width_ratios=[1, 1, 1, 1, 1, 1.6], left=0.035, right=0.99, top=0.91,
                                  bottom=0.05, wspace=0.34, hspace=0.42)
            ax = {}
            letters = "abcdefghijklmno"
            for k in range(15):
                ax[letters[k]] = fig.add_subplot(gs[k // 5, k % 5])
            ax_v = fig.add_subplot(gs[:, 5])
            fig.suptitle(rl("{}  ·  {}".format(s["title"], s["filter"])), x=0.035, ha="left", fontsize=14,
                         fontweight="bold", color=TEXT)
            fig.text(0.035, 0.94, rl("Paris et al. (2017) metrics. Top: r80 loci. Middle: within individuals (each "
                                     "sample's ustacks output) and coverage. Bottom: the over-merging "
                                     "signs Paris names, and each individual's SNPs compared with the catalog (fixed "
                                     "SNPs). Error bars and step labels are ± 2 standard errors; * = beyond 2 SE."),
                     fontsize=9.5, color=TEXT2)

            # (a) r80 loci
            a = ax["a"]
            line(a, xs, [r["loci"] for r in runs])
            a.set_title(rl("(a) r80 loci"))
            a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
            xaxis(a, s)
            headroom(a, 0.1)

            # (b) polymorphic r80 loci: the r80 maximum and the Paris pick
            a = ax["b"]
            P = s["P"]
            line(a, xs, P)
            a.set_title(rl("(b) r80 polymorphic loci"))
            a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
            xaxis(a, s)
            headroom(a, 0.3, 0.12)
            b, p = s["best"], s["pick"]
            if p is not None and p == b:
                mark(a, xs[b], P[b], rl("r80 max = Paris pick\n{:,}".format(P[b])))
            else:
                mark(a, xs[b], P[b], rl("r80 max\n{:,}".format(P[b])))
                if p is not None:
                    mark(a, xs[p], P[p], "Paris pick\n{:,}".format(P[p]))

            # (c) SNPs (Paris Fig. 3c)
            a = ax["c"]
            line(a, xs, [r["snps"] for r in runs])
            a.set_title(rl("(c) SNPs in r80 loci"))
            a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
            xaxis(a, s)
            headroom(a, 0.1)

            # (d) new polymorphic loci per step (Paris Fig. 4)
            a = ax["d"]
            sd = [st["d"] for st in s["steps"]]
            a.bar(mx, sd, width=0.45 * step, color=BLUE, zorder=2)
            a.axhline(0, color=TEXT2, lw=0.8, zorder=3)
            a.set_title(rl("(d) New r80 polymorphic loci per step (Fig. 4)"))
            a.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "{:+,.0f}".format(v) if v else "0"))
            xaxis(a, s)
            lo, hi = min(sd + [0]), max(sd + [0])
            span = (hi - lo) or 1
            a.set_ylim(lo - 0.22 * span, hi + 0.22 * span)
            for st, x in zip(s["steps"], mx):
                a.annotate("{:+,}\n{}".format(st["d"], pct(st["rel"])), (x, st["d"]),
                           xytext=(0, 4 if st["d"] >= 0 else -4), textcoords="offset points", ha="center",
                           va="bottom" if st["d"] >= 0 else "top", fontsize=7.5, color=TEXT2)

            # (f) loci per individual
            a = ax["f"]
            a.set_title("(f) Loci per individual (ustacks)")
            xaxis(a, s)
            if ind["ind_loci"]["ok"]:
                boxes(a, s, [sorted(r["ind_loci"].values()) for r in runs])
                a.plot(xs, ind["ind_loci"]["means"], "D", color=BLUE_DARK, ms=5, zorder=4, label="mean")
                a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
                a.set_ylabel("one point per individual")
                a.legend(fontsize=7.5, loc="upper right")
                xaxis(a, s)
            else:
                missing_panel(a, "no per-sample ustacks lines\nin denovo_map.log")

            # (g) allele collapse (Paris Table S6a)
            a = ax["g"]
            a.set_title("(g) Allele collapse: stacks merged into loci")
            xaxis(a, s)
            if ind["ind_collapse"]["ok"]:
                cm = [v * 100 for v in ind["ind_collapse"]["means"]]
                line(a, xs, cm)
                headroom(a, 0.25, 0.1)
                for x, v in zip(xs, cm):
                    a.annotate("{:.1f}%".format(v), (x, v), xytext=(0, 8), textcoords="offset points", ha="center",
                               va="bottom", fontsize=7.5, color=TEXT2)
                a.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "{:.0f}%".format(v)))
                a.set_ylabel("mean share of each individual's stacks")
            else:
                missing_panel(a, "no per-sample ustacks lines\nin denovo_map.log")

            # (h), (i) polymorphism within individuals, with paired steps and where it first levels off
            for key, letter, title, digits, scale, ylab in (
                    ("ind_het", "h", "(h) Heterozygous loci per individual", 0, 1.0, "mean over individuals (ustacks)"),
                    ("ind_snps", "i", "(i) Heterozygous SNPs per individual", 0, 1.0,
                     "mean over individuals (ustacks)")):
                a = ax[letter]
                a.set_title(title)
                xaxis(a, s)
                if ind[key]["ok"]:
                    line(a, xs, [v * scale for v in ind[key]["means"]])
                    step_labels(a, s, ind[key]["steps"], digits, scale)
                    room(a)
                    if scale != 1.0:
                        a.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "{:.1f}%".format(v)))
                    else:
                        a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
                    a.set_ylabel(ylab)
                    k = s["wpl"].get(key)
                    if k is not None:
                        vpick(a, xs[k], "first levels off", top=False)
                else:
                    missing_panel(a, "not available\n(see Warnings)")

            # (e) SNPs-per-locus distribution
            a = ax["e"]
            cl = s["class_labels"]
            bx = list(range(len(cl)))
            for sh, c, lb in zip(s["shares"], ramp(len(runs)), lab):
                a.plot(bx[:-1], sh[:-1], color=c, lw=2, label=lb, solid_capstyle="round")
                a.plot(bx[-1:], sh[-1:], "o", color=c, ms=6, mec=SURFACE, mew=1)
            a.axvline(len(cl) - 1.5, color=GRID, lw=0.8)
            a.set_yscale("log")
            a.yaxis.set_major_locator(FixedLocator([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]))
            a.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "{:g}%".format(round(v * 100, 3))))
            a.yaxis.set_minor_locator(NullLocator())
            a.yaxis.set_minor_formatter(NullFormatter())
            a.set_xticks(bx)
            a.set_xticklabels(cl)
            a.set_title(rl("(e) SNPs per polymorphic r80 locus"))
            a.set_xlabel("SNPs in locus ({} pooled, not on the line)".format(cl[-1]))
            a.set_ylabel(rl("share of polymorphic r80 loci (log scale)"))
            a.legend(fontsize=7.5, ncol=2, loc="upper right", handlelength=1.5)

            # (k) SNP-dense loci (Paris's over-merging sign)
            a = ax["k"]
            a.set_title("(k) Loci with {}".format(s["tail_what"]))
            xaxis(a, s)
            if all(v is not None and v == v for v in s["tail_vals"]):
                share_panel(a, s, s["tail_vals"], [r["poly"] for r in runs])
                a.set_ylabel(rl("share of polymorphic r80 loci"))
            else:
                missing_panel(a, "no SNP counts")

            # (l) loci heterozygous in most individuals (Paris's paralog sign)
            a = ax["l"]
            a.set_title("(l) Loci with a ≥{:.0%} heterozygous site".format(args.het_site))
            xaxis(a, s)
            if all(r["het_frac"] is not None for r in runs):
                share_panel(a, s, [r["het_frac"] for r in runs], [r["poly"] for r in runs])
                a.set_ylabel("all samples pooled")
            else:
                missing_panel(a, "populations.sumstats.tsv\nmissing for some runs")

            # (j) coverage (Paris Fig. 1: ustacks coverage)
            a = ax["j"]
            a.set_title("(j) Coverage per sample (ustacks; Fig. 1)")
            xaxis(a, s)
            if all(r["ustacks_cov"] for r in runs):
                boxes(a, s, [sorted(r["ustacks_cov"].values()) for r in runs])
                for yref, txt in ((LOW_SAMPLE_COV, "10x"), (LOW_MEAN_COV, "25x")):
                    a.axhline(yref, color=MUTED, lw=0.8, zorder=1)
                    a.annotate(txt, (1, yref), xycoords=("axes fraction", "data"), xytext=(3, 0),
                               textcoords="offset points", ha="left", va="center", fontsize=7.5, color=TEXT2)
                a.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "{:.0f}x".format(v)))
                a.set_ylabel("after merging alleles")
                xaxis(a, s)
                a.set_ylim(bottom=0)
            else:
                missing_panel(a, "no per-sample coverage\nin denovo_map.log")

            # (m)-(o) each individual's SNPs compared with the catalog (Paris's count_fixed_catalog_snps.py)
            fx = [r["fixed"] for r in runs]
            a = ax["m"]
            a.set_title("(m) Het. SNPs per individual placed on catalog")
            xaxis(a, s)
            if all(r["ind_snps_cat_mean"] is not None for r in runs):
                line(a, xs, [r["ind_snps_cat_mean"] for r in runs])
                headroom(a, 0.2, 0.1)
                a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
                a.set_ylabel("mean over individuals")
            else:
                missing_panel(a, "no SAMPLE.snps / .matches\nor catalog.snps files")
            a = ax["n"]
            a.set_title("(n) Catalog SNPs")
            xaxis(a, s)
            if all(fx):
                line(a, xs, [f["cat_snps"] for f in fx], color=BLUE_DARK, label="all catalog SNPs")
                line(a, xs, [f["variable"] for f in fx], color=RAMP[0], label="heterozygous in ≥1 individual")
                a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
                headroom(a, 0.3, 0.1)
                a.legend(fontsize=7.5, loc="upper left")
            else:
                missing_panel(a, "no SAMPLE.snps / .matches\nor catalog.snps files")
            a = ax["o"]
            a.set_title("(o) Fixed SNPs in the catalog (Fig. S3)")
            xaxis(a, s)
            if all(fx):
                bottom = [0] * len(fx)
                for key, color, lb in (("between", BLUE_DARK, "2+ bases among matched individuals"),
                                       ("nomatch", BLUE, "in loci no sample matched"),
                                       ("unplaced", BLUE_LIGHT, "≤1 base among matched individuals")):
                    v = [f[key] for f in fx]
                    a.bar(xs, v, bottom=bottom, width=0.55 * step, color=color, label=lb, zorder=2)
                    bottom = [b2 + x for b2, x in zip(bottom, v)]
                for x, f in zip(xs, fx):
                    a.annotate("{:,}\n{:.1%}".format(f["fixed"], f["fixed_frac"]), (x, f["fixed"]), xytext=(0, 3),
                               textcoords="offset points", ha="center", va="bottom", fontsize=7.5, color=TEXT2)
                a.set_ylim(0, max(f["fixed"] for f in fx) * 1.45)
                a.yaxis.set_major_formatter(FuncFormatter(count_fmt))
                a.set_ylabel("catalog SNPs with no heterozygous individual")
                a.legend(fontsize=7, loc="upper left", ncol=1)
                if s["kind"] == "M=n":
                    a.set_xlabel("M and n change together here")
            else:
                missing_panel(a, "no SAMPLE.snps / .matches\nor catalog.snps files")

            # verdict text, full height on the right
            ax_v.axis("off")
            width = 78
            body = [textwrap.fill(ln.strip(), width=width, subsequent_indent="    " if ln.startswith("  ") else "  ",
                                  initial_indent="  " if ln.startswith("  ") else "") for ln in s["verdict"]]
            budget = 92 - sum(x.count("\n") + 1 for x in body) - 1
            shown = 0
            for nt in s["notes"]:
                blk = textwrap.fill(nt, width=width, initial_indent="• ", subsequent_indent="  ")
                if blk.count("\n") + 1 > budget - 2:
                    break
                if shown == 0:
                    body.append("")
                body.append(blk)
                budget -= blk.count("\n") + 1
                shown += 1
            if shown < len(s["notes"]):
                body.append("")
                body.append("({} more note(s) in paris_report.txt)".format(len(s["notes"]) - shown))
            ax_v.text(0, 1, "\n".join(body), va="top", ha="left", fontsize=8, color=TEXT, transform=ax_v.transAxes,
                      linespacing=1.3)
            ax_v.set_title("Verdict")

            pdf.savefig(fig)
            if image:
                ext = "jpg" if image == "jpeg" else "png"
                fn = os.path.join(image_dir, "paris_plot_{}.{}".format(si + 1, ext))
                if ext == "jpg":
                    try:
                        fig.savefig(fn, dpi=100, format="jpeg", pil_kwargs={"quality": 92})
                    except TypeError:       # matplotlib < 3.3
                        fig.savefig(fn, dpi=100, format="jpeg")
                else:
                    fig.savefig(fn, dpi=100)
            plt.close(fig)


# --------------------------------------------------------------------------- main

def available_cpus():
    try:
        return len(os.sched_getaffinity(0)) or 1
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Choose Stacks 2 de novo parameters (m, M, n) following Paris et al. (2017), from "
                    "`populations -R 0.8` output. Writes paris_report.txt (recommendation and main evidence), "
                    "paris_plots.pdf and paris_summary.xlsx (supporting tables). The screen shows only the "
                    "recommendation and any warnings; --verbose prints the full report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  paris_optimize.py M1_R80 M2_R80 M3_R80 M4_R80 -o paris_results\n"
               "  paris_optimize.py --runs runs.tsv --image png\n\n"
               "A --runs table is tab-separated with a header line: populations_dir (required), and optionally\n"
               "denovo_dir, m, M, n. Relative paths are relative to the table. See README_paris_optimize.md.")
    ap.add_argument("dirs", nargs="*", help="populations output directories (each run with -R 0.8)")
    ap.add_argument("--runs", help="tab-separated table of runs (instead of, or as well as, directories)")
    ap.add_argument("-o", "--out", default="paris_results", help="output directory (default: paris_results)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="print the full report to the screen (it is always saved as paris_report.txt)")
    ap.add_argument("--image", choices=["png", "jpeg", "jpg"],
                    help="also save each plot page as an image file: png or jpeg")
    ap.add_argument("--plateau", type=float, default=0.01,
                    help="Paris Fig. 4 reading: a step adding fewer new polymorphic r80 loci than this share of "
                         "the previous count adds 'very few' (default 0.01). Used only when polymorphic r80 loci "
                         "are still rising at the highest value and the within-individual metrics are missing.")
    ap.add_argument("--dist-max-snps", type=int, default=6,
                    help="SNPs-per-locus plot: classes 1..N; loci with more SNPs are pooled (default 6)")
    ap.add_argument("--tail-snps", type=int, default=10,
                    help="a locus with at least this many SNPs counts as SNP-dense (default 10)")
    ap.add_argument("--fwd-tail-snps", type=int, default=None,
                    help="paired-end data: SNP-dense threshold for the forward read (default: --tail-snps scaled "
                         "by forward-read length / locus length)")
    ap.add_argument("--het-site", type=float, default=0.7,
                    help="Paris's paralog sign: a site where at least this share of all samples is "
                         "heterozygous ('heterozygous in a large majority of individuals') is flagged (default 0.7)")
    ap.add_argument("-t", "--threads", type=int, default=min(8, available_cpus()),
                    help="worker processes for reading the per-sample files of the fixed-SNP count (default: the "
                         "CPUs available, up to 8)")
    ap.add_argument("--max-snp-col", type=int, default=12,
                    help="last SNP-count column in the workbook; higher counts are pooled (default 12)")
    args = ap.parse_args(argv)
    if not args.dirs and not args.runs:
        ap.error("give populations directories or --runs")
    if args.image == "jpg":
        args.image = "jpeg"

    runs = []
    try:
        specs = [{"popdir": d} for d in args.dirs]
        if args.runs:
            specs += read_runs_table(args.runs)
        seen = set()
        for spec in specs:
            real = os.path.realpath(spec["popdir"])
            if real in seen:
                warn("{} was given more than once; using it once".format(spec["popdir"]))
                continue
            seen.add(real)
            runs.append(load_run(spec, args))
    except ParseError as e:
        sys.exit("error: {}".format(e))
    flush_missing(len(runs))

    for key, label in (("n_samples", "sample count"), ("R", "-R value"), ("other_flags", "populations options")):
        vals = sorted(set(str(r[key]) for r in runs))
        if len(vals) > 1:
            warn("runs differ in {}: {}. Runs are only compared within the same sample count and -R value."
                 .format(label, "; ".join("'{}'".format(v) for v in vals)))
    for R in sorted(set(r["R"] for r in runs)):
        if abs(R - 0.8) > 1e-9:
            warn("populations was run with -R {:g}; Paris et al.'s rule uses -R 0.8. Loci are labelled {} "
                 "instead of r80.".format(R, rlabel(R)))

    ti = tail_settings(runs, args)
    apply_tails(runs, ti)
    series = [analyse(s, args, ti) for s in find_series(runs)]
    recs = recommend(series, runs)
    in_series = set(id(r) for s in series for r in s["runs"])
    for r in runs:
        if id(r) not in in_series:
            warn("{} ({}) is not part of any series, so it does not enter a decision.".format(
                r["popdir"], r["label"]))

    try:
        import openpyxl  # noqa: F401
        have_xlsx = True
    except ImportError:
        have_xlsx = False
        sys.stderr.write("warning: openpyxl not available, so no paris_summary.xlsx\n")
    try:
        import matplotlib  # noqa: F401
        have_plots = bool(series)
    except ImportError:
        have_plots = False
        sys.stderr.write("warning: matplotlib not available, so no paris_plots.pdf\n")
    extra = []
    if have_plots:
        extra.append("plots in paris_plots.pdf")
    if have_xlsx:
        extra.append("supporting tables (per run, per series value, per individual, SNPs per locus, fixed SNPs, "
                     "high-heterozygosity loci, coverage per sample) in paris_summary.xlsx")
    footer = "See also: {}.".format(join_list(extra)) if extra else None

    os.makedirs(args.out, exist_ok=True)
    report = build_report(runs, series, recs, args, ti, footer)
    say(report if args.verbose else build_brief(recs, args.out))
    with open(os.path.join(args.out, "paris_report.txt"), "w", encoding="utf-8") as fh:
        fh.write(report)
    written = ["paris_report.txt"]
    if have_plots:
        write_plots(os.path.join(args.out, "paris_plots.pdf"), series, args, image=args.image, image_dir=args.out)
        written.append("paris_plots.pdf")
        if args.image:
            written.append("paris_plot_1–{}.{}".format(len(series), "jpg" if args.image == "jpeg" else "png")
                           if len(series) > 1 else "paris_plot_1.{}".format("jpg" if args.image == "jpeg" else "png"))
    elif args.image:
        sys.stderr.write("warning: --image needs matplotlib and at least one series; no images written\n")
    if have_xlsx:
        write_xlsx(os.path.join(args.out, "paris_summary.xlsx"), runs, series, args, ti)
        written.append("paris_summary.xlsx")
    say("{}Wrote {} in {}/\n".format("\n" if args.verbose else "", ", ".join(written), args.out))


if __name__ == "__main__":
    main()
