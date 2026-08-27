# Sweep commands — read grouping

Copy-paste. One parameter changes per run; combinations are not swept.

Build `ec21b6e8`. Argument order differs by script kind, because the apply-group scripts carry the
vote arguments the grouping scripts have no use for:

```
grouping scripts (…_group.sh that call `metabuli grouping`)
hash  min-overlap-ratio  weak-band-ratio  max-kmer-quantile  common-kmer-span  threads  partitions
 1           2                  3                4                  5             6         7

Metabuli apply-group scripts (…_metabuli_…, …_clinical_viral_…)
hash  min-overlap-ratio  weak-band-ratio  max-kmer-quantile  common-kmer-span  weight-mode  min-vote-score  group-threads  group-partitions
 1           2                  3                4                  5              6              7               8               9

Kraken2 / Centrifuge / SPAdes apply-group scripts
hash  min-overlap-ratio  weak-band-ratio  max-kmer-quantile  common-kmer-span  group-threads  group-partitions
 1           2                  3                4                  5                6               7
```

On the apply-group scripts none of arguments 1–5 and the trailing two run anything: they name the
grouping run tag, which is how the script finds the grouping directory to read. They must match the
grouping run exactly or the directory will not exist.

Run tag, identical on both kinds:

```
<hash>_<rho>_<weak-band-ratio>[_cs<span>]_q<quantile>[_p<partitions>][_t<threads>]
```

`_p` appears only when partitions ≠ 0, `_t` only when threads ≠ 64.

Baseline: `rho 0.40, weak-band-ratio 0.3333, quantile 0.995, span 0`, which is what
`setGroupGenerationDefaults` also uses. `--min-edge` existed only between `93b8e45a` and
`0b44e30d`; it was removed and `--weak-band-ratio` restored, because an absolute 10 is 10/33 = 0.30
of the core on strain-madness but 10/26 = 0.38 on species-exclusion — the absolute form is the one
that means different things on different data.

**`--cpus-per-task` must match the thread count.** The scripts request 64; override it on the
sbatch line whenever argument 6 is not 64, or the run is packed onto fewer cores than it spawns
threads for.

All jobs pinned to `super001` with `-w`.

## 0. Build

```bash
cd /home/lunajang/src/Metabuli_work && git pull origin master
mkdir -p build_ec21b6e8 && cd build_ec21b6e8
cmake -DCMAKE_BUILD_TYPE=Release .. && make -j64
cd /home/lunajang/src/read-grouping-analysis && git pull origin luna
```

## 1. Baseline grouping — 6 datasets [x]

```bash
sbatch -w super001 benchmark/benchmark_exclusion_group.sh ec21b6e8 0.40 0.3333 0.995 0
sbatch -w super001 benchmark/benchmark_inclusion_group.sh ec21b6e8 0.40 0.3333 0.995 0
sbatch -w super001 cami2/cami2_strain_madness_group.sh    ec21b6e8 0.40 0.3333 0.995 0
sbatch -w super001 cami2/cami2_marine_group.sh            ec21b6e8 0.40 0.3333 0.995 0
sbatch -w super001 cami2/cami2_plant_associated_group.sh  ec21b6e8 0.40 0.3333 0.995 0
sbatch -w super001 cami2/cami2_clinical_pathogen_group.sh ec21b6e8 0.40 0.3333 0.995 0
```

Submitted batch job 550749
Submitted batch job 550750
Submitted batch job 550751
Submitted batch job 550752
Submitted batch job 550753
Submitted batch job 550754

## 1b. Thread independence — strain-madness, partitions pinned

The point of the run: with `--partitions` fixed the routing no longer follows `--threads`, so the
three outputs must be identical. Leaving partitions at 0 would change the routing along with the
thread count and prove nothing.

```bash
sbatch -w super001 --cpus-per-task=4  cami2/cami2_strain_madness_group.sh ec21b6e8 0.40 0.3333 0.995 0 4  16
sbatch -w super001 --cpus-per-task=16 cami2/cami2_strain_madness_group.sh ec21b6e8 0.40 0.3333 0.995 0 16 16
sbatch -w super001 --cpus-per-task=64 cami2/cami2_strain_madness_group.sh ec21b6e8 0.40 0.3333 0.995 0 64 16
```

Tags: `ec21b6e8_0.40_0.3333_q0.995_p16_t4`, `_p16_t16`, `_p16` (64 threads is the default, so no
`_t`). Compare with `cmp` on `groupMap`, `groupRep` and `groups`.

## 2. Sweep grouping — species-exclusion + strain-madness

Two datasets, not six: exclusion has a per-read key and controlled coverage, strain-madness has
several genomes per species. 16 runs.

```bash
G1=benchmark/benchmark_exclusion_group.sh
G2=cami2/cami2_strain_madness_group.sh

# --min-overlap-ratio
for r in 0.3 0.5;       do sbatch -w super001 $G1 ec21b6e8 $r   0.3333 0.995 0; sbatch -w super001 $G2 ec21b6e8 $r   0.3333 0.995 0; done
# --weak-band-ratio
for b in 0.20 0.50;     do sbatch -w super001 $G1 ec21b6e8 0.40 $b     0.995 0; sbatch -w super001 $G2 ec21b6e8 0.40 $b     0.995 0; done
# --max-kmer-quantile
for q in 0.990 0.999;   do sbatch -w super001 $G1 ec21b6e8 0.40 0.3333 $q    0; sbatch -w super001 $G2 ec21b6e8 0.40 0.3333 $q    0; done
# --common-kmer-span
for s in 1 2;           do sbatch -w super001 $G1 ec21b6e8 0.40 0.3333 0.995 $s; sbatch -w super001 $G2 ec21b6e8 0.40 0.3333 0.995 $s; done
```

## 3. Label propagation — after step 1

Reads the baseline grouping result; grouping is not recomputed. Kraken2 and Centrifuge output
carries no per-read score, so those scripts are fixed at `--weight-mode 0` and take no vote
arguments — the weight sweep is Metabuli only.

### 3a. weight-mode 2 — baseline

Output goes to `metabuli_group_<tag>` with no vote suffix.

```bash
sbatch -w super001 benchmark/benchmark_exclusion_metabuli_group.sh ec21b6e8 0.40 0.3333 0.995 0 2 0.15
sbatch -w super001 benchmark/benchmark_inclusion_metabuli_group.sh ec21b6e8 0.40 0.3333 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_strain_madness_group.sh    ec21b6e8 0.40 0.3333 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_marine_group.sh            ec21b6e8 0.40 0.3333 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_plant_associated_group.sh  ec21b6e8 0.40 0.3333 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_clinical_pathogen_group.sh ec21b6e8 0.40 0.3333 0.995 0 2 0.15
```

To read a non-default grouping run, append its threads and partitions, e.g. the `_p16_t4` run from
step 1b:

```bash
sbatch -w super001 cami2/cami2_metabuli_strain_madness_group.sh ec21b6e8 0.40 0.3333 0.995 0 2 0.15 4 16
```

### 3b. weight-mode 0 and 1

Output goes to `metabuli_group_<tag>_w0` / `_w1`, so 3a is not overwritten.

```bash
for w in 0 1; do
  sbatch -w super001 benchmark/benchmark_exclusion_metabuli_group.sh ec21b6e8 0.40 0.3333 0.995 0 $w 0.15
  sbatch -w super001 benchmark/benchmark_inclusion_metabuli_group.sh ec21b6e8 0.40 0.3333 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_strain_madness_group.sh    ec21b6e8 0.40 0.3333 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_marine_group.sh            ec21b6e8 0.40 0.3333 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_plant_associated_group.sh  ec21b6e8 0.40 0.3333 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_clinical_pathogen_group.sh ec21b6e8 0.40 0.3333 0.995 0 $w 0.15
done
```

### 3c. Kraken2 / Centrifuge — baseline only

```bash
for s in benchmark/benchmark_exclusion_kraken2_group.sh \
         benchmark/benchmark_inclusion_kraken2_group.sh \
         benchmark/benchmark_exclusion_centrifuge_group.sh \
         benchmark/benchmark_inclusion_centrifuge_group.sh \
         cami2/cami2_kraken2_strain_madness_group.sh \
         cami2/cami2_kraken2_marine_group.sh \
         cami2/cami2_kraken2_plant_associated_group.sh \
         cami2/cami2_kraken2_clinical_pathogen_group.sh \
         cami2/cami2_centrifuge_strain_madness_group.sh \
         cami2/cami2_centrifuge_marine_group.sh \
         cami2/cami2_centrifuge_plant_associated_group.sh ; do
  sbatch -w super001 $s ec21b6e8 0.40 0.3333 0.995 0
done
```

## 4. Propagation on the sweep points — optional

Only worth running where a grouping point is a candidate operating point. Same argument list as
step 2 plus the vote arguments, e.g. span 1 on exclusion:

```bash
sbatch -w super001 benchmark/benchmark_exclusion_metabuli_group.sh ec21b6e8 0.40 0.3333 0.995 1 2 0.15
```

## Notes

`--weak-band-ratio` near 0 is not in the sweep, but if it is ever run: the band becomes
`(~0, core]`, nearly the whole weight range, so Phase 2 and Phase 3 merge far more. Disk and wall
time can jump and purity can fall hard. Run it alone first. metabuli refuses values outside
`(0, 1)` outright.

`clinical-pathogen` has no per-read answer sheet (`9_prepare_clinical_pathogen.sh` does not build
one), so its `gradeGroup` numbers are not meaningful. Use that dataset for group-structure
statistics and the pathogen read-out only.

strain-madness took 3h 40m on 16 threads. These run on 64, but the disk-fold and sort stages are
not compute-bound, so do not expect the full 4x. The first baseline run gives the real figure.

Phases, for reading the logs: Phase 1 core group formation, Phase 2 support-based merging (was
"Phase 1.5"), Phase 3 singleton linking (was "Phase 2").
