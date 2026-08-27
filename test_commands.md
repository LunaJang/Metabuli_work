# Sweep commands — read grouping

Copy-paste. One parameter changes per run; combinations are not swept.

Build `ec21b6e8`. Argument order for every script here:

```
hash  min-overlap-ratio  min-edge  max-kmer-quantile  common-kmer-span  [weight-mode  min-vote-score]
 1           2               3             4                 5              6              7
```

The last two exist only on the Metabuli apply-group scripts. Run tag is
`<hash>_<rho>_e<min-edge>[_cs<span>]_q<quantile>`, so both absolute thresholds are in the
directory name and nothing collides with results from earlier commits.

Baseline: `rho 0.40, min-edge 10, quantile 0.995, span 0`. min-edge 10 because that is what the
most recent strain-madness run actually used — its weak band was `(11, 33]` at span 0 and
`(10, 29]` at span 1, and `--weak-band-ratio 0.3333` is gone.

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
sbatch -w super001 benchmark/benchmark_exclusion_group.sh ec21b6e8 0.40 10 0.995 0
sbatch -w super001 benchmark/benchmark_inclusion_group.sh ec21b6e8 0.40 10 0.995 0
sbatch -w super001 cami2/cami2_strain_madness_group.sh    ec21b6e8 0.40 10 0.995 0
sbatch -w super001 cami2/cami2_marine_group.sh            ec21b6e8 0.40 10 0.995 0
sbatch -w super001 cami2/cami2_plant_associated_group.sh  ec21b6e8 0.40 10 0.995 0
sbatch -w super001 cami2/cami2_clinical_pathogen_group.sh ec21b6e8 0.40 10 0.995 0
```

Submitted batch job 550749
Submitted batch job 550750
Submitted batch job 550751
Submitted batch job 550752
Submitted batch job 550753
Submitted batch job 550754

## 2. Sweep grouping — species-exclusion + strain-madness

Two datasets, not six: exclusion has a per-read key and controlled coverage, strain-madness has
several genomes per species. 16 runs.

```bash
G1=benchmark/benchmark_exclusion_group.sh
G2=cami2/cami2_strain_madness_group.sh

# --min-overlap-ratio
for r in 0.3 0.5;     do sbatch -w super001 $G1 ec21b6e8 $r  10 0.995 0; sbatch -w super001 $G2 ec21b6e8 $r  10 0.995 0; done
# --min-edge
for e in 5 15;        do sbatch -w super001 $G1 ec21b6e8 0.40 $e 0.995 0; sbatch -w super001 $G2 ec21b6e8 0.40 $e 0.995 0; done
# --max-kmer-quantile
for q in 0.990 0.999; do sbatch -w super001 $G1 ec21b6e8 0.40 10 $q    0; sbatch -w super001 $G2 ec21b6e8 0.40 10 $q    0; done
# --common-kmer-span
for s in 1 2;         do sbatch -w super001 $G1 ec21b6e8 0.40 10 0.995 $s; sbatch -w super001 $G2 ec21b6e8 0.40 10 0.995 $s; done
```

## 3. Label propagation — after step 1

Reads the baseline grouping result; grouping is not recomputed. Kraken2 and Centrifuge output
carries no per-read score, so those scripts are fixed at `--weight-mode 0` and take no vote
arguments — the weight sweep is Metabuli only.

### 3a. weight-mode 2 — baseline

Output goes to `metabuli_group_<tag>` with no vote suffix.

```bash
sbatch -w super001 benchmark/benchmark_exclusion_metabuli_group.sh ec21b6e8 0.40 10 0.995 0 2 0.15
sbatch -w super001 benchmark/benchmark_inclusion_metabuli_group.sh ec21b6e8 0.40 10 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_strain_madness_group.sh    ec21b6e8 0.40 10 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_marine_group.sh            ec21b6e8 0.40 10 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_plant_associated_group.sh  ec21b6e8 0.40 10 0.995 0 2 0.15
sbatch -w super001 cami2/cami2_metabuli_clinical_pathogen_group.sh ec21b6e8 0.40 10 0.995 0 2 0.15
```

### 3b. weight-mode 0 and 1

Output goes to `metabuli_group_<tag>_w0` / `_w1`, so 3a is not overwritten.

```bash
for w in 0 1; do
  sbatch -w super001 benchmark/benchmark_exclusion_metabuli_group.sh ec21b6e8 0.40 10 0.995 0 $w 0.15
  sbatch -w super001 benchmark/benchmark_inclusion_metabuli_group.sh ec21b6e8 0.40 10 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_strain_madness_group.sh    ec21b6e8 0.40 10 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_marine_group.sh            ec21b6e8 0.40 10 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_plant_associated_group.sh  ec21b6e8 0.40 10 0.995 0 $w 0.15
  sbatch -w super001 cami2/cami2_metabuli_clinical_pathogen_group.sh ec21b6e8 0.40 10 0.995 0 $w 0.15
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
  sbatch -w super001 $s ec21b6e8 0.40 10 0.995 0
done
```

## 4. Propagation on the sweep points — optional

Only worth running where a grouping point is a candidate operating point. Same argument list as
step 2 plus the vote arguments, e.g. span 1 on exclusion:

```bash
sbatch -w super001 benchmark/benchmark_exclusion_metabuli_group.sh ec21b6e8 0.40 10 0.995 1 2 0.15
```

## Notes

`--min-edge 1` is not in the sweep, but if it is ever run: the band becomes `(1, core]`, nearly the
whole weight range, so Phase 2 and Phase 3 merge far more. Disk and wall time can jump and purity
can fall hard. Run it alone first.

`clinical-pathogen` has no per-read answer sheet (`9_prepare_clinical_pathogen.sh` does not build
one), so its `gradeGroup` numbers are not meaningful. Use that dataset for group-structure
statistics and the pathogen read-out only.

strain-madness took 3h 40m on 16 threads. These run on 64, but the disk-fold and sort stages are
not compute-bound, so do not expect the full 4x. The first baseline run gives the real figure.

Phases, for reading the logs: Phase 1 core group formation, Phase 2 support-based merging (was
"Phase 1.5"), Phase 3 singleton linking (was "Phase 2").
