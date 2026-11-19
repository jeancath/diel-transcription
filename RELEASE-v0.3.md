# Release archive — v0.3

Diel transcriptional programs across cyanobacteria.

> **Addendum.** This file is the manifest of **v0.3 as released**, and its checksums are not
> updated. It lists some files that this public snapshot does not carry, because the snapshot
> is an allowlisted export rather than a copy of the working repository. For a manifest of
> exactly what is here, see `MANIFEST.sha256`.

## What this file is

The index of the archived release: every file that belongs to v0.3, with a
SHA-256 of the exact bytes. A file not listed here is not part of the release.
Checksums are of the tagged tree, which differs from the outputs commit only
by this index file, so `sha256sum -c` run from `diel-conservation/` after checking out the
tag reproduces every line below.

## Identity

| | |
|---|---|
| Source freeze | `b9c978b` |
| Outputs commit | `9f5a089` |
| Tagged commit | this file's commit, which adds only this index |
| Annotated tag object | `277db66` |
| Build | `make clean && make check && make all` in a fresh clone of the source freeze |

`git diff <source freeze> <outputs commit> -- src Makefile MANUSCRIPT.template.md`
is empty: the tagged outputs were produced by exactly the frozen source, and
the outputs commit changed no code.

## Verdicts

| Domain | Verdict |
|---|---|
| Pipeline reproducibility, content level | Pass |
| Pipeline reproducibility, byte level | Incomplete: one gzip header |
| Scientific validity | Qualified pass |

The three verdicts are separate on purpose. Reproducing a number is not the
same as the number answering the question it is presented as answering; the
reasoning behind each verdict is in `CLOSURE-REPORT.md`.

## Reproducing it

```bash
git clone <repo> && cd <repo>/diel-conservation
git checkout v0.3
pip install -r requirements.txt
make clean && make check && make all
python3 src/compare_builds.py <a results dir> results
```

`make check` runs the semantic assertions, the mutation tests, the renderer
audit and the hand-typed-statistics scan, and fails the build on any of them.
The full run recomputes all 93 simulation cells; on the machine this release
was built on that is a few hours, single-threaded by design, because the
per-cell RNG is seeded from the cell coordinates rather than drawn from one
advancing stream.

## Known non-reproducible byte

`results/rhythms_per_gene.csv.gz` is not byte-reproducible: gzip writes the
current time into its header. Its contents are identical across builds, and
`compare_builds.py` compares it by reading through the compression. The fix is
one argument to the writer and is left for a later release so that this one is
tagged exactly as it was verified.

## Contents

### 1. Manuscript

```
fadd4198fe81c94b2c1e9908ee126fc9563897672441810d29a90ee95862bb21  MANUSCRIPT.md
1c200edbbc51871f2c91857d82667238aa3cb03571f9337d5883fa8fd11358e8  MANUSCRIPT.template.md
```

### 2. Rendered tables and figures

```
d6bbdd58f893eda6ea024d8c7e2707ee1990955c8731fd1346a625a8b81b72de  figures/fig1_effect_size.png
3fea74bc9003db777e21b2d8e0dc5173f21cc5c3c20f19c7c8e1ebc77fb17c7b  figures/fig2_phase_difference_ecdf.png
b03ae3db898c1c126dd481c36d84f1b7389dadb970cccda03ccea31f4b7c9d23  figures/fig3_phase_marginals.png
```

### 3. Result files

```
b31fba2018cf3a71e39071ed78e7f8bf11a33343147ae6b4f1dad0ef87f5a8de  results/ASSERTION_COVERAGE.csv
e0e973e04737224946abca230ea3ff2060b34dd684f0f947096213bedfc7b44b  results/NUMBER_PROVENANCE.csv
5306b06edbbb3e7f49d02e671ccbea7b29878650641953acb851645dd2456b62  results/clock_gene_control.csv
729969035f5e5e12c96401dc8b15231a099a96674ad952fb6f7a80587268c60c  results/dataset_summary.csv
c70548135d1562782427b9357ca5157f1657e2b628d74ea1ab3f59214cf74e2c  results/lit_ferrari_citing_screen.csv
115716ec19f8950eea553cc00d5b2dd0e2dc0dd040bb69fb16b5efbddb6a6dfe  results/lit_query_log.csv
14c513045aa32d0aee5ddf22233dce24ed9ed55c49cb4db7f197740a5d418d7d  results/lit_screen.csv
4e309a16e687bed82ac1677285e0076a33ca6093a8faf005de5cefc6cea5b809  results/n10_leave_one_out.csv
06fcd8e91b53cb1a2f9676c0315ef62fa0b209e8007f52968458a0564947663f  results/n10_pairs_bh_corrected.csv
1afff6491cdfe2d2204af40c230cf535a3fd2e6a5fb4c7b530cfd612c5f2765e  results/n11_effect_vs_pileup.csv
6cc86f0fd8a583692f8d57241b1f54cb8cead605a9dca9c0c18fa89acd3a2f63  results/n12_power_grid.csv
a5c1940d175ec08e7561ea6aac830b3ab8988b2bd42e626c51f7d50cdbd90264  results/n13_three_nulls_1000.csv
897c8eae810c02ed408c82995f938a2fd4e7994780a0b410c1a5db52cb36ab88  results/n14_power_surface.csv
3d2e19c5ffe7883f3298001e3f61dfc582a1a106ef792900d763b06e5a52aafe  results/n15_stratified_null.csv
3b38ab1c3881041f7e321c4cb4a5bf10a74d1537f2f043a0cf3f797ea8f91178  results/n16_design_matched_nulls.csv
8eb23ff871ce6f1625576a3af3288485665ee21bd10fd426293d2360ae05d33c  results/n17_power_small_fractions.csv
d86ad890ac0f3e562c63c238d1f97af62fb16c4cbf52ced41fc1f9f9db48632c  results/n18_stratified_closure_A.csv
1c933a9f4d7089f45a2ee07e33d95f6b1d592639428033ba9d3b0520873167ba  results/n18_stratified_closure_A_last.csv
ded45125e3743ac8c1f703adf9f73603c53cf09388e991fe497d90ff1ce6a362  results/n18_stratified_closure_B.csv
95f8a620988318b775ed714af1ec85612cb06dfca5787c608d008ef7c242de53  results/n18_stratified_closure_B_last.csv
60e2fb8735d8ba706e1f27318980c97c19646183bc0f66d7159ba00ebad5fc56  results/n19_brite_multiassign.csv
446190f9f90ec64b112e08de95d93b63717d84aa052a192035eed4f203eaf0fa  results/n1_two_nulls.csv
e7b9e1b896416bcfe7a2c9844d67fe9e6440a09edadac2dcbb629a21a909984a  results/n2_threshold_sensitivity.csv
262da046d6ea57dd297f2c4dab969e366e989ff3c1418e95f107ef80c49988c6  results/n3_rbh_orthologs.csv
9906e6a1b1adf3d7e1e3f82d568f04e9866f63bbf2a17841fffa0069578f0441  results/n4_null_simulation.csv
95bef1e857dd18d104fc842e249f05e87c7e5da97b6c30ee3a07414073682ba0  results/n5_conserved_core.csv
b7332866248cd215d440548846b9b62aa7832956261d763fa3df65dd0a51d174  results/n5_core_pathway_enrichment.csv
6723284394ee189519844f0d84d43cfc4d83967da030a7954476703c1ff95b1e  results/n5_core_sd_null.csv
9c0058a0483a83c10b001d172624e5c71f5ceb537904aeda0bd5c64aa6fea09e  results/n6_power_analysis.csv
023ae4188f2194413ecdc75262303aa58c452f90b34ccc481e0955283245972e  results/n7_independence_check.csv
f5724a399967d86c4fa1eb0a32690f90adec433b28b0878b8aaba915de7c87b1  results/n8_block_null.csv
39ee8de61a7d38fe48e79ef592a3a9f859f43da43dc9bbc9a7205f7d77c3287e  results/q1_within_species_reproducibility.csv
c555511b364c95012a9dbf2c705e3cdd4c4f81bf6418e8729c5d5d417357c622  results/q2_expected_vs_observed.csv
fca0ebbfbd270fe3d68ebba288a76cf1af7e3b3a951d5b29420b5a29ffde6bc5  results/q2_ko_conservation.csv
cb9237db4f27df79831fa8974bdb9ec2808e19486237a31a77745779ebf5499f  results/q3_phase_sd_by_conservation.csv
261ae4adcb216d54daa75d50f45467c4b6a506178617c4e6ab17998d0c3c7cfe  results/r1_rank_sweep.csv
65e28bd0bd10d6d3598d6028ee93610716896fbd3e56ce499fd8a224301c0f8f  results/r2_phase_test.csv
eba961195c9520353eb7ba7f05712b3cd5a6108c86e418b8fac6d15fc44bcd23  results/r3_phase_marginals.csv
b28f8b25817baeae61077066e2ca81872e1c9f06d1e1a7510930439d198ef780  results/PROVENANCE.json
91019694c0575ae90c328f3cd5a8c6bd80d92f452bf3941180a4670ff30481fb  results/determinism_report.json
a54152f0672821f1a018a3e3bf04588dcd7dfbe96a7cef6abc8984fbd2b91693  results/n13_three_nulls_1000.manifest.json
459b73c6e67185660045b133294905fd608d63f8e08151ba4f93ba602894c43a  results/n14_power_surface.manifest.json
4c9feb3f75b4e34fbfbaf8484c0609b3a99382be65726f73b8aa6ea6e378f885  results/n16_design_matched_nulls.manifest.json
a388b36847d642d6946d46c426df1e7bd4e18f9d6e66b0146c47e3ddda4b8fe6  results/n17_power_small_fractions.manifest.json
ebea18c98bdc2039028767ba3f77544f2b06fd921bdfc6d90c8d6b55bfe7b641  results/n3_rbh_params.json
1ff01fbbe882aedddc0a9e0450ba32e61c20ec696c16941c93d049cde63c12c4  results/numbers.json
2443b55fdd42950f10dafd8aea01b75909b676de6ab65fd86f2263f3ead5b196  results/rhythms_per_gene.csv.gz
```

### 4. Analysis source

```
e9fbff7059f7df23dc94b598eec292b1828f0aa04b350fd4f155f5d667e9316f  src/analyze.py
a4a7132cc36777e2866e448b543b540fbb9491e4252fe2c910cbbd88b8511d17  src/audit_renderer.py
49548e535003187999e6513ad07f032b4e9e3b078d72811960cb501d4aa07a24  src/block_null.py
a27794b4f3cec887c724a631ce6ad4e218bd06821862d86d9ab32d00a7dcf9fb  src/check_manuscript.py
8cbc18c911f46000c67892b5d5df548fdb4a7cc46128fc5976e77b1b298ed00a  src/compare_builds.py
ed3f500b3c304941cf2599bf83eb6d06c9f082f5f7931d7d0701471fb73e4673  src/conserved_core.py
77484492102d5bc815b5615030f2b2fbd3fe7bb7dbae3d5a1deb56e3c93dc761  src/coverage_report.py
757990bb91b679ebe68b3604f7175c324d1ee17cb7ee90786bd65083424a25af  src/datasets.py
d67db345b5f49c1af77a5225048aed06381fabb6841368a12af5bfdc05e9530a  src/dependence_of_pairs.py
6628cbe5a015ed50c28cf59b6376a79b926a9b914c205ff3bb851ed9cdc415b5  src/design_matched_nulls.py
b14c0530b92a8fb5e52240f4b3515cae1483ead539969532b9180a286a66bb98  src/effect_vs_pileup.py
cce6d508838e1b3729f273079f9048f55a69bffb93ab20656c1deeee3817da63  src/fetch_data.py
b1504a61239e15aa057a507b7747e78aa10de09c5b0be30d5d8afaead8bd4261  src/figures.py
d4e4c42621332e81e4257ced3d7ab078756d39342236a39f86e3f57d59139e33  src/independence_check.py
33fd53e8010798f40abe058910012fd56b90966620998c91d6fb28fa222a257a  src/literature_screen.py
7dff4488d4d4366635aea3faa50b4554717947474cfff0355502a03ae131f95c  src/load.py
cbede23811b4f93cb5bb9645ea46ebe19545d7a7e113ebd32c62c3d9357b2bd9  src/manuscript_numbers.py
3ac21a7d43623cf60eaac9fff10881104a85229e07f429c0627e9276cef45bbe  src/null_comparison.py
da6a3631611a9a5370d18bdf435fc20643f848f9ffe7398610427df8ea941f1e  src/null_simulation.py
5025e3a582c262342730c57c4ec752a408c402c70534822fcf3a7aefbc630e80  src/number_provenance.py
4cb03726488edd6d53ddfdb7fca8bd26c7e3b3552f5558f252fc16953a786191  src/orthology.py
c0a7ff7dab308121139ed08c0348afd99d70519551068e19d42732767ad44f67  src/paths.py
2dd5e4cf60e19ff6486918833838dbd9ee01643fc6083d2b6d2ac99756d6d741  src/power_analysis.py
8fb0ede3182ddfdbc9f9a6ad0c873b613fa00cc6332d0b01d623040b850872b4  src/power_grid.py
5206c7b4bcbf05dc57eca2d8b9405180e8db9ed6faa70557229bfb006e240ff7  src/power_grid_full.py
eeba469a2dbfbe244d865d90bab4587aafdd7dbe2434f05c72daa313b521b2a8  src/power_small_fractions.py
61dd9ce21731e49c42e569dc2114ae42386c27d0642f02b8b307a2bc54f9e6d5  src/provenance.py
dbadd51b00e29fc02c12f18ba8976eecb362a2259d6c9cb6b70212e66a115c12  src/rbh.py
75578215e4e280beda6fdadaf6a08082bcc7257066c6ef6bf8dd3c408ee72331  src/rbh_check.py
b000b0813638091380645a1a61239514a033db0573843ce602d950ac3f6352db  src/render_manuscript.py
575cc5b4249d5ce3995346696abf55eaf4f4782ece27c2d9c9afa9267bb7afe3  src/rhythm.py
c7f7a6c2e8bfff4c376e60ecf0891b8f049526e30b84abd422a9c423dd6f9bdf  src/robustness.py
38e23a98d76cc54f4c98f026c01e139b14c0c0ce7c0169f5dd68378f59c28fe7  src/run_rhythms.py
125357401ef5b338bca2de6e7347fc72241508d6c0ba586965ad0f84e4f3eb3f  src/simcheck.py
45750337d405875710692ff3e065a926d314e844965ef0fccaf8bb57ccf8b43b  src/stratified_null.py
aa328ea40e6695689d0c3d25b4d43de7b582c382b96203291a7a87feae9c6508  src/test_cosinor_equivalence.py
8a4f076096338ed0276eb1cd1f70f209a053576b14859bd5e039354e47971665  src/test_determinism.py
36dae544a7f26b31e117ff7a551727ae3478f68be416b8b4eba7e256edd8c009  src/test_semantic_mutations.py
8c389e76545155b3b9fa261162b651ae1c67147a434de7785214bf60b1c12d65  src/three_nulls_fast.py
48f65f913163e804e068e3c756ef485c64c051f1ab42d8e58587118e43b09619  src/threshold_sensitivity.py
5f0f98c009ed245d493c1aa21a43df73632a43722f0c5062bef9bae06522ed88  src/validate_semantics.py
3f51211c9ae61a912cf3fb5d974b127342be960ce63aee000a4dce440b762af7  Makefile
```

### 5. Environment and dependencies

```
bb7fdd5ae75689d33e49af4f228b20a28b0dacb0505403ebec19e6970a26f187  requirements.txt
```

### 6. Provenance manifest

```
b28f8b25817baeae61077066e2ca81872e1c9f06d1e1a7510930439d198ef780  results/PROVENANCE.json
e0e973e04737224946abca230ea3ff2060b34dd684f0f947096213bedfc7b44b  results/NUMBER_PROVENANCE.csv
```

### 7. Simulation configuration and seeds

```
a54152f0672821f1a018a3e3bf04588dcd7dfbe96a7cef6abc8984fbd2b91693  results/n13_three_nulls_1000.manifest.json
459b73c6e67185660045b133294905fd608d63f8e08151ba4f93ba602894c43a  results/n14_power_surface.manifest.json
4c9feb3f75b4e34fbfbaf8484c0609b3a99382be65726f73b8aa6ea6e378f885  results/n16_design_matched_nulls.manifest.json
a388b36847d642d6946d46c426df1e7bd4e18f9d6e66b0146c47e3ddda4b8fe6  results/n17_power_small_fractions.manifest.json
```

### 8. Validation evidence

```
b31fba2018cf3a71e39071ed78e7f8bf11a33343147ae6b4f1dad0ef87f5a8de  results/ASSERTION_COVERAGE.csv
91019694c0575ae90c328f3cd5a8c6bd80d92f452bf3941180a4670ff30481fb  results/determinism_report.json
```

### 9. Deviations log

```
184f32d290094f74d058c19b133141d293e9af4b17c7063d3690f9a2ac5cc961  DEVIATIONS.md
```

### 10. Closure report

```
6600b79b5914723724e1c1953a47decaaf9c3659bf04b719da2915686f87181c  CLOSURE-REPORT.md
2ba3d10b35b6a7bea9786a9dc9bef247a672851964d42c68c8b56affd056a4cc  README.md
```

106 files.
