# supermatrix #
scripts to add new proteins to an existing alignment using hmms, and simple diagnostics/visualizations

The software here is intended to extend existing supermatrix alignments with new taxa, although a number of general purpose programs are also given here. If something does not work, please email me.

### Jump to: ###
* [add_taxa_to_align.py](https://github.com/wrf/supermatrix#add_taxa_to_align) main script to add new taxa to an existing alignment
* [check_supermatrix_alignments.py](https://github.com/wrf/supermatrix#check_supermatrix_alignments) diagnostic for supermatrix occupancy, to be plotted with the Rscript [draw_matrix_occupancy.R](https://github.com/wrf/supermatrix#draw_matrix_occupancy)
* [a long list of examples](https://github.com/wrf/supermatrix#test-alignments-and-occupancy-matrices) from various published datasets

Please also note that some similar diagnostic/manipulation tools exist in other packages by various people, including [AMAS](https://github.com/marekborowiec/AMAS), [SCaFoS](http://megasun.bch.umontreal.ca/Software/scafos/scafos.html) and [BuddySuite](https://github.com/biologyguy/BuddySuite)

**Note that most Python scripts below require** [BioPython](http://biopython.org/wiki/Download).

## [add_taxa_to_align](https://github.com/wrf/supermatrix/blob/master/add_taxa_to_align.py) ##
Script to add new taxa to an existing protein supermatrix alignment. Proteins from new taxa are untrimmed, though a trimming step may be implemented. Input alignments could be in a number of formats, as a supermatrix with a separate partition file, or individual alignments as separate files.

* Multiple new taxa can be added one of two ways, with `-t`, as space-separate list of protein files (could be gene models or translated from transcriptomes) with `-T` as species names. Fasta files from `-t` and species names from `-T` must be in the same order. Alternatively, a tabular input file may be used (with `-X`) to specify both fasta file and species name for each new species.
* By default, only the single best hit is taken (changed with `-m`), and is renamed to the corresponding species given in `-T`.
* Several folders with many intermediate files are generated (`-d`, `-E`, `-I`, and `-S`), in case it is needed to re-examine later. These are automatically named with each run (so they cannot overwrite each other), and probably do not need to be changed.
* For alignment format (`-f`), most cases *phylip* format is actually *phylip-relaxed*.
* By default, the e-value cutoff is determined uniquely for each gene based on the lower limit of that hmm against the original gene set. This is to reduce the chance of finding out-paralogs. However, a static e-value cutoff for `hmmsearch` can be given using `--ev-threshold`, though this is not advised. See [below for an explanation](https://github.com/wrf/supermatrix#determination-of-evalues-for-each-partition).
* To generate the new supermatrix with the added taxa, specify the name of the new file with `-U`. A new partition file will be automatically generated if `-U` is specified.

For example, to specify new species with `-t` and `-T`:

`add_taxa_to_align.py -a philippe2009_FullAlignment.phy -i philippe2009_partitions.txt -t ~/genomes/apis_mellifera/amel_OGSv3.2_pep.fa -T Apis_mellifera -f phylip-relaxed -U philippe2009_w_amel.aln`

To instead use a tabular input to specify multiple input files and species, specify the text file with `-X` and **do not use options** `-t` **and** `-T`.

`add_taxa_to_align.py -a philippe2009_FullAlignment.phy -i philippe2009_partitions.txt -X new_taxa.txt -f phylip-relaxed -U philippe2009_w_amel.aln`

The tabular input file is a text file where each line has one fasta file and then the corresponding taxon name, separated by a tab.

```
~/genomes/apis_mellifera/amel_OGSv3.2_pep.fa	Apis_mellifera
~/genomes/trichinella_spiralis/t_spiralis.WS248.protein.fa	Trichinella_spiralis
...
```

**Requires** [BioPython](http://biopython.org/wiki/Download), [hmmsearch and hmmbuild](http://hmmer.org/), and [mafft v7.3.10](http://mafft.cbrc.jp/alignment/software/source.html), though could be modified to use any aligner. Older versions fo mafft are compatible if using the `-r` option (current script requires an option in v7.3), which skips the mafft alignment-trimming step.

Binaries are assumed to be in the user's PATH. This can be changed with the options `--mafft`, `--hmmbin`, and `--fasttree`. Both `--mafft` and `--fasttree` should point to binaries, but `--hmmbin` points to a folder containing both hmmsearch and hmmbuild.

`add_taxa_to_align.py -a philippe2009_FullAlignment.phy -i philippe2009_partitions.txt -t ~/genomes/apis_mellifera/amel_OGSv3.2_pep.fa -T Apis_mellifera -f phylip-relaxed -U philippe2009_w_amel.aln --mafft ~/programs/mafft --hmmbin ~/programs/`

All messages and reports can be captured as standard error (using `2>`), such as `2> philippe2009_w_new_taxa.log`. It is recommended to do this (possibly in verbose mode `-v`) as the log contains information as to why every hmmsearch hit was selected or rejected, and it may be necessary to manually correct some sequences.

## [join_alignments](https://github.com/wrf/supermatrix/blob/master/join_alignments.py) ##
Join multiple individual alignment files into a supermatrix, allowing for only one occurence of any taxa in each alignment. Names must be the same, though can have unique identifiers (like gene names or numbers) as long as they can be systematically split from the taxon names (using `-d`).

`join_alignments.py -a hehenberger2017_alignments/* -d "@" -u hehenberger2017_supermatrix.fasta`

This can also be used to rejoin alignments produced by `add_taxa_to_align.py` that are manually edited. Use the `-A` option to detect the order from the default output naming scheme of the alignments.

This will automatically generate a partition file for the supermatrix, adding `.partition.txt` to the name from `-u`.

## [split_supermatrix_to_genes](https://github.com/wrf/supermatrix/blob/master/split_supermatrix_to_genes.py) ##
Split a supermatrix back into individual alignment files in fasta format, one for gene defined by the partition file. An optional output directory can be given with `-d`, otherwise files are placed in the present working directory. This is the reverse operation of `join_alignments.py`.

`./split_supermatrix_to_genes.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -d simion_genes -p simion2017_partitions.txt`

## [split_supermatrix_to_taxa](https://github.com/wrf/supermatrix/blob/master/split_supermatrix_to_taxa.py) ##
Split a supermatrix into fasta files, one for each taxa where individual proteins are defined by the partition file. Empty proteins are ignored, but gaps are retained, meaning gaps may need to be removed later depending on the next step. This **is NOT** the reverse operation of `join_alignments.py`, which joins multiple alignment files into a supermatrix.

`./split_supermatrix_to_taxa.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -d simion_taxa -p simion2017_partitions.txt`

## [check_supermatrix_alignments](https://github.com/wrf/supermatrix/blob/master/check_supermatrix_alignments.py) ##
Quick diagnostic script to check matrix occupancy. Adjust format accordingly based on the alignment using the `-f` option. As above, in most cases *phylip* format is probably *phylip-relaxed*.

`check_supermatrix_alignments.py -a philippe2009_FullAlignment.phy -p philippe2009_partitions.txt -f phylip-relaxed`

To generate a chart of matrix occupancy, add the `-m` option with the name of the new output file. This matrix can be visualized using the R script [draw_matrix_occupancy.R](https://github.com/wrf/supermatrix/blob/master/draw_matrix_occupancy.R).

`check_supermatrix_alignments.py -a philippe2009_FullAlignment.phy -p philippe2009_partitions.txt -f phylip-relaxed -m philippe2009_occupancy_matrix.tab`

To reorder the matrix based on taxa in a rooted tree, use the `-T` option with a nexus-format tree. For instance, open any tree in [figtree](http://tree.bio.ed.ac.uk/software/figtree/), rotate branches, etc, then copy and paste the tree to a text file, and this is in nexus format. Then run `draw_matrix_occupancy.R` below.

`check_supermatrix_alignments.py -p Metazoa_full_Models_short.txt -a Metazoa_full-fix.phy -f phylip-relaxed -m Metazoa_full_matrix.tab -T Metazoa_full.nex`

Because `Present` is coded as 50-100% of the full length, this may hide taxa that have mostly partial sequences. To instead output the occupancy matrix *as the percentage* of the full length gene, use the `--percent` option. The R script `draw_matrix_occupancy.R` can be used on this matrix as well.

`check_supermatrix_alignments.py -p Metazoa_full_Models_short.txt -a Metazoa_full-fix.phy -f phylip-relaxed --percent -m Metazoa_full_percent_matrix.tab -T Metazoa_full.nex`

## [compare_supermatrix_alignments](https://github.com/wrf/supermatrix/blob/master/compare_supermatrix_alignments.py) ##
Directly compare two output supermatrices, say from two different runs of `add_taxa_to_align.py` using slightly different parameters. This will show genes that missing in one or the other, or are different between the two, perhaps due to finding incorrect genes or different splice variants.

Note that partitions must be the same, meaning only vary by presence or absence.

## [filter_supermatrix](https://github.com/wrf/supermatrix/blob/master/filter_supermatrix.py) ##
Filter supermatrices based on coverage for each gene (not removing individual sites). Minimum coverage is given by `-c` for values between 0 and 1. A new partition file is automatically generated based on the output name `-o`. 

`filter_supermatrix.py -a simion2017_97sp_401632pos_1719genes.fasta -c 0.75 -p simion2017_partitions.txt -o simion2017_97sp_75cov.fasta`

## [draw_matrix_occupancy](https://github.com/wrf/supermatrix/blob/master/draw_matrix_occupancy.R) ##
A graph of matrix occupancy can be generated with the accompanied R script. Genes that are present are colored blue, partial genes are red, and absent genes are white, for each taxa. Taxa with 100% occupancy are colored green, while those with under 50% are colored purple. By default, the taxon order in the graph is the same order as given in the supermatrix alignment, which is then preserved in the matrix occupancy chart above. Taxa can be reordered arbitrarily in either file.

The script can be run in the terminal as:

`Rscript draw_matrix_occupancy.R philippe2009_occupancy_matrix.tab`

![philippe2009_occupancy_matrix.png](https://github.com/wrf/supermatrix/blob/master/philippe2009_occupancy_matrix.png)

## [constant_breaker_plot](https://github.com/wrf/supermatrix/blob/master/constant_breaker_plot.R) ##
Again, making use of [check_supermatrix_alignments.py](https://github.com/wrf/supermatrix/blob/master/check_supermatrix_alignments.py) with the `-b` mode, this plot shows probably assembly errors (for genome/txome) or annotation errors (for genomes, bad exon/CDS calls). The plot displays the maximum "constant breaker" score for each gene for each taxon. The "constant breaker" score is an accumulated penalty for a long run of amino acids that differ from an overall conserved region. The principle is that sites that are constant in nearly all species are improbably changed in any given taxon, **including fast-evolving taxa**, and this is particularly true [for Glycine or Proline sites](https://github.com/wrf/heteropecilly/tree/master/aa_counts), where a series of "constant-breaking" sites are very likely indicative of technical errors.

The penalty score is proportional to the frequency of the "breaking" amino acid, with a maximum of 1. Meaning, out of 100 taxa, if 99 of them have `G` at some site, and one species has `N`, then the penalty is either `1/freq-G`, which is `1/99`, or `1/freq-N`, which is `1/1`. This means that a totally unique AA is high-scoring, but that even a frequency of 2 or 3 substantially reduces the score to `0.5` or `0.333`.

![constant_breaker_scoring_v2.png](https://github.com/wrf/supermatrix/blob/master/misc_analyses/constant_breaker_scoring_v2.png)

Below shows an example plot using an expanded alignment from the Philippe 2009 matrix. Dark purple spots indicate a gene with a long constant-breaking stretch. The color depends on the length of the longest run for that gene, typically only 1 per gene if any. Unsurprisingly, the human genome, being well-annotated and assembled, appears to have 0 of these errors.

![metazoa_v11_select_const_break_matrix.png](https://github.com/wrf/supermatrix/blob/master/misc_analyses/metazoa_v11_select_const_break_matrix.png)

Example alignments from gene 2 and gene 4 are shown below, showing the regions that do not align correctly, and are almost certainly errors.

![metazoa_const_breaker_ex_g2.png](https://github.com/wrf/supermatrix/blob/master/misc_analyses/metazoa_const_breaker_ex_g2.png)

![metazoa_const_breaker_ex_g4.png](https://github.com/wrf/supermatrix/blob/master/misc_analyses/metazoa_const_breaker_ex_g4.png)


## [coverage_by_site](https://github.com/wrf/supermatrix/blob/master/coverage_by_site.py) ##
Simple diagnostic for checking coverage by site, providing a histogram of frequency of number of gaps. Change format using `-f`. The script can accept gzipped files.

`coverage_by_site.py -a alignments/philippe2009_FullAlignment.phy -f phylip-relaxed`

Some results are summarized below. With the exception of the Borowiec 2015 study, which made use of entirely genomic data, none of them have *a single site* that is covered by all taxa. Even for the Borowiec study, this is only 180 sites out of 384k, so not even close to 1%.

| Dataset         | Taxa | Genes | Occupancy % | Total sites | Max cov (gaps) | Sites w/ max | Avg gaps % |
|-----------------|------|-------|-------------|-------------|----------------|--------------|------------|
| Dunn 2008       |  77  | 65    | 48.07 (2.6) |    9918     |      70 (7)    |          21  |   42.60    |
| Hejnol 2009     |  94  | 1486  | 17.61 (1.2) |   270392    |      77 (17)   |          53  |   79.35    |
| Philippe 2009   |  55  | 128   | 79.32 (2.6) |   30257     |      54 (1)    |          91  |   14.79    |
| Ryan 2013 EST   |  70  | 406   | 46.18 (4.1) |   88384     |      58 (12)   |          52  |   41.52    |
| Nosenko 2013 R  |  71  | 87    | 76.15 (2.2) |   14612     |      63 (8)    |         105  |   19.24    |
| Nosenko 2013 nR |  71  | 35    | 70.14 (4.8) |    9187     |      67 (4)    |           3  |   24.98    |
| Whelan 2015 D10 |  70  | 210   | 68.10 (7.8) |   59733     |      68 (2)    |          99  |   26.65    |
| Whelan 2015 D16 |  70  | 87    | 71.48 (10.6)|   23680     |      68 (2)    |          70  |   24.51    |
| Borowiec 2015   |  36  | 1080  | 79.84 (6.7) |   384981    |      36 (0)    |         180  |    8.69    |
| Cannon 2016     |  78  | 212   | 80.25 (0.2) |   44896     |      75 (3)    |          59  |   24.25    |
| Simion 2017     |  97  | 1719  | 65.88 (8.2) |   401632    |      95 (2)    |          28  |   38.07    |

## test alignments and occupancy matrices ##
Some alignments can be found in the [alignments folder](https://github.com/wrf/supermatrix/tree/master/alignments), and PDFs of occupancy matrices can be found in the [matrix folder](https://github.com/wrf/supermatrix/tree/master/matrix). Datasets used are:

* [Dunn 2008](https://treebase.org/treebase-web/search/study/summary.html?id=2020) dataset of 65 genes, from [Broad phylogenomic sampling improves resolution of the animal tree of life](https://www.nature.com/articles/nature06614)
* Philippe et al (2009) dataset of 128 genes and 30k positions, from [Phylogenomics Revives Traditional Views on Deep Animal Relationships](https://www.sciencedirect.com/science/article/pii/S0960982209008057)
* [Hejnol 2009](https://bitbucket.org/caseywdunn/hejnol_etal_2009) dataset of 1486 genes for 270k positions, from [Assessing the root of bilaterian animals with scalable phylogenomic methods](http://rspb.royalsocietypublishing.org/content/276/1677/4261)
* [Schierwater 2009](https://doi.org/10.1371/journal.pbio.1000020.st004) datasets of 24 and 73 species for 17k sites, from [Concatenated Analysis Sheds Light on Early Metazoan Evolution and Fuels a Modern Urmetazoon Hypothesis](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000020)
* [Erwin 2011](http://science.sciencemag.org/highwire/filestream/593073/field_highwire_adjunct_files/0/Databases_S1-S4.zip) dataset of 10 genes (7 protein + 3 rRNA, 6k sites), from [The Cambrian Conundrum: Early Divergence and Later Ecological Success in the Early History of Animals](http://science.sciencemag.org/content/334/6059/1091)
* [Parfrey 2011](https://treebase.org/treebase-web/search/study/summary.html?id=10562) dataset of 16 genes, from [Broadly Sampled Multigene Analyses Yield a Well-resolved Eukaryotic Tree of Life](https://academic.oup.com/sysbio/article/59/5/518/1645425) and [Estimating the timing of early eukaryotic diversification with multigene molecular clocks](http://www.pnas.org/content/108/33/13624.full)
* [Philippe 2011](https://media.nature.com/original/nature-assets/nature/journal/v470/n7333/extref/nature09676-s3.txt) dataset of 197 genes, from [Acoelomorph flatworms are deuterostomes related to Xenoturbella](https://www.nature.com/articles/nature09676)
* [Ryan 2013](https://research.nhgri.nih.gov/manuscripts/Baxevanis/science2013_supplement/) dataset of 406 genes ("EST set"), from [The Genome of the Ctenophore Mnemiopsis leidyi and Its Implications for Cell Type Evolution](http://science.sciencemag.org/content/342/6164/1242592)
* [Nosenko 2013](https://data.ub.uni-muenchen.de/55/) datasets of 35 genes/9k sites and 87 genes/14k sites from 71 taxa, from [Deep metazoan phylogeny: When different genes tell different stories](http://www.sciencedirect.com/science/article/pii/S1055790313000298)
* [Cannon 2014](http://datadryad.org/resource/doi:10.5061/dryad.20s7c) datasets of 299 and 185 genes, from [Phylogenomic Resolution of the Hemichordate and Echinoderm Clade](http://www.sciencedirect.com/science/article/pii/S0960982214012925)
* [Misof 2014](http://datadryad.org/resource/doi:10.5061/dryad.3c0f1) dataset of 1478 genes for 584k sites, from [Phylogenomics resolves the timing and pattern of insect evolution](http://science.sciencemag.org/content/346/6210/763)
* [Weigert 2014](https://datadryad.org/resource/doi:10.5061/dryad.g2qp5) dataset of 421 genes for 104k sites, from [Illuminating the base of the annelid tree using transcriptomics](https://academic.oup.com/mbe/article/31/6/1391/1009370)
* [Dos Reis 2015](https://figshare.com/articles/Uncertainty_in_the_timing_of_origin_of_animals_and_the_limits_of_precision_in_molecular_timescales/1525089) of 203 genes, mostly taken from the Philippe 2009 set, from [Uncertainty in the Timing of Origin of Animals and the Limits of Precision in Molecular Timescales](https://www.sciencedirect.com/science/article/pii/S096098221501177X)
* [Whelan 2015](https://figshare.com/articles/Error_signal_and_the_placement_of_Ctenophora_sister_to_all_other_animals/1334306) Dataset-10 (for Fig 3), of 210 genes and 59k sites, from [Error, signal, and the placement of Ctenophora sister to all other animals](http://www.pnas.org/content/112/18/5773.full)
* [Zapata 2015](https://bitbucket.org/caseywdunn/cnidaria2014) dataset of 1262 genes for 365k positions, from [Phylogenomic Analyses Support Traditional Relationships within Cnidaria](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139068)
* [Borowiec 2015](http://datadryad.org/resource/doi:10.5061/dryad.k6tq2) dataset of 1080 genes for 384k positions from 36 taxa with genomes, from [Extracting phylogenetic signal and accounting for bias in whole-genome data sets supports the Ctenophora as sister to remaining Metazoa](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2146-4)
* [Cannon 2016](https://datadryad.org/resource/doi:10.5061/dryad.493b7) dataset of 212 genes, from [Xenacoelomorpha is the sister group to Nephrozoa](http://www.nature.com/nature/journal/v530/n7588/full/nature16520.html)
* [Kocot 2016](http://datadryad.org/resource/doi:10.5061/dryad.30k4v.2) dataset of 638 genes, from [Phylogenomics of Lophotrochozoa with Consideration of Systematic Error](https://academic.oup.com/sysbio/article/66/2/256/2449704)
* [Simion 2017](https://github.com/psimion/SuppData_Metazoa_2017) dataset of 1719 genes, where partition file has been reduced to only the numbers and ? in the supermatrix are replaced with gaps, from [A Large and Consistent Phylogenomic Dataset Supports Sponges as the Sister Group to All Other Animals](http://www.sciencedirect.com/science/article/pii/S0960982217301999)
* [Tanner 2017](https://datadryad.org/resource/doi:10.5061/dryad.180nh), dataset of 180 genes for 30k sites, from [Molecular clocks indicate turnover and diversification of modern coleoid cephalopods during the Mesozoic Marine Revolution](http://rspb.royalsocietypublishing.org/content/284/1850/20162818)
* [Hehenberger 2017](http://datadryad.org/resource/doi:10.5061/dryad.26bv4) dataset of 255 genes from 38 taxa, from [Novel Predators Reshape Holozoan Phylogeny and Reveal the Presence of a Two-Component Signaling System in the Ancestor of Animals](http://www.sciencedirect.com/science/article/pii/S0960982217307078)
* [Whelan 2017](https://figshare.com/articles/Ctenophora_Phylogeny_Datasets_and_Core_Orthologs/4484138) datasets of 350, 212, and 117 genes, from [Ctenophore relationships and their placement as the sister group to all other animals](https://www.nature.com/articles/s41559-017-0331-3)
* [Betts 2018](https://bitbucket.org/bzxdp/betts_et_al_2017) dataset of 29 genes across all kingdoms, from [Integrated genomic and fossil evidence illuminates life’s early evolution and eukaryote origin](https://doi.org/10.1038/s41559-018-0644-x)
* [Schwentner 2018](https://datadryad.org/resource/doi:10.5061/dryad.sn35910) dataset of 455 genes (519 partitions guessed) from 96 taxa, from [Tetraconatan phylogeny with special focus on Malacostraca and Branchiopoda—Highlighting the strength of taxon-specific matrices in phylogenomics](https://doi.org/10.1098/rspb.2018.1524)
* [Marletaz 2019](https://zenodo.org/record/1403005) dataset of 1174 genes from 103 taxa, from [A New Spiralian Phylogeny Places the Enigmatic Arrow Worms among Gnathiferans](https://doi.org/10.1016/j.cub.2018.11.042)

## determination of evalues for each partition ##
For programs like BLAST or HMMSEARCH, the bitscore, and ultimately the E-value, is dependent on the length of the matched portion. Thus, a good match of a short protein will never have a bitscore as high as a good match for a long protein. For this reason, a static E-value cutoff is not suitable for identifying orthologs in new species.

Considering the chart below, each point represents a self-hit from the HMM profile against the original dataset that was used to make the profile. The longest proteins are dark blue (~700AAs) and the shortest ones are red (~100AAs), with a gradient in between. Proteins that are 80% of the length of the alignment are indicated by the dark circles.

![philippe2009_w_coral_selfhits.png](https://github.com/wrf/supermatrix/blob/master/philippe2009_w_coral_selfhits.png)

It is immediately evident that the highest E-values belong to the longest proteins. Thus, it is clear that a single value cannot be used as a filter for all partitions in a supermatrix. However, even for a long protein, it is clear that many proteins in the original set have E-values substantially lower than the max. These are partial proteins that are kept in the matrix. Thus, the threshold for each partition must be determined primarily by the long proteins, not the lowest value.

A single threshold by E-value based on only long sequences would work, but this would never allow partial matches, as the target proteins would probably have to be the same length as most of the complete sequences. For this, an additional heuristic is used, based on the [bitscore-per-length (BpL)](https://github.com/wrf/supermatrix/blob/master/misc_analyses/philippe2009_w_coral_selfhits_normalized.pdf). This measurement can effectively sort out out-paralogs, but can also help to identify partial sequences. This is necessary because a closely related protein may be full length, and ultimately get a higher bitscore, than a real protein that is only partially complete (say in a transcriptome). Thus, any hits that have a higher BpL than the mean for that partition are kept anyway, even if they are too short, and this takes precedent over a longer hit with a much lower BpL.

## How good is BLAST bitscore for finding true homologs ##
I investigated this by blasting the [Parra 2007](https://doi.org/10.1093/bioinformatics/btm071) KOG set against itself, allowing for 16 hits each. This is 6 for each KOG (that is, 6 species) and then up to 10 more.
This produced the graph below, there self hits are the black hollow circles, hits to the same KOG are within on color (colors are random) and hits to other KOGs are black. The dotted lines are the cutoff at bitscore 200 and the shortest protein with a bitscore of 200 (which is 100AAs).

All sequences hit themselves best, followed by all other sequences within the same KOG; that was not surprising. Again, there is a clear pattern of max bitscore as a function of length, with an average around `2bits/AA`. Here is a zoom of the shorter proteins, which is most of them.

![parra_blast_v_bitscore_chart_v2.jpg](https://github.com/wrf/supermatrix/blob/master/misc_analyses/bitscore_length/parra_blast_v_bitscore_chart_v2.jpg)

You can see the dark points at the bottom of the graph, which for the most part stay below 200 bits. That is fine for large proteins, but for under 200AAs, that is the bulk of the true hits as well. There are a few strange zones between 400 and 600. Those are probably due to certain proteins have good hits to the next best KOG which are better than some real hits of similar sized proteins. To say another way, real homologs will have a relatively better bitscore than off-target hits, but not an absolutely better one.

Another dotted line is added for the function of `l/2`. This type seems to better capture the off target hits.

Next I try to normalize to spread them out and better view the length-bitscore relationship. The graph below divides the bitscore by length. Self-hits are nearly all at `2l`.

![parra_blast_v_bitscore_chart_v3.jpg](https://github.com/wrf/supermatrix/blob/master/misc_analyses/bitscore_length/parra_blast_v_bitscore_chart_v3.jpg)

In the previous graphs, the off target hits were at the bottom, but if a cutoff of 0.4 (bitscore/query length) is used, this excludes basically only off targets and is not restrictive on the size of proteins.

**A static cutoff is not an appropriate filter for gene clustering, but this is still commonly used!** 

Unfortunately a few off targets slip through that filter too. That would suggest that e-value or bitscore cutoffs basically have to be protein specific, or have some other way of correctly sorting out the clusters. In this set, those off hits will have true hits that are still higher in bitscore, so maybe this could alternatively use some threshold of taking no more sequences than number of species used in the clade or set. 


