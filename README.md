# supermatrix #
scripts to add new proteins to an existing alignment using hmms, and simple diagnostics/visualizations

The software here is intended to extend existing supermatrix alignments with new taxa, although a number of general purpose programs are also given here. If something does not work, please email me.

Please also note that some similar diagnostic/manipulation tools exist in other packages by various people, including [AMAS](https://github.com/marekborowiec/AMAS), [SCaFoS](http://megasun.bch.umontreal.ca/Software/scafos/scafos.html) and [BuddySuite](https://github.com/biologyguy/BuddySuite)

## add_taxa_to_align.py ##
Script to add new taxa to an existing protein supermatrix alignment. Proteins from new taxa are untrimmed, though a trimming step may be implemented. Input alignments could be in a number of formats, as a supermatrix with a separate partition file, or individual alignments as separate files.
* Multiple new taxa can be added with `-t`, as space-separate protein files (could be gene models or translated from transcriptomes). 
* By default, only the single best hit is taken (changed with `-m`), and is renamed to the corresponding species given in `-T`. 
* Species names from `-T` and files from `-t` must be in the same order. 
* Several folders with many intermediate files are generated (`-d`, `-E`, `-I`, and `-S`), in case it is needed to re-examine later.
* For alignment format (`-f`), most cases *phylip* format is actually *phylip-relaxed*.
* By default, the e-value cutoff is determined uniquely for each gene based on the lower limit of that hmm against the original gene set. This is to reduce the chance of finding out-paralogs. However, a static e-value cutoff for `hmmsearch` can be given using `--ev-threshold`, though this is not advised. See [below](https://github.com/wrf/supermatrix#determination-of-evalues-for-each-partition) for an explanation.
* To generate the new supermatrix with the added taxa, specify the name of the new file with `-U`. A new partition file will be automatically generated if `-U` is specified.

`add_taxa_to_align.py -a philippe2009_FullAlignment.phy -i philippe2009_partitions.txt -t ~/genomes/apis_mellifera/amel_OGSv3.2_pep.fa -T Apis_mellifera -f phylip-relaxed -U philippe2009_w_amel.aln`

Requires [BioPython](http://biopython.org/wiki/Download), [hmmsearch and hmmbuild](http://hmmer.org/), and [mafft v7.3.10](http://mafft.cbrc.jp/alignment/software/source.html), though could be modified to use any aligner. Older versions fo mafft are compatible if using the `-r` option  (current script requires an option in v7.3), which skips the mafft alignment-trimming step.

Binaries are assumed to be in the user's PATH. This can be changed with the options `--mafft`, `--hmmbin`, and `--fasttree`. Both `--mafft` and `--fasttree` should point to binaries, but `--hmmbin` points to a folder containing both hmmsearch and hmmbuild.

`add_taxa_to_align.py -a philippe2009_FullAlignment.phy -i philippe2009_partitions.txt -t ~/genomes/apis_mellifera/amel_OGSv3.2_pep.fa -T Apis_mellifera -f phylip-relaxed -U philippe2009_w_amel.aln --mafft ~/programs/mafft --hmmbin ~/programs/`

All messages and reports can be captured as standard error (using `2>`), such as `2> philippe2009_w_new_taxa.log`. It is recommended to do this.

## join_alignments.py ##
Join multiple individual alignment files into a supermatrix, allowing for only one occurence of any taxa in each alignment. Names must be the same, though can have unique identifiers (like gene names or numbers) as long as they can be systematically split from the taxon names (using `-d`).

`join_alignments.py -a hehenberger2017_alignments/* -d "@" -u hehenberger2017_supermatrix.fasta`

This can also be used to rejoin alignments produced by `add_taxa_to_align.py` that are manually edited. Use the `-A` option to detect the order from the default output naming scheme of the alignments.

## split_supermatrix_to_taxa.py ##
Split a supermatrix into fasta files, one for each taxa where individual proteins are defined by the partition file. Empty proteins are ignored, but gaps are retained. This is NOT the reverse operation of `join_alignments.py`, which joins multiple alignment files into a supermatrix.

`./split_supermatrix_to_taxa.py -a simion2017_97sp_401632pos_1719genes.fasta.gz -d simion_taxa -p simion2017_partitions.txt`

## check_supermatrix_alignments.py ##
Quick diagnostic script to check matrix occupancy. Adjust format accordingly based on the alignment using the `-f` option. As above, in most cases *phylip* format is probably *phylip-relaxed*.

`check_supermatrix_alignments.py -a philippe2009_FullAlignment.phy -p philippe2009_partitions.txt -f phylip-relaxed`

Requires [BioPython](http://biopython.org/wiki/Download)

To generate a chart of matrix occupancy, add the `-m` option with the name of the new output file. This matrix can be visualized using the R script `draw_matrix_occupancy.R`.

`check_supermatrix_alignments.py -a philippe2009_FullAlignment.phy -p philippe2009_partitions.txt -f phylip-relaxed -m philippe2009_occupancy_matrix.tab`

To reorder the matrix based on taxa in a rooted tree, use the `-T` option with a nexus-format tree. For instance, open any tree in figtree, rotate branches, etc, then copy and paste the tree to a text file, and this is in nexus format. Then run `draw_matrix_occupancy.R` below.

`check_supermatrix_alignments.py -p Metazoa_full_Models_short.txt -a Metazoa_full-fix.phy -f phylip-relaxed -m Metazoa_full_matrix.tab -T Metazoa_full.nex`

## compare_supermatrix_alignments.py ##
Directly compare two output supermatrices, say from two different runs of `add_taxa_to_align.py` using slightly different parameters. This will show genes that missing in one or the other, or are different between the two, perhaps due to finding incorrect genes or different splice variants.

Note that partitions must be the same, meaning only vary by presence or absence.

## filter_supermatrix.py ##
Filter supermatrices based on coverage for each gene. Minimum coverage is given by `-c` for values between 0 and 1. A new partition file is automatically generated based on the output name `-o`. 

`filter_supermatrix.py -a simion2017_97sp_401632pos_1719genes.fasta -c 0.75 -p simion2017_partitions.txt -o simion2017_97sp_75cov.fasta`

Requires [BioPython](http://biopython.org/wiki/Download)

## draw_matrix_occupancy.R ##
A graph of matrix occupancy can be generated with the accompanied R script. Genes that are present are colored blue, partial genes are red, and absent genes are white, for each taxa. Taxa with 100% occupancy are colored green, while those with under 50% are colored purple. By default, the taxon order in the graph is the same order as given in the supermatrix alignment, which is then preserved in the matrix occupancy chart above. Taxa can be reordered arbitrarily in either file.

The script can be run in the terminal as:

`Rscript draw_matrix_occupancy.R philippe2009_occupancy_matrix.tab`

![philippe2009_occupancy_matrix.png](https://github.com/wrf/supermatrix/blob/master/philippe2009_occupancy_matrix.png)

## coverage_by_site.py ##
Simple diagnostic for checking coverage by site, providing a histogram of frequency of number of gaps. Change format using `-f`. The script can accept gzipped files.

`coverage_by_site.py -a alignments/philippe2009_FullAlignment.phy -f phylip-relaxed`

Some results are summarized below. With the exception of the Borowiec 2015 study, which made use of entirely genomic data, none of them have a single site that is covered by all taxa. Even for the Borowiec study, this is only 180 sites out of 384k, so not even close to 1%.

| Dataset         | Taxa |  Total sites | Max cov (gaps) | Sites w/ max | Avg gaps % |
|-----------------|------|--------------|----------------|--------------|------------|
| Dunn 2008       |  77  |     9918     |      70 (7)    |          21  |   42.60    |
| Philippe 2009   |  55  |    30257     |      54 (1)    |          91  |   14.79    |
| Ryan 2013 EST   |  70  |    88384     |      58 (12)   |          52  |   41.52    |
| Nosenko 2013 R  |  71  |    14612     |      63 (8)    |         105  |   19.24    |
| Nosenko 2013 nR |  71  |     9187     |      67 (4)    |           3  |   24.98    |
| Whelan 2015 D10 |  70  |    59733     |      68 (2)    |          99  |   26.65    |
| Whelan 2015 D16 |  70  |    23680     |      68 (2)    |          70  |   24.51    |
| Borowiec 2015   |  36  |    384981    |      36 (0)    |         180  |    8.69    |
| Cannon 2016     |  78  |    44896     |      75 (3)    |          59  |   24.25    |
| Simion 2017     |  97  |    401632    |      95 (2)    |          28  |   38.07    |

## test alignments and occupancy matrices ##
* Dunn 2008 dataset of 65 genes, from [Broad phylogenomic sampling improves resolution of the animal tree of life](https://www.nature.com/articles/nature06614)
* Philippe et al (2009) dataset of 128 genes and 30k positions, from [Phylogenomics Revives Traditional Views on Deep Animal Relationships](https://www.sciencedirect.com/science/article/pii/S0960982209008057)
* [Hejnol 2009](https://bitbucket.org/caseywdunn/hejnol_etal_2009) dataset of 1486 genes for 270k positions, from [Assessing the root of bilaterian animals with scalable phylogenomic methods](http://rspb.royalsocietypublishing.org/content/276/1677/4261)
* [Schierwater 2009](https://doi.org/10.1371/journal.pbio.1000020.st004) datasets of 24 and 73 species for 17k sites, from [Concatenated Analysis Sheds Light on Early Metazoan Evolution and Fuels a Modern Urmetazoon Hypothesis](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000020)
* [Erwin 2011](http://science.sciencemag.org/highwire/filestream/593073/field_highwire_adjunct_files/0/Databases_S1-S4.zip) dataset of 10 genes (7 protein + 3 rRNA, 6k sites), from [The Cambrian Conundrum: Early Divergence and Later Ecological Success in the Early History of Animals](http://science.sciencemag.org/content/334/6059/1091)
* [Parfrey 2011](https://treebase.org/treebase-web/search/study/summary.html?id=10562) dataset of 16 genes, from [Broadly Sampled Multigene Analyses Yield a Well-resolved Eukaryotic Tree of Life](https://academic.oup.com/sysbio/article/59/5/518/1645425) and [Estimating the timing of early eukaryotic diversification with multigene molecular clocks](http://www.pnas.org/content/108/33/13624.full)
* [Ryan 2013](https://research.nhgri.nih.gov/manuscripts/Baxevanis/science2013_supplement/) dataset of 406 genes ("EST set"), from [The Genome of the Ctenophore Mnemiopsis leidyi and Its Implications for Cell Type Evolution](http://science.sciencemag.org/content/342/6164/1242592)
* [Nosenko 2013](https://data.ub.uni-muenchen.de/55/) datasets of 35 genes/9k sites and 87 genes/14k sites from 71 taxa, from [Deep metazoan phylogeny: When different genes tell different stories](http://www.sciencedirect.com/science/article/pii/S1055790313000298)
* [Cannon 2014](http://datadryad.org/resource/doi:10.5061/dryad.20s7c) datasets of 299 and 185 genes, from [http://www.sciencedirect.com/science/article/pii/S0960982214012925](Phylogenomic Resolution of the Hemichordate and Echinoderm Clade)
* [Misof 2014](http://datadryad.org/resource/doi:10.5061/dryad.3c0f1) dataset of 1478 genes for 584k sites, from [http://science.sciencemag.org/content/346/6210/763](Phylogenomics resolves the timing and pattern of insect evolution)
* [Weigert 2014](https://datadryad.org/resource/doi:10.5061/dryad.g2qp5) dataset of 421 genes (the 428 partitions were guessed from the matrix, so 7 are wrong), from [Illuminating the base of the annelid tree using transcriptomics](https://academic.oup.com/mbe/article/31/6/1391/1009370)
* [Whelan 2015](https://figshare.com/articles/Error_signal_and_the_placement_of_Ctenophora_sister_to_all_other_animals/1334306) Dataset-10 (for Fig 3), of 210 genes and 59k sites, from [Error, signal, and the placement of Ctenophora sister to all other animals](http://www.pnas.org/content/112/18/5773.full)
* [Zapata 2015](https://bitbucket.org/caseywdunn/cnidaria2014) dataset of 1262 genes for 365k positions, from [Phylogenomic Analyses Support Traditional Relationships within Cnidaria](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139068)
* [Borowiec 2015](http://datadryad.org/resource/doi:10.5061/dryad.k6tq2) dataset of 1080 genes for 384k positions from 36 taxa with genomes, from [Extracting phylogenetic signal and accounting for bias in whole-genome data sets supports the Ctenophora as sister to remaining Metazoa](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2146-4)
* Cannon et al (2016) dataset of 212 genes, from [Xenacoelomorpha is the sister group to Nephrozoa](http://www.nature.com/nature/journal/v530/n7588/full/nature16520.html) 
* [Simion et al 2017](https://github.com/psimion/SuppData_Metazoa_2017) dataset of 1719 genes, where partition file has been reduced to only the numbers and ? in the supermatrix are replaced with gaps, from [A Large and Consistent Phylogenomic Dataset Supports Sponges as the Sister Group to All Other Animals](http://www.sciencedirect.com/science/article/pii/S0960982217301999)
* [Hehenberger 2017](http://datadryad.org/resource/doi:10.5061/dryad.26bv4) dataset of 255 genes from 38 taxa, from [Novel Predators Reshape Holozoan Phylogeny and Reveal the Presence of a Two-Component Signaling System in the Ancestor of Animals](http://www.sciencedirect.com/science/article/pii/S0960982217307078)
* [Whelan 2017](https://figshare.com/articles/Ctenophora_Phylogeny_Datasets_and_Core_Orthologs/4484138) datasets of 350, 212, and 117 genes, from [Ctenophore relationships and their placement as the sister group to all other animals](https://www.nature.com/articles/s41559-017-0331-3)

## determination of evalues for each partition ##
For programs like BLAST or HMMSEARCH, the bitscore, and ultimately the E-value, is dependent on the length of the matched portion. Thus, a good match of a short protein will never have a bitscore as high as a good match for a long protein. For this reason, a static E-value cutoff is not suitable for identifying orthologs in new species.

Considering the chart below, each point represents a self-hit from the HMM profile against the original dataset that was used to make the profile. The longest proteins are dark blue (~700AAs) and the shortest ones are red (~100AAs), with a gradient in between. Proteins that are 80% of the length of the alignment are indicated by the dark circles.

![philippe2009_w_coral_selfhits.png](https://github.com/wrf/supermatrix/blob/master/philippe2009_w_coral_selfhits.png)

It is immediately evident that the highest E-values belong to the longest proteins. Thus, it is clear that a single value cannot be used as a filter for all partitions in a supermatrix. However, even for a long protein, it is clear that many proteins in the original set have E-values substantially lower than the max. These are partial proteins that are kept in the matrix. Thus, the threshold for each partition must be determined primarily by the long proteins, not the lowest value.

A single threshold by E-value based on only long sequences would work, but this would never allow partial matches, as the target proteins would probably have to be the same length as most of the complete sequences. For this, an additional heuristic is used, based on the bitscore-per-length (BpL). This measurement can effectively sort out out-paralogs, but can also help to identify partial sequences. This is necessary because a closely related protein may be full length, and ultimately get a higher bitscore, than a real protein that is only partially complete (say in a transcriptome). Thus, any hits that have a higher BpL than the mean for that partition are kept anyway, even if they are too short, and this takes precedent over a longer hit with a much lower BpL.
