# REtarget

REtarget is a computational tool to find suitable genomic regions and design gRNAs for Recursive Editing, a strategy to enhance the efficiency of Homology-directed repair (HDR) by retargeting undesired editing outcomes.

Recursive Editing and REtarget are described in the following paper:

[Moeller, L., Aird, E. J.\*, Schroeder, M. S., Kobel, L., Kissling, L., van de Venn, L., Corn, J. E.\* Recursive Editing improves homology-directed repair through retargeting of undesired outcomes, *Nat. Commun.*, 2022]()

Please cite this paper when using REtarget. By using REtarget you accept the terms of use.

---

## Overview
- [Online tool]()
- [Setup for local usage](#setup-for-local-usage)
- [Local usage](#local-usage)
- [Search parameters](#search-parameters)
- [Format of results](#format-of-results)
- [Database of pre-computed REtarget results](#database-of-pre-computed-REtarget-results)

---

## Online tool
A web-based version of REtarget can be accessed under https://recursive-editing.herokuapp.com/. Note that on-the-fly off-target analysis is only available for local usage. We recommend to run the online tool using the Chrome or Firefox browser.

---

## Setup for local usage
Note that the all functions were tested with python 3.7.10 and a Linux operating system.

### 1. Download and setup required software:
#### **Lindel** 
Chen, W. *et al.* Nucleic Acids Res. 47, 7989–8003 (2019)
```
git clone https://github.com/shendurelab/Lindel.git 
```
#### **inDelphi**
Shen, M. W. *et al.* Nature 563, 646–651 (2018)
```
git clone https://github.com/maxwshen/inDelphi-model
mv inDelphi-model/ indelphi_local/
```
Other edit outcome prediction tools are not supported by this tool.

#### **FlashFry for off-target prediction**
McKenna, A. *et al.* BMC Biol 16, 74 (2018)\
Make sure to run on Java 8, then download FlashFry-implementation:
```
mkdir flash_fry
cd flash_fry
wget https://github.com/mckennalab/FlashFry/releases/download/1.12/FlashFry-assembly-1.12.jar
mkdir tmp
cd ..
```
#### **REtarget**
Please ask authors for access to the REtarget repository.
```
git clone https://github.com/lkmoeller/retarget.git
```
### 2. Setup conda environment
```
cd recursive_editing
conda env create -f environment.yml
conda activate retarget
```
### 3. Add cloned folders to python path
Either permanently add directories to the python path by adding the following line to `~/.bashrc` (Linux, replace directory by the path that contains the cloned directories)
```
export PYTHONPATH="/path/to/Lindel:/path/to/indelphi_local:/path/to/retarget:$PYTHONPATH"
```
or add the following line to the top your script
```
import sys
sys.path.append('/path/to/Lindel/')
sys.path.append('/path/to/indelphi_local/')
sys.path.append('/path/to/retarget/')
```
### 4. Download genome of interest for off-target predictions
REtarget uses FlashFry to predict the off-target activity of gRNAs. In order to use this functionality, it is required to generate a off-target database (see FlashFry-documentation for more information) before running REtarget. Depending on the genome used, this can take several hours. The step only needs to be executed once. Make sure to distribute sufficient memory. This step can be omitted if you do not want REtarget to check off-targets.
```
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O hg38.fa.gz
java -Xmx4g -jar ../flash_fry/FlashFry-assembly-1.12.jar \
 index \
 --tmpLocation ../flash_fry/tmp \
 --database hg38_cas9ngg_database \
 --reference path/to/data/hg38.fa.gz
```
It is necessary to specify the path to the corresponding java executable file in `config.py` (`JDK_PATH`).

---

## Local usage
### 1. Check specific genomic position for suitability for Recursive Editing
Example usage:
```
python recursive_editing/hdr_optimization.py -s gaatttagtctcccagcaggagaaataagagtaaataaatgtccagtggaatccacagcagacccaccctcaccttcagttttattgttttgctccaaacataactctgctgcttccactgctctggggctggtaaaaatgagtcccccgtaatcttcaggatgagaaagctgcacaccaaaaagcaataaagacatttt -c 100 -o rev
```
Possible arguments:
- Input sequence (`-s` or `--seq`, required): specify forward strand of sequence to be searched
- Initial cutsite of first gRNA (`-c` or `--cutsite`, required): specify cutsite (integer, 0-based)
- Orientation of input sequence and first gRNA (`-o` or `--orientation`, required): specify `fwd` if gRNA and search sequence have same orientation and `rev` if gRNA and search sequence have different orientation
- Name of output directory (`-n` or `--name`, not required): name of directory where results of analysis will be stored
- Name of FleshFry database (`-f` or `--database`, not required)
- Sequence of HDR template (`-t` or `--template`, not required)
- gRNA blacklist (`-b` or `--blacklist`, not required): list of gRNAs that should not be used by REtarget


### 2. Search larger sequence region for possible Recursive Editing sites
Example usage:
```
python recursive_editing/sequence_search.py -q gaatttagtctcccagcaggagaaataagagtaaataaatgtccagtggaatccacagcagacccaccctcaccttcagttttattgttttgctccaaacataactctgctgcttccactgctctggggctggtaaaaatgagtcccccgtaatcttcaggatgagaaagctgcacaccaaaaagcaataaagacatttt
```
Possible arguments:
- Input sequence (`-q` or `--seq`, required): specify forward strand of sequence to be searched
- Start position of search (`-s` or `--start`, not required): specify integer
- End position of search (`-e` or `--end`, not required): specify integer
- Name of output directory (`-n` or `--name`, not required): name of directory where results of analysis will be stored
- Name of FleshFry database (`-f` or `--database`, not required)
- Sequence of HDR template (`-t` or `--template`, not required)
- gRNA blacklist (`-b` or `--blacklist`, not required): list of gRNAs that should not be used by REtarget

---

## Search parameters
Parameters for running REtarget are specified in the file `recursive_editing/config.py`. Users can adjust the following parameters (variable names used in the python implementation referenced in brackets):
- **Maximal number of gRNAs** (``MAX_GUIDE_NUM``): the maximal number of gRNAs that REtarget will generate as output. Once the maximal number of gRNAs is reached, REtarget will compare the REtarget scores of newly designed and previously selected gRNAs and keep those gRNAs resulting in the highest total REtarget score. *Default: 5*
- **Maximal number of levels** (``MAX_LEVEL_NUM``): gRNAs in the same level require the same number of editing rounds for their target genotype to be generated. Level 1 comprises, by definition, only the first gRNA targeting the WT sequence (= 1st round of Editing). All editing outcomes resulting from gRNA1 give rise to potential level 2 gRNAs (= 2nd round of Editing). The maximal number of levels serves as a stopping criterion. REtarget will terminate after ``MAX_LEVEL_NUM`` levels if the optimization process has not been stopped before due to other criteria. *Default: 3*
- **Number of predicted editing outcomes considered per gRNA** (``MAX_GUIDES``): Specifies how many predicted indels resulting from a single editing event are considered in the optimization process. The value should be lower or equal to the maximal number of gRNAs. *Default: 3*
- **Number of NGG-PAMs considered to find new candidate gRNAs** (``MAX_PAM_NUM``): Each sufficiently abundant editing outcome will serve as a seed sequence for the generation of new candidate gRNAs. For this, a particular sequence window around the previous cut site will be searched for forward and reverse NGG-PAMs. This parameter specifies how many of the corresponding gRNAs should be further evaluated. Before this filter is applied, corresponding gRNAs are sorted according to their on-target cutting efficiency. Thus, only the most promising gRNAs are handed over to the next step. *Default: 5*
- **Minimal frequency of editing outcomes** (``SINGLE_FREQ_CUT``): All editing outcomes with lower predicted frequencies will be disregarded. Note that meaningful cutoffs can be markedly different depending on the prediction tool utilized. *Default: 5%*
- **Minimal level REtarget score** (``LEVEL_FREQ_CUT``): This parameter serves as the primary stopping criterion for the algorithm. Once the level REtarget score falls below this user-defined threshold, the optimization process will terminate and output the results. Note that REtarget scores for higher layers are usually very small. Thus, setting a high threshold may unintentionally result in early termination of the optimization process. *Default: 0.05*
- **Factor to down weight cutting scores of higher levels** (``LEVEL_WEIGHT_FACTOR``): A factor smaller than 1 will downweigh the cutting scores of higher levels (*cf.* Supplementary Notes). *Default: 1* (no down weighting)
- **Factor to adjust REtarget scores of higher layers** (``LEVEL_SCORE_FACTOR``): A factor smaller than 1 will down weight the Level REtarget scores of higher levels (*cf.* Supplementary Notes). *Default: 1* (no down weighting)
- **Model used to predict editing outcomes** (``MODEL``): 0 = Lindel (*default*), 1 = inDelphi. We recommend to use Lindel for running REtarget.
- **Cell type to initialize inDelphi** (``CELLTYPE_MODEL``): value will be disregarded if Lindel is selected as a model. *Default: mESC*
- **Minimal overlap of candidate gRNAs with the right sequence arm of the previous cut site** (``MIN_OVERLAP``): Ensures that gRNAs target unique sequences that are generated in the process of Recursive Editing, not lower-level editing outcomes or the WT-sequence. A value of 4, for instance, specifies that there is an overlap of at least 3 bp between the candidate gRNA and the right sequence arm of the previous cut site. Accordingly, a value of -1 means that the sequence "GG" of a potential PAM is centered around the previous cut site (corresponding to an overlap of -2 between candidate gRNA and the right sequence arm of the previous cut site. Specificity due to novel PAM-generation). *Default: -1*
- **Maximal overlap of candidate gRNAs with the right sequence arm of the previous cut site** (``MAX_OVERLAP``): A value of 10, for instance, means that no candidate gRNA would be selected with an overlap between candidate gRNA and the right sequence arm of the previous cut site larger than 10 bp. *Default: 10*
- **Minimal Level REtarget Score of Level 1** (``L1_CUT``). The search will terminate without returning results if the criterion is not met. *Default: 0.4 (recommended value for Lindel, please adjust value when using inDelphi)*
Minimal Level REtarget Score of Level 2 (``L2_CUT``). The search will terminate without returning results if the criterion is not met. *Default: 0.2 (recommended value for Lindel, please adjust value when using inDelphi)*
- **Minimal number of Levels** (``L_MIN_SAVE``): Results will not be stored if the number of resulting levels is smaller than this parameter. Not relevant for the online version of REtarget. *Default: 2*
- **Check for endogenous recognition sites in proximity to the initial cut site** (``CHECK_OFFT_PROX``): Defines if gRNAs in levels > 1 targeting the WT sequence in proximity to the initial cut site should be sorted out. gRNAs that target lower-level editing outcomes but not the WT sequence will not be sorted out as this would also reduce undesired editing outcomes and potentially increase HDR. *Default: True*
- **Check for genome-wide off-targets** (``CHECK_OFFT_GENOME``): defines if FlashFry should be used to search for off-targets (up to 4 mismatches). Not relevant for the online version of REtarget. *Default: True.* The name of the previously generated FlashFry database can be specified using the parameter ``FF_DATABASE_NAME``, the path to the required java executable using ``JDK_PATH``. FlashFry will be started utilizing the ``taskset`` command if ``USE_TASKSET_JKD`` is set to True (*default*).
- **Maximal number of off-targets with 0-4 mismatches** (``MAX_OFF_TAR``): threshold for off-target based gRNA filtering. Not relevant for the online version of REtarget. *Default: [0, 2, 20, 200, 2000]*
- **Filter gRNAs based on off-target scores** (``CHECK_OFFT_SCORES``). Available scores: MIT specificity score, Doench CFD score. For description of scores, please consult the FlashFry wiki. Not relevant for the online version of REtarget. *Default: True*
- **Minimal Doench 2014 on-target efficiency score** (``D2014ON_MIN``): gRNAs with a Doench score lower than this threshold will be disregarded. *Default: 0.05*
- **Maximal Doench CFD score for gRNA with highest off-target activity**(``DCFD_MAXOT_MAX``): gRNAs with a score higher than this threshold will be disregarded. *Default: 0.75*
- **Minimal Doench CFD off-target specificity score** (``DCFD_SPEC_MIN``): gRNAs with a Doench CFD score lower than this threshold will be disregarded. *Default: 0.5*
- **Minimal MIT specificity score** (``HSU2013_MIN``): gRNAs with an MIT specificity score lower than this threshold will be disregarded. *Default: 0.5*
- **Filter gRNAs based on GC-content** (``CHECK_GC``): *Default: False*
- **Minimal GC-content of gRNA** (``GC_LOW``): gRNAs with a GC-content lower than this threshold will be disregarded. *Default: 0.1*
- **Maximal GC-content of gRNA** (``GC_HIGH``): gRNAs with a GC-content higher than this threshold will be disregarded. *Default: 0.9*
- **Filter gRNAs based on homopolymers** (``CHECK_POLY_N``). *Default: False*
- **Check if gRNAs target given HDR template** (``CHECK_TEMPLATE``). *Default: False.* HDR template can be defined using the parameter ``HDR_TEMPLATE``.
- **Filter gRNAs based user-defined blacklist** (``CHECK_BLACKLIST``). *Default: False*. gRNA blacklist can be defined using the parameter ``BLACKLIST``.
- Specify if gRNAs of final level should be designed according to recursive optimization scheme (highest REtarget score) or based on their on-target efficiency (Doench 2014 score, ``FINAL_BY_D2014ON``): *Default: False (final layer not based on Doench 2014 score) - Online version default: True*
- **Logging mode** (`VERBOSE`): *Default: True*

---

## Format of results
REtarget generates three types of result files:
- **res_df.csv**: will be generated for each search and contains a single entry for each genomic site amenable to Recursive Editing. The columns contain the following information for each site: Id (identifier for the genomic site in the format id_[position relative to start position of input sequence]_[orientation of gRNA in relation of input sequence]), Number of levels, Number of guides, Guide RNAs (list of all gRNAs), Total REtarget score, REtarget score level 1, REtarget score level 2.
- **[id]_guide_df.csv**: will be generated for each site listed in **res_df.csv** and contains a single entry for each gRNA involved in Recursive Editing. The columns contain the following information for each gRNA: Guide_ID, Guide_RNA, Level, REtarget_Score, MGA_REtarget_Score (accumulated REtarget scores if single gRNA appears multiple times), Cutting_Score, Editing_Score, Doench2014_Score, DoenchCFD_MaxOT (OT: off-target), DoenchCFD_SpecScore (off-target score), Hsu2013_Score, Number_Off_Targets, Mismatches_Closest_OT,  Number_Closest_OT, Guide_Orientation, PAM, Target, Target_Type, Cut_Position, Target_Origin, Multi (specifies if gRNA appears multiple times). The corresponding files will be stored in a folder denoting the number of generated levels.
- **[id]_edit_df.csv**: will be generated for each site listed in **res_df.csv** and contains all predicted editing outcomes for each gRNA. The columns contain the following information for each retargetable indel: Category (insertion or deletion), Genotype position (start position of indel in relation to cut site), Inserted Bases, Length (length of indel), Predicted frequency (relative indel frequency in %), Genotype, REtarget_Score_Contribution, Seq_Type (orientation of genotype in relation to input sequence), Cutsite, Edit_ID, Guide_ID, Level. The corresponding files will be stored in a folder denoting the number of generated levels.
- In the online version of REtarget, users can also download a file containing the parameters applied for the search. Parameters are listed comma-separated in the following order:
``MAX_GUIDE_NUM``, ``MAX_LEVEL_NUM``, ``MAX_GUIDES``, ``MAX_PAM_NUM``, ``SINGLE_FREQ_CUT``, ``LEVEL_FREQ_CUT``, ``LEVEL_WEIGHT_FACTOR``, ``LEVEL_SCORE_FACTOR``, ``MIN_DSCORE``, ``GC_LOW``, ``GC_HIGH``, ``BLACKLIST``, ``MIN_OVERLAP``, ``MAX_OVERLAP``, ``L1_CUT``, ``L2_CUT``, ``MODEL``, ``CELLTYPE_MODEL``, ``CHECK_OFFT_PROX``, ``CHECK_GC``, ``CHECK_POLY_N``, ``CHECK_TTT``, ``CHECK_TEMPLATE``, ``CHECK_BLACKLIST``, ``FINAL_BY_D2014ON``

---

## Database of pre-computed REtarget results
### 1. Genome-wide search for loci amenable to Recursive Editing
We used REtarget to search the human genome (GRCh38, downloaded from Genbank) for sites amenable for Recursive Editing. A summary of all results with initial target sequences and genomic positions can be found in the file ``data/genome_search_summary.csv``. Due to the selected search parameter, this list is not necessarily exhaustive. Genomic sites not present in the dataset could be suitable for Recursive Editing, as well.

Search parameters applied for genome-wide search:
````
MAX_GUIDE_NUM = 10
MAX_LEVEL_NUM = 3
MAX_GUIDES = 3
MAX_PAM_NUM = 10
SINGLE_FREQ_CUT = 10
LEVEL_FREQ_CUT = 0.1
LEVEL_WEIGHT_FACTOR = 1
LEVEL_SCORE_FACTOR = 1
MODEL = 0
CELLTYPE_MODEL = 'mESC'
MIN_OVERLAP = -1
MAX_OVERLAP = 10
L1_CUT = 0.6
L2_CUT = 0.3
L_MIN_SAVE = 2
CHECK_OFFT_GENOME = True
CHECK_OFFT_PROX = True
MAX_OFF_TAR = [0, 2, 20, 200, 2000]
CHECK_OFFT_SCORES = False
D2014ON_MIN = 0.25
DCFD_MAXOT_MAX = 0.75
DCFD_SPEC_MIN = 0.5
HSU2013_MIN = 0.5
CHECK_GC = True
GC_LOW = 0.1
GC_HIGH = 0.9
CHECK_POLY_N = False
CHECK_TTT = False
FINAL_BY_D2014ON = False
````


### 2. Search of ClinVar database for loci amenable to Recursive Editing
We used REtarget to find the best Recursive Editing gRNA set for each of the 94,000+ annotated pathogenic mutations in ClinVar (version as of 11/2021), excluding indels >50 bp that are less ideal for ssODN donors and applying looser parameters than in the previous genome-wide search for globally optimal reagents. A summary of all results with initial target sequences and genomic positions can be found in the file ``data/clinvar_search_summary.csv``.

Search parameters applied for ClinVar search:
````
MAX_GUIDE_NUM=5
MAX_LEVEL_NUM=3
MAX_GUIDES=3
MAX_PAM_NUM=10
SINGLE_FREQ_CUT=10
LEVEL_FREQ_CUT=0.1
LEVEL_WEIGHT_FACTOR=1
LEVEL_SCORE_FACTOR=1
MODEL=0
CELLTYPE_MODEL='mESC'
MIN_OVERLAP = -1
MAX_OVERLAP = 10
L1_CUT = 0.35
L2_CUT = 0.15
L_MIN_SAVE = 2
CHECK_OFFT_GENOME = True
CHECK_OFFT_PROX = True
MAX_OFF_TAR = [0, 2, 20, 200, 2000]
CHECK_OFFT_SCORES = False
D2014ON_MIN = 0.05
DCFD_MAXOT_MAX = 0.75
DCFD_SPEC_MIN = 0.5
HSU2013_MIN = 0.5
CHECK_GC = False
GC_LOW = 0.1
GC_HIGH = 0.9
CHECK_POLY_N = False
CHECK_TTT = False
FINAL_BY_D2014ON = False
````


### 3. Genome-wide start and stop codon search for loci amenable to Recursive Editing
We used REtarget to search all start and stop codons of the human genome (GRCh38) for sites amenable for Recursive Editing. A summary of all results with initial target sequences and genomic positions can be found in the file ``data/codon_search_summary.csv``.

Search parameters applied for start and stop codon search:
````
MAX_GUIDE_NUM=5
MAX_LEVEL_NUM=3
MAX_GUIDES=3
MAX_PAM_NUM=10
SINGLE_FREQ_CUT=10
LEVEL_FREQ_CUT=0.1
LEVEL_WEIGHT_FACTOR=1
LEVEL_SCORE_FACTOR=1
MODEL=0
CELLTYPE_MODEL='mESC'
MIN_OVERLAP = -1
MAX_OVERLAP = 10
L1_CUT = 0.35
L2_CUT = 0.15
L_MIN_SAVE = 2
CHECK_OFFT_GENOME = True
CHECK_OFFT_PROX = True
MAX_OFF_TAR = [0, 2, 20, 200, 2000]
CHECK_OFFT_SCORES = False
D2014ON_MIN = 0.05
DCFD_MAXOT_MAX = 0.75
DCFD_SPEC_MIN = 0.5
HSU2013_MIN = 0.5
CHECK_GC = False
GC_LOW = 0.1
GC_HIGH = 0.9
CHECK_POLY_N = False
CHECK_TTT = False
FINAL_BY_D2014ON = False
````

### 4. Database generation method
Note that all databases were initially generated without using FlashFry for off-target prediction to limit the computational complexity. All WT-targeting gRNAs were subsequently compiled and filtered according to their off-targets, yielding a reduced dataset, from which we compiled the summary files listed above.

---