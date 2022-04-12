#######################################################################################################################
#
# Set parameters for Recursive Editing search
#
#######################################################################################################################

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
L1_CUT = 0.3
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
CHECK_TEMPLATE = False
CHECK_BLACKLIST = False
BLACKLIST = []
HDR_TEMPLATE = ''
FINAL_BY_D2014ON = False
FF_DATABASE_NAME = 'hg38_cas9_database'
JDK_PATH = '/home/cornlab/miniconda2/envs/jdk8/bin/java'
USE_TASKSET_JKD = True
VERBOSE = False