#include "AminoAcidProbs.h"
#include "PrmGraph.h"
#include "PeptideRankScorer.h"
#include "DeNovoSolutions.h"
#include "PepNovo_auxfun.h"

const char * ScoreModelFields_SAA_names[]={
"SAA_CONST",	"SAA_IND_SEQ_RANK_1",	"SAA_IND_SEQ_RANK_2",	"SAA_IND_SEQ_RANK_3",	"SAA_IND_SEQ_RANK_4",	"SAA_IND_SEQ_RANK_5",	"SAA_IND_SEQ_RANK_6",	"SAA_IND_SEQ_RANK_8",	"SAA_IND_SEQ_RANK_10",	"SAA_IND_SEQ_RANK_12",	"SAA_IND_SEQ_RANK_14",	"SAA_IND_SEQ_RANK_16",	"SAA_IND_SEQ_RANK_20",	"SAA_IND_SEQ_RANK_MORE",	"SAA_IND_AA_POS_20",	"SAA_IND_AA_POS_40",	"SAA_IND_AA_POS_60",	"SAA_IND_AA_POS_80",	
"SAA_IND_AA_POS_100",	"SAA_IND_AA_REL_SCORE_20",	"SAA_IND_AA_REL_SCORE_40",	"SAA_IND_AA_REL_SCORE_60",	"SAA_IND_AA_REL_SCORE_80",	"SAA_IND_AA_REL_SCORE_100",	"SAA_MAX_NODE_SCORE",	"SAA_MIN_NODE_SCORE",	"SAA_N_SCORE_ABOVE_TWO_THIRDS",	"SAA_C_SCORE_ABOVE_TWO_THIRDS",	"SAA_IND_TWO_THIRDS_BOTH_ABOVE",	"SAA_N_SCORE_ABOVE_THIRD",	"SAA_C_SCORE_ABOVE_THIRD",	"SAA_IND_THIRD_BOTH_ABOVE",	
"SAA_N_SCORE_ABOVE_ZERO",	"SAA_C_SCORE_ABOVE_ZERO",	"SAA_IND_ZERO_BOTH_ABOVE",	"SAA_MAX_SCORE_RANK",	"SAA_MIN_SCORE_RANK",	"SAA_N_SCORE_RANK",	"SAA_C_SCORE_RANK",	"SAA_SCORE_RANK_SUM",	"SAA_SCORE_RANK_DIFF", "SAA_SCORE_RANK_ABS_DIFF",	"SAA_N_NUM_FRAGS",	"SAA_C_NUM_FRAGS",	"SAA_NUM_FRAG_DIFF",	"SAA_ABS_NUM_FRAG_DIFF",	"SAA_IND_N_STRONGER_INTEN",	"SAA_N_STRONGER_LOG_NC_INTEN_RATIO",	
"SAA_IND_C_STRONGER_INTEN",	"SAA_C_STRONGER_LOG_NC_INTEN_RATIO",	"SAA_IND_N_IS_MAX_IDX_TO_C",	"SAA_IND_C_IS_MAX_IDX_FROM_N",	"SAA_DIFF_N_MAX_IN_C_SCORE_RANKS",	"SAA_IND_N_NOT_MAX_IDX_TO_C",	"SAA_IND_C_NOT_MAX_IDX_FROM_N",	"SAA_DIFF_C_MAX_OUT_N_SCORE_RANKS",	"SAA_IND_BOTH_CONNECT_TO_MAX",	"SAA_NODE_MASS_DIFF",	"SAA_NODE_SQR_MASS_DIFF",	"SAA_NUM_FRAG_PAIRS",	"SAA_AVG_PEAK_DIFF",	
"SAA_AVG_PEAK_SQR_DIFF",	"SAA_BEST_PEAK_DIFF",	"SAA_BEST_PEAK_SQR_DIFF",	"SAA_IND_NO_PEAK_DIFF",	"SAA_AVG_PEAK_DIFF_TIMES_SCORE_RANK_SUM",	"SAA_AVG_PEAK_DIFF_TIMES_SCORE_ABS_DIFF",	"SAA_AVG_PEAK_DIFF_DIV_NUM_FRAG_PAIRS",	"SAA_IND_BOTH_CONNECT_TO_MAX_TIMES_AVG_DIFF",	"SAA_IND_N_SCORE_MORE_THAN_MIRROR",	"SAA_IND_C_SCORE_MORE_THAN_MIRROR",	"SAA_IND_BOTH_SCORE_MORE_THAN_MIRROR",	
"SAA_LOG_DIFF_MORE_THAN_MIRROR",	"SAA_LOG_DIFF_LESS_THAN_MIRROR",	"SAA_NO_MIRROR",	"SAA_IND_AA_N_TERM",	"SAA_IND_AA_C_TERM",	"SAA_IND_AA_Gap",	"SAA_IND_AA_Xle",	"SAA_IND_AA_Ala",	"SAA_IND_AA_Arg",	"SAA_IND_AA_Asn",	"SAA_IND_AA_Asp",	"SAA_IND_AA_Cys",	"SAA_IND_AA_Gln",	"SAA_IND_AA_Glu",	"SAA_IND_AA_Gly",	"SAA_IND_AA_His",	"SAA_IND_AA_Ile",	"SAA_IND_AA_Leu",	"SAA_IND_AA_Lys",	"SAA_IND_AA_Met",	
"SAA_IND_AA_Phe",	"SAA_IND_AA_Pro",	"SAA_IND_AA_Ser",	"SAA_IND_AA_Thr",	"SAA_IND_AA_Trp",	"SAA_IND_AA_Tyr",	"SAA_IND_AA_Val",	"SAA_SCORE_RANK_SUM_N_TERM",	"SAA_SCORE_RANK_SUM_C_TERM",	"SAA_SCORE_RANK_SUM_Gap",	"SAA_SCORE_RANK_SUM_Xle",	"SAA_SCORE_RANK_SUM_Ala",	"SAA_SCORE_RANK_SUM_Arg",	"SAA_SCORE_RANK_SUM_Asn",	"SAA_SCORE_RANK_SUM_Asp",	"SAA_SCORE_RANK_SUM_Cys",	"SAA_SCORE_RANK_SUM_Gln",	
"SAA_SCORE_RANK_SUM_Glu",	"SAA_SCORE_RANK_SUM_Gly",	"SAA_SCORE_RANK_SUM_His",	"SAA_SCORE_RANK_SUM_Ile",	"SAA_SCORE_RANK_SUM_Leu",	"SAA_SCORE_RANK_SUM_Lys",	"SAA_SCORE_RANK_SUM_Met",	"SAA_SCORE_RANK_SUM_Phe",	"SAA_SCORE_RANK_SUM_Pro",	"SAA_SCORE_RANK_SUM_Ser",	"SAA_SCORE_RANK_SUM_Thr",	"SAA_SCORE_RANK_SUM_Trp",	"SAA_SCORE_RANK_SUM_Tyr",	"SAA_SCORE_RANK_SUM_Val",	"SAA_SCORE_RANK_DIFF_N_TERM",	
"SAA_SCORE_RANK_DIFF_C_TERM",	"SAA_SCORE_RANK_DIFF_Gap",	"SAA_SCORE_RANK_DIFF_Xle",	"SAA_SCORE_RANK_DIFF_Ala",	"SAA_SCORE_RANK_DIFF_Arg",	"SAA_SCORE_RANK_DIFF_Asn",	"SAA_SCORE_RANK_DIFF_Asp",	"SAA_SCORE_RANK_DIFF_Cys",	"SAA_SCORE_RANK_DIFF_Gln",	"SAA_SCORE_RANK_DIFF_Glu",	"SAA_SCORE_RANK_DIFF_Gly",	"SAA_SCORE_RANK_DIFF_His",	"SAA_SCORE_RANK_DIFF_Ile",	"SAA_SCORE_RANK_DIFF_Leu",	
"SAA_SCORE_RANK_DIFF_Lys",	"SAA_SCORE_RANK_DIFF_Met",	"SAA_SCORE_RANK_DIFF_Phe",	"SAA_SCORE_RANK_DIFF_Pro",	"SAA_SCORE_RANK_DIFF_Ser",	"SAA_SCORE_RANK_DIFF_Thr",	"SAA_SCORE_RANK_DIFF_Trp",	"SAA_SCORE_RANK_DIFF_Tyr",	"SAA_SCORE_RANK_DIFF_Val"};

const char * ScoreModelFields_SAANCD_names[]={
"SAANCD_CONST",	"SAANCD_IND_SEQ_RANK_1",	"SAANCD_IND_SEQ_RANK_2",	"SAANCD_IND_SEQ_RANK_3",	"SAANCD_IND_SEQ_RANK_4",	"SAANCD_IND_SEQ_RANK_5",	"SAANCD_IND_SEQ_RANK_6",	"SAANCD_IND_SEQ_RANK_8",	"SAANCD_IND_SEQ_RANK_10",	"SAANCD_IND_SEQ_RANK_12",	"SAANCD_IND_SEQ_RANK_14",	"SAANCD_IND_SEQ_RANK_16",	"SAANCD_IND_SEQ_RANK_20",	"SAANCD_IND_SEQ_RANK_MORE",	"SAANCD_IND_CONNECTS_TO_NODE_WITH_NO_FRAGS",	
"SAANCD_IND_CONNECTS_TO_N_TERMINAL",	"SAANCD_IND_CONNECTS_TO_C_TERMINAL",	"SAANCD_IND_CONNECTS_TO_DIGEST",	"SAANCD_IND_HAS_MAX_SCORE_FROM_N",	"SAANCD_SCORE_FROM_N",	"SAANCD_IND_HAS_MAX_SCORE_TO_C",	"SAANCD_SCORE_TO_C",	"SAANCD_IND_HAS_MAX_SCORE_TO_DIGEST",	"SAANCD_SCORE_TO_DIGEST",	"SAANCD_NODE_MASS_DIFF_FROM_N",	"SAANCD_NODE_SQR_MASS_DIFF_FROM_N",	"SAANCD_NODE_MASS_DIFF_TO_C",	
"SAANCD_NODE_SQR_MASS_DIFF_TO_C",	"SAANCD_NODE_MASS_DIFF_TO_DIGEST",	"SAANCD_NODE_SQR_MASS_DIFF_TO_DIGEST",	"SAANCD_IND_TERM_TO_DIGEST",	"SAANCD_IND_TERM_TO_DIGEST_WITH_FRAGS",	"SAANCD_IND_TERM_TO_DIGEST_W_POS_SCORE",	"SAANCD_IND_NT_AA_N_TERM",	"SAANCD_IND_NT_AA_C_TERM",	"SAANCD_IND_NT_AA_Gap",	"SAANCD_IND_NT_AA_Xle",	"SAANCD_IND_NT_AA_Ala",	"SAANCD_IND_NT_AA_Arg",	"SAANCD_IND_NT_AA_Asn",	
"SAANCD_IND_NT_AA_Asp",	"SAANCD_IND_NT_AA_Cys",	"SAANCD_IND_NT_AA_Gln",	"SAANCD_IND_NT_AA_Glu",	"SAANCD_IND_NT_AA_Gly",	"SAANCD_IND_NT_AA_His",	"SAANCD_IND_NT_AA_Ile",	"SAANCD_IND_NT_AA_Leu",	"SAANCD_IND_NT_AA_Lys",	"SAANCD_IND_NT_AA_Met",	"SAANCD_IND_NT_AA_Phe",	"SAANCD_IND_NT_AA_Pro",	"SAANCD_IND_NT_AA_Ser",	"SAANCD_IND_NT_AA_Thr",	"SAANCD_IND_NT_AA_Trp",	"SAANCD_IND_NT_AA_Tyr",	
"SAANCD_IND_NT_AA_Val"};

const char * ScoreModelFields_DAA_names[]={
"DAA_CONST",	"DAA_IND_SEQ_RANK_1",	"DAA_IND_SEQ_RANK_2",	"DAA_IND_SEQ_RANK_3",	"DAA_IND_SEQ_RANK_4",	"DAA_IND_SEQ_RANK_5",	"DAA_IND_SEQ_RANK_6",	"DAA_IND_SEQ_RANK_8",	"DAA_IND_SEQ_RANK_10",	"DAA_IND_SEQ_RANK_12",	"DAA_IND_SEQ_RANK_14",	"DAA_IND_SEQ_RANK_16",	"DAA_IND_SEQ_RANK_20",	"DAA_IND_SEQ_RANK_MORE",	"DAA_IND_AA_POS_20",	"DAA_IND_AA_POS_40",	"DAA_IND_AA_POS_60",	"DAA_IND_AA_POS_80",	
"DAA_IND_AA_POS_100",	"DAA_IND_AA_REL_SCORE_20",	"DAA_IND_AA_REL_SCORE_40",	"DAA_IND_AA_REL_SCORE_60",	"DAA_IND_AA_REL_SCORE_80",	"DAA_IND_AA_REL_SCORE_100",	"DAA_MAX_NODE_SCORE",	"DAA_MIN_NODE_SCORE",	"DAA_N_SCORE_ABOVE_TWO_THIRDS",	"DAA_C_SCORE_ABOVE_TWO_THIRDS",	"DAA_IND_TWO_THIRDS_BOTH_ABOVE",	"DAA_N_SCORE_ABOVE_THIRD",	"DAA_C_SCORE_ABOVE_THIRD",	"DAA_IND_THIRD_BOTH_ABOVE",	
"DAA_N_SCORE_ABOVE_ZERO",	"DAA_C_SCORE_ABOVE_ZERO",	"DAA_IND_ZERO_BOTH_ABOVE",	"DAA_IND_BOTH_ABOVE_FIVE",	"DAA_IND_N_ABOVE_FIVE",	"DAA_IND_C_ABOVE_FIVE",	"DAA_MAX_SCORE_RANK",	"DAA_MIN_SCORE_RANK",	"DAA_N_SCORE_RANK",	"DAA_C_SCORE_RANK",	"DAA_SCORE_RANK_SUM",	"DAA_SCORE_RANK_DIFF",	"DAA_SCORE_RANK_ABS_DIFF",	"DAA_N_NUM_FRAGS",	"DAA_C_NUM_FRAGS",	"DAA_NUM_FRAG_DIFF",	"DAA_ABS_NUM_FRAG_DIFF",	
"DAA_IND_N_IS_MAX_IDX_TO_C",	"DAA_IND_C_IS_MAX_IDX_FROM_N",	"DAA_DIFF_N_MAX_IN_C_SCORE_RANKS",	"DAA_IND_N_NOT_MAX_IDX_TO_C",	"DAA_IND_C_NOT_MAX_IDX_FROM_N",	"DAA_DIFF_C_MAX_OUT_N_SCORE_RANKS",	"DAA_IND_BOTH_CONNECT_TO_MAX",	"DAA_NODE_MASS_DIFF",	"DAA_NODE_SQR_MASS_DIFF",	"DAA_NUM_FRAG_PAIRS",	"DAA_AVG_PEAK_DIFF",	"DAA_AVG_PEAK_SQR_DIFF",	"DAA_BEST_PEAK_DIFF",	"DAA_BEST_PEAK_SQR_DIFF",	
"DAA_IND_NO_PEAK_DIFF",	"DAA_AVG_PEAK_DIFF_TIMES_SCORE_RANK_SUM",	"DAA_AVG_PEAK_DIFF_TIMES_SCORE_ABS_DIFF",	"DAA_AVG_PEAK_DIFF_DIV_NUM_FRAG_PAIRS",	"DAA_IND_BOTH_CONNECT_TO_MAX_TIMES_AVG_DIFF",	"DAA_NUM_DOUBLE_EDGE_ROUTES",	"DAA_IND_NO_DOUBLE_EDGE_ROUTES",	"DAA_SCORE_RANK_DOUBLE_EDGE_ROUTES",	"DAA_IND_IND_MAX_SCORE_ALTERNATE_MORE_THAN_ZERO",	"DAA_IND_MAX_ALTERNATE_IS_MAX_OUT_N",	
"DAA_IND_MAX_ALTERNATE_IS_MAX_IN_C",	"DAA_NODE_OFFSETS_ALTERNTE",	"DAA_SQR_NODE_ODFFSETS_ALTERNATE",	"DAA_IND_N_SCORE_MORE_THAN_MIRROR",	"DAA_IND_C_SCORE_MORE_THAN_MIRROR",	"DAA_IND_BOTH_SCORE_MORE_THAN_MIRROR",	"DAA_LOG_DIFF_MORE_THAN_MIRROR",	"DAA_LOG_DIFF_LESS_THAN_MIRROR",	"DAA_NO_MIRROR",	"DAA_IND_PROBLEMATIC_PAIR_OF_AAS",	"DAA_IND_N_AA_N_TERM",	"DAA_IND_N_AA_C_TERM",	"DAA_IND_N_AA_Gap",	
"DAA_IND_N_AA_Xle",	"DAA_IND_N_AA_Ala",	"DAA_IND_N_AA_Arg",	"DAA_IND_N_AA_Asn",	"DAA_IND_N_AA_Asp",	"DAA_IND_N_AA_Cys",	"DAA_IND_N_AA_Gln",	"DAA_IND_N_AA_Glu",	"DAA_IND_N_AA_Gly",	"DAA_IND_N_AA_His",	"DAA_IND_N_AA_Ile",	"DAA_IND_N_AA_Leu",	"DAA_IND_N_AA_Lys",	"DAA_IND_N_AA_Met",	"DAA_IND_N_AA_Phe",	"DAA_IND_N_AA_Pro",	"DAA_IND_N_AA_Ser",	"DAA_IND_N_AA_Thr",	"DAA_IND_N_AA_Trp",	
"DAA_IND_N_AA_Tyr",	"DAA_IND_N_AA_Val",	"DAA_SCORE_RANK_SUM_N_N_TERM",	"DAA_SCORE_RANK_SUM_N_C_TERM",	"DAA_SCORE_RANK_SUM_N_Gap",	"DAA_SCORE_RANK_SUM_N_Xle",	"DAA_SCORE_RANK_SUM_N_Ala",	"DAA_SCORE_RANK_SUM_N_Arg",	"DAA_SCORE_RANK_SUM_N_Asn",	"DAA_SCORE_RANK_SUM_N_Asp",	"DAA_SCORE_RANK_SUM_N_Cys",	"DAA_SCORE_RANK_SUM_N_Gln",	"DAA_SCORE_RANK_SUM_N_Glu",	"DAA_SCORE_RANK_SUM_N_Gly",	
"DAA_SCORE_RANK_SUM_N_His",	"DAA_SCORE_RANK_SUM_N_Ile",	"DAA_SCORE_RANK_SUM_N_Leu",	"DAA_SCORE_RANK_SUM_N_Lys",	"DAA_SCORE_RANK_SUM_N_Met",	"DAA_SCORE_RANK_SUM_N_Phe",	"DAA_SCORE_RANK_SUM_N_Pro",	"DAA_SCORE_RANK_SUM_N_Ser",	"DAA_SCORE_RANK_SUM_N_Thr",	"DAA_SCORE_RANK_SUM_N_Trp",	"DAA_SCORE_RANK_SUM_N_Tyr",	"DAA_SCORE_RANK_SUM_N_Val",	"DAA_SCORE_RANK_DIFF_N_N_TERM",	"DAA_SCORE_RANK_DIFF_N_C_TERM",	
"DAA_SCORE_RANK_DIFF_N_Gap",	"DAA_SCORE_RANK_DIFF_N_Xle",	"DAA_SCORE_RANK_DIFF_N_Ala",	"DAA_SCORE_RANK_DIFF_N_Arg",	"DAA_SCORE_RANK_DIFF_N_Asn",	"DAA_SCORE_RANK_DIFF_N_Asp",	"DAA_SCORE_RANK_DIFF_N_Cys",	"DAA_SCORE_RANK_DIFF_N_Gln",	"DAA_SCORE_RANK_DIFF_N_Glu",	"DAA_SCORE_RANK_DIFF_N_Gly",	"DAA_SCORE_RANK_DIFF_N_His",	"DAA_SCORE_RANK_DIFF_N_Ile",	"DAA_SCORE_RANK_DIFF_N_Leu",	"DAA_SCORE_RANK_DIFF_N_Lys",	
"DAA_SCORE_RANK_DIFF_N_Met",	"DAA_SCORE_RANK_DIFF_N_Phe",	"DAA_SCORE_RANK_DIFF_N_Pro",	"DAA_SCORE_RANK_DIFF_N_Ser",	"DAA_SCORE_RANK_DIFF_N_Thr",	"DAA_SCORE_RANK_DIFF_N_Trp",	"DAA_SCORE_RANK_DIFF_N_Tyr",	"DAA_SCORE_RANK_DIFF_N_Val",	"DAA_IND_C_AA_N_TERM",	"DAA_IND_C_AA_C_TERM",	"DAA_IND_C_AA_Gap",	"DAA_IND_C_AA_Xle",	"DAA_IND_C_AA_Ala",	"DAA_IND_C_AA_Arg",	"DAA_IND_C_AA_Asn",	
"DAA_IND_C_AA_Asp",	"DAA_IND_C_AA_Cys",	"DAA_IND_C_AA_Gln",	"DAA_IND_C_AA_Glu",	"DAA_IND_C_AA_Gly",	"DAA_IND_C_AA_His",	"DAA_IND_C_AA_Ile",	"DAA_IND_C_AA_Leu",	"DAA_IND_C_AA_Lys",	"DAA_IND_C_AA_Met",	"DAA_IND_C_AA_Phe",	"DAA_IND_C_AA_Pro",	"DAA_IND_C_AA_Ser",	"DAA_IND_C_AA_Thr",	"DAA_IND_C_AA_Trp",	"DAA_IND_C_AA_Tyr",	"DAA_IND_C_AA_Val",	"DAA_SCORE_RANK_SUM_C_N_TERM",	"DAA_SCORE_RANK_SUM_C_C_TERM",	
"DAA_SCORE_RANK_SUM_C_Gap",	"DAA_SCORE_RANK_SUM_C_Xle",	"DAA_SCORE_RANK_SUM_C_Ala",	"DAA_SCORE_RANK_SUM_C_Arg",	"DAA_SCORE_RANK_SUM_C_Asn",	"DAA_SCORE_RANK_SUM_C_Asp",	"DAA_SCORE_RANK_SUM_C_Cys",	"DAA_SCORE_RANK_SUM_C_Gln",	"DAA_SCORE_RANK_SUM_C_Glu",	"DAA_SCORE_RANK_SUM_C_Gly",	"DAA_SCORE_RANK_SUM_C_His",	"DAA_SCORE_RANK_SUM_C_Ile",	"DAA_SCORE_RANK_SUM_C_Leu",	"DAA_SCORE_RANK_SUM_C_Lys",	
"DAA_SCORE_RANK_SUM_C_Met",	"DAA_SCORE_RANK_SUM_C_Phe",	"DAA_SCORE_RANK_SUM_C_Pro",	"DAA_SCORE_RANK_SUM_C_Ser",	"DAA_SCORE_RANK_SUM_C_Thr",	"DAA_SCORE_RANK_SUM_C_Trp",	"DAA_SCORE_RANK_SUM_C_Tyr",	"DAA_SCORE_RANK_SUM_C_Val",	"DAA_SCORE_RANK_DIFF_C_N_TERM",	"DAA_SCORE_RANK_DIFF_C_C_TERM",	"DAA_SCORE_RANK_DIFF_C_Gap",	"DAA_SCORE_RANK_DIFF_C_Xle",	"DAA_SCORE_RANK_DIFF_C_Ala",	"DAA_SCORE_RANK_DIFF_C_Arg",	
"DAA_SCORE_RANK_DIFF_C_Asn",	"DAA_SCORE_RANK_DIFF_C_Asp",	"DAA_SCORE_RANK_DIFF_C_Cys",	"DAA_SCORE_RANK_DIFF_C_Gln",	"DAA_SCORE_RANK_DIFF_C_Glu",	"DAA_SCORE_RANK_DIFF_C_Gly",	"DAA_SCORE_RANK_DIFF_C_His",	"DAA_SCORE_RANK_DIFF_C_Ile",	"DAA_SCORE_RANK_DIFF_C_Leu",	"DAA_SCORE_RANK_DIFF_C_Lys",	"DAA_SCORE_RANK_DIFF_C_Met",	"DAA_SCORE_RANK_DIFF_C_Phe",	"DAA_SCORE_RANK_DIFF_C_Pro",	"DAA_SCORE_RANK_DIFF_C_Ser",	
"DAA_SCORE_RANK_DIFF_C_Thr",	"DAA_SCORE_RANK_DIFF_C_Trp",	"DAA_SCORE_RANK_DIFF_C_Tyr",	"DAA_SCORE_RANK_DIFF_C_Val"};

const char * ScoreModelFields_DAANCD_names[]={
"DAANCD_CONST",	"DAANCD_IND_SEQ_RANK_1",	"DAANCD_IND_SEQ_RANK_2",	"DAANCD_IND_SEQ_RANK_3",	"DAANCD_IND_SEQ_RANK_4",	"DAANCD_IND_SEQ_RANK_5",	"DAANCD_IND_SEQ_RANK_6",	"DAANCD_IND_SEQ_RANK_8",	"DAANCD_IND_SEQ_RANK_10",	"DAANCD_IND_SEQ_RANK_12",	"DAANCD_IND_SEQ_RANK_14",	"DAANCD_IND_SEQ_RANK_16",	"DAANCD_IND_SEQ_RANK_20",	"DAANCD_IND_SEQ_RANK_MORE",	"DAANCD_IND_CONNECTS_TO_NODE_WITH_NO_FRAGS",	
"DAANCD_IND_CONNECTS_TO_N_TERMINAL",	"DAANCD_IND_CONNECTS_TO_C_TERMINAL",	"DAANCD_IND_CONNECTS_TO_DIGEST",	"DAANCD_IND_HAS_MAX_SCORE_FROM_N",	"DAANCD_SCORE_FROM_N",	"DAANCD_IND_HAS_MAX_SCORE_TO_C",	"DAANCD_SCORE_TO_C",	"DAANCD_IND_HAS_MAX_SCORE_TO_DIGEST",	"DAANCD_SCORE_TO_DIGEST",	"DAANCD_NODE_MASS_DIFF_FROM_N",	"DAANCD_NODE_SQR_MASS_DIFF_FROM_N",	"DAANCD_NODE_MASS_DIFF_TO_C",	
"DAANCD_NODE_SQR_MASS_DIFF_TO_C",	"DAANCD_NODE_MASS_DIFF_TO_DIGEST",	"DAANCD_NODE_SQR_MASS_DIFF_TO_DIGEST",	"DAANCD_NUM_DOUBLE_EDGE_ROUTES",	"DAANCD_IND_NO_DOUBLE_EDGE_ROUTES",	"DAANCD_SCORE_RANK_DOUBLE_EDGE_ROUTES",	"DAANCD_IND_IND_MAX_SCORE_ALTERNATE_MORE_THAN_ZERO",	"DAANCD_IND_MAX_ALTERNATE_IS_MAX_OUT_N",	"DAANCD_IND_MAX_ALTERNATE_IS_MAX_IN_C",	"DAANCD_NODE_OFFSETS_ALTERNTE",	
"DAANCD_SQR_NODE_ODFFSETS_ALTERNATE",	"DAANCD_IND_DIGEST_AA_IS_GOOD",	"DAANCD_IND_NT_N_AA_N_TERM",	"DAANCD_IND_NT_N_AA_C_TERM",	"DAANCD_IND_NT_N_AA_Gap",	"DAANCD_IND_NT_N_AA_Xle",	"DAANCD_IND_NT_N_AA_Ala",	"DAANCD_IND_NT_N_AA_Arg",	"DAANCD_IND_NT_N_AA_Asn",	"DAANCD_IND_NT_N_AA_Asp",	"DAANCD_IND_NT_N_AA_Cys",	"DAANCD_IND_NT_N_AA_Gln",	"DAANCD_IND_NT_N_AA_Glu",	"DAANCD_IND_NT_N_AA_Gly",	
"DAANCD_IND_NT_N_AA_His",	"DAANCD_IND_NT_N_AA_Ile",	"DAANCD_IND_NT_N_AA_Leu",	"DAANCD_IND_NT_N_AA_Lys",	"DAANCD_IND_NT_N_AA_Met",	"DAANCD_IND_NT_N_AA_Phe",	"DAANCD_IND_NT_N_AA_Pro",	"DAANCD_IND_NT_N_AA_Ser",	"DAANCD_IND_NT_N_AA_Thr",	"DAANCD_IND_NT_N_AA_Trp",	"DAANCD_IND_NT_N_AA_Tyr",	"DAANCD_IND_NT_N_AA_Val",	"DAANCD_IND_NT_C_AA_N_TERM",	"DAANCD_IND_NT_C_AA_C_TERM",	"DAANCD_IND_NT_C_AA_Gap",	
"DAANCD_IND_NT_C_AA_Xle",	"DAANCD_IND_NT_C_AA_Ala",	"DAANCD_IND_NT_C_AA_Arg",	"DAANCD_IND_NT_C_AA_Asn",	"DAANCD_IND_NT_C_AA_Asp",	"DAANCD_IND_NT_C_AA_Cys",	"DAANCD_IND_NT_C_AA_Gln",	"DAANCD_IND_NT_C_AA_Glu",	"DAANCD_IND_NT_C_AA_Gly",	"DAANCD_IND_NT_C_AA_His",	"DAANCD_IND_NT_C_AA_Ile",	"DAANCD_IND_NT_C_AA_Leu",	"DAANCD_IND_NT_C_AA_Lys",	"DAANCD_IND_NT_C_AA_Met",	"DAANCD_IND_NT_C_AA_Phe",	
"DAANCD_IND_NT_C_AA_Pro",	"DAANCD_IND_NT_C_AA_Ser",	"DAANCD_IND_NT_C_AA_Thr",	"DAANCD_IND_NT_C_AA_Trp",	"DAANCD_IND_NT_C_AA_Tyr",	"DAANCD_IND_NT_C_AA_Val"};

const char * AAProbTypes_names[]={
"AAP_SAA",	"AAP_SAANCD",	"AAP_DAA",	"AAP_DAANCD",	"AAP_NUM_TYPES",	"AAProbTypes"};



const int rank_threshes[]={1,2,3,4,5,6,8,10,12,14,16,20,POS_INF};
const int num_ranks = sizeof(rank_threshes)/sizeof(int);
int calc_rank_idx(int r)
{
	int i;
	for (i=0; i<num_ranks; i++)
		if (r<rank_threshes[i])
			break;
	return i;
}


void fill_fval_vector_for_saa(const PrmGraph* prm_ptr, 
							  const int me_idx, 
							  const int* variant_ptr, 
							  const int  seq_rank,
							  ME_Regression_Sample& sam)
{
	const Config *config = prm_ptr->get_config();
	const vector<mass_t>& aa2mass =config->get_aa2mass();
	const vector<int>& org_aa = config->get_org_aa();
	const vector<int>& forbidden_node_idxs = prm_ptr->get_forbidden_node_idxs();
	const vector<score_t>& cummulative_node_scores = prm_ptr->get_cummulative_scores();
	const score_t total_positive_score = cummulative_node_scores[cummulative_node_scores.size()-1];

	if (*variant_ptr != 1)
	{
		cout << "Error: using a single aa score for a variant with " << *variant_ptr << " aas." << endl;
		exit(1);
	}

	const int var_aa = *(variant_ptr+1); // this is the variant aa (the amino acid of the edge)

	const mass_t exp_mass = aa2mass[var_aa];
	const score_t third_score = (prm_ptr->get_max_node_score()<0) ? 0 : prm_ptr->get_max_node_score() * 0.3333;
	const score_t two_thirds_score = (prm_ptr->get_max_node_score()<0) ? 0 : prm_ptr->get_max_node_score() * 0.66666;

	
	const MultiEdge& me = prm_ptr->get_multi_edge(me_idx);
	const Node& n_node  = prm_ptr->get_node(me.n_idx);
	const Node& c_node  = prm_ptr->get_node(me.c_idx);
	const int rank_offset = calc_rank_idx(seq_rank);

	const int rel_mass_pos  = (int)((n_node.mass / prm_ptr->get_pm_with_19()) * 5.0);
	const int rel_score_pos =  (int)((cummulative_node_scores[me.n_idx]/total_positive_score) * 5.0);

	score_t n_node_score = n_node.score;
	score_t c_node_score = c_node.score;
	if (n_node_score<-30)
		n_node_score = -30;
	if (n_node_score>30)
		n_node_score=30;
	if (c_node_score<-30)
		c_node_score = -30;
	if (c_node_score>30)
		c_node_score=30;

	const bool n_above_two_thirds = (n_node.score >= two_thirds_score);
	const bool c_above_two_thirds = (c_node.score >= two_thirds_score);
	const bool n_above_third = (n_node.score >= third_score);
	const bool c_above_third = (c_node.score >= third_score);
	const bool n_above_zero = (n_node.score >= 0);
	const bool c_above_zero = (c_node.score >= 0);

	vector<fval>& fvals = sam.f_vals;
	fvals.clear();

	sam.add_feature(SAA_CONST,1.0);
	sam.add_feature(SAA_IND_SEQ_RANK_1+rank_offset,1.0);
	sam.add_feature(SAA_IND_AA_POS_20+rel_mass_pos,1.0);
	sam.add_feature(SAA_IND_AA_REL_SCORE_20+rel_score_pos,1.0);

	if (n_node_score>c_node_score)
	{
		sam.add_feature(SAA_MAX_NODE_SCORE,n_node_score);
		sam.add_feature(SAA_MIN_NODE_SCORE,c_node_score);
	}
	else
	{
		sam.add_feature(SAA_MAX_NODE_SCORE,c_node_score);
		sam.add_feature(SAA_MIN_NODE_SCORE,n_node_score);
	}

	if (n_above_two_thirds)
		fvals.push_back(fval(SAA_N_SCORE_ABOVE_TWO_THIRDS,1));

	if (c_above_two_thirds)
		fvals.push_back(fval(SAA_C_SCORE_ABOVE_TWO_THIRDS,1));

	if (n_above_two_thirds && c_above_two_thirds)
		fvals.push_back(fval(SAA_IND_TWO_THIRDS_BOTH_ABOVE,1));

	if (! n_above_two_thirds && n_above_third)
		fvals.push_back(fval(SAA_N_SCORE_ABOVE_THIRD, 1 ));

	if (! c_above_two_thirds && c_above_third)
		fvals.push_back(fval(SAA_C_SCORE_ABOVE_THIRD, 1 ));

	if (! (n_above_two_thirds && c_above_two_thirds) && n_above_third && c_above_third)
		fvals.push_back(fval(SAA_IND_THIRD_BOTH_ABOVE, 1 ));

	if ( ! n_above_third && n_above_zero)
		fvals.push_back(fval(SAA_N_SCORE_ABOVE_ZERO, 1));

	if ( ! c_above_third && c_above_zero)
		fvals.push_back(fval(SAA_C_SCORE_ABOVE_ZERO, 1));

	if (! (n_above_third && c_above_third) && n_above_zero && c_above_zero)
		fvals.push_back(fval(SAA_IND_ZERO_BOTH_ABOVE, 1 ));

	float max_score_rank=-1,min_score_rank=-1;
	if (n_node_score > c_node_score)
	{
		max_score_rank = n_node.log_rank;
		min_score_rank = c_node.log_rank;
	}
	else
	{
		max_score_rank = c_node.log_rank;
		min_score_rank = n_node.log_rank;
	}

	fvals.push_back(fval(SAA_MAX_SCORE_RANK,max_score_rank));
	fvals.push_back(fval(SAA_MIN_SCORE_RANK,min_score_rank));

	fvals.push_back(fval(SAA_N_SCORE_RANK, n_node.log_rank)); 
	fvals.push_back(fval(SAA_C_SCORE_RANK, c_node.log_rank));
	fvals.push_back(fval(SAA_SCORE_RANK_SUM, n_node.log_rank+c_node.log_rank));
	fvals.push_back(fval(SAA_SCORE_RANK_DIFF, n_node.log_rank-c_node.log_rank));
	fvals.push_back(fval(SAA_SCORE_RANK_ABS_DIFF, fabs(n_node.log_rank-c_node.log_rank) ));

	const int num_n_frags = n_node.breakage.fragments.size();
	const int num_c_frags = c_node.breakage.fragments.size();

	fvals.push_back(fval(SAA_N_NUM_FRAGS,num_n_frags ));
	fvals.push_back(fval(SAA_C_NUM_FRAGS,num_c_frags ));
	fvals.push_back(fval(SAA_NUM_FRAG_DIFF, num_n_frags - num_c_frags ));
	fvals.push_back(fval(SAA_ABS_NUM_FRAG_DIFF, abs(num_n_frags - num_c_frags) ));

	// intensity ratio

	const intensity_t& n_inten = n_node.breakage.total_intensity;
	const intensity_t& c_inten = c_node.breakage.total_intensity;

	if (n_inten<=0 && c_inten<=0)
	{
		
	}
	else if (c_inten<=0)
	{
		fvals.push_back(fval(SAA_IND_N_STRONGER_INTEN,1));
		fvals.push_back(fval(SAA_N_STRONGER_LOG_NC_INTEN_RATIO,5.0));
	}
	else if (n_inten<=0)
	{
		fvals.push_back(fval(SAA_IND_C_STRONGER_INTEN,1));
		fvals.push_back(fval(SAA_C_STRONGER_LOG_NC_INTEN_RATIO,5.0));
	}
	else
	{
		if (n_inten>=c_inten)
		{
			fvals.push_back(fval(SAA_IND_N_STRONGER_INTEN,1));
			fvals.push_back(fval(SAA_N_STRONGER_LOG_NC_INTEN_RATIO,log(n_inten/c_inten) ));
		}
		else
		{
			fvals.push_back(fval(SAA_IND_C_STRONGER_INTEN,1));
			fvals.push_back(fval(SAA_C_STRONGER_LOG_NC_INTEN_RATIO,log(c_inten/n_inten) ));
		}
	}

	// connection features
	int max_in_c_idx  = c_node.idx_max_in_score_node;
	int max_out_n_idx = n_node.idx_max_out_score_node;

	if (max_in_c_idx<0 || max_out_n_idx<0)
	{
		cout << "Error: max_in and max_out idxs not filled correctly!" << endl;
		exit(1);
	}

	if (max_in_c_idx == me.n_idx)
	{
		fvals.push_back(fval(SAA_IND_N_IS_MAX_IDX_TO_C,1));
	}
	else
	{
		fvals.push_back(fval(SAA_IND_N_NOT_MAX_IDX_TO_C,1));
		fvals.push_back(fval(SAA_DIFF_N_MAX_IN_C_SCORE_RANKS,
			n_node.log_rank - prm_ptr->get_node(max_in_c_idx).log_rank));
	}

	if (max_out_n_idx == me.c_idx)
	{
		fvals.push_back(fval(SAA_IND_C_IS_MAX_IDX_FROM_N,1));
	}
	else
	{
		fvals.push_back(fval(SAA_IND_C_NOT_MAX_IDX_FROM_N,1));
		fvals.push_back(fval(SAA_DIFF_C_MAX_OUT_N_SCORE_RANKS,
			c_node.log_rank - prm_ptr->get_node(max_out_n_idx).log_rank));
	}

	bool both_connect_to_max =  (max_in_c_idx == me.n_idx && max_out_n_idx == me.c_idx);
	if (both_connect_to_max)
		fvals.push_back(fval(SAA_IND_BOTH_CONNECT_TO_MAX,1));

	
	// node mass diff

	mass_t node_mass_diff = fabs(c_node.mass - n_node.mass - exp_mass);

	if (node_mass_diff>3.0)
	{
		cout << "Error: large node mass diff: " << node_mass_diff << endl;
		exit(1);
	}

	fvals.push_back(fval(SAA_NODE_MASS_DIFF,node_mass_diff));
	fvals.push_back(fval(SAA_NODE_SQR_MASS_DIFF,node_mass_diff*node_mass_diff));

	// peak mass diff

	mass_t best_mass_diff = 1000.0;
	int i;
	int num_pairs=0;
	mass_t avg_diff=0;
	const vector<BreakageFragment> & n_fragments = n_node.breakage.fragments;
	const vector<BreakageFragment> & c_fragments = c_node.breakage.fragments;
	for (i=0; i<n_fragments.size(); i++)
	{
		const int& frag_type_idx = n_fragments[i].frag_type_idx;
		const int pos = c_node.breakage.get_position_of_frag_idx(frag_type_idx);

		if (pos<0)
			continue;

		num_pairs++;

		const int charge=config->get_fragment(frag_type_idx).charge;

		mass_t mass_diff = fabs(n_fragments[i].mass - c_fragments[pos].mass);
		mass_diff *= charge;
		mass_diff -= exp_mass;
		mass_diff = fabs(mass_diff);

		avg_diff+=mass_diff;

		if (mass_diff<best_mass_diff)
			best_mass_diff = mass_diff;
	}


	if (best_mass_diff<100.0)
	{
		avg_diff /= num_pairs;
		fvals.push_back(fval(SAA_NUM_FRAG_PAIRS,num_pairs));
		fvals.push_back(fval(SAA_AVG_PEAK_DIFF,avg_diff));
		fvals.push_back(fval(SAA_AVG_PEAK_SQR_DIFF,avg_diff*avg_diff));
		fvals.push_back(fval(SAA_BEST_PEAK_DIFF,best_mass_diff));
		fvals.push_back(fval(SAA_BEST_PEAK_SQR_DIFF,best_mass_diff*best_mass_diff));

		fvals.push_back(fval(SAA_AVG_PEAK_DIFF_TIMES_SCORE_RANK_SUM,
			avg_diff * (n_node.log_rank+c_node.log_rank) ));
		fvals.push_back(fval(SAA_AVG_PEAK_DIFF_TIMES_SCORE_ABS_DIFF,
			avg_diff * fabs(n_node.log_rank-c_node.log_rank) ));

		fvals.push_back(fval(SAA_AVG_PEAK_DIFF_DIV_NUM_FRAG_PAIRS,
			avg_diff / num_pairs));

		if (both_connect_to_max)
			fvals.push_back(fval(SAA_IND_BOTH_CONNECT_TO_MAX_TIMES_AVG_DIFF,avg_diff));

	}
	else
		fvals.push_back(fval(SAA_IND_NO_PEAK_DIFF,1.0));


	
	// mirror edge features
	const int n_mirror_idx = forbidden_node_idxs[me.n_idx];
	const int c_mirror_idx = forbidden_node_idxs[me.c_idx];
	const float log_sum = n_node.log_rank + c_node.log_rank;
	float mirror_log_sum=0;
	bool n_more=false;
	if (n_mirror_idx>=0)
	{
		const Node& n_mirror_node = prm_ptr->get_node(n_mirror_idx);
		n_more = (n_mirror_node.score<n_node_score);
		if (n_more)
			sam.add_feature(SAA_IND_N_SCORE_MORE_THAN_MIRROR,1.0);

		mirror_log_sum +=n_mirror_node.log_rank;
	}

	if (c_mirror_idx>=0)
	{
		const Node& c_mirror_node = prm_ptr->get_node(c_mirror_idx);	
		mirror_log_sum += c_mirror_node.log_rank;

		if (c_mirror_node.score<c_node_score)
		{
			sam.add_feature(SAA_IND_C_SCORE_MORE_THAN_MIRROR,1.0);

			if (n_more)
				sam.add_feature(SAA_IND_BOTH_SCORE_MORE_THAN_MIRROR,1.0);
		}


		if (n_mirror_idx>=0)
		{
			if (log_sum>mirror_log_sum)
			{
				sam.add_feature(SAA_LOG_DIFF_MORE_THAN_MIRROR,log_sum-mirror_log_sum);
			}
			else
				sam.add_feature(SAA_LOG_DIFF_LESS_THAN_MIRROR,log_sum-mirror_log_sum);
		}
	}
	else
	{
		if (n_mirror_idx<0)
			sam.add_feature(SAA_NO_MIRROR,1.0);
	}

	// amino acid features
	int aa = org_aa[var_aa];
	if (aa == Ile)
		aa = Leu;

	sam.add_feature(SAA_IND_AA_N_TERM+aa,1.0);
	sam.add_feature(SAA_SCORE_RANK_SUM_N_TERM+aa,log_sum);
	sam.add_feature(SAA_SCORE_RANK_DIFF_N_TERM+aa,n_node.log_rank - c_node.log_rank);

	sort(fvals.begin(),fvals.end()); 

/*	if (fvals[fvals.size()-1].f_idx>=SAA_NUM_FIELDS)
	{
		cout << "SAA_SCORE_RANK_DIFF_N_TERM: " << SAA_SCORE_RANK_DIFF_N_TERM << endl;
		cout << "aa: " << aa << endl;
		int i;
		for (i=0; i<fvals.size(); i++)
			cout << i << "\t" << fvals[i].f_idx << "\t" << fvals[i].val << endl;
		exit(0);
	}*/
}

void fill_fval_vector_for_saancd(const PrmGraph* prm_ptr, 
								 const int me_idx, 
								 const int* variant_ptr, 
								 const int  seq_rank,
								 ME_Regression_Sample& sam)
{
	const Config *config = prm_ptr->get_config();
	const vector<mass_t>& aa2mass =config->get_aa2mass();
	const vector<int>& org_aa = config->get_org_aa();
	const vector<int>& forbidden_node_idxs = prm_ptr->get_forbidden_node_idxs();
	const vector<score_t>& cummulative_node_scores = prm_ptr->get_cummulative_scores();
	const score_t total_positive_score = cummulative_node_scores[cummulative_node_scores.size()-1];

	if (*variant_ptr != 1)
	{
		cout << "Error: using a single aa score for a variant with " << *variant_ptr << " aas." << endl;
		exit(1);
	}

	int var_aa = *(variant_ptr+1); // this is the variant aa (the amino acid of the edge)
	const mass_t exp_mass = aa2mass[var_aa];

	var_aa = org_aa[var_aa];
	if (var_aa==Ile)
		var_aa=Leu;

	
	const MultiEdge& edge = prm_ptr->get_multi_edge(me_idx);
	const Node& n_node  = prm_ptr->get_node(edge.n_idx);
	const Node& c_node  = prm_ptr->get_node(edge.c_idx);
	const int rank_offset = calc_rank_idx(seq_rank);

	const mass_t mass_diff = fabs(exp_mass - c_node.mass + n_node.mass);
	const mass_t sqr_diff  = mass_diff * mass_diff;

	score_t n_node_score = n_node.score;
	score_t c_node_score = c_node.score;
	if (n_node_score<-30)
		n_node_score = -30;
	if (n_node_score>30)
		n_node_score=30;
	if (c_node_score<-30)
		c_node_score = -30;
	if (c_node_score>30)
		c_node_score=30;

	vector<fval>& fvals = sam.f_vals;
	fvals.clear();

	fvals.push_back(fval(SAANCD_CONST,1.0));

	sam.add_feature(SAANCD_IND_SEQ_RANK_1+rank_offset,1.0);

	if (n_node.breakage.fragments.size() == 0 &&
		c_node.breakage.fragments.size() == 0)
		fvals.push_back(fval(SAANCD_IND_CONNECTS_TO_NODE_WITH_NO_FRAGS,1.0));

	
	if (n_node.type == NODE_N_TERM)
	{
		fvals.push_back(fval(SAANCD_IND_CONNECTS_TO_N_TERMINAL,1.0));
		if (n_node.idx_max_out_score_node == edge.c_idx)
			fvals.push_back(fval(SAANCD_IND_HAS_MAX_SCORE_FROM_N,1.0));

		fvals.push_back(fval(SAANCD_SCORE_FROM_N,c_node_score));
		fvals.push_back(fval(SAANCD_NODE_MASS_DIFF_FROM_N,mass_diff));
		fvals.push_back(fval(SAANCD_NODE_SQR_MASS_DIFF_FROM_N,sqr_diff));


		if (c_node.type == NODE_DIGEST)
		{
			fvals.push_back(fval(SAANCD_IND_TERM_TO_DIGEST,1.0));
			if (c_node.breakage.fragments.size()>0)
				fvals.push_back(fval(SAANCD_IND_TERM_TO_DIGEST_WITH_FRAGS,1.0));
			if (c_node.score>0)
				fvals.push_back(fval(SAANCD_IND_TERM_TO_DIGEST_W_POS_SCORE,1.0));
		}
	}
	else
	if (c_node.type == NODE_C_TERM)
	{
		fvals.push_back(fval(SAANCD_IND_CONNECTS_TO_C_TERMINAL,1.0));
		if (c_node.idx_max_in_score_node == edge.n_idx)
			fvals.push_back(fval(SAANCD_IND_HAS_MAX_SCORE_TO_C,1.0));

		fvals.push_back(fval(SAANCD_SCORE_TO_C,n_node_score));
		fvals.push_back(fval(SAANCD_NODE_MASS_DIFF_TO_C,mass_diff));
		fvals.push_back(fval(SAANCD_NODE_SQR_MASS_DIFF_TO_C,sqr_diff));

		if (n_node.type == NODE_DIGEST)
		{
			fvals.push_back(fval(SAANCD_IND_TERM_TO_DIGEST,1.0));
			if (n_node.breakage.fragments.size()>0)
				fvals.push_back(fval(SAANCD_IND_TERM_TO_DIGEST_WITH_FRAGS,1.0));
			if (n_node.score>0)
				fvals.push_back(fval(SAANCD_IND_TERM_TO_DIGEST_W_POS_SCORE,1.0));
		}
	}

	if (n_node.type == NODE_DIGEST && c_node.type != NODE_C_TERM)
	{
		fvals.push_back(fval(SAANCD_IND_CONNECTS_TO_DIGEST,1.0));
		if (c_node.idx_max_in_score_node == edge.n_idx)
			fvals.push_back(fval(SAANCD_IND_HAS_MAX_SCORE_TO_DIGEST,1.0));

		fvals.push_back(fval(SAANCD_SCORE_TO_DIGEST,c_node_score));
		fvals.push_back(fval(SAANCD_NODE_MASS_DIFF_TO_DIGEST,mass_diff));
		fvals.push_back(fval(SAANCD_NODE_SQR_MASS_DIFF_TO_DIGEST,sqr_diff));
	}

	if (c_node.type == NODE_DIGEST && n_node.type != NODE_N_TERM)
	{
		fvals.push_back(fval(SAANCD_IND_CONNECTS_TO_DIGEST,1.0));
		if (c_node.idx_max_in_score_node == edge.n_idx)
			fvals.push_back(fval(SAANCD_IND_HAS_MAX_SCORE_TO_DIGEST,1.0));

		fvals.push_back(fval(SAANCD_SCORE_TO_DIGEST,n_node_score));
		fvals.push_back(fval(SAANCD_NODE_MASS_DIFF_TO_DIGEST,mass_diff));
		fvals.push_back(fval(SAANCD_NODE_SQR_MASS_DIFF_TO_DIGEST,sqr_diff));
	}


	if (n_node.type == NODE_N_TERM)
	{
		fvals.push_back(fval(SAANCD_IND_NT_AA_N_TERM+var_aa,1.0));
	}


	vector<fval> ff = fvals;

	sort(fvals.begin(),fvals.end());

	if (fvals[fvals.size()-1].f_idx>=SAANCD_NUM_FIELDS)
	{
		int i;
		for (i=0; i<SAANCD_NUM_FIELDS; i++)
		{
			cout << i << "\t" << ScoreModelFields_SAANCD_names[i] << endl;
		}
		cout << endl;

		for (i=0; i<fvals.size(); i++)
			cout << i << "\t" << ff[i].f_idx << "\t" << ff[i].val << endl;
		exit(0);
	}
}

void fill_fval_vector_for_daa(const PrmGraph* prm_ptr, 
							  const int me_idx, 
							  const int* variant_ptr, 
							  const int  seq_rank,
							  ME_Regression_Sample& sam)
{

	const Config *config = prm_ptr->get_config();
	const vector<mass_t>& aa2mass =config->get_aa2mass();
	const vector<int>& org_aa = config->get_org_aa();
	const vector<int>& forbidden_node_idxs = prm_ptr->get_forbidden_node_idxs();
	const vector<score_t>& cummulative_node_scores = prm_ptr->get_cummulative_scores();
	const score_t total_positive_score = cummulative_node_scores[cummulative_node_scores.size()-1];

	const int num_aa = *variant_ptr;
	if (num_aa < 2)
	{
		cout << "Error: using a double aa score for a variant with " << num_aa << " aas." << endl;
		exit(1);
	}

	const vector<string>& aa2label= config->get_aa2label();

//	cout << "]] " << me_idx << "\t" << *variant_ptr << " " << aa2label[*(variant_ptr+1)] << " " << aa2label[*(variant_ptr+2)] << endl;

	const score_t third_score = (prm_ptr->get_max_node_score()<0) ? 0 : prm_ptr->get_max_node_score() * 0.3333;
	const score_t two_thirds_score = (prm_ptr->get_max_node_score()<0) ? 0 : prm_ptr->get_max_node_score() * 0.66666;

	const MultiEdge& me = prm_ptr->get_multi_edge(me_idx);
	const Node& n_node  = prm_ptr->get_node(me.n_idx);
	const Node& c_node  = prm_ptr->get_node(me.c_idx);
	const int rank_offset = calc_rank_idx(seq_rank);

	const int rel_mass_pos  = (int)((n_node.mass / prm_ptr->get_pm_with_19()) * 5.0);
	const int rel_score_pos =  (int)((cummulative_node_scores[me.n_idx]/total_positive_score) * 5.0);

	const int n_aa = variant_ptr[1];
	const int c_aa = variant_ptr[num_aa];

	score_t n_node_score = n_node.score;
	score_t c_node_score = c_node.score;
	if (n_node_score<-30)
		n_node_score = -30;
	if (n_node_score>30)
		n_node_score=30;
	if (c_node_score<-30)
		c_node_score = -30;
	if (c_node_score>30)
		c_node_score=30;

	const bool n_above_two_thirds = (n_node.score >= two_thirds_score);
	const bool c_above_two_thirds = (c_node.score >= two_thirds_score);
	const bool n_above_third = (n_node.score >= third_score);
	const bool c_above_third = (c_node.score >= third_score);
	const bool n_above_zero = (n_node.score >= 0);
	const bool c_above_zero = (c_node.score >= 0);
	const bool n_above_five = (n_node.score >= 5.0);
	const bool c_above_five = (c_node.score >= 5.0);

	mass_t exp_mass=0;
	int i;
	for (i=0; i<num_aa; i++)
		exp_mass += aa2mass[variant_ptr[i+1]];


	vector<fval>& fvals = sam.f_vals;
	fvals.clear();

	sam.add_feature(DAA_CONST,1.0);
	sam.add_feature(DAA_IND_SEQ_RANK_1+rank_offset,1.0);
	sam.add_feature(DAA_IND_AA_POS_20+rel_mass_pos,1.0);
	sam.add_feature(DAA_IND_AA_REL_SCORE_20+rel_score_pos,1.0);

	if (n_node_score>c_node_score)
	{
		sam.add_feature(DAA_MAX_NODE_SCORE,n_node_score);
		sam.add_feature(DAA_MIN_NODE_SCORE,c_node_score);
	}
	else
	{
		sam.add_feature(DAA_MAX_NODE_SCORE,c_node_score);
		sam.add_feature(DAA_MIN_NODE_SCORE,n_node_score);
	}

	if (n_above_two_thirds)
		fvals.push_back(fval(DAA_N_SCORE_ABOVE_TWO_THIRDS,1));

	if (c_above_two_thirds)
		fvals.push_back(fval(DAA_C_SCORE_ABOVE_TWO_THIRDS,1));

	if (n_above_two_thirds && c_above_two_thirds)
		fvals.push_back(fval(DAA_IND_TWO_THIRDS_BOTH_ABOVE,1));

	if (! n_above_two_thirds && n_above_third)
		fvals.push_back(fval(DAA_N_SCORE_ABOVE_THIRD, 1 ));

	if (! c_above_two_thirds && c_above_third)
		fvals.push_back(fval(DAA_C_SCORE_ABOVE_THIRD, 1 ));

	if (! (n_above_two_thirds && c_above_two_thirds) && n_above_third && c_above_third)
		fvals.push_back(fval(DAA_IND_THIRD_BOTH_ABOVE, 1 ));

	if ( ! n_above_third && n_above_zero)
		fvals.push_back(fval(DAA_N_SCORE_ABOVE_ZERO, 1));

	if ( ! c_above_third && c_above_zero)
		fvals.push_back(fval(DAA_C_SCORE_ABOVE_ZERO, 1));

	if (! (n_above_third && c_above_third) && n_above_zero && c_above_zero)
		fvals.push_back(fval(DAA_IND_ZERO_BOTH_ABOVE, 1 ));

	if (n_above_five && c_above_five)
	{
		fvals.push_back(fval(DAA_IND_BOTH_ABOVE_FIVE, 1 ));
	}
	else if (n_above_five)
	{
		fvals.push_back(fval(DAA_IND_N_ABOVE_FIVE, 1 ));
	}
	else if (c_above_five)
	{
		fvals.push_back(fval(DAA_IND_C_ABOVE_FIVE, 1 ));
	}

	float max_score_rank=-1,min_score_rank=-1;

	if (n_node_score > c_node_score)
	{
		max_score_rank = n_node.log_rank;
		min_score_rank = c_node.log_rank;
	}
	else
	{
		max_score_rank = c_node.log_rank;
		min_score_rank = n_node.log_rank;
	}

	fvals.push_back(fval(DAA_MAX_SCORE_RANK,max_score_rank));
	fvals.push_back(fval(DAA_MIN_SCORE_RANK,min_score_rank));


	fvals.push_back(fval(DAA_N_SCORE_RANK, n_node.log_rank)); 
	fvals.push_back(fval(DAA_C_SCORE_RANK, c_node.log_rank));
	fvals.push_back(fval(DAA_SCORE_RANK_SUM, n_node.log_rank+c_node.log_rank));
	fvals.push_back(fval(DAA_SCORE_RANK_DIFF, n_node.log_rank-c_node.log_rank));
	fvals.push_back(fval(DAA_SCORE_RANK_ABS_DIFF, fabs(n_node.log_rank-c_node.log_rank) ));

	const int num_n_frags = n_node.breakage.fragments.size();
	const int num_c_frags = c_node.breakage.fragments.size();

	fvals.push_back(fval(DAA_N_NUM_FRAGS,num_n_frags ));
	fvals.push_back(fval(DAA_C_NUM_FRAGS,num_c_frags ));
	fvals.push_back(fval(DAA_NUM_FRAG_DIFF, num_n_frags - num_c_frags ));
	fvals.push_back(fval(DAA_ABS_NUM_FRAG_DIFF, abs(num_n_frags - num_c_frags) ));



	// connection features
	int max_in_c_idx  = c_node.idx_max_in_score_node;
	int max_out_n_idx = n_node.idx_max_out_score_node;

	if (max_in_c_idx<0 || max_out_n_idx<0)
	{
		cout << "Error: max_in and max_out idxs not filled correctly!" << endl;
		exit(1);
	}

	if (max_in_c_idx == me.n_idx)
	{
		fvals.push_back(fval(DAA_IND_N_IS_MAX_IDX_TO_C,1));
	}
	else
	{
		fvals.push_back(fval(DAA_IND_N_NOT_MAX_IDX_TO_C,1));
		fvals.push_back(fval(DAA_DIFF_N_MAX_IN_C_SCORE_RANKS,
			n_node.log_rank - prm_ptr->get_node(max_in_c_idx).log_rank));
	}

	if (max_out_n_idx == me.c_idx)
	{
		fvals.push_back(fval(DAA_IND_C_IS_MAX_IDX_FROM_N,1));
	}
	else
	{
		fvals.push_back(fval(DAA_IND_C_NOT_MAX_IDX_FROM_N,1));
		fvals.push_back(fval(DAA_DIFF_C_MAX_OUT_N_SCORE_RANKS,
			c_node.log_rank - prm_ptr->get_node(max_out_n_idx).log_rank));
	}

	bool both_connect_to_max =  (max_in_c_idx == me.n_idx && max_out_n_idx == me.c_idx);
	if (both_connect_to_max)
		fvals.push_back(fval(DAA_IND_BOTH_CONNECT_TO_MAX,1));

	
	// node mass diff

	mass_t node_mass_diff = fabs(c_node.mass - n_node.mass - exp_mass);

	if (node_mass_diff>3.0)
	{
		cout << "ME: " << me_idx << endl;
		cout << "num aa: " << *variant_ptr << endl;
		cout << "aas: " << *(variant_ptr+1) <<" " << *(variant_ptr+2) << endl;
		cout << "Error: large node mass diff: " << node_mass_diff << endl;
		cout << "N NODE: " << n_node.mass << endl;
		cout << "C NODE: " << c_node.mass << endl;
		cout << endl;
		prm_ptr->print();
		exit(1);
	}

	fvals.push_back(fval(DAA_NODE_MASS_DIFF,node_mass_diff));
	fvals.push_back(fval(DAA_NODE_SQR_MASS_DIFF,node_mass_diff*node_mass_diff));

	// peak mass diff

	mass_t best_mass_diff = 1000.0;
	int num_pairs=0;
	mass_t avg_diff=0;
	const vector<BreakageFragment> & n_fragments = n_node.breakage.fragments;
	const vector<BreakageFragment> & c_fragments = c_node.breakage.fragments;
	for (i=0; i<n_fragments.size(); i++)
	{
		const int& frag_type_idx = n_fragments[i].frag_type_idx;
		const int pos = c_node.breakage.get_position_of_frag_idx(frag_type_idx);

		if (pos<0)
			continue;

		num_pairs++;

		const int charge=config->get_fragment(frag_type_idx).charge;

		mass_t mass_diff = fabs(n_fragments[i].mass - c_fragments[pos].mass);
		mass_diff *= charge;
		mass_diff -= exp_mass;
		mass_diff = fabs(mass_diff);

		avg_diff+=mass_diff;

		if (mass_diff<best_mass_diff)
			best_mass_diff = mass_diff;
	}


	if (best_mass_diff<100.0)
	{
		avg_diff /= num_pairs;
		fvals.push_back(fval(DAA_NUM_FRAG_PAIRS,num_pairs));
		fvals.push_back(fval(DAA_AVG_PEAK_DIFF,avg_diff));
		fvals.push_back(fval(DAA_AVG_PEAK_SQR_DIFF,avg_diff*avg_diff));
		fvals.push_back(fval(DAA_BEST_PEAK_DIFF,best_mass_diff));
		fvals.push_back(fval(DAA_BEST_PEAK_SQR_DIFF,best_mass_diff*best_mass_diff));

		fvals.push_back(fval(DAA_AVG_PEAK_DIFF_TIMES_SCORE_RANK_SUM,
			avg_diff * (n_node.log_rank+c_node.log_rank) ));
		fvals.push_back(fval(DAA_AVG_PEAK_DIFF_TIMES_SCORE_ABS_DIFF,
			avg_diff * fabs(n_node.log_rank-c_node.log_rank) ));

		fvals.push_back(fval(DAA_AVG_PEAK_DIFF_DIV_NUM_FRAG_PAIRS,
			avg_diff / num_pairs));

		if (both_connect_to_max)
			fvals.push_back(fval(DAA_IND_BOTH_CONNECT_TO_MAX_TIMES_AVG_DIFF,avg_diff));

	}
	else
		fvals.push_back(fval(DAA_IND_NO_PEAK_DIFF,1.0));


	// features for the alternate route that consists of two amino acid edges (instead
	// of this double edge).

	int      num_alternate_routes=0;
	score_t best_alternate_route_score=NEG_INF;
	int      best_n_edge_idx=-1, best_c_edge_idx=-1, best_node_idx=-1;

	for (i=0; i<n_node.out_edge_idxs.size(); i++)
	{
		const MultiEdge& n_edge = prm_ptr->get_multi_edge(n_node.out_edge_idxs[i]);
		if (n_edge.num_aa>1)
			continue;

		int j;
		for (j=0; j<c_node.in_edge_idxs.size(); j++)
		{
			const MultiEdge& c_edge = prm_ptr->get_multi_edge(c_node.in_edge_idxs[j]);
			if (c_edge.num_aa>1)
				continue;

			if (n_edge.c_idx == c_edge.n_idx)
			{
				num_alternate_routes++;
				score_t node_score = prm_ptr->get_node(n_edge.c_idx).score;
				if (node_score>best_alternate_route_score)
				{
					best_alternate_route_score = node_score;
					best_n_edge_idx = n_node.out_edge_idxs[i];
					best_c_edge_idx = c_node.in_edge_idxs[j];
					best_node_idx = c_edge.n_idx;
				}
			}
		}
	}

	fvals.push_back(fval(DAA_NUM_DOUBLE_EDGE_ROUTES,(float)num_alternate_routes));

	if (best_alternate_route_score>NEG_INF && best_alternate_route_score<200)
	{
		const Node& node = prm_ptr->get_node(best_node_idx);
		
		fvals.push_back(fval(DAA_SCORE_RANK_DOUBLE_EDGE_ROUTES,node.log_rank));
		if (node.score>0)
			fvals.push_back(fval(DAA_IND_IND_MAX_SCORE_ALTERNATE_MORE_THAN_ZERO,1.0));

		if (best_node_idx == n_node.idx_max_out_score_node)
			fvals.push_back(fval(DAA_IND_MAX_ALTERNATE_IS_MAX_OUT_N,1.0));

		if (best_node_idx == c_node.idx_max_in_score_node)
			fvals.push_back(fval(DAA_IND_MAX_ALTERNATE_IS_MAX_IN_C,1.0));

		mass_t e1_mass = aa2mass[*(prm_ptr->get_multi_edge(best_n_edge_idx).variant_ptrs[0]+1)];
		mass_t e2_mass = aa2mass[*(prm_ptr->get_multi_edge(best_c_edge_idx).variant_ptrs[0]+1)];

		mass_t off1 = fabs(e1_mass - node.mass + n_node.mass);
		mass_t off2 = fabs(e2_mass - c_node.mass + node.mass);

		if (off1>5 || off2>5)
		{
			cout << "Error: mismatches in the offset feature of double edges!" << endl;
			exit(1);
		}

		fvals.push_back(fval(DAA_NODE_OFFSETS_ALTERNTE, off1+off2));
		fvals.push_back(fval(DAA_SQR_NODE_ODFFSETS_ALTERNATE, off1*off1 + off2*off2));

	}
	else
		fvals.push_back(fval(DAA_IND_NO_DOUBLE_EDGE_ROUTES,1.0));


	// mirror edge features
	const int n_mirror_idx = forbidden_node_idxs[me.n_idx];
	const int c_mirror_idx = forbidden_node_idxs[me.c_idx];
	const float log_sum = n_node.log_rank + c_node.log_rank;
	float mirror_log_sum=0;
	bool n_more=false;
	if (n_mirror_idx>=0)
	{
		const Node& n_mirror_node = prm_ptr->get_node(n_mirror_idx);
		n_more = (n_mirror_node.score<n_node_score);
		if (n_more)
			sam.add_feature(DAA_IND_N_SCORE_MORE_THAN_MIRROR,1.0);

		mirror_log_sum +=n_mirror_node.log_rank;
	}

	if (c_mirror_idx>=0)
	{
		const Node& c_mirror_node = prm_ptr->get_node(c_mirror_idx);	
		mirror_log_sum += c_mirror_node.log_rank;

		if (c_mirror_node.score<c_node_score)
		{
			sam.add_feature(DAA_IND_C_SCORE_MORE_THAN_MIRROR,1.0);

			if (n_more)
				sam.add_feature(DAA_IND_BOTH_SCORE_MORE_THAN_MIRROR,1.0);
		}


		if (n_mirror_idx>=0)
		{
			if (log_sum>mirror_log_sum)
			{
				sam.add_feature(DAA_LOG_DIFF_MORE_THAN_MIRROR,log_sum-mirror_log_sum);
			}
			else
				sam.add_feature(DAA_LOG_DIFF_LESS_THAN_MIRROR,log_sum-mirror_log_sum);
		}
	}
	else
	{
		if (n_mirror_idx<0)
			sam.add_feature(DAA_NO_MIRROR,1.0);
	}



	const int problem_aa1[]={Gly,Gly,Glu,Val,Ser,Ala,Asp,Ala,Gly};
	const int problem_aa2[]={Gly,Glu,Gly,Ser,Val,Asp,Ala,Gly,Ala};
	const int num_problem_aas = sizeof(problem_aa1)/sizeof(int);

	for (i=0; i<num_problem_aas; i++)
		if (n_aa == problem_aa1[i] && c_aa == problem_aa2[i])
		{
			fvals.push_back(fval(DAA_IND_PROBLEMATIC_PAIR_OF_AAS,1.0));
			break;
		}


	int org_n_aa = org_aa[n_aa];
	int org_c_aa = org_aa[c_aa];
	if (org_n_aa == Ile)
		org_n_aa = Leu;
	if (org_c_aa == Ile)
		org_c_aa = Leu;


	const float log_diff = n_node.log_rank - c_node.log_rank;

	sam.add_feature(DAA_IND_N_AA_N_TERM+org_n_aa,1.0);
	sam.add_feature(DAA_SCORE_RANK_SUM_N_N_TERM+org_n_aa,log_sum);
	sam.add_feature(DAA_SCORE_RANK_DIFF_N_N_TERM+org_n_aa,log_diff);

	sam.add_feature(DAA_IND_C_AA_N_TERM+org_c_aa,1.0);
	sam.add_feature(DAA_SCORE_RANK_SUM_C_N_TERM+org_c_aa,log_sum);
	sam.add_feature(DAA_SCORE_RANK_DIFF_C_N_TERM+org_c_aa,log_diff);
}

void fill_fval_vector_for_daancd(const PrmGraph* prm_ptr, 
								 const int me_idx, 
								 const int* variant_ptr, 
								 const int  seq_rank,
								 ME_Regression_Sample& sam)
{

	const Config *config = prm_ptr->get_config();
	const vector<mass_t>& aa2mass =config->get_aa2mass();
	const vector<int>& org_aa = config->get_org_aa();
	const vector<int>& forbidden_node_idxs = prm_ptr->get_forbidden_node_idxs();
	const vector<score_t>& cummulative_node_scores = prm_ptr->get_cummulative_scores();
	const score_t total_positive_score = cummulative_node_scores[cummulative_node_scores.size()-1];

	const int num_aa = *variant_ptr;
	if (num_aa < 2)
	{
		cout << "Error: using a double aa score for a variant with " << num_aa << " aas." << endl;
		exit(1);
	}

	const MultiEdge& edge = prm_ptr->get_multi_edge(me_idx);
	const Node& n_node  = prm_ptr->get_node(edge.n_idx);
	const Node& c_node  = prm_ptr->get_node(edge.c_idx);
	const int rank_offset = calc_rank_idx(seq_rank);

	score_t n_node_score = n_node.score;
	score_t c_node_score = c_node.score;
	if (n_node_score<-30)
		n_node_score = -30;
	if (n_node_score>30)
		n_node_score=30;
	if (c_node_score<-30)
		c_node_score = -30;
	if (c_node_score>30)
		c_node_score=30;

	mass_t exp_mass=0;
	int i;
	for (i=0; i<num_aa; i++)
		exp_mass += aa2mass[variant_ptr[i+1]];

	const mass_t mass_diff = fabs(exp_mass - c_node.mass + n_node.mass);
	const mass_t sqr_diff  = mass_diff * mass_diff;


	const int n_aa = variant_ptr[1];
	const int c_aa = variant_ptr[num_aa];


	vector<fval>& fvals = sam.f_vals;
	fvals.clear();

	fvals.push_back(fval(DAANCD_CONST,1.0));
	fvals.push_back(fval(DAANCD_IND_SEQ_RANK_1+rank_offset,1.0));

	if (n_node.breakage.fragments.size() == 0 &&
		c_node.breakage.fragments.size() == 0)
		fvals.push_back(fval(DAANCD_IND_CONNECTS_TO_NODE_WITH_NO_FRAGS,1.0));

	
	if (n_node.type == NODE_N_TERM)
	{
		fvals.push_back(fval(DAANCD_IND_CONNECTS_TO_N_TERMINAL,1.0));
		if (n_node.idx_max_out_score_node == edge.c_idx)
			fvals.push_back(fval(DAANCD_IND_HAS_MAX_SCORE_FROM_N,1.0));

		fvals.push_back(fval(DAANCD_SCORE_FROM_N,c_node_score));
		fvals.push_back(fval(DAANCD_NODE_MASS_DIFF_FROM_N,mass_diff));
		fvals.push_back(fval(DAANCD_NODE_SQR_MASS_DIFF_FROM_N,sqr_diff));
	}

	if (c_node.type == NODE_C_TERM)
	{
		fvals.push_back(fval(DAANCD_IND_CONNECTS_TO_C_TERMINAL,1.0));
		if (c_node.idx_max_in_score_node == edge.n_idx)
			fvals.push_back(fval(DAANCD_IND_HAS_MAX_SCORE_TO_C,1.0));

		fvals.push_back(fval(DAANCD_SCORE_TO_C,n_node_score));
		fvals.push_back(fval(DAANCD_NODE_MASS_DIFF_TO_C,mass_diff));
		fvals.push_back(fval(DAANCD_NODE_SQR_MASS_DIFF_TO_C,sqr_diff));
	}

	if (n_node.type == NODE_DIGEST && c_node.type != NODE_C_TERM)
	{
		fvals.push_back(fval(DAANCD_IND_CONNECTS_TO_DIGEST,1.0));
		if (c_node.idx_max_in_score_node == edge.n_idx)
			fvals.push_back(fval(DAANCD_IND_HAS_MAX_SCORE_TO_DIGEST,1.0));

		fvals.push_back(fval(DAANCD_SCORE_TO_DIGEST,c_node_score));
		fvals.push_back(fval(DAANCD_NODE_MASS_DIFF_TO_DIGEST,mass_diff));
		fvals.push_back(fval(DAANCD_NODE_SQR_MASS_DIFF_TO_DIGEST,sqr_diff));
	}

	if (c_node.type == NODE_DIGEST && n_node.type != NODE_N_TERM)
	{
		fvals.push_back(fval(DAANCD_IND_CONNECTS_TO_DIGEST,1.0));
		if (c_node.idx_max_in_score_node == edge.n_idx)
			fvals.push_back(fval(DAANCD_IND_HAS_MAX_SCORE_TO_DIGEST,1.0));

		fvals.push_back(fval(DAANCD_SCORE_TO_DIGEST,n_node_score));
		fvals.push_back(fval(DAANCD_NODE_MASS_DIFF_TO_DIGEST,mass_diff));
		fvals.push_back(fval(DAANCD_NODE_SQR_MASS_DIFF_TO_DIGEST,sqr_diff));
	}


	// features for the alternate route that consists of two amino acid edges (instead
	// of this double edge).

	int      num_alternate_routes=0;
	score_t best_alternate_route_score=NEG_INF;
	int      best_n_edge_idx=-1, best_c_edge_idx=-1, best_node_idx=-1;

	for (i=0; i<n_node.out_edge_idxs.size(); i++)
	{
		const MultiEdge& n_edge = prm_ptr->get_multi_edge(n_node.out_edge_idxs[i]);
		if (n_edge.num_aa>1)
			continue;

		int j;
		for (j=0; j<c_node.in_edge_idxs.size(); j++)
		{
			const MultiEdge& c_edge = prm_ptr->get_multi_edge(c_node.in_edge_idxs[j]);
			if (c_edge.num_aa>1)
				continue;

			if (n_edge.c_idx == c_edge.n_idx)
			{
				num_alternate_routes++;
				score_t node_score = prm_ptr->get_node(n_edge.c_idx).score;
				if (node_score>best_alternate_route_score)
				{
					best_alternate_route_score = node_score;
					best_n_edge_idx = n_node.out_edge_idxs[i];
					best_c_edge_idx = c_node.in_edge_idxs[j];
					best_node_idx = c_edge.n_idx;
				}
			}
		}
	}

	fvals.push_back(fval(DAANCD_NUM_DOUBLE_EDGE_ROUTES,(float)num_alternate_routes));

	if (best_alternate_route_score>-100 && best_alternate_route_score<200)
	{
		const Node& node = prm_ptr->get_node(best_node_idx);
		
		fvals.push_back(fval(DAANCD_SCORE_RANK_DOUBLE_EDGE_ROUTES,node.log_rank));
		if (node.score>0)
			fvals.push_back(fval(DAANCD_IND_IND_MAX_SCORE_ALTERNATE_MORE_THAN_ZERO,1.0));

		if (best_node_idx == n_node.idx_max_out_score_node)
			fvals.push_back(fval(DAANCD_IND_MAX_ALTERNATE_IS_MAX_OUT_N,1.0));

		if (best_node_idx == c_node.idx_max_in_score_node)
			fvals.push_back(fval(DAANCD_IND_MAX_ALTERNATE_IS_MAX_IN_C,1.0));

		const mass_t e1_mass = aa2mass[*(prm_ptr->get_multi_edge(best_n_edge_idx).variant_ptrs[0]+1)];
		const mass_t e2_mass = aa2mass[*(prm_ptr->get_multi_edge(best_c_edge_idx).variant_ptrs[0]+1)];
		const mass_t off1 = fabs(e1_mass - node.mass + n_node.mass);
		const mass_t off2 = fabs(e2_mass - c_node.mass + node.mass);

		if (off1>5 || off2>5)
		{
			cout << "Error: mismatches in the offset feature of double edges!" << endl;
			exit(1);
		}

		fvals.push_back(fval(DAANCD_NODE_OFFSETS_ALTERNTE, off1+off2));
		fvals.push_back(fval(DAANCD_SQR_NODE_ODFFSETS_ALTERNATE, off1*off1 + off2*off2));

	}
	else
		fvals.push_back(fval(DAANCD_IND_NO_DOUBLE_EDGE_ROUTES,1.0));

	if (n_node.type == NODE_N_TERM)
	{
		int n_aa = org_aa[variant_ptr[1]];
		int c_aa = org_aa[variant_ptr[num_aa]];

		if (n_aa == Ile)
			n_aa = Leu;
		if (c_aa == Ile)
			c_aa = Leu;

		fvals.push_back(fval(DAANCD_IND_NT_N_AA_N_TERM+n_aa,1.0));
		fvals.push_back(fval(DAANCD_IND_NT_C_AA_N_TERM+c_aa,1.0));
	}

	sort(fvals.begin(),fvals.end());


}

int  AminoAcidProbs::get_edge_type(const void* prm_void_ptr, int me_idx) const
{
	const PrmGraph * prm_ptr = (const PrmGraph *)prm_void_ptr;
	const MultiEdge& edge = prm_ptr->get_multi_edge(me_idx);
	const Node& n_node = prm_ptr->get_node(edge.n_idx);
	const Node& c_node = prm_ptr->get_node(edge.c_idx);


	if (edge.num_aa == 1)
	{
		if (n_node.type == NODE_N_TERM || c_node.type == NODE_C_TERM)
			return AAP_SAANCD;

		if (n_node.breakage.fragments.size()>0 && c_node.breakage.fragments.size()>0)
			return AAP_SAA;

		return AAP_SAANCD;

	}
	
	
	if (n_node.type == NODE_N_TERM || c_node.type == NODE_C_TERM)
		return AAP_DAANCD;

	if (n_node.breakage.fragments.size()>0 && c_node.breakage.fragments.size()>0)
		return AAP_DAA;

	return AAP_DAANCD;
}


void AminoAcidProbs::fill_aa_prob_fval_vector(int type,
						  const void* prm_void_ptr, 
						  const int me_idx, 
						  const int* variant_ptr, 
						  const int  seq_rank,
						  ME_Regression_Sample& sam) const
{
	const PrmGraph * prm_ptr = (const PrmGraph *)prm_void_ptr;
	if (type ==  AAP_SAA)
	{
		fill_fval_vector_for_saa(prm_ptr,me_idx,variant_ptr,seq_rank,sam);
		return;
	}
		
	if (type == AAP_SAANCD)
	{
		fill_fval_vector_for_saancd(prm_ptr,me_idx,variant_ptr,seq_rank,sam);
		return;
	}
		
	if (type == AAP_DAA)
	{
		fill_fval_vector_for_daa(prm_ptr,me_idx,variant_ptr,seq_rank,sam);
		return;
	}
		
	if (type == AAP_DAANCD)
	{
		fill_fval_vector_for_daancd(prm_ptr,me_idx,variant_ptr,seq_rank,sam);
		return;
	}

	cout << "Error: illegel amino acid prob edge type: " << type << endl;
	exit(1);
}



float AminoAcidProbs::calc_variant_prob(const void *prm_void_ptr, 
										int me_idx, 
										int* variant_ptr, 
										int seq_rank) const
{
	const PrmGraph *prm = (const PrmGraph *)prm_void_ptr;
	const int edge_type = get_edge_type(prm_void_ptr,me_idx);
	const int charge    = prm->get_charge();
	const int size_idx  = prm->get_size_idx();

	if (aa_prob_me_models.size()<=charge ||
		aa_prob_me_models[charge].size() <= size_idx ||
		aa_prob_me_models[charge][size_idx].size() <= edge_type)
	{
		return -2;
	}

	ME_Regression_Sample sam;

	fill_aa_prob_fval_vector(edge_type,prm_void_ptr,me_idx,variant_ptr,seq_rank,sam);

	return aa_prob_me_models[charge][size_idx][edge_type]->p_y_given_x(0,sam);
}



/*
void AminoAcidProbs::train_amino_acid_prob_models(const FileManager& fm, 
												  void *model_ptr, 
												  int specific_charge,
												  int specific_size_idx)
{	
	// These are th idxs of sample to be used in training
	const int num_rerank_sols = 300;
	const int sol_idxs[] ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
						   22,24,26,28,30,36,51,77,100,145,195,202,245,281};
	const int num_sols = sizeof(sol_idxs)/sizeof(int);

	const int num_features[]={SAA_NUM_FIELDS, SAANCD_NUM_FIELDS, DAA_NUM_FIELDS, DAANCD_NUM_FIELDS};
	const char ** feature_names[]={ScoreModelFields_SAA_names,ScoreModelFields_SAANCD_names,
										   ScoreModelFields_DAA_names,ScoreModelFields_DAANCD_names};


	AllScoreModels *model = (AllScoreModels *)model_ptr;
	Config *config = model->get_config();

	const vector< vector<mass_t> >& size_threshes =  config->get_size_thresholds();
	
	PeptideRankScorer *rank_model = (PeptideRankScorer *)model->get_rank_model_ptr(1);
	static vector<PrmGraph *> prm_ptrs;
	prm_ptrs.resize(8,NULL);

	// resize model arrays
	aa_prob_me_models.resize(size_threshes.size());
	int charge;
	for (charge=1; charge<size_threshes.size(); charge++)
	{
		aa_prob_me_models[charge].resize(size_threshes[charge].size());
		int size_idx;
		for (size_idx=0; size_idx<size_threshes[charge].size(); size_idx++)
			aa_prob_me_models[charge][size_idx].resize(AAP_NUM_TYPES,NULL);
	}
	
	config->set_use_spectrum_charge(1);

	// train models
	for (charge=1; charge<size_threshes.size(); charge++)
	{
		if (specific_charge>0 && charge != specific_charge)
			continue;

		int size_idx;
		for (size_idx=0; size_idx<size_threshes[charge].size(); size_idx++)
		{

			if (specific_size_idx>=0 && size_idx != specific_size_idx)
				continue;

			const mass_t min_mass = (size_idx == 0 ? 0 : size_threshes[charge][size_idx-1]);
			const mass_t max_mass = size_threshes[charge][size_idx];


			FileSet fs;

			fs.select_files_in_mz_range(fm,min_mass/charge,max_mass/charge,charge);
			fs.randomly_reduce_ssfs(10000);

			cout << "Training charge " << charge << " size " << size_idx << ", found " << 
				fs.get_total_spectra() << " spectra.." << endl;

			bool perform_rerank = true;
			if (! rank_model->get_ind_part_model_was_initialized(charge,size_idx))
			{
				cout << "Warning: no re rank model for de novo sequences of these models..." << endl;
				perform_rerank = false;
			}
	
			vector<ME_Regression_DataSet> sample_sets;
			sample_sets.resize(AAP_NUM_TYPES);

			BasicSpecReader bsr;
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			int sc;
			for (sc=0; sc<all_ssf.size(); sc++)
			{
				static vector<QCPeak> peaks;
				SingleSpectrumFile *ssf = all_ssf[sc];
				if (peaks.size()<ssf->num_peaks)
				{
					int new_size = ssf->num_peaks*2;
					if (new_size<2500)
						new_size=2500;
					peaks.resize(new_size); 
				}	

				if (sc>0 && sc % 100 == 0)
				{
					cout << "Done " << sc << "/" << all_ssf.size() << endl;
				}

				const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
				Spectrum s;
				s.init_from_QCPeaks(config,&peaks[0],num_peaks,ssf);

				vector<SeqPath> solutions;
				solutions.clear();

				vector<mass_t> pms_with_19;
				vector<int>    charges;
				pms_with_19.clear();
				charges.clear();		
				BasicSpectrum bs;
				bs.ssf = ssf;
				bs.peaks = &peaks[0];
				bs.num_peaks = num_peaks;

				// output m/z and prob values for the different charge states
				vector<PmcSqsChargeRes> pmcsqs_res;
				//model->select_pms_and_charges(config,bs,pms_with_19,charges,&pmcsqs_res);
				
				if (pms_with_19[0]<100)
				{
					continue;
				}
					
				generate_denovo_solutions_from_several_pms(
						prm_ptrs,
						model,
						&s,
						true, 
						300,
						7,
						16,
						pms_with_19,
						charges,
						solutions,
						false);
				
				vector<score_pair> scores;
				if (perform_rerank)
				{
				//	rank_model->score_denovo_sequences(solutions,ssf,&peaks[0],num_peaks,scores,size_idx);
					sort(scores.begin(),scores.end());
				}
				else
				{
					scores.resize(solutions.size());
					int i;
					for (i=0; i<solutions.size(); i++)
						scores[i].idx=i;
				}

				vector<int> correct_aas = bs.ssf->peptide.get_amino_acids();
				vector<mass_t> exp_cuts;
				bs.ssf->peptide.calc_expected_breakage_masses(config,exp_cuts);

				int i;
				for (i=0; i<correct_aas.size(); i++)
					if (correct_aas[i] == Ile)
						correct_aas[i] = Leu;
				
				// extract training samples from solutions
				for (i=0; i<num_sols; i++)
				{
					if (sol_idxs[i]>=solutions.size())
						break;

					const int sol_idx =sol_idxs[i];
					const SeqPath& sol = solutions[sol_idx];

					int j;
					for (j=0; j<sol.positions.size(); j++)
					{
						PathPos pos = sol.positions[j];
						
						const int edge_idx = sol.positions[j].edge_idx;
						const int *variant_ptr = sol.positions[j].variant_ptr;

						if (edge_idx>=0 && variant_ptr)
						{
							const int edge_type = get_edge_type(sol.prm_ptr,edge_idx);
							const int num_aa = sol.prm_ptr->get_multi_edge(edge_idx).num_aa;

						//	if (edge_type == AAP_SAA && myRandom()<0.75)
						//		continue;

						//	if (edge_type == AAP_SAANCD && myRandom()<0.5)
						//		continue;

							ME_Regression_Sample sam;
							fill_aa_prob_fval_vector(edge_type,sol.prm_ptr,edge_idx,variant_ptr,sol_idx,sam);

							sam.weight=1;
						
							// check if edge is correct
							int cut_idx;
							for (cut_idx=0; cut_idx<exp_cuts.size(); cut_idx++)
								if (fabs(pos.mass-exp_cuts[cut_idx])<2.0)
									break;

							if (cut_idx == exp_cuts.size())
							{
								sam.label = 1;
							}
							else
							{
							
								bool skip = false;
								int k;
								for (k=0; k<num_aa; k++)
								{
									const int corr_aa = correct_aas[cut_idx+k];
									const int edge_aa = sol.positions[j+k].aa;
									if ( corr_aa == Lys && edge_aa == Gln)
									{
										skip = true;
										break;
									}

									if (corr_aa != edge_aa)
										break;
								}

								if (skip && k == num_aa)
									continue;

								sam.label = ( k == num_aa ? 0 : 1);
							}
									
							sample_sets[edge_type].add_sample(sam);
						}
					}
				}
			}

			// normalize dataset and train me model


			int type;
			for (type=0; type<sample_sets.size(); type++)
			{
				sample_sets[type].num_classes=2;
				sample_sets[type].num_features = num_features[type];
				if (! aa_prob_me_models[charge][size_idx][type])
					aa_prob_me_models[charge][size_idx][type] = new ME_Regression_Model;

				cout << endl << "Training ME model for type " << type << ", #features = " << num_features[type] << endl;
				sample_sets[type].tally_samples();
				sample_sets[type].print_summary();
				sample_sets[type].print_feature_summary(cout,feature_names[type]);
				aa_prob_me_models[charge][size_idx][type]->train_cg(sample_sets[type],300,1E-5);

				aa_prob_me_models[charge][size_idx][type]->print_ds_probs(sample_sets[type]);
				cout << endl;
			}

			// write models
			char fname[128];
			sprintf(fname,"%s_AAP_%d_%d.txt",config->getModelName().c_str(),charge,size_idx);

			ofstream ofs(fname);
			ofs << "#MODELS\t" << charge << "\t" << size_idx << endl;
			for (type=0; type<sample_sets.size(); type++)
				aa_prob_me_models[charge][size_idx][type]->write_regression_model(ofs);
			ofs.close();
		}
	}

	ind_model_was_initialized = true;
}
*/

void AminoAcidProbs::read_amino_acid_prob_models(Config *config, char *file_name)
{
	string path = config->get_resource_dir() + "/" + file_name;
	ifstream ifs(path.c_str());
	if (! ifs.good()  || ! ifs.is_open())
	{
		cout << "Error: could't open aa prob model file for reading: " << path << endl;
		exit(1);
	}


	const vector< vector<mass_t> >& size_threshes =  config->get_size_thresholds();
	// resize model arrays
	aa_prob_me_models.resize(size_threshes.size());
	int charge;
	for (charge=1; charge<size_threshes.size(); charge++)
	{
		aa_prob_me_models[charge].resize(size_threshes[charge].size());
		int size_idx;
		for (size_idx=0; size_idx<size_threshes[charge].size(); size_idx++)
			aa_prob_me_models[charge][size_idx].resize(AAP_NUM_TYPES,NULL);
	}
	
	char buff[128];
	while (! ifs.eof())
	{
		ifs.getline(buff,128);


		if (strncmp("#MODELS",buff,7))
			continue;	

		int charge=-1,size_idx=-1;
		if (sscanf(buff,"#MODELS\t%d\t%d",&charge,&size_idx) != 2)
		{
			cout << "Error: bad line in " << path << endl << "LINE:" << buff << endl;
			exit(1);
		}

		int type;
		for (type=0; type<AAP_NUM_TYPES; type++)
		{
			if (! aa_prob_me_models[charge][size_idx][type])
				aa_prob_me_models[charge][size_idx][type] = new ME_Regression_Model;

			aa_prob_me_models[charge][size_idx][type]->read_regression_model(ifs);
		}
	}

	ifs.close();
	ind_model_was_initialized = true;
}


void AminoAcidProbs::write_amino_acid_prob_models(const char *path) const
{

	if (! ind_model_was_initialized)
		return;

	ofstream ofs(path);
	if (! ofs.good()  || ! ofs.is_open())
	{
		cout << "Error: could't open aa prob model file for writing: " << path << endl;
		exit(1);
	}

	int charge;
	for (charge=1; charge<this->aa_prob_me_models.size(); charge++)
	{
		int size_idx;
		for (size_idx=0; size_idx<aa_prob_me_models[charge].size(); size_idx++)
		{
			if (aa_prob_me_models[charge][size_idx].size() != AAP_NUM_TYPES)
				continue;

			int type;
			for (type=0; type<AAP_NUM_TYPES; type++)
				if (! aa_prob_me_models[charge][size_idx][type] )
					break;

			if (type<AAP_NUM_TYPES)
				continue;

			ofs << "#MODELS\t" << charge << "\t" << size_idx << endl;
			for (type=0; type<AAP_NUM_TYPES; type++)
				aa_prob_me_models[charge][size_idx][type]->write_regression_model(ofs);
		}
	}
	ofs.close();
}
