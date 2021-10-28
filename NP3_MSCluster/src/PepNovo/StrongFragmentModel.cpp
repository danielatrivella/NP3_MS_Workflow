#include "StrongFragmentModel.h"

const char * ScoreModelFields_SI_names[]={
"SI_CONST",	"SI_IND_MIRROR1_VIZ",	"SI_IND_MIRROR1_NOT_VIZ",	"SI_IND_HAS_MIRROR1_INTEN",	"SI_MIRROR1_ISO_LEVEL",	"SI_IND_MIRROR1_HAS_MINUS_1",	"SI_MIRROR1_MINUS_1_INTEN_DIFF",	"SI_IND_MIRROR1_HAS_MINUS_2",	"SI_MIRROR1_MINUS_2_INTEN_DIFF",	"SI_IND_MIRROR1_HAS_PLUS_1",	"SI_MIRROR1_PLUS_1_INTEN_DIFF",	"SI_IND_MIRROR1_HAS_PLUS_2",	"SI_MIRROR1_PLUS_2_INTEN_DIFF",	"SI_MIRROR1_MASS_DIFF25",	
"SI_MIRROR1_MASS_DIFF75",	"SI_MIRROR1_MASS_DIFF_LARGE",	"SI_IND_MIRROR2_VIZ",	"SI_IND_MIRROR2_NOT_VIZ",	"SI_IND_HAS_MIRROR2_INTEN",	"SI_MIRROR2_ISO_LEVEL",	"SI_IND_MIRROR2_HAS_MINUS_1",	"SI_MIRROR2_MINUS_1_INTEN_DIFF",	"SI_IND_MIRROR2_HAS_MINUS_2",	"SI_MIRROR2_MINUS_2_INTEN_DIFF",	"SI_IND_MIRROR2_HAS_PLUS_1",	"SI_MIRROR2_PLUS_1_INTEN_DIFF",	"SI_IND_MIRROR2_HAS_PLUS_2",	"SI_MIRROR2_PLUS_2_INTEN_DIFF",	
"SI_MIRROR2_MASS_DIFF25",	"SI_MIRROR2_MASS_DIFF75",	"SI_MIRROR2_MASS_DIFF_LARGE",	"SI_IND_PARENT1_VIZ",	"SI_IND_PARENT1_NOT_VIZ",	"SI_IND_PARENT1_NO_INTEN",	"SI_PARENT1_ISO_LEVEL",	"SI_IND_PARENT1_LESS_THAN_100_MIN_MAX",	"SI_IND_PARENT1_LESS_THAN_200_MIN_MAX",	"SI_IND_PARENT1_INTEN_MORE",	"SI_PARENT1_INTEN_DIFF_MORE",	"SI_IND_PARENT1_INTEN_LESS",	"SI_PARENT1_INTEN_DIFF_LESS",	
"SI_IND_PARENT2_VIZ",	"SI_IND_PARENT2_NOT_VIZ",	"SI_IND_PARENT2_NO_INTEN",	"SI_PARENT2_ISO_LEVEL",	"SI_IND_PARENT2_LESS_THAN_100_MIN_MAX",	"SI_IND_PARENT2_LESS_THAN_200_MIN_MAX",	"SI_IND_PARENT2_INTEN_MORE",	"SI_PARENT2_INTEN_DIFF_MORE",	"SI_IND_PARENT2_INTEN_LESS",	"SI_PARENT2_INTEN_DIFF_LESS",	"SI_LOG_LOCAL_RANK",	"SI_LOG_GLOBAL_RANK",	"SI_ISO_LEVEL",	"SI_IND_HAS_MINUS_1",	
"SI_IND_HAS_MINUS_2",	"SI_IND_HAS_PLUS_1",	"SI_IND_HAS_PLUS_2",	"SI_IND_LOG_INTEN_LESS1",	"SI_LOG_INTEN_LESS1",	"SI_IND_LOG_INTEN_LESS2",	"SI_LOG_INTEN_LESS2",	"SI_IND_LOG_INTEN_LESS3",	"SI_LOG_INTEN_LESS3",	"SI_IND_LOG_INTEN_LESS4",	"SI_LOG_INTEN_LESS4",	"SI_IND_LOG_INTEN_MORE",	"SI_LOG_INTEN_MORE",	"SI_IND_DIS_FROM_MINMAX_LESS_50",	"SI_DIS_FROM_MINMAX0",	"SI_LOG_INTEN_DIS_50",	
"SI_IND_DIS_FROM_MINMAX_LESS_150",	"SI_DIS_FROM_MINMAX50",	"SI_LOG_INTEN_DIS_150",	"SI_IND_DIS_FROM_MINMAX_LESS_250",	"SI_DIS_FROM_MINMAX150",	"SI_LOG_INTEN_DIS_250",	"SI_IND_DIS_FROM_MINMAX_MORE",	"SI_DIS_FROM_MINMAX250",	"SI_LOG_INTEN_DIS_MORE",	"SI_REL_POS0",	"SI_REL_POS1",	"SI_REL_POS2",	"SI_REL_POS3",	"SI_REL_POS4",	"SI_REL_POS5",	"SI_REL_POS6",	"SI_REL_POS7",	"SI_REL_POS8",	
"SI_REL_POS9",	"SI_IND_HAS_PLUS_NH3",	"SI_PLUS_NH3_INTEN_DIFF",	"SI_IND_HAS_PLUS_H2O",	"SI_PLUS_H2O_INTEN_DIFF",	"SI_IND_HAS_PLUS_CO",	"SI_PLUS_CO_INTEN_DIFF",	"SI_IND_HAS_PLUS_H2ONH3",	"SI_PLUS_H2ONH3_INTEN_DIFF",	"SI_IND_HAS_PLUS_H2OH2O",	"SI_PLUS_H2OH2O_INTEN_DIFF",	"SI_IND_HAS_CHARGE_PLUS1",	"SI_CHARGE_PLUS1_INTEN_DIFF",	"SI_IND_HAS_CHARGE_MINUS1",	"SI_CHARGE_MINUS1_INTEN_DIFF",	
"SI_IND_N_IS_GAP",	"SI_IND_C_IS_GAP",	"SI_N_TERM_CAT20",	"SI_N_TERM_CAT18",	"SI_N_TERM_CAT12",	"SI_N_TERM_CAT8",	"SI_N_TERM_CAT4",	"SI_N_TERM_CAT2",	"SI_N_EDGE_CAT20",	"SI_N_EDGE_CAT18",	"SI_N_EDGE_CAT12",	"SI_N_EDGE_CAT8",	"SI_N_EDGE_CAT4",	"SI_N_EDGE_CAT2",	"SI_C_TERM_CAT20",	"SI_C_TERM_CAT18",	"SI_C_TERM_CAT12",	"SI_C_TERM_CAT8",	"SI_C_TERM_CAT4",	"SI_C_TERM_CAT2",	"SI_C_EDGE_CAT20",	
"SI_C_EDGE_CAT18",	"SI_C_EDGE_CAT12",	"SI_C_EDGE_CAT8",	"SI_C_EDGE_CAT4",	"SI_C_EDGE_CAT2",	"SI_SPAN_CAT20",	"SI_SPAN_CAT18",	"SI_SPAN_CAT12",	"SI_SPAN_CAT8",	"SI_SPAN_CAT4",	"SI_SPAN_CAT2",	"SI_ND_SPAN_CAT20",	"SI_ND_SPAN_CAT18",	"SI_ND_SPAN_CAT12",	"SI_ND_SPAN_CAT8",	"SI_ND_SPAN_CAT4",	"SI_ND_SPAN_CAT2",	"SI_CD_SPAN_CAT20",	"SI_CD_SPAN_CAT18",	"SI_CD_SPAN_CAT12",	"SI_CD_SPAN_CAT8",	
"SI_CD_SPAN_CAT4",	"SI_CD_SPAN_CAT2",	"SI_IND_CONNECTS_TO_N_TERM",	"SI_IND_PREFERRED_DIGEST_AA_N_TERM",	"SI_PREFERRED_DIGEST_AA_N_TERM_INTEN",	"SI_NON_PREFERRED_DIGEST_AA_N_TERM_INTEN",	"SI_IND_CONNECTS_TO_C_TERM",	"SI_IND_PREFERRED_DIGEST_AA_C_TERM",	"SI_PREFERRED_DIGEST_AA_C_TERM_INTEN",	"SI_NON_PREFERRED_DIGEST_AA_C_TERM_INTEN",	"SI_IND_NOT_CONNECTED_TO_TERMS",	"SI_IND_MISSED_CLEAVAGE",	
"SI_IND_N_FRAG_NOT_VIZ",	"SI_IND_C_FRAG_NOT_VIZ",	"SI_IND_N_INTEN",	"SI_IND_C_INTEN",	"SI_IND_N_NO_INTEN",	"SI_IND_C_NO_INTEN",	"SI_IND_N_N_TERM",	"SI_IND_N_C_TERM",	"SI_IND_N_Gap",	"SI_IND_N_Xle",	"SI_IND_N_Ala",	"SI_IND_N_Arg",	"SI_IND_N_Asn",	"SI_IND_N_Asp",	"SI_IND_N_Cys",	"SI_IND_N_Gln",	"SI_IND_N_Glu",	"SI_IND_N_Gly",	"SI_IND_N_His",	"SI_IND_N_Ile",	"SI_IND_N_Leu",	"SI_IND_N_Lys",	
"SI_IND_N_Met",	"SI_IND_N_Phe",	"SI_IND_N_Pro",	"SI_IND_N_Ser",	"SI_IND_N_Thr",	"SI_IND_N_Trp",	"SI_IND_N_Tyr",	"SI_IND_N_Val",	"SI_IND_N_N_TERM_INTEN",	"SI_IND_N_C_TERM_INTEN",	"SI_IND_N_Gap_INTEN",	"SI_IND_N_Xle_INTEN",	"SI_IND_N_Ala_INTEN",	"SI_IND_N_Arg_INTEN",	"SI_IND_N_Asn_INTEN",	"SI_IND_N_Asp_INTEN",	"SI_IND_N_Cys_INTEN",	"SI_IND_N_Gln_INTEN",	"SI_IND_N_Glu_INTEN",	"SI_IND_N_Gly_INTEN",	
"SI_IND_N_His_INTEN",	"SI_IND_N_Ile_INTEN",	"SI_IND_N_Leu_INTEN",	"SI_IND_N_Lys_INTEN",	"SI_IND_N_Met_INTEN",	"SI_IND_N_Phe_INTEN",	"SI_IND_N_Pro_INTEN",	"SI_IND_N_Ser_INTEN",	"SI_IND_N_Thr_INTEN",	"SI_IND_N_Trp_INTEN",	"SI_IND_N_Tyr_INTEN",	"SI_IND_N_Val_INTEN",	"SI_IND_N_SE",	"SI_SE_IND_HAS_N_FRAG_INTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	
"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_SE_IND_N_FRAG_DIFF_01",	"SI_SE_IND_N_FRAG_DIFF_05",	"SI_SE_IND_N_FRAG_DIFF_LARGE",	"SI_SE_IND_N_N_TERM_DIFF_INTEN",	"SI_SE_IND_N_C_TERM_DIFF_INTEN",	"SI_SE_IND_N_Gap_DIFF_INTEN",	"SI_SE_IND_N_Xle_DIFF_INTEN",	"SI_SE_IND_N_Ala_DIFF_INTEN",	"SI_SE_IND_N_Arg_DIFF_INTEN",	"SI_SE_IND_N_Asn_DIFF_INTEN",	"SI_SE_IND_N_Asp_DIFF_INTEN",	
"SI_SE_IND_N_Cys_DIFF_INTEN",	"SI_SE_IND_N_Gln_DIFF_INTEN",	"SI_SE_IND_N_Glu_DIFF_INTEN",	"SI_SE_IND_N_Gly_DIFF_INTEN",	"SI_SE_IND_N_His_DIFF_INTEN",	"SI_SE_IND_N_Ile_DIFF_INTEN",	"SI_SE_IND_N_Leu_DIFF_INTEN",	"SI_SE_IND_N_Lys_DIFF_INTEN",	"SI_SE_IND_N_Met_DIFF_INTEN",	"SI_SE_IND_N_Phe_DIFF_INTEN",	"SI_SE_IND_N_Pro_DIFF_INTEN",	"SI_SE_IND_N_Ser_DIFF_INTEN",	"SI_SE_IND_N_Thr_DIFF_INTEN",	
"SI_SE_IND_N_Trp_DIFF_INTEN",	"SI_SE_IND_N_Tyr_DIFF_INTEN",	"SI_SE_IND_N_Val_DIFF_INTEN",	"SI_SE_IND_HAS_NO_N_FRAG_INTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_IND_N_ME",	"SI_ME_IND_HAS_N_FRAG_INTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	
"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_ME_IND_N_FRAG_DIFF_01",	"SI_ME_IND_N_FRAG_DIFF_05",	"SI_ME_IND_N_FRAG_DIFF_LARGE",	"SI_ME_IND_HAS_NO_N_FRAG_INTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_IND_C_N_TERM",	"SI_IND_C_C_TERM",	"SI_IND_C_Gap",	"SI_IND_C_Xle",	
"SI_IND_C_Ala",	"SI_IND_C_Arg",	"SI_IND_C_Asn",	"SI_IND_C_Asp",	"SI_IND_C_Cys",	"SI_IND_C_Gln",	"SI_IND_C_Glu",	"SI_IND_C_Gly",	"SI_IND_C_His",	"SI_IND_C_Ile",	"SI_IND_C_Leu",	"SI_IND_C_Lys",	"SI_IND_C_Met",	"SI_IND_C_Phe",	"SI_IND_C_Pro",	"SI_IND_C_Ser",	"SI_IND_C_Thr",	"SI_IND_C_Trp",	"SI_IND_C_Tyr",	"SI_IND_C_Val",	"SI_IND_C_N_TERM_INTEN",	"SI_IND_C_C_TERM_INTEN",	"SI_IND_C_Gap_INTEN",	
"SI_IND_C_Xle_INTEN",	"SI_IND_C_Ala_INTEN",	"SI_IND_C_Arg_INTEN",	"SI_IND_C_Asn_INTEN",	"SI_IND_C_Asp_INTEN",	"SI_IND_C_Cys_INTEN",	"SI_IND_C_Gln_INTEN",	"SI_IND_C_Glu_INTEN",	"SI_IND_C_Gly_INTEN",	"SI_IND_C_His_INTEN",	"SI_IND_C_Ile_INTEN",	"SI_IND_C_Leu_INTEN",	"SI_IND_C_Lys_INTEN",	"SI_IND_C_Met_INTEN",	"SI_IND_C_Phe_INTEN",	"SI_IND_C_Pro_INTEN",	"SI_IND_C_Ser_INTEN",	"SI_IND_C_Thr_INTEN",	
"SI_IND_C_Trp_INTEN",	"SI_IND_C_Tyr_INTEN",	"SI_IND_C_Val_INTEN",	"SI_IND_C_SE",	"SI_SE_IND_HAS_C_FRAG_INTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_SE_IND_C_FRAG_DIFF_01",	"SI_SE_IND_C_FRAG_DIFF_05",	"SI_SE_IND_C_FRAG_DIFF_LARGE",	"SI_SE_IND_C_N_TERM_DIFF_INTEN",	"SI_SE_IND_C_C_TERM_DIFF_INTEN",	
"SI_SE_IND_C_Gap_DIFF_INTEN",	"SI_SE_IND_C_Xle_DIFF_INTEN",	"SI_SE_IND_C_Ala_DIFF_INTEN",	"SI_SE_IND_C_Arg_DIFF_INTEN",	"SI_SE_IND_C_Asn_DIFF_INTEN",	"SI_SE_IND_C_Asp_DIFF_INTEN",	"SI_SE_IND_C_Cys_DIFF_INTEN",	"SI_SE_IND_C_Gln_DIFF_INTEN",	"SI_SE_IND_C_Glu_DIFF_INTEN",	"SI_SE_IND_C_Gly_DIFF_INTEN",	"SI_SE_IND_C_His_DIFF_INTEN",	"SI_SE_IND_C_Ile_DIFF_INTEN",	"SI_SE_IND_C_Leu_DIFF_INTEN",	
"SI_SE_IND_C_Lys_DIFF_INTEN",	"SI_SE_IND_C_Met_DIFF_INTEN",	"SI_SE_IND_C_Phe_DIFF_INTEN",	"SI_SE_IND_C_Pro_DIFF_INTEN",	"SI_SE_IND_C_Ser_DIFF_INTEN",	"SI_SE_IND_C_Thr_DIFF_INTEN",	"SI_SE_IND_C_Trp_DIFF_INTEN",	"SI_SE_IND_C_Tyr_DIFF_INTEN",	"SI_SE_IND_C_Val_DIFF_INTEN",	"SI_SE_IND_HAS_NO_C_FRAG_INTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	
"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_IND_C_ME",	"SI_ME_IND_HAS_C_FRAG_INTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_ME_IND_C_FRAG_DIFF_01",	"SI_ME_IND_C_FRAG_DIFF_05",	"SI_ME_IND_C_FRAG_DIFF_LARGE",	"SI_ME_IND_HAS_NO_C_FRAG_INTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	
"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_NUM_FEATURES",	"ScoreModelFields_SI"};

const char * ScoreModelFields_SNI_names[]={
"SNI_CONST",	"SNI_IND_MIRROR1_VIZ",	"SNI_IND_MIRROR1_NOT_VIZ",	"SNI_IND_HAS_MIRROR1_INTEN",	"SNI_IND_MIRROR1_NO_INTEN",	"SNI_MIRROR1_ISO_LEVEL",	"SNI_IND_MIRROR1_HAS_MINUS_1",	"SNI_MIRROR1_MINUS_1_INTEN_DIFF",	"SNI_IND_MIRROR1_HAS_MINUS_2",	"SNI_MIRROR1_MINUS_2_INTEN_DIFF",	"SNI_IND_MIRROR1_HAS_PLUS_1",	"SNI_MIRROR1_PLUS_1_INTEN_DIFF",	"SNI_IND_MIRROR1_HAS_PLUS_2",	"SNI_MIRROR1_PLUS_2_INTEN_DIFF",	
"SNI_IND_MIRROR2_VIZ",	"SNI_IND_MIRROR2_NOT_VIZ",	"SNI_IND_HAS_MIRROR2_INTEN",	"SNI_IND_MIRROR2_NO_INTEN",	"SNI_MIRROR2_ISO_LEVEL",	"SNI_IND_MIRROR2_HAS_MINUS_1",	"SNI_MIRROR2_MINUS_1_INTEN_DIFF",	"SNI_IND_MIRROR2_HAS_MINUS_2",	"SNI_MIRROR2_MINUS_2_INTEN_DIFF",	"SNI_IND_MIRROR2_HAS_PLUS_1",	"SNI_MIRROR2_PLUS_1_INTEN_DIFF",	"SNI_IND_MIRROR2_HAS_PLUS_2",	"SNI_MIRROR2_PLUS_2_INTEN_DIFF",	
"SNI_IND_PARENT1_VIZ",	"SNI_IND_PARENT1_NOT_VIZ",	"SNI_IND_PARENT1_INTEN",	"SNI_IND_PARENT1_NO_INTEN",	"SNI_PARENT1_ISO_LEVEL",	"SNI_PARENT1_LOG_INTEN",	"SNI_PARENT1_LOG_GLOBAL_RANK",	"SNI_IND_PARENT2_VIZ",	"SNI_IND_PARENT2_NOT_VIZ",	"SNI_IND_PARENT2_INTEN",	"SNI_IND_PARENT2_NO_INTEN",	"SNI_PARENT2_ISO_LEVEL",	"SNI_PARENT2_LOG_INTEN",	"SNI_PARENT2_LOG_GLOBAL_RANK",	"SNI_IND_DIS_FROM_MINMAX_LESS_50",	
"SNI_DIS_FROM_MINMAX0",	"SNI_IND_DIS_FROM_MINMAX_LESS_150",	"SNI_DIS_FROM_MINMAX50",	"SNI_IND_DIS_FROM_MINMAX_LESS_250",	"SNI_DIS_FROM_MINMAX150",	"SNI_IND_DIS_FROM_MINMAX_MORE",	"SNI_DIS_FROM_MINMAX250",	"SNI_REL_POS0",	"SNI_REL_POS1",	"SNI_REL_POS2",	"SNI_REL_POS3",	"SNI_REL_POS4",	"SNI_REL_POS5",	"SNI_REL_POS6",	"SNI_REL_POS7",	"SNI_REL_POS8",	"SNI_REL_POS9",	"SNI_IND_N_IS_GAP",	
"SNI_IND_C_IS_GAP",	"SNI_N_TERM_CAT20",	"SNI_N_TERM_CAT18",	"SNI_N_TERM_CAT12",	"SNI_N_TERM_CAT8",	"SNI_N_TERM_CAT4",	"SNI_N_TERM_CAT2",	"SNI_N_EDGE_CAT20",	"SNI_N_EDGE_CAT18",	"SNI_N_EDGE_CAT12",	"SNI_N_EDGE_CAT8",	"SNI_N_EDGE_CAT4",	"SNI_N_EDGE_CAT2",	"SNI_C_TERM_CAT20",	"SNI_C_TERM_CAT18",	"SNI_C_TERM_CAT12",	"SNI_C_TERM_CAT8",	"SNI_C_TERM_CAT4",	"SNI_C_TERM_CAT2",	"SNI_C_EDGE_CAT20",	
"SNI_C_EDGE_CAT18",	"SNI_C_EDGE_CAT12",	"SNI_C_EDGE_CAT8",	"SNI_C_EDGE_CAT4",	"SNI_C_EDGE_CAT2",	"SNI_SPAN_CAT20",	"SNI_SPAN_CAT18",	"SNI_SPAN_CAT12",	"SNI_SPAN_CAT8",	"SNI_SPAN_CAT4",	"SNI_SPAN_CAT2",	"SNI_ND_SPAN_CAT20",	"SNI_ND_SPAN_CAT18",	"SNI_ND_SPAN_CAT12",	"SNI_ND_SPAN_CAT8",	"SNI_ND_SPAN_CAT4",	"SNI_ND_SPAN_CAT2",	"SNI_CD_SPAN_CAT20",	"SNI_CD_SPAN_CAT18",	"SNI_CD_SPAN_CAT12",	
"SNI_CD_SPAN_CAT8",	"SNI_CD_SPAN_CAT4",	"SNI_CD_SPAN_CAT2",	"SNI_IND_CONNECTS_TO_N_TERM",	"SNI_IND_CONNECTS_TO_C_TERM",	"SNI_IND_PREFERRED_DIGEST_AA_C_TERM",	"SNI_IND_PREFERRED_DIGEST_AA_N_TERM",	"SNI_IND_NOT_CONNECTED_TO_TERMS",	"SNI_IND_MISSED_CLEAVAGE",	"SNI_IND_N_NOT_VIZ",	"SNI_IND_C_NOT_VIZ",	"SNI_IND_N_INTEN",	"SNI_IND_N_NO_INTEN",	"SNI_IND_C_INTEN",	"SNI_IND_C_NO_INTEN",	
"SNI_IND_N_N_TERM",	"SNI_IND_N_C_TERM",	"SNI_IND_N_Gap",	"SNI_IND_N_Xle",	"SNI_IND_N_Ala",	"SNI_IND_N_Arg",	"SNI_IND_N_Asn",	"SNI_IND_N_Asp",	"SNI_IND_N_Cys",	"SNI_IND_N_Gln",	"SNI_IND_N_Glu",	"SNI_IND_N_Gly",	"SNI_IND_N_His",	"SNI_IND_N_Ile",	"SNI_IND_N_Leu",	"SNI_IND_N_Lys",	"SNI_IND_N_Met",	"SNI_IND_N_Phe",	"SNI_IND_N_Pro",	"SNI_IND_N_Ser",	"SNI_IND_N_Thr",	"SNI_IND_N_Trp",	
"SNI_IND_N_Tyr",	"SNI_IND_N_Val",	"SNI_IND_N_SE",	"SNI_SE_IND_HAS_N_FRAG_INTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_SE_IND_N_N_TERM_DIFF_INTEN",	"SNI_SE_IND_N_C_TERM_DIFF_INTEN",	"SNI_SE_IND_N_Gap_DIFF_INTEN",	"SNI_SE_IND_N_Xle_DIFF_INTEN",	"SNI_SE_IND_N_Ala_DIFF_INTEN",	
"SNI_SE_IND_N_Arg_DIFF_INTEN",	"SNI_SE_IND_N_Asn_DIFF_INTEN",	"SNI_SE_IND_N_Asp_DIFF_INTEN",	"SNI_SE_IND_N_Cys_DIFF_INTEN",	"SNI_SE_IND_N_Gln_DIFF_INTEN",	"SNI_SE_IND_N_Glu_DIFF_INTEN",	"SNI_SE_IND_N_Gly_DIFF_INTEN",	"SNI_SE_IND_N_His_DIFF_INTEN",	"SNI_SE_IND_N_Ile_DIFF_INTEN",	"SNI_SE_IND_N_Leu_DIFF_INTEN",	"SNI_SE_IND_N_Lys_DIFF_INTEN",	"SNI_SE_IND_N_Met_DIFF_INTEN",	"SNI_SE_IND_N_Phe_DIFF_INTEN",	
"SNI_SE_IND_N_Pro_DIFF_INTEN",	"SNI_SE_IND_N_Ser_DIFF_INTEN",	"SNI_SE_IND_N_Thr_DIFF_INTEN",	"SNI_SE_IND_N_Trp_DIFF_INTEN",	"SNI_SE_IND_N_Tyr_DIFF_INTEN",	"SNI_SE_IND_N_Val_DIFF_INTEN",	"SNI_SE_IND_HAS_NO_N_FRAG_INTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_IND_N_ME",	
"SNI_ME_IND_HAS_N_FRAG_INTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_ME_IND_HAS_NO_N_FRAG_INTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_IND_C_N_TERM",	"SNI_IND_C_C_TERM",	
"SNI_IND_C_Gap",	"SNI_IND_C_Xle",	"SNI_IND_C_Ala",	"SNI_IND_C_Arg",	"SNI_IND_C_Asn",	"SNI_IND_C_Asp",	"SNI_IND_C_Cys",	"SNI_IND_C_Gln",	"SNI_IND_C_Glu",	"SNI_IND_C_Gly",	"SNI_IND_C_His",	"SNI_IND_C_Ile",	"SNI_IND_C_Leu",	"SNI_IND_C_Lys",	"SNI_IND_C_Met",	"SNI_IND_C_Phe",	"SNI_IND_C_Pro",	"SNI_IND_C_Ser",	"SNI_IND_C_Thr",	"SNI_IND_C_Trp",	"SNI_IND_C_Tyr",	"SNI_IND_C_Val",	"SNI_IND_C_SE",	
"SNI_SE_IND_HAS_C_FRAG_INTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_SE_IND_C_N_TERM_DIFF_INTEN",	"SNI_SE_IND_C_C_TERM_DIFF_INTEN",	"SNI_SE_IND_C_Gap_DIFF_INTEN",	"SNI_SE_IND_C_Xle_DIFF_INTEN",	"SNI_SE_IND_C_Ala_DIFF_INTEN",	"SNI_SE_IND_C_Arg_DIFF_INTEN",	"SNI_SE_IND_C_Asn_DIFF_INTEN",	
"SNI_SE_IND_C_Asp_DIFF_INTEN",	"SNI_SE_IND_C_Cys_DIFF_INTEN",	"SNI_SE_IND_C_Gln_DIFF_INTEN",	"SNI_SE_IND_C_Glu_DIFF_INTEN",	"SNI_SE_IND_C_Gly_DIFF_INTEN",	"SNI_SE_IND_C_His_DIFF_INTEN",	"SNI_SE_IND_C_Ile_DIFF_INTEN",	"SNI_SE_IND_C_Leu_DIFF_INTEN",	"SNI_SE_IND_C_Lys_DIFF_INTEN",	"SNI_SE_IND_C_Met_DIFF_INTEN",	"SNI_SE_IND_C_Phe_DIFF_INTEN",	"SNI_SE_IND_C_Pro_DIFF_INTEN",	"SNI_SE_IND_C_Ser_DIFF_INTEN",	
"SNI_SE_IND_C_Thr_DIFF_INTEN",	"SNI_SE_IND_C_Trp_DIFF_INTEN",	"SNI_SE_IND_C_Tyr_DIFF_INTEN",	"SNI_SE_IND_C_Val_DIFF_INTEN",	"SNI_SE_IND_HAS_NO_C_FRAG_INTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_IND_C_ME",	"SNI_ME_IND_HAS_C_FRAG_INTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	
"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_ME_IND_HAS_NO_C_FRAG_INTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_NUM_FEATURES",	"ScoreModelFields_SNI"};


bool StrongFragmentModel::write_model(ostream &os) const
{
	if (! ind_has_models)
		return false;

	os << model_frag_idx << " " << mirror1_idx << " " << mirror2_idx << " " <<
		parent1_idx << " " << parent2_idx << endl;

	os << model_frag_charge << " " << mirror1_charge << " " <<  mirror2_charge << " " << 
		parent1_charge << " " << parent2_charge << endl;

	os << setprecision(4) << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;

	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		cout << "Model for frag " << config_->get_fragment(model_frag_idx).label << " had scaling problems: " << endl;
		cout << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;
		cout << "NOT writing the model, fix it!" << endl;
		exit(1);
	}

	inten_model.write_regression_model(os);
	no_inten_model.write_regression_model(os);

	return true;
}


bool StrongFragmentModel::read_model(const Config* config, istream& is, bool silent_ind)
{
	config_ = config;

	char buff[256];
	is.getline(buff,256);
	istringstream iss(buff);
	iss >> model_frag_idx >> mirror1_idx >> mirror2_idx >> parent1_idx >> parent2_idx;

	is.getline(buff,256);
	iss.str(buff);
	iss >> model_frag_charge >> mirror1_charge >> mirror2_charge >> parent1_charge >> parent2_charge;

	is.getline(buff,256);
	sscanf(buff,"%f %f",&inten_log_scaling_factor,&no_inten_log_scaling_factor);

	inten_model.read_regression_model(is);
	no_inten_model.read_regression_model(is);

	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		cout << "Model for frag " << config->get_fragment(model_frag_idx).label << " had scaling problems: " << endl;
		cout << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;
		exit(1);
	}

	if (! inten_model.get_has_weights() || ! no_inten_model.get_has_weights())
	{
		
		cout << "Model for frag " << config->get_fragment(model_frag_idx).label << " has no weights!" << endl;
		exit(1);
	}

	ind_has_models = true;

	return true;
}



void StrongFragmentModel::fill_constant_vals(
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage, 
							   vector<fval>& f_vals) const
{
	const mass_t iso_tolerance = spec->getConfig()->getTolerance()*0.5;
	FragStats frag_stats;
	FragStats parent1_stats, parent2_stats;
	FragStats mirror1_stats, mirror2_stats;

	frag_stats.fill_from_breakage(breakage,spec,model_frag_idx);
	if (parent1_idx>=0)
		parent1_stats.fill_from_breakage(breakage,spec,parent1_idx);
	if (parent2_idx>=0)
		parent2_stats.fill_from_breakage(breakage,spec,parent2_idx);
	if (mirror1_idx>=0)
		mirror1_stats.fill_from_breakage(breakage,spec,mirror1_idx);
	if (mirror2_idx>=0)
		mirror2_stats.fill_from_breakage(breakage,spec,mirror2_idx);

	f_vals.clear();
	if (frag_stats.has_intensity) // fill features for visible fragment
	{
		f_vals.push_back(fval(SI_CONST,1.0));

		const float log_inten = frag_stats.log_intensity;

		// Mirror1 features
		if (mirror1_idx>=0)
		{
			if (mirror1_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_MIRROR1_VIZ,1.0));
				if (mirror1_stats.has_intensity)
				{
					f_vals.push_back(fval(SI_IND_HAS_MIRROR1_INTEN,1.0));
					if (mirror1_stats.iso_level!=0)
						f_vals.push_back(fval(SI_MIRROR1_ISO_LEVEL,mirror1_stats.iso_level));
				
					vector<float> iso_intens;
					spec->findIsotopicEnvelope(mirror1_stats.peak_idx, iso_intens, iso_tolerance, mirror1_charge);	
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR1_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR1_MINUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR1_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR1_MINUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR1_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR1_PLUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR1_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR1_PLUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[3]));	
						}
					}
						

					const mass_t sum_masses = frag_stats.mass * model_frag_charge + 
										mirror1_stats.mass * mirror1_charge;

					const mass_t offset = fabs(pm_with_19 - sum_masses + (model_frag_charge + mirror1_charge -1)*MASS_PROTON);
					const float offset_level = offset * one_over_tolerance;

					if (offset_level<0.25)
					{
						f_vals.push_back(fval(SI_MIRROR1_MASS_DIFF25,1.0));
					} 
					else if (offset_level<0.75)
					{
						f_vals.push_back(fval(SI_MIRROR1_MASS_DIFF75,1.0));
					}
					else
						f_vals.push_back(fval(SI_MIRROR1_MASS_DIFF_LARGE,1.0));
				}
			}
			else
				f_vals.push_back(fval(SI_IND_MIRROR1_NOT_VIZ,1.0));
		}

		// Mirror2 features
		if (mirror2_idx>=0)
		{
			if (mirror2_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_MIRROR2_VIZ,1.0));
				if (mirror2_stats.has_intensity)
				{
					f_vals.push_back(fval(SI_IND_HAS_MIRROR2_INTEN,1.0));
					if (mirror2_stats.iso_level != 0)
						f_vals.push_back(fval(SI_MIRROR2_ISO_LEVEL,mirror2_stats.iso_level));

					vector<float> iso_intens;
					spec->findIsotopicEnvelope(mirror2_stats.peak_idx, iso_intens, iso_tolerance, mirror2_charge);
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR2_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR2_MINUS_1_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR2_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR2_MINUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR2_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR2_PLUS_1_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR2_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR2_PLUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[3]));	
						}
					}



					const mass_t sum_masses = frag_stats.mass * model_frag_charge + 
										mirror2_stats.mass * mirror2_charge;

					const mass_t offset = fabs(pm_with_19 - sum_masses + (model_frag_charge + mirror2_charge -1)*MASS_PROTON);
					const float offset_level = offset * one_over_tolerance;

					if (offset_level<0.25)
					{
						f_vals.push_back(fval(SI_MIRROR2_MASS_DIFF25,1.0));
					} 
					else if (offset_level<0.75)
					{
						f_vals.push_back(fval(SI_MIRROR2_MASS_DIFF75,1.0));
					}
					else
						f_vals.push_back(fval(SI_MIRROR2_MASS_DIFF_LARGE,1.0));
				}
			}
			else
				f_vals.push_back(fval(SI_IND_MIRROR2_NOT_VIZ,1.0));
		}

		// Parent 1 features
		if (parent1_idx>=0)
		{
			if (parent1_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_PARENT1_VIZ,1.0));
				if (parent1_stats.has_intensity)
				{
					const float inten_diff = parent1_stats.log_intensity - log_inten;
					const mass_t dis_min = parent1_stats.mass - spec->get_min_peak_mass();
					const mass_t dis_max = spec->get_max_peak_mass() - parent1_stats.mass;
					const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);
					
					if (parent1_stats.iso_level != 0)
						f_vals.push_back(fval(SI_PARENT1_ISO_LEVEL,parent1_stats.iso_level));
					if (dis<100)
					{
						f_vals.push_back(fval(SI_IND_PARENT1_LESS_THAN_100_MIN_MAX,1.0));
					} 
					else if (dis<200)
						f_vals.push_back(fval(SI_IND_PARENT1_LESS_THAN_200_MIN_MAX,1.0));

					if (inten_diff>0)
					{
						f_vals.push_back(fval(SI_IND_PARENT1_INTEN_MORE,1.0));
						f_vals.push_back(fval(SI_PARENT1_INTEN_DIFF_MORE,inten_diff));
					}
					else
					{
						f_vals.push_back(fval(SI_IND_PARENT1_INTEN_LESS,1.0));
						f_vals.push_back(fval(SI_PARENT1_INTEN_DIFF_LESS,inten_diff));
					}
				}
				else
					f_vals.push_back(fval(SI_IND_PARENT1_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SI_IND_PARENT1_NOT_VIZ,1.0));
		}

		// Parent 2 features
		if (parent2_idx>=0)
		{
			if (parent2_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_PARENT2_VIZ,1.0));
				if (parent2_stats.has_intensity)
				{
					const float inten_diff = parent2_stats.log_intensity - log_inten;
					const mass_t dis_min = parent2_stats.mass - spec->get_min_peak_mass();
					const mass_t dis_max = spec->get_max_peak_mass() - parent2_stats.mass;
					const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);
					
					if (parent2_stats.iso_level != 0)
						f_vals.push_back(fval(SI_PARENT2_ISO_LEVEL,parent2_stats.iso_level));
					if (dis<100)
					{
						f_vals.push_back(fval(SI_IND_PARENT2_LESS_THAN_100_MIN_MAX,1.0));
					} 
					else if (dis<200)
						f_vals.push_back(fval(SI_IND_PARENT2_LESS_THAN_200_MIN_MAX,1.0));

					if (inten_diff>0)
					{
						f_vals.push_back(fval(SI_IND_PARENT2_INTEN_MORE,1.0));
						f_vals.push_back(fval(SI_PARENT2_INTEN_DIFF_MORE,inten_diff));
					}
					else
					{
						f_vals.push_back(fval(SI_IND_PARENT2_INTEN_LESS,1.0));
						f_vals.push_back(fval(SI_PARENT2_INTEN_DIFF_LESS,inten_diff));
					}
				}
				else
					f_vals.push_back(fval(SI_IND_PARENT2_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SI_IND_PARENT2_NOT_VIZ,1.0));
		}

		// self intensity
		f_vals.push_back(fval(SI_LOG_LOCAL_RANK,frag_stats.log_local_rank));
		f_vals.push_back(fval(SI_LOG_GLOBAL_RANK,frag_stats.log_global_rank));
		if (frag_stats.iso_level != 0)
			f_vals.push_back(fval(SI_ISO_LEVEL,frag_stats.iso_level));
		
		if (log_inten<1.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS1,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS1,log_inten));
		} 
		else if (log_inten<2.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS2,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS2,log_inten-1.0));
		}
		else if (log_inten<3.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS3,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS3,log_inten-2.0));
		}
		else if (log_inten<4.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS4,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS4,log_inten-3.0));
		}
		else
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_MORE,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_MORE,log_inten-4.0));
		}

		// self distance
		const mass_t dis_min = frag_stats.mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - frag_stats.mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		if (dis<50)
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX0,dis));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_50,log_inten));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX50,dis-50.0));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_150,log_inten));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX150,dis-150.0));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_250,log_inten));
		}
		else
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX250,dis-250.0));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_MORE,log_inten));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(SI_REL_POS0+rel_pos,1.0));

		const float one_over_model_frag_charge = 1.0 / model_frag_charge;
		const mass_t forward_offsets[]={MASS_NH3,MASS_H2O,MASS_CO,MASS_H2ONH3,MASS_H2OH2O};
		const int num_forward_offsets = sizeof(forward_offsets)/sizeof(mass_t);
		const mass_t peak_mass = frag_stats.mass;

		int t;
		for (t=0; t<num_forward_offsets; t++)
		{
			const int forward_idx = spec->findPeakWithMaxIntensity(peak_mass + forward_offsets[t]*one_over_model_frag_charge,frag_tolerance);
			if (forward_idx>0)
			{
				f_vals.push_back(fval(SI_IND_HAS_PLUS_NH3+2*t,1.0));
				f_vals.push_back(fval(SI_IND_HAS_PLUS_NH3+2*t+1,spec->get_peak_log_intensity(forward_idx)-log_inten));
			}
		}

		
		const int plus_idx = spec->findPeakWithMaxIntensity((peak_mass * model_frag_charge + MASS_PROTON)/(model_frag_charge+1.0),frag_tolerance);
		if (plus_idx>=0)
		{
			f_vals.push_back(fval(SI_IND_HAS_CHARGE_PLUS1,1.0));
			f_vals.push_back(fval(SI_CHARGE_PLUS1_INTEN_DIFF,spec->get_peak_log_intensity(plus_idx)-log_inten));
		}

		if (model_frag_charge>1)
		{
			const int minus_idx = spec->findPeakWithMaxIntensity((peak_mass * model_frag_charge - MASS_PROTON)/(model_frag_charge-1.0),frag_tolerance);
			if (minus_idx>=0)
			{
				f_vals.push_back(fval(SI_IND_HAS_CHARGE_MINUS1,1.0));
				f_vals.push_back(fval(SI_CHARGE_MINUS1_INTEN_DIFF,spec->get_peak_log_intensity(minus_idx)-log_inten));
			}
		}
	}
	else // Fill features for non-visible
	{
		f_vals.push_back(fval(SNI_CONST,1.0));

		// Mirror1 features
		if (mirror1_idx>=0)
		{
			if (mirror1_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_MIRROR1_VIZ,1.0));
				if (mirror1_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_HAS_MIRROR1_INTEN,1.0));
					if (mirror1_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_MIRROR1_ISO_LEVEL,mirror1_stats.iso_level));

					vector<float> iso_intens;
					spec->findIsotopicEnvelope(mirror1_stats.peak_idx, iso_intens, iso_tolerance, mirror1_charge);
					
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR1_MINUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR1_MINUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR1_PLUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR1_PLUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[3]));	
						}
					}
				}
				else
					f_vals.push_back(fval(SNI_IND_MIRROR1_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_MIRROR1_NOT_VIZ,1.0));
		}

		// Mirror2 features
		if (mirror2_idx>=0)
		{
			if (mirror2_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_MIRROR2_VIZ,1.0));
				if (mirror2_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_HAS_MIRROR2_INTEN,1.0));
					if (mirror2_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_MIRROR2_ISO_LEVEL,mirror2_stats.iso_level));

					vector<float> iso_intens;
					spec->findIsotopicEnvelope(mirror2_stats.peak_idx, iso_intens, iso_tolerance, mirror2_charge);
					
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR1_MINUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR2_MINUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR2_PLUS_1_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR2_PLUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[3]));	
						}
					}
				}
				else
					f_vals.push_back(fval(SNI_IND_MIRROR2_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_MIRROR2_NOT_VIZ,1.0));
		}

		// Parent 1 features
		if (parent1_idx>=0)
		{
			if (parent1_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_PARENT1_VIZ,1.0));
				if (parent1_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_PARENT1_INTEN,1.0));
					if (parent1_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_PARENT1_ISO_LEVEL,parent1_stats.iso_level));
					f_vals.push_back(fval(SNI_PARENT1_LOG_INTEN,parent1_stats.log_intensity));
					f_vals.push_back(fval(SNI_PARENT1_LOG_GLOBAL_RANK,parent1_stats.log_global_rank));
				}
				else
					f_vals.push_back(fval(SNI_IND_PARENT1_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_PARENT1_NOT_VIZ,1.0));
		}

		// Parent 2 features
		if (parent2_idx>=0)
		{
			if (parent2_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_PARENT2_VIZ,1.0));
				if (parent2_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_PARENT2_INTEN,1.0));
					if (SNI_PARENT2_ISO_LEVEL,parent2_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_PARENT2_ISO_LEVEL,parent2_stats.iso_level));
					f_vals.push_back(fval(SNI_PARENT2_LOG_INTEN,parent2_stats.log_intensity));
					f_vals.push_back(fval(SNI_PARENT2_LOG_GLOBAL_RANK,parent2_stats.log_global_rank));
				}
				else
					f_vals.push_back(fval(SNI_IND_PARENT2_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_PARENT2_NOT_VIZ,1.0));
		}

		// self distance
		const mass_t expected_mass = config_->get_fragment(model_frag_idx).calc_expected_mass(breakage->mass,pm_with_19);
		const mass_t dis_min = expected_mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - expected_mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		if (dis<50)
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX0,dis));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX50,dis-50.0));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX150,dis-150.0));
		}
		else
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX250,dis-250.0));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(SNI_REL_POS0+rel_pos,1.0));
	}
}



/**********************************************************************************
***********************************************************************************/
void StrongFragmentModel::fill_aa_variable_vals(
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage,
							   const BreakageInfo* info,
							   vector<fval>& f_vals) const
{
	const vector<int>& org_aas = config_->get_org_aa();
	const int n_aa = (info->n_aa>=0 ? org_aas[info->n_aa] : Gap);
	const int c_aa = (info->c_aa>=0 ? org_aas[info->c_aa] : Gap);
	const int pos = breakage->get_position_of_frag_idx(model_frag_idx);
	
	const bool do_n_features = (n_aa != Gap);
	const bool do_c_features = (c_aa != Gap);

	const FragmentType& fragment = config_->get_fragment(model_frag_idx);
	const mass_t exp_peak_mass = (info->breakage ? 
		fragment.calc_expected_mass(info->breakage->mass,pm_with_19) : NEG_INF);
	const mass_t exp_n_peak_mass = (do_n_features && info->n_break ? 
		fragment.calc_expected_mass(info->n_break->mass,pm_with_19) : NEG_INF);
	const mass_t exp_c_peak_mass = (do_c_features && info->c_break ? 
		fragment.calc_expected_mass(info->c_break->mass,pm_with_19) : NEG_INF);	
	
	static vector<int> threshes;
	if (threshes.size()==0)
	{
		threshes.push_back(18);
		threshes.push_back(12);
		threshes.push_back(8);
		threshes.push_back(4);
		threshes.push_back(2);
		threshes.push_back(NEG_INF);
	}

	// fill intensity
	if (pos>=0)
	{
		if (! do_n_features)
			f_vals.push_back(fval(SI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(SI_IND_C_IS_GAP,1.0));

		// aa category feature
		if (do_n_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->n_side_cat>threshes[k])
					break;
			if (info->connects_to_N_term)
			{
				f_vals.push_back(fval(SI_N_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SI_N_EDGE_CAT20+k,1.0));

			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;
		}

		if (do_c_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;

			if (info->connects_to_C_term)
			{
				f_vals.push_back(fval(SI_C_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SI_C_EDGE_CAT20+k,1.0));

			if (do_n_features)
			{
				int k;
				for (k=0; k<threshes.size(); k++)
					if (info->span_cat>threshes[k])
						break;
				f_vals.push_back(fval(SI_SPAN_CAT20+k,1.0));

				if (info->n_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->n_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SI_ND_SPAN_CAT20+k,1.0));
				}

				if (info->c_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->c_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SI_CD_SPAN_CAT20+k,1.0));
				}
			}
		}

		const int peakIdx = breakage->fragments[pos].peak_idx;
		const Peak& peak = spec->getPeak(peakIdx);
		const float peakLogIntensity = spec->get_peak_log_intensity(peakIdx);

		if (do_n_features)
		{
			if (info->connects_to_N_term)
				f_vals.push_back(fval(SI_IND_CONNECTS_TO_N_TERM,1.0));

			if (info->preferred_digest_aa_N_term)
				f_vals.push_back(fval(SI_IND_PREFERRED_DIGEST_AA_N_TERM,1.0));
		}

		if (do_c_features)
		{
			if (info->connects_to_C_term)
				f_vals.push_back(fval(SI_IND_CONNECTS_TO_C_TERM,1.0));

			if (info->preferred_digest_aa_C_term)
				f_vals.push_back(fval(SI_IND_PREFERRED_DIGEST_AA_C_TERM,1.0));
		}

		if (!info->connects_to_N_term && ! info->connects_to_C_term)
			f_vals.push_back(fval(SI_IND_NOT_CONNECTED_TO_TERMS,1.0));
		
		if (info->missed_cleavage)
			f_vals.push_back(fval(SI_IND_MISSED_CLEAVAGE,1.0));
	

		if (do_n_features)
		{
			if (! info->n_break || ! info->n_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SI_IND_N_FRAG_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SI_IND_N_N_TERM,1.0));
				f_vals.push_back(fval(SI_IND_N_N_TERM+n_aa,1.0));
				
				const int n_pos=info->n_break->get_position_of_frag_idx(model_frag_idx);
				if (n_pos>=0)
				{
					f_vals.push_back(fval(SI_IND_N_INTEN,1.0));
					f_vals.push_back(fval(SI_IND_N_N_TERM_INTEN,peakLogIntensity));
					f_vals.push_back(fval(SI_IND_N_N_TERM_INTEN+n_aa,peakLogIntensity));
				}
				else
					f_vals.push_back(fval(SI_IND_N_NO_INTEN,1.0));

				const mass_t dis_min =exp_n_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_n_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->n_edge_is_single)
				{
					f_vals.push_back(fval(SI_IND_N_SE,1.0));
					if (n_pos>=0)
					{
						const Peak& n_peak = spec->getPeak(info->n_break->fragments[n_pos].peak_idx);
						const float nPosPeakLogIntensity = spec->get_peak_log_intensity(info->n_break->fragments[n_pos].peak_idx);
						const mass_t n_diff = fabs(fabs((peak.mass - n_peak.mass)* model_frag_charge) - info->exp_n_edge_mass);

						f_vals.push_back(fval(SI_SE_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (n_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_N_FRAG_DIFF_01,1.0));
						}
						else if (n_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_N_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_SE_IND_N_FRAG_DIFF_LARGE,1.0));

						const float diff_inten = peakLogIntensity - nPosPeakLogIntensity;
						f_vals.push_back(fval(SI_SE_IND_N_N_TERM_DIFF_INTEN,diff_inten));
						f_vals.push_back(fval(SI_SE_IND_N_N_TERM_DIFF_INTEN+n_aa,diff_inten));
					}
					else
					{
					
						f_vals.push_back(fval(SI_SE_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SI_IND_N_ME,1.0));
					if (n_pos>=0)
					{
						const Peak& n_peak = spec->getPeak(info->n_break->fragments[n_pos].peak_idx);
						const mass_t n_diff = fabs(fabs((peak.mass - n_peak.mass)* model_frag_charge) - info->exp_n_edge_mass);

						f_vals.push_back(fval(SI_ME_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (n_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_N_FRAG_DIFF_01,1.0));
						}
						else if (n_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_N_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_ME_IND_N_FRAG_DIFF_LARGE,1.0));
					}
					else
					{
						f_vals.push_back(fval(SI_ME_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}

		if (do_c_features)
		{
			if (! info->c_break || ! info->c_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SI_IND_C_FRAG_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SI_IND_C_N_TERM,1.0));
				f_vals.push_back(fval(SI_IND_C_N_TERM+c_aa,1.0));
				
				const int c_pos=info->c_break->get_position_of_frag_idx(model_frag_idx);
				if (c_pos>=0)
				{
					f_vals.push_back(fval(SI_IND_C_INTEN,1.0));
					f_vals.push_back(fval(SI_IND_C_N_TERM_INTEN,peakLogIntensity));
					f_vals.push_back(fval(SI_IND_C_N_TERM_INTEN+c_aa,peakLogIntensity));
				}
				else
					f_vals.push_back(fval(SI_IND_C_NO_INTEN,1.0));

				const mass_t dis_min =exp_c_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_c_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->c_edge_is_single)
				{
					f_vals.push_back(fval(SI_IND_C_SE,1.0));
					if (c_pos>=0)
					{
						const Peak& c_peak = spec->getPeak(info->c_break->fragments[c_pos].peak_idx);
						const float cPosPeakLogIntensity = spec->get_peak_log_intensity(info->c_break->fragments[c_pos].peak_idx);
						const mass_t c_diff = fabs(fabs((peak.mass - c_peak.mass)* model_frag_charge) - info->exp_c_edge_mass);

						f_vals.push_back(fval(SI_SE_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (c_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_C_FRAG_DIFF_01,1.0));
						}
						else if (c_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_C_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_SE_IND_C_FRAG_DIFF_LARGE,1.0));

						const float diff_inten = peakLogIntensity - cPosPeakLogIntensity;
						f_vals.push_back(fval(SI_SE_IND_C_N_TERM_DIFF_INTEN,diff_inten));
						f_vals.push_back(fval(SI_SE_IND_C_N_TERM_DIFF_INTEN+c_aa,diff_inten));
					}
					else
					{
						f_vals.push_back(fval(SI_SE_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SI_IND_C_ME,1.0));
					if (c_pos>=0)
					{
						const Peak& c_peak = spec->getPeak(info->c_break->fragments[c_pos].peak_idx);
						const mass_t c_diff = fabs(fabs((peak.mass - c_peak.mass)* model_frag_charge) - info->exp_c_edge_mass);

						f_vals.push_back(fval(SI_ME_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (c_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_C_FRAG_DIFF_01,1.0));
						}
						else if (c_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_C_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_ME_IND_C_FRAG_DIFF_LARGE,1.0));
					}
					else
					{
						f_vals.push_back(fval(SI_ME_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}
	}
	else  // fill features for no intensity
	{
		if (! do_n_features)
			f_vals.push_back(fval(SNI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(SNI_IND_C_IS_GAP,1.0));

		// aa category feature
		if (do_n_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->n_side_cat>threshes[k])
					break;
			if (info->connects_to_N_term)
			{
				f_vals.push_back(fval(SNI_N_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SNI_N_EDGE_CAT20+k,1.0));

			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;
		}

		if (do_c_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;

			if (info->connects_to_C_term)
			{
				f_vals.push_back(fval(SNI_C_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SNI_C_EDGE_CAT20+k,1.0));

			for (k=0; k<threshes.size(); k++)
				if (info->span_cat>threshes[k])
					break;
		
			if (do_n_features)
			{
				f_vals.push_back(fval(SNI_SPAN_CAT20+k,1.0));

				if (info->n_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->n_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SNI_ND_SPAN_CAT20+k,1.0));
				}

				if (info->c_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->c_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SNI_CD_SPAN_CAT20+k,1.0));
				}
			}
		}

		if (do_n_features)
		{
			if (info->connects_to_N_term)
				f_vals.push_back(fval(SNI_IND_CONNECTS_TO_N_TERM,1.0));

			if (info->preferred_digest_aa_N_term)
				f_vals.push_back(fval(SNI_IND_PREFERRED_DIGEST_AA_N_TERM,1.0));
		}

		if (do_c_features)
		{
			if (info->connects_to_C_term)
				f_vals.push_back(fval(SNI_IND_CONNECTS_TO_C_TERM,1.0));

			if (info->preferred_digest_aa_C_term)
				f_vals.push_back(fval(SNI_IND_PREFERRED_DIGEST_AA_C_TERM,1.0));
		}


		if (! info->connects_to_N_term && ! info->connects_to_C_term)
			f_vals.push_back(fval(SNI_IND_NOT_CONNECTED_TO_TERMS,1.0));

		if (info->missed_cleavage)
			f_vals.push_back(fval(SNI_IND_MISSED_CLEAVAGE,1.0));
	

		if (do_n_features)
		{
			if (! info->n_break || ! info->n_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SNI_IND_N_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SNI_IND_N_N_TERM,1.0));
				f_vals.push_back(fval(SNI_IND_N_N_TERM+n_aa,1.0));
			
				const int n_pos=info->n_break->get_position_of_frag_idx(model_frag_idx);
				if (n_pos>=0)
				{
					f_vals.push_back(fval(SNI_IND_N_INTEN,1.0));
				}
				else
					f_vals.push_back(fval(SNI_IND_N_NO_INTEN,1.0));

				const mass_t dis_min = exp_n_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_n_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->n_edge_is_single)
				{
					f_vals.push_back(fval(SNI_IND_N_SE,1.0));
					if (n_pos>=0)
					{
						const Peak& n_peak = spec->getPeak(info->n_break->fragments[n_pos].peak_idx);
						const float nPosPeakLogIntensity = spec->get_peak_log_intensity(info->n_break->fragments[n_pos].peak_idx);

						f_vals.push_back(fval(SNI_SE_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));

						f_vals.push_back(fval(SNI_SE_IND_N_N_TERM_DIFF_INTEN,nPosPeakLogIntensity));
						f_vals.push_back(fval(SNI_SE_IND_N_N_TERM_DIFF_INTEN+n_aa,nPosPeakLogIntensity));
					}
					else
					{
					
						f_vals.push_back(fval(SNI_SE_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SNI_IND_N_ME,1.0));
					if (n_pos>=0)
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
					}
					else
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}


		if (do_c_features)
		{
			if (! info->c_break || ! info->c_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SNI_IND_C_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SNI_IND_C_N_TERM,1.0));
				f_vals.push_back(fval(SNI_IND_C_N_TERM+c_aa,1.0));
			
				const int c_pos=info->c_break->get_position_of_frag_idx(model_frag_idx);
				if (c_pos>=0)
				{
					f_vals.push_back(fval(SNI_IND_C_INTEN,1.0));
				}
				else
					f_vals.push_back(fval(SNI_IND_C_NO_INTEN,1.0));

				const mass_t dis_min = exp_c_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_c_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->c_edge_is_single)
				{
					f_vals.push_back(fval(SNI_IND_C_SE,1.0));
					if (c_pos>=0)
					{
						const Peak& c_peak = spec->getPeak(info->c_break->fragments[c_pos].peak_idx);
						const float cPosPeakLogIntensity = spec->get_peak_log_intensity(info->c_break->fragments[c_pos].peak_idx);
					
						f_vals.push_back(fval(SNI_SE_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));

						f_vals.push_back(fval(SNI_SE_IND_C_N_TERM_DIFF_INTEN,cPosPeakLogIntensity));
						f_vals.push_back(fval(SNI_SE_IND_C_N_TERM_DIFF_INTEN+c_aa,cPosPeakLogIntensity));
					}
					else
					{
					
						f_vals.push_back(fval(SNI_SE_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SNI_IND_C_ME,1.0));
					if (c_pos>=0)
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
					}
					else
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}
	}
}






void StrongFragmentModel::fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const
{
	vector<fval> const_vals;

	fill_constant_vals(spec,pm_with_19,breakage,const_vals);
	samples.resize(infos.size());
	int i;
	for (i=0; i<infos.size(); i++)
	{
		vector< fval > var_vals;
		fill_aa_variable_vals(spec,pm_with_19,breakage,&infos[i],var_vals);

		samples[i].f_vals = const_vals;
		int j;
		for (j=0; j<var_vals.size(); j++)
			samples[i].f_vals.push_back(var_vals[j]);
	}
}


void StrongFragmentModel::fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const
{

	fill_constant_vals(spec,pm_with_19,breakage,f_vals);
	fill_aa_variable_vals(spec,pm_with_19,breakage,&info,f_vals);
	sort(f_vals.begin(),f_vals.end());
}




