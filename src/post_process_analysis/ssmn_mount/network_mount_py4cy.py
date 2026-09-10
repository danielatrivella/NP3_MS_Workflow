#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 17:24:52 2021

@author: Cris e Luiz e Rafa

Script para montar automaticamente a rede SSMN e SSMN [M+H]+ no cytoscape
Caso existam colunas de COR_, Cria o Ranking de Correlação (Spearman) em formato decimal (Rank + Repeticoes/1000) - Codigo Enzo

Argumentos preliminares:
    selfloop SSMN
    selfloop SSMN [M+H]+
    clean table csv
    numero de bioscores - numero de linhas para pular, caso correlacao
    # TODO add estilo cytoscape
    
Passos:
  1. Ler a tabela clean e aplicar ranking
  2. Importa a SSMN
  3. Importa a SSMN [M+H]+
  4. Importa a tabela clean
  5. (Manual) Salvar projeto do cytoscape
"""

import py4cytoscape as py4
import pandas as pd
import argparse
import numpy as np
import requests
from pathlib import Path
import sys
from openpyxl.workbook import Workbook

def apply_corr_ranking(df):
    """
    Função que aplica o ranking decrescente codificado em formato decimal:
    Valor = Rank + (Contagem_de_Repeticoes / 1000.0)

    Exemplo:
    - Rank 4 com 10 repeticoes -> 4.010
    - Rank 4 com 1 repeticao   -> 4.001

    Retorna o dataframe modificado caso foram encontradas colunas de cor_
    """
    # Detecta colunas COR_
    cols_cor_ = [col for col in df.columns if col.startswith("COR_")]
    
    if not cols_cor_:
        print("AVISO: Nenhuma coluna iniciando com 'COR_' foi encontrada. O ranking nao sera gerado.")
        return df
    
    print("* Calculando rankings decrescentes com contagem decimal (3 casas) para empates *")
    
    # Processamento do Ranking Decimal
    for col in cols_cor_:
        col_ranked_name = col + "_RANKED"
        
        # se nenhum valor not NaN, skip
        if not df[col].notna().any():
            continue
        
        # 1. Calcula o ranking decrescente (metodo denso)
        ranks = df[col].rank(ascending=False, method="dense")
        # 2. Conta quantas vezes cada ranking se repete
        counts = ranks.map(ranks.value_counts())
        # 3. Seta nova coluna com codificacao em formato decimal (Rank + Contagem / 1000.0)
        df.loc[:,col_ranked_name] = ranks + (counts / 1000.0)
        
    print("  * Done! Rankings calculados com sucesso! *")
    return df

def create_collaborators_tables(df, clean_table_file, filename_suffix="protonated_collaborators"):
    print("* Creating table for the collaborators in CSV and Excel *")
    # Criar coluna 'curated_superclass' caso nao exista, preenchendo com valor padrao vazio
    if 'best_origin_curated_superclass_grouping' not in df.columns:
        df['best_origin_curated_superclass_grouping'] = ''
    
    # Passo 1: remove blanks e beds
    # Filtro: manter apenas linhas com BLANKS_TOTAL == 0
    if 'BLANKS_TOTAL' in df.columns:
        df = df.loc[df['BLANKS_TOTAL'] == 0,:].copy()
    # Se a coluna BEDS_TOTAL existir, aplicar o filtro adicional
    if 'BEDS_TOTAL' in df.columns:
        df = df.loc[df['BEDS_TOTAL'] == 0,:].copy()
    
    # Passo 2a
    # Remove todas as barras '/' da string
    df['best_origin_SMILES'] = df['best_origin_SMILES'].str.replace('/', '', regex=False)

    # Passo 2b - Remover linhas com best_origin_SMILES vazio ou nulo
    df = df[df['best_origin_SMILES'].notna() & (df['best_origin_SMILES'] != '')].copy()
    
    if df.shape[0] == 0:
        print("  - WARNING: No valid 'best_origin_SMILES' left, collaborators table will not be created. ")
        return
    
    # %%
    # Passo 2c - Novo DataFrame com colunas especificas
    count_cols = df.columns.str.endswith("_area")
    if not count_cols.any():
        count_cols = df.columns.str.endswith("_spectra")
    count_cols = df.columns[count_cols].values
    colunas_selection = np.asarray([
        'msclusterID',
        'mzConsensus',
        'rtMean_minutes',
        'sumInts',
        'basePeakInt',
        'curated_identification_best_origin']+ list(count_cols) + ['best_origin_curated_superclass_grouping',
        'best_origin_SMILES',
        'gnps_Smiles',
        'gnps_SpectrumID',
        'gnps_Compound_Name',
        'gnps_MZErrorPPM',
        'gnps_MQScore',
        'gnps_SharedPeaks',
        'gnps_Organism',
        'gnps_npclassifier_superclass',
        'gnps_npclassifier_class',
        'gnps_npclassifier_pathway',
        'gnps_category',
        'gnps_score',
        'tremolo_SMILES_best',
        'tremolo_UNPD_IDs_best',
        'tremolo_chemicalNames_best',
        'tremolo_mzErrorPPM_best',
        'tremolo_MQScore_best',
        'tremolo_numSharedPeaks_best',
        'tremolo_NPClassifier_superclass_best',
        'tremolo_NPClassifier_class_best',
        'tremolo_NPClassifier_pathway_best',
        'tremolo_UNPD_category_best',
        'tremolo_UNPD_score_best',
        'tanimoto_unpd_gnps'])
    # filtrar apenas as colunas selecionadas que aparecem no df
    df_smiles = df.loc[:,colunas_selection[np.isin(colunas_selection, df.columns)]].copy().reset_index(drop=True)
    # %%
    # Ordenar pelo valor decrescente da coluna "sumInts"
    df_smiles = df_smiles.sort_values(by='sumInts', ascending=False).reset_index(drop=True)
    #print(df_smiles)
    # %%
    # Saida final
    df_smiles.to_csv(clean_table_file.as_posix().replace(".csv", "_"+filename_suffix+".csv"), index=False, float_format="%.4f")
    # Saida final em Excel
    df_smiles.to_excel(clean_table_file.as_posix().replace(".csv", "_"+filename_suffix+".xlsx"), index=False, float_format="%.4f")
    print("  * Done! *")

def build_ssmn_protonated_cytoscape(ssmn_file, ssmn_protonated_file, clean_table_file, cy_style_file,
                                    number_bioactivities_scores):
    if not clean_table_file.exists() or not clean_table_file.is_file():
        sys.exit("The provided path to the clean table does not exists: " +clean_table_file.as_posix() +
                 ". Provide the correct path and try again.")
    if not ssmn_file.exists() or not ssmn_file.is_file():
        sys.exit("The provided path to the filtered SSMN file does not exists: " +ssmn_file.as_posix() +
                 ". Provide the correct path and try again.")
    if not ssmn_protonated_file.exists() or not ssmn_protonated_file.is_file():
        sys.exit("The provided path to the protonated SSMN file does not exists: " +ssmn_protonated_file.as_posix() +
                 ". Provide the correct path and try again.")
    if (not cy_style_file.is_file() and cy_style_file.as_posix() != ".") or not cy_style_file.exists():
        print("\nWARNING: The provided path to the Cytoscape visual style file does not exists. The networks will be created without style.\n")
        cy_style_file = Path("")
    
    clean_table = pd.read_csv(clean_table_file, escapechar="\\", skiprows=number_bioactivities_scores,
                              converters={'msclusterID':str})
    if clean_table.columns.str.startswith("COR_").any():
        # calcula o ranking das biocorrelacoes existentes
        clean_table = apply_corr_ranking(clean_table)
    
    print("* Loading SSMN and SSMN [M+H]+ to Cytoscape *")
    
    # Ensure Cytoscape is running and connected
    print("  - Check if Cytoscape is running...")
    try:
        py4.cytoscape_ping()
    except requests.exceptions.RequestException as e:
        sys.exit("  - Cytoscape is not running... =( Open Cytoscape application and try again. Error:",e)
    
    print("  - Import SSMN filtered and SSMN [M+H]+ filtered and check isolated nodes (selfloops)...")
    # import the informed network files
    ssmn_name = "SSMN filtered"
    ssmn_cy = py4.import_network_from_tabular_file(ssmn_file, first_row_as_column_names=True,
                                                   start_load_row=1,
                                                   column_type_list="s,t,ea,ea,ea,ea,ea,ea", delimiters="\\,")
    py4.rename_network(ssmn_name, ssmn_cy['networks'][0])
    ssmn = pd.read_csv(ssmn_file,converters={'msclusterID_source':str,'msclusterID_target':str})
    ssmn_isolated_nodes = ssmn.loc[(ssmn.msclusterID_source == ssmn.msclusterID_target),"msclusterID_source"].values
    clean_table["ssmn_filtered_isolated"] = False
    clean_table.loc[clean_table.msclusterID.isin(ssmn_isolated_nodes), "ssmn_filtered_isolated"] = True
    del ssmn
    ssmn_protonated_name = "SSMN [M+H]+ filtered"
    ssmn_protonated_cy = py4.import_network_from_tabular_file(ssmn_protonated_file, first_row_as_column_names=True,
                                                              start_load_row=1,
                                                              column_type_list="s,t,ea,ea,ea,ea,ea,ea",
                                                              delimiters="\\,")
    py4.rename_network(ssmn_protonated_name, ssmn_protonated_cy['networks'][0])
    ssmn_protonated = pd.read_csv(ssmn_protonated_file,converters={'msclusterID_source':str,'msclusterID_target':str})
    ssmn_protonated_isolated_nodes = ssmn_protonated.loc[ssmn_protonated.msclusterID_source == ssmn_protonated.msclusterID_target, "msclusterID_source"].values
    del ssmn_protonated
    # filter the protonated nodes only
    protonated_nodes = py4.get_all_nodes(ssmn_protonated_name)
    clean_table["protonated_selected"] = False
    clean_table.loc[clean_table.msclusterID.isin(protonated_nodes), "protonated_selected"] = True
    clean_table_protonated = clean_table.loc[clean_table.protonated_selected, :]
    clean_table_protonated = clean_table_protonated.reset_index(drop=True)
    clean_table_protonated["ssmn_protonated_filtered_isolated"] = False
    clean_table_protonated.loc[
        clean_table_protonated.msclusterID.isin(ssmn_protonated_isolated_nodes), "ssmn_protonated_filtered_isolated"] = True
    del ssmn_isolated_nodes, ssmn_protonated_isolated_nodes
    
    print("  - Load the clean table data to Cytoscape node table...")
    # load the clean table with the new data to the networks
    py4.load_table_data(clean_table_protonated, data_key_column='msclusterID', table='node', table_key_column='name', network=ssmn_protonated_cy['networks'][0])
    py4.load_table_data(clean_table, data_key_column='msclusterID', table='node', table_key_column='name', network=ssmn_cy['networks'][0])
    
    print("  - Save the clean table for the SSMN and SSMN protonated with the ranked results: ")
    # Salve o novo DataFrame em um arquivo CSV
    clean_table_ssmn_file = clean_table_file.as_posix().replace(".csv", "_ssmn.csv")
    clean_table_ssmn_protonated_file = clean_table_file.as_posix().replace(".csv", "_ssmn_protonated.csv")
    clean_table.to_csv(clean_table_ssmn_file, index=False, float_format="%.4f")
    clean_table_protonated.to_csv(clean_table_ssmn_protonated_file, index=False,float_format="%.4f")
    print("    - ",clean_table_ssmn_file, "\n    - ",clean_table_ssmn_protonated_file)
    
    # Import the visual style from an XML file
    if cy_style_file.as_posix() != ".":
        print("  - Add visual style to the cytoscape session")
        loaded_styles = py4.import_visual_styles(cy_style_file.as_posix())
        print("    - Loaded styles:", loaded_styles)
        # Apply the newly loaded style to your ssmn networks
        py4.set_visual_style(loaded_styles[0], network=ssmn_name)
        py4.set_visual_style(loaded_styles[0],network=ssmn_protonated_name)
        
    print("  * Done! *")
    # criar e salvar tabela para os colaboradores
    create_collaborators_tables(clean_table_protonated, clean_table_file, filename_suffix="protonated_collaborators")

if __name__ == "__main__":
    # Parser for arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--ssmn_file",
                        help="Path to the filtered SSMN .selfloop file. "
                             "In the molecular_networking folder from the NP^3 output folder.",
                        type=str, required=True)
    parser.add_argument("-p", "--ssmn_protonated_file",
                        help="Path to the [M+H]+ filtered SSMN .selfloop file. "
                             "In the molecular_networking folder from the NP^3 output folder.",
                        type=str, required=True)
    parser.add_argument("-c", "--clean_table_file",
                        help="Path to the clean table with the final list of consensus spectra - "
                             "the nodes information of the networks (.csv). A table with the correlation scores may "
                             "be used here, and properly set the number_bioactivities_scores if needed.",
                        type=str, required=True)
    parser.add_argument("-y", "--cytoscape_style_file",
                        help="Path to the cytoscape visual style file in XML to be imported in the networks view.",
                        type=str, default="")
    parser.add_argument("-b", "--number_bioactivities_scores",
                        help="The number of bioactivities scores present in the clean_table_file, default to 0. "
                             "This value will be used to skip the first rows from the clean_table_file "
                             "in case it contains correlation scores.",
                        type=int, default=0)

    # Parse arguments
    args = parser.parse_args()

    ssmn_protonated_file = Path(args.ssmn_protonated_file)
    ssmn_file = Path(args.ssmn_file)
    clean_table_file = Path(args.clean_table_file)
    number_bioactivities_scores = args.number_bioactivities_scores
    cy_style_file = Path(args.cytoscape_style_file)

    build_ssmn_protonated_cytoscape(ssmn_file, ssmn_protonated_file, clean_table_file, cy_style_file,
                                    number_bioactivities_scores)