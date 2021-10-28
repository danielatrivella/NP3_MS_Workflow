# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 13:31:34 2021

@author: np3
"""

import streamlit as st
from src.pages import PAGE_MAP, PAGE_COMP
from src.state import provide_state
from src.data import initial_positive, initial_negative, rand_pep, rand_pepB
from src.utils import frags_to_text, choice_bar_at_side, text_to_frags, text_to_frags_spec, float_or_none


@provide_state()
def single_spectra_choices(state=None):
    current_page = st.sidebar.selectbox("Choose your input",
                                        list(PAGE_MAP),
                                        index=1)
    
    if current_page == "MGF File":
        mgf_file_pos = st.sidebar.file_uploader("Upload your file",
                                                type=["mgf"],
                                                accept_multiple_files=False)
        state.client_config["mgf_file_pos"] = mgf_file_pos
                                          
    elif current_page == "Peak List":
        pepmass_pos = st.sidebar.text_input('Parent Mass',
                                            value=rand_pep,
                                            max_chars=10)
        pks_pos = st.sidebar.text_area(label="Peak List",
                                       value=frags_to_text(initial_positive))
        
        state.client_config["pepmass_pos"] = float_or_none(pepmass_pos)
        state.client_config["peak_list_pos"] = pks_pos

    elif current_page == "UNPD-ISDB":
        st.sidebar.write("UNPD")
    PAGE_MAP[current_page](state=state).write()
    

@provide_state()
def compare_two_choices(state=None):
    
    st.sidebar.markdown("## Spectra A")

    # TODO Uncomment when UNPD-ISDB is implemented
    # inp_list = ["MGF File", "Peak List", "UNPD-ISDB"]
    inp_list = ["MGF File", "Peak List"]
    
    input_pos = st.sidebar.selectbox("Choose your positive input", 
                                     inp_list,
                                     index=1)
    
    if input_pos == "MGF File":
        mgf_file_pos = st.sidebar.file_uploader("Upload your file\nfor positive plot",
                                                type=["mgf"],
                                                accept_multiple_files=False)
        
        if mgf_file_pos is None:
            state.client_config["spec_ms_pos"] = text_to_frags_spec("0 0")
        else:
            # TODO abort try catch
            try:
                state.client_config["spec_ms_pos"] = choice_bar_at_side(mgf_file_pos, key="mgf_pos")
            except Exception:
                st.warning("MGF A without titles")
                state.client_config["spec_ms_pos"] = text_to_frags_spec("0 0")
                                          
    elif input_pos == "Peak List":
        pepmass_pos = st.sidebar.text_input('Parent Mass',
                                            value=rand_pep,
                                            max_chars=10,
                                            key="pep_pos")

        pks_pos = st.sidebar.text_area(label="Peak List",
                                       value=frags_to_text(initial_positive),
                                       key="pl_pos")

        state.client_config["spec_ms_pos"] = text_to_frags_spec(pks_pos,
                                                                id_spec='customList',
                                                                title='UpperSpectra',
                                                                pepmass=float_or_none(pepmass_pos))

    elif input_pos == "UNPD-ISDB":
        st.sidebar.write("UNPD")
    
    # TODO colocar um separador melhor, talvez inserir cores
    st.sidebar.markdown("## Spectra B")
    
    input_neg = st.sidebar.selectbox("Choose your negative input", 
                                     inp_list,
                                     index=1)
    if input_neg == "MGF File":
        mgf_file_neg = st.sidebar.file_uploader("Upload your file\nfor negative plot",
                                                type=["mgf"],
                                                accept_multiple_files=False)
        
        if mgf_file_neg is None:
            state.client_config["spec_ms_neg"] = text_to_frags_spec("0 0", title="spectrumB")
        else:
            # TODO abort this try catch
            try:
                state.client_config["spec_ms_neg"] = choice_bar_at_side(mgf_file_neg, key="mgf_neg")
            except Exception:
                st.warning("MGF B without titles")
                state.client_config["spec_ms_neg"] = text_to_frags_spec("0 0", title="spectrumB")

    elif input_neg == "Peak List":
        pepmass_neg = st.sidebar.text_input('Parent Mass',
                                            value=rand_pepB,
                                            max_chars=10,
                                            key="pep_neg")

        pks_neg = st.sidebar.text_area(label="Peak List",
                                       value=frags_to_text(initial_negative),
                                       key="pl_neg")
        
        state.client_config["spec_ms_neg"] = text_to_frags_spec(pks_neg,
                                                                id_spec='customList',
                                                                title='LowerSpectra',
                                                                pepmass=float_or_none(pepmass_neg))

    elif input_neg == "UNPD-ISDB":
        st.sidebar.write("UNPD")

    PAGE_COMP["Compare"](state=state).write()
