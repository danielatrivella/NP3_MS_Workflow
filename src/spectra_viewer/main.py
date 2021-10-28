# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:59:03 2021

@author: np3
"""

import streamlit as st

from src.sidebar import single_spectra_choices, compare_two_choices

# from src.pages import PAGE_MAP
# from src.state import provide_state
# from src.data import initial_positive
# from src.utils import frags_to_text


#@provide_state()
#def main(state=None):
def main():
    st.sidebar.title("NP³ Spectra Viewer")
    
    #Check box to compare spectra
    if st.sidebar.checkbox("Compare two spectra", value=False):
        compare_two_choices()
    else:
        single_spectra_choices()

if __name__ == "__main__":
    main()
