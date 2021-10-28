import streamlit as st
from pyteomics import mgf

from ..utils import Page, get_all_titles, pyteomics_2_matchms, show_meta_details
from ..charts import complete_simple_chart_sp


class mgf_page(Page):
    def __init__(self, state):
        self.state = state

    def write(self):
        st.title("MGF View")
        file = self.state.client_config["mgf_file_pos"]

        if file is not None:
            st.write("Select your spectrum")
            spectra = mgf.IndexedMGF(file)
            # List all spectra to choose
            all_title = get_all_titles(spectra)
            if len(all_title) != 0:
                idx_sel = st.selectbox("Select your spectrum", all_title)

                # Plot Title
                plt_title = st.text_input(label="Rename the plot below:",
                                          value="Peak list view of " + idx_sel,
                                          max_chars=80)

                # Plot Spectra
                sp_to_show = pyteomics_2_matchms(spectra, idx_sel)
                complete_simple_chart_sp(sp_to_show,
                                         title=plt_title)
                show_meta_details(sp_to_show)
            else:
                st.warning("Your mgf file does not have the TITLE parameter")

        else:
            st.info("Please upload a mgf file at the sidebar")
