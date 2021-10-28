import streamlit as st
from ..utils import Page
from ..utils import text_to_frags, text_to_frags_spec
from ..charts import complete_simple_chart_sp


class peak_list(Page):
    def __init__(self, state):
        self.state = state

    def write(self):
        txt = self.state.client_config["peak_list_pos"]
        pep_mass = self.state.client_config["pepmass_pos"]
        st.title("Peak List plot")
        plt_title = st.text_input(label="Rename the plot below:",
                                  value="Peak list view")
        
        spec = text_to_frags_spec(txt,
                                  id_spec='customList',
                                  pepmass=pep_mass)
        complete_simple_chart_sp(spec, plt_title)