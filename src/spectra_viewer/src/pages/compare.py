import streamlit as st
from ..utils import Page, show_meta_details
from ..charts import complete_double_chart_sp, double_chart_sp



class compare(Page):
    def __init__(self, state):
        self.state = state

    def write(self):
        st.title("Comparing Spectra")
        st.markdown("## Choose the spectra at the left sidebar")

        # Get spectra objects
        sp_pos = self.state.client_config["spec_ms_pos"]
        sp_neg = self.state.client_config["spec_ms_neg"]

        # Use Spectrum titles as plot title
        ttl_p = sp_pos.get('title')
        ttl_p_ttl = ttl_p[:]
        ttl_n = sp_neg.get('title')
        ttl_n_ttl = ttl_n[:]
        if ttl_p is None:
            ttl_p = 'A'
            ttl_p_ttl = 'A'
        else:
            if len(ttl_p) > 20:
                ttl_p_ttl = ttl_p_ttl[:3]+'...'+ttl_p_ttl[-10:]
        if ttl_n is None:
            ttl_n = 'B'
            ttl_n_ttl = 'B'
        else:
            if len(ttl_n) > 20:
                ttl_n_ttl = ttl_n_ttl[:3]+'...'+ttl_n_ttl[-10:]

        plt_title = st.text_input(label="Rename the plot below:",
                                  value=f"{ttl_p_ttl} vs {ttl_n_ttl}")

        
        complete_double_chart_sp(sp_pos,
                                 sp_neg,
                                 title=plt_title)
        
        cont_details = st.container()
        with cont_details:
            spA_col, spB_col = st.columns(2)

            # Needs key name to not conflict using two check boxes
            with spA_col:
                st.text(ttl_p+"\nParameters")
                show_meta_details(sp_pos,
                                  check_key='meta_key_p')
            with spB_col:
                st.text(ttl_n+"\nParameters")
                show_meta_details(sp_neg,
                                  check_key='meta_key_n')

