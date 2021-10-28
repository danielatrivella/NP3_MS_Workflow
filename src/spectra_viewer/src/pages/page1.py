import streamlit as st
from ..utils import Page


class Page1(Page):
    def __init__(self, state):
        self.state = state

    def write(self):
        st.title("Work in progress")

        #slider_value = st.slider(
        #    "Set Value from here See it on Page 2",
        #    value=self.state.client_config["slider_value"],
        #)
        #self.state.client_config["slider_value"] = slider_value
        
        st.write("This webpage is a work in progress, wait for futher updates")
        st.balloons()