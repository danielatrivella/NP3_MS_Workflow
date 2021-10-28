from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Span, Range1d
import streamlit as st
import numpy as np

from matchms.similarity.spectrum_similarity_functions import collect_peak_pairs
from matchms.filtering import normalize_intensities, default_filters, reduce_to_number_of_peaks, \
    select_by_relative_intensity, select_by_mz
from src.utils import dot_product_sp, dot_product_shift_sp, log_sp, sqrt_sp, try_get_pep, sp_2_np, remove_not_shifted_common_pks


def simple_barchart_sp(spectra,
                       title,
                       x_grid=True,
                       y_grid=True,
                       w_bar=1.5,
                       out_format='svg',
                       color="#4c78a8",
                       title_pt=13,
                       tick_pt=10,
                       tick_lab_pt=12):
    """
    Barchart bokeh composer for peak lists from spectra
    """
    # Get data
    source = ColumnDataSource({'mz': spectra.peaks.mz,
                               'int': spectra.peaks.intensities})
    hoover_tool_tips = [("File", title),
                        ("m/z", "@mz"),
                        ("int", "@int")]
    p = figure(
        title=title,
        plot_height=300,
        x_axis_label="m/z",
        y_axis_label="Intensity",
        tools="xpan,ypan,xwheel_zoom,ywheel_zoom,box_zoom,reset,hover,save",
        tooltips=hoover_tool_tips,
        toolbar_location="above")

    # To save as svg, output_backend must be definided as svg
    # To save as png, output_backend must be definided as webgl (webgl also accelerate interface)
    if out_format == 'svg':
        p.output_backend = 'svg'
    else:
        p.output_backend = 'webgl'

    p.vbar(x="mz", top="int", width=w_bar, source=source, color=color)
    p.xgrid.visible = y_grid
    p.ygrid.visible = x_grid
    p.title.text_font_size = f"{title_pt}pt"
    p.xaxis.major_label_text_font_size = f"{tick_pt}pt"
    p.yaxis.major_label_text_font_size = f"{tick_pt}pt"
    p.xaxis.axis_label_text_font_size = f"{tick_lab_pt}pt"
    p.yaxis.axis_label_text_font_size = f"{tick_lab_pt}pt"

    return p

def complete_simple_chart_sp(spectra, title):
    final_sp = spectra.clone()
    color_grid = st.expander(label="Color, grid and output type", expanded=False)
    sizes = st.expander(label="Labels sizes and bar width", expanded=False)
    scales_boxes = st.expander(label="Scales", expanded=False)
    filters = st.expander(label="Filters", expanded=False)

    with color_grid:
        col_color, col_x, col_y, out_type = st.columns([1, 1, 1, 1])
    with col_color:
        # st.write("Pick a color")
        color = st.color_picker("Pick a color", "#4c78a8")
    with col_x:
        st.write("H grid", )
        grid_x = st.checkbox("", value=True, key="x")
    with col_y:
        st.write("V grid")
        grid_y = st.checkbox("", value=False, key="y")
    with out_type:
        out_t = st.radio("Output type",
                         ["svg", "png"],
                         key='out_radio')

    with sizes:
        large_bar, title_size, labels_size, ticks_size = st.columns([1, 1, 1, 1])
    with large_bar:
        # st.write("Bar width")
        w_bar = st.number_input("Bar width",
                                min_value=0.0,
                                max_value=50.0,
                                step=0.1,
                                value=1.5,
                                format="%.1f",
                                key="bar_len")
    with title_size:
        t_size = st.number_input("Title size",
                                 min_value=1,
                                 max_value=50,
                                 step=1,
                                 value=13,
                                 key="t_size")
    with labels_size:
        l_size = st.number_input("Labels size",
                                 min_value=1,
                                 max_value=50,
                                 step=1,
                                 value=12,
                                 key="l_size")
    with ticks_size:
        tk_size = st.number_input("Axis labels size",
                                  min_value=1,
                                  max_value=50,
                                  step=1,
                                  value=10,
                                  key="tk_size")

    with scales_boxes:
        cont01 = st.container()
        with cont01:
            # st.write("Root square or Natural log at intensities")
            do_witch = st.selectbox("Root square or Natural log at intensities",
                                    ["None", "sqrt", "log"],
                                    index=0)
        do_norm = st.checkbox("Normalize intensities of peaks to unit height.", value=False, key="norm")

    if do_witch == "None":
        pass
    elif do_witch == "sqrt":
        final_sp = sqrt_sp(final_sp)
    elif do_witch == "log":
        final_sp = log_sp(final_sp)

    if do_norm:
        final_sp = normalize_intensities(final_sp)

    with filters:
        sel_mz = st.checkbox("Select by m/z", value=False, key="sel_mz")
        if sel_mz:
            st.write("Keep only peaks between mz_from and mz_to (keep if mz_from >= m/z >= mz_to).\n")
            cont02 = st.container()
            with cont02:
                mz_from_col, mz_to_col = st.columns(2)
                with mz_from_col:
                    mz_from = st.number_input("Lower threshold",
                                              value=0.0,
                                              min_value=0.0,
                                              step=0.01,
                                              format="%.2f",
                                              key='mz_from')
                with mz_to_col:
                    mz_to = st.number_input("Upper threshold",
                                            value=10000.0,
                                            min_value=0.0,
                                            step=0.01,
                                            format="%.2f",
                                            key='mz_to')
            # Correction if user set greater value to from
            # swipe values
            mz_from = min(mz_from, mz_to)
            mz_to = max(mz_from, mz_to)

        sel_rel = st.checkbox("Select by relative intensity", value=False, key="sel_rel")
        if sel_rel:
            st.write(
                "Keep only peaks within set relative intensity range (keep if intensity_from >= intensity >= intensity_to).\n")
            cont03 = st.container()
            with cont03:
                int_from_col, int_to_col = st.columns(2)
                with int_from_col:
                    int_from = st.number_input("Lower threshold",
                                               value=0.0,
                                               max_value=1.0,
                                               min_value=0.0,
                                               step=0.01,
                                               format="%.2f",
                                               key='int_from')
                with int_to_col:
                    int_to = st.number_input("Upper threshold",
                                             value=1.0,
                                             max_value=1.0,
                                             min_value=0.0,
                                             step=0.01,
                                             format="%.2f",
                                             key='int_to')
            # Correction if user set greater value to from
            # swipe values
            int_from = min(int_from, int_to)
            int_to = max(int_from, int_to)

        max_pks = st.checkbox("Reduce to n peaks", value=False, key="max_pks")
        if max_pks:
            # st.write("Lowest intensity peaks will be removed when it has more peaks than desired.")
            max_pks_val = st.number_input("Lowest intensity peaks will be removed when it has more peaks than desired.",
                                          min_value=1,
                                          value=10,
                                          key='max_pks_val')
    # Aplling filters
    if sel_mz:
        final_sp = select_by_mz(final_sp,
                                mz_from=mz_from,
                                mz_to=mz_to)
    if sel_rel:
        final_sp = select_by_relative_intensity(final_sp,
                                                intensity_from=int_from,
                                                intensity_to=int_to)
    if max_pks:
        final_sp = reduce_to_number_of_peaks(final_sp,
                                             n_max=max_pks_val)

    plot = simple_barchart_sp(final_sp,
                              title=title,
                              x_grid=grid_x,
                              y_grid=grid_y,
                              tick_pt=tk_size,
                              tick_lab_pt=l_size,
                              title_pt=t_size,
                              w_bar=w_bar,
                              color=color,
                              out_format=out_t)
    st.bokeh_chart(plot,
                   use_container_width=True)


def double_chart_sp(spP,
                    spN,
                    title,
                    colorP="#4c78a8",
                    colorN="#d62e2e",
                    cos_comp=None,
                    x_grid=True,
                    y_grid=True,
                    w_bar=1.5,
                    out_format='svg',
                    title_pt=13,
                    tick_pt=10,
                    tick_lab_pt=12):
    sourceP = ColumnDataSource({'mz': spP.peaks.mz,
                                'int': spP.peaks.intensities})

    sourceN = ColumnDataSource({'mz': spN.peaks.mz,
                                'int': -1 * spN.peaks.intensities})

    hoover_tool_tips = [("Plot", title),
                        ("m/z", "@mz"),
                        ("int", "@int")]
    if cos_comp != None:
        title += " | cosine: " + str(cos_comp)

    p = figure(
        title=title,
        plot_height=400,
        x_axis_label="m/z",
        y_axis_label="Intensity",
        tools="xpan,ypan,xwheel_zoom,ywheel_zoom,box_zoom,reset,hover,save",
        tooltips=hoover_tool_tips,
        toolbar_location="above")

    # To save as svg, output_backend must be definided as svg
    # To save as png, output_backend must be definided as webgl (webgl also accelerate interface)
    if out_format == 'svg':
        p.output_backend = 'svg'
    else:
        p.output_backend = 'webgl'

    p.vbar(x="mz", top="int", width=w_bar, source=sourceP, color=colorP)
    p.vbar(x="mz", top="int", width=w_bar, source=sourceN, color=colorN)
    p.xgrid.visible = y_grid
    p.ygrid.visible = x_grid
    p.title.text_font_size = f"{title_pt}pt"
    p.xaxis.major_label_text_font_size = f"{tick_pt}pt"
    p.yaxis.major_label_text_font_size = f"{tick_pt}pt"
    p.xaxis.axis_label_text_font_size = f"{tick_lab_pt}pt"
    p.yaxis.axis_label_text_font_size = f"{tick_lab_pt}pt"
    return p


def double_chart_shared_sp(spP,
                           spN,
                           sha_array,
                           title,
                           colorP="#4c78a8",
                           colorN="#d62e2e",
                           colorP_sha="#152e4b",
                           colorN_sha="#8b1b1b",
                           cos_comp=None,
                           x_grid=True,
                           y_grid=True,
                           w_bar=1.5,
                           out_format='svg',
                           title_pt=13,
                           tick_pt=10,
                           tick_lab_pt=12):
    sourceP = ColumnDataSource({'mz': spP.peaks.mz,
                                'int': spP.peaks.intensities})

    sourceN = ColumnDataSource({'mz': spN.peaks.mz,
                                'int': -1 * spN.peaks.intensities})

    sourceP_sha = ColumnDataSource({'mz': np.take(spP.peaks.mz, list(sha_array[:, 0])),
                                    'int': np.take(spP.peaks.intensities, list(sha_array[:, 0])), })

    sourceN_sha = ColumnDataSource({'mz': np.take(spN.peaks.mz, list(sha_array[:, 1])),
                                    'int': -1 * np.take(spN.peaks.intensities, list(sha_array[:, 1])), })

    hoover_tool_tips = [("Plot", title),
                        ("m/z", "@mz"),
                        ("int", "@int")]
    if cos_comp != None:
        title += " | cosine: " + str(cos_comp)

    p = figure(
        title=title,
        plot_height=400,
        x_axis_label="m/z",
        y_axis_label="Intensity",
        tools="xpan,ypan,xwheel_zoom,ywheel_zoom,box_zoom,reset,hover,save",
        tooltips=hoover_tool_tips,
        toolbar_location="above")

    # To save as svg, output_backend must be definided as svg
    # To save as png, output_backend must be definided as webgl (webgl also accelerate interface)
    if out_format == 'svg':
        p.output_backend = 'svg'
    else:
        p.output_backend = 'webgl'

    p.vbar(x="mz", top="int", width=w_bar, source=sourceP, color=colorP)
    p.vbar(x="mz", top="int", width=w_bar, source=sourceN, color=colorN)
    p.vbar(x="mz", top="int", width=w_bar, source=sourceP_sha, color=colorP_sha)
    p.vbar(x="mz", top="int", width=w_bar, source=sourceN_sha, color=colorN_sha)
    p.xgrid.visible = y_grid
    p.ygrid.visible = x_grid
    p.title.text_font_size = f"{title_pt}pt"
    p.xaxis.major_label_text_font_size = f"{tick_pt}pt"
    p.yaxis.major_label_text_font_size = f"{tick_pt}pt"
    p.xaxis.axis_label_text_font_size = f"{tick_lab_pt}pt"
    p.yaxis.axis_label_text_font_size = f"{tick_lab_pt}pt"
    return p


def complete_double_chart_sp(spP, spN, title):
    final_spP = spP.clone()
    final_spN = spN.clone()

    color_grid = st.expander(label="Color, grid and output type", expanded=False)
    sizes = st.expander(label="Bar width, labels and ticks sizes", expanded=False)
    scales_boxes = st.expander(label="Scales", expanded=False)
    filters = st.expander(label="Filters", expanded=False)
    similarity = st.expander(label="Similarity", expanded=False)

    with color_grid:
        col_color, col_x, col_y, out_type = st.columns([1, 1, 1, 1])
    with col_color:
        # st.write("Pick a color")
        color_p = st.color_picker("Upper spectra color", "#4c78a8")
        color_n = st.color_picker("Lower spectra color", "#d62e2e")
    with col_x:
        st.text("H grid", )
        grid_x = st.checkbox("", value=True, key="x")
    with col_y:
        st.text("V grid")
        grid_y = st.checkbox("", value=False, key="y")
    with out_type:
        out_t = st.radio("Output type",
                         ["svg", "png"],
                         key='out_radio')

    with sizes:
        large_bar, title_size, labels_size, ticks_size = st.columns([1, 1, 1, 1])
    with large_bar:
        # st.write("Bar width")
        w_bar = st.number_input("Bar width",
                                min_value=0.0,
                                max_value=50.0,
                                step=0.1,
                                value=1.5,
                                format="%.1f",
                                key="bar_len")
    with title_size:
        t_size = st.number_input("Title size",
                                 min_value=1,
                                 max_value=50,
                                 step=1,
                                 value=13,
                                 key="t_size")
    with labels_size:
        l_size = st.number_input("Axis title size",
                                 min_value=1,
                                 max_value=50,
                                 step=1,
                                 value=12,
                                 key="l_size")
    with ticks_size:
        tk_size = st.number_input("Ticks label size",
                                  min_value=1,
                                  max_value=50,
                                  step=1,
                                  value=10,
                                  key="tk_size")

    with scales_boxes:
        cont01 = st.container()
        with cont01:
            # st.write("Root square or Natural log at intensities")
            do_witch = st.selectbox("Root square or Natural log at intensities",
                                    ["None", "sqrt", "log"],
                                    index=0)
        do_norm = st.checkbox("Normalize intensities of peaks to unit height.", value=False, key="norm")

    if do_witch == "None":
        pass
    elif do_witch == "sqrt":
        final_spP = sqrt_sp(final_spP)
        final_spN = sqrt_sp(final_spN)
    elif do_witch == "log":
        final_spP = log_sp(final_spP)
        final_spN = log_sp(final_spN)

    if do_norm:
        final_spP = normalize_intensities(final_spP)
        final_spN = normalize_intensities(final_spN)

    with filters:
        sel_mz = st.checkbox("Select by m/z", value=False, key="sel_mz")
        if sel_mz:
            st.write("Keep only peaks between mz_from and mz_to (keep if mz_from >= m/z >= mz_to).\n")
            cont02 = st.container()
            with cont02:
                mz_from_col, mz_to_col = st.columns(2)
                with mz_from_col:
                    mz_from = st.number_input("Lower threshold",
                                              value=0.0,
                                              min_value=0.0,
                                              step=0.01,
                                              format="%.2f",
                                              key='mz_from')
                with mz_to_col:
                    mz_to = st.number_input("Upper threshold",
                                            value=10000.0,
                                            min_value=0.0,
                                            step=0.01,
                                            format="%.2f",
                                            key='mz_to')
            # Correction if user set greater value to from
            # swipe values
            mz_from = min(mz_from, mz_to)
            mz_to = max(mz_from, mz_to)

        sel_rel = st.checkbox("Select by relative intensity", value=False, key="sel_rel")
        if sel_rel:
            st.write(
                "Keep only peaks within set relative intensity range (keep if intensity_from >= intensity >= "
                "intensity_to).\n")
            cont03 = st.container()
            with cont03:
                int_from_col, int_to_col = st.columns(2)
                with int_from_col:
                    int_from = st.number_input("Lower threshold",
                                               value=0.0,
                                               max_value=1.0,
                                               min_value=0.0,
                                               step=0.01,
                                               format="%.2f",
                                               key='int_from')
                with int_to_col:
                    int_to = st.number_input("Upper threshold",
                                             value=1.0,
                                             max_value=1.0,
                                             min_value=0.0,
                                             step=0.01,
                                             format="%.2f",
                                             key='int_to')
            # Correction if user set greater value to from
            # swipe values
            int_from = min(int_from, int_to)
            int_to = max(int_from, int_to)

        max_pks = st.checkbox("Reduce to n peaks", value=False, key="max_pks")
        if max_pks:
            # st.write("Lowest intensity peaks will be removed when it has more peaks than desired.")
            max_pks_val = st.number_input("Lowest intensity peaks will be removed when it has more peaks than desired.",
                                          min_value=1,
                                          value=10,
                                          key='max_pks_val')
    # Applying filters
    if sel_mz:
        final_spP = select_by_mz(final_spP,
                                 mz_from=mz_from,
                                 mz_to=mz_to)

        final_spN = select_by_mz(final_spN,
                                 mz_from=mz_from,
                                 mz_to=mz_to)
    if sel_rel:
        final_spP = select_by_relative_intensity(final_spP,
                                                 intensity_from=int_from,
                                                 intensity_to=int_to)

        final_spN = select_by_relative_intensity(final_spN,
                                                 intensity_from=int_from,
                                                 intensity_to=int_to)
    if max_pks:
        final_spP = reduce_to_number_of_peaks(final_spP,
                                              n_max=max_pks_val)

        final_spN = reduce_to_number_of_peaks(final_spN,
                                              n_max=max_pks_val)

    with similarity:
        # Get cosines
        pepP = try_get_pep(spP)
        pepN = try_get_pep(spN)

        st.info("To reproduce default NP³ MS Workflow cosine, make sure that **\'sqrt\'** scale is activated at the "
                "**\'Scales\'** menu")

        bin_size_numb = st.number_input("m/z tolerance",
                                        min_value=0.0,
                                        max_value=5.0,
                                        step=0.001,
                                        value=0.05,
                                        format="%.3f",
                                        key="bin_size_numb")

        # Regular cosine without pep mass
        # TODO Join peakids checkbox
        if pepP is None or pepN is None:
            cos_org = dot_product_sp(spP, spN,
                                     bin_size=bin_size_numb)
            cos_flt = dot_product_sp(final_spP, final_spN,
                                     bin_size=bin_size_numb)
            shared_raw = collect_peak_pairs(sp_2_np(spP),
                                            sp_2_np(spN),
                                            tolerance=bin_size_numb)

            shared_flt = collect_peak_pairs(sp_2_np(final_spP),
                                            sp_2_np(final_spN),
                                            tolerance=bin_size_numb)
        else:
            cos_org = dot_product_shift_sp(spP, spN,
                                           bin_size=bin_size_numb)
            cos_flt = dot_product_shift_sp(final_spP, final_spN,
                                           bin_size=bin_size_numb)

            #TODO Put it in a function

            # Get peaks shared without shift
            shared_raw_no_shift = collect_peak_pairs(sp_2_np(spP),
                                                     sp_2_np(spN),
                                                     tolerance=bin_size_numb)
            if shared_raw_no_shift is not None:
                # From original remove marked peaks on previous function
                remained_pks_raw_pos = remove_not_shifted_common_pks(sp_2_np(spP),
                                                                     shared_raw_no_shift,
                                                                     col=0)
                remained_pks_raw_neg = remove_not_shifted_common_pks(sp_2_np(spN),
                                                                     shared_raw_no_shift,
                                                                     col=1)
                # Get peaks shared with the shift (from remaining peaks)
                shared_raw_with_shift = collect_peak_pairs(remained_pks_raw_pos,
                                                           remained_pks_raw_neg,
                                                           tolerance=bin_size_numb,
                                                           shift=pepP - pepN)
                # Append both
                if shared_raw_with_shift is not None:
                    shared_raw = np.append(shared_raw_no_shift,
                                           shared_raw_with_shift,
                                           axis=0)
                else:
                    shared_raw = shared_raw_no_shift

            else:
                shared_raw = None

            # Get peaks shared without shift
            shared_flt_no_shift = collect_peak_pairs(sp_2_np(final_spP),
                                                     sp_2_np(final_spN),
                                                     tolerance=bin_size_numb)
            if shared_flt_no_shift is not None:
                # From original remove marked peaks on previous function
                remained_pks_flt_pos = remove_not_shifted_common_pks(sp_2_np(final_spP),
                                                                     shared_flt_no_shift,
                                                                     col=0)
                remained_pks_flt_neg = remove_not_shifted_common_pks(sp_2_np(final_spN),
                                                                     shared_flt_no_shift,
                                                                     col=1)

                # Get peaks shared with the shift (from remaining peaks)
                shared_flt_with_shift = collect_peak_pairs(remained_pks_flt_pos,
                                                           remained_pks_flt_neg,
                                                           tolerance=bin_size_numb,
                                                           shift=pepP - pepN)

                # Append both
                if shared_flt_with_shift is not None:
                    shared_flt = np.append(shared_flt_no_shift,
                                           shared_flt_with_shift,
                                           axis=0)
                else:
                    shared_flt = shared_flt_no_shift
            else:
                shared_flt = None

        col_cos, col_peaks, col_shared = st.columns(3)
        with col_cos:
            st.text(f'Cosine\nraw spectra:\n\t{round(cos_org, 4)}')
            st.text(f'Cosine\nfiltred spectra:\n\t{round(cos_flt, 4)}')
        with col_peaks:
            if shared_raw is None:
                st.text(f'Shared peaks\nraw spectra:\n\t0')
            else:
                st.text(f'Shared peaks\nraw spectra:\n\t{shared_raw.shape[0]}')

            if shared_flt is None:
                st.text(f'Shared peaks\nraw spectra:\n\t0')
            else:
                st.text(f'Shared peaks\nraw spectra:\n\t{shared_flt.shape[0]}')
        with col_shared:
            show_sha = st.checkbox("Show shared peaks", value=False)
            if show_sha:
                color_p_s = st.color_picker("Upper shared color", "#152e4b")
                color_n_s = st.color_picker("Lower shared color", "#8b1b1b")

    if show_sha and shared_flt is not None:
        plot = double_chart_shared_sp(final_spP,
                                      final_spN,
                                      shared_flt,
                                      colorP=color_p,
                                      colorN=color_n,
                                      colorP_sha=color_p_s,
                                      colorN_sha=color_n_s,
                                      x_grid=grid_x,
                                      y_grid=grid_y,
                                      title=title,
                                      w_bar=w_bar,
                                      out_format=out_t,
                                      title_pt=t_size,
                                      tick_pt=tk_size,
                                      tick_lab_pt=l_size)

    else:
        plot = double_chart_sp(final_spP,
                               final_spN,
                               colorP=color_p,
                               colorN=color_n,
                               x_grid=grid_x,
                               y_grid=grid_y,
                               title=title,
                               w_bar=w_bar,
                               out_format=out_t,
                               title_pt=t_size,
                               tick_pt=tk_size,
                               tick_lab_pt=l_size)

    st.bokeh_chart(plot,
                   use_container_width=True)
