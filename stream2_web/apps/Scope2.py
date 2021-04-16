# -*- coding: utf-8 -*-
import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

import numpy as np
import stream as st

import matplotlib

matplotlib.use('Agg')

from app import app

adata = st.read(file_name='./SampleData/SCoPE2_2020/stream_result_var.pkl', workdir='./stream_result')


available_samples = ['Nestorowa, S. et al. 2016', 'Harrison, S. et al. 2021', 'Trapnell, C. et al. 2014',
                     ' Tang, Q. et al. 2017']
available_projections = ['dimension_reduction', 'visualization_2D', 'flat_tree', 'branches']
available_colors = adata.obs.columns
available_stream = ['single_cell_stream', 'stream']

layout = html.Div([
    dbc.Container([

        dbc.Col(html.H1("Visualization", className="text-center"), className="mb-5 mt-5"),
        dbc.Col(html.H5("Click the following tabs for step-by-step visualization.", className="text-center"),
                className="mb-5 mt-5"),
        dbc.Col(html.Hr(), className="mb-5 mt-5"),

        dbc.Card([
            dbc.CardHeader(
                dbc.Tabs(
                    [
                        dbc.Tab(label='1. Data Selection', tab_id="tab-ds", label_style={"font-size": "20px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='2. Quality Control', tab_id="tab-qc", label_style={"font-size": "20px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='3. Dimension Reduction', tab_id="tab-dr", label_style={"font-size": "20px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='4. Stream Plot', tab_id="tab-sp", label_style={"font-size": "20px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"})
                    ],
                    id="card-tabs", card=True, active_tab='tab-ds')
            ),
            dbc.CardHeader(id="card-header"),
            dbc.CardBody(id="card-content", style={"height": "60rem"}),
            dbc.CardFooter("This is the footer")
        ], color="dark", outline=True),

    ])
])

card_ds = dbc.FormGroup([
    dbc.Row([dbc.Col(
        dbc.Alert(
            "You can simply choose from our sample datasets or upload your own STREAM object (currently only *.pkl and *.h5ad formats are supported)",
            color="success"),
        className="mb-5 mt-5")]),

    dcc.RadioItems(
        id='Data',
        options=[
            {'label': 'Sample Datasets', 'value': 'sample'},
            {'label': 'Personal Datasets', 'value': 'personal'}
        ],
        value='sample',
        labelStyle={'display': 'inline-block'}
    ),

    html.Div(id='data-selection')
])

card_qc = dbc.FormGroup([
    dbc.Row([dbc.Col(
        dbc.Alert("You can explore QC metrics here!", color="success"),
        className="mb-5 mt-5")])
])

card_dr = dbc.FormGroup([
    dbc.Row([dbc.Col(
        dbc.Alert(
            "You can modify detailed dimension reduction figures here!",
            color="success"),
        className="mb-5 mt-5")]),

    dcc.Markdown("""
 **Color Cells by**
 """),
    dcc.Dropdown(
        id='color_cells_by',
        options=[{'label': i, 'value': i} for i in available_colors],
        value='n_counts'
    ),
    dcc.Markdown("""
 **Alpha**
 """),
    dcc.Slider(
        id='Alpha',
        min=0,
        max=1,
        value=0.8,
        marks={'%2.1f' % a: '%2.1f' % a for a in np.linspace(0, 1, 11)},
        step=None
    )
])

card_sp = dbc.FormGroup([
    dbc.Row([dbc.Col(
        dbc.Alert("You can explore STREAM plot here!", color="success"),
        className="mb-5 mt-5")]),

    dcc.Markdown("""
        **Color Cells by**
    """),
    dcc.Dropdown(
        id='color_cells_by_st',
        options=[{'label': i, 'value': i} for i in available_colors],
        value='n_counts'
    ),
    dcc.Markdown("""
**Alpha**
"""),
    dcc.Slider(
        id='Alpha_sp',
        min=0,
        max=1,
        value=0.8,
        marks={'%2.1f' % a: '%2.1f' % a for a in np.linspace(0, 1, 11)},
        step=None
    )
])


### tab activate
@app.callback(
    Output("card-content", "children"), Output("card-header", "children"), [Input("card-tabs", "active_tab")]
)
def tab_content(active_tab):
    if active_tab == "tab-qc":

        return html.Div(id= 'qc-graphic'), card_qc
    elif active_tab == "tab-dr":
        return html.Div(id="dr-graphic"), card_dr
    elif active_tab == "tab-sp":
        return html.Div(id="st-graphic"), card_sp
    else:
        fig = st.plot_stream(adata, root='S1', return_svg=True)
        return html.Img(src="data:image/svg+xml;base64,{}".format(fig)), card_ds


### Data Selection
@app.callback(
    Output('data-selection', 'children'),
    Input('Data', 'value'))
def update_dr(Data):
    if Data == "sample":
        return dcc.Dropdown(id='ds_sample', options=[{'label': i, 'value': i} for i in available_samples],
                            value='Harrison, S. et al. 2021')
    else:
        return dcc.Upload([
            'Drag and Drop or ',
            html.A('Select a File')
        ], style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center'
        })


### Demension Reduction
@app.callback(
    Output('dr-graphic', 'children'),
    Input('color_cells_by', 'value'),
    Input('Alpha', 'value'))
def update_dr(colors, alpha):
    fig1 = st.plot_dimension_reduction(adata, color=[colors], n_components=3, alpha=alpha, show_graph=True,
                                       show_text=False,
                                       plotly=True, return_fig=True)
    fig2 = st.plot_visualization_2D(adata, method='umap', n_neighbors=50, alpha=alpha,
                                    color=[colors], use_precomputed=True, plotly=True, return_fig=True)
    fig3 = st.plot_flat_tree(adata, color=[colors], alpha=alpha,
                             dist_scale=0.5, show_graph=True, show_text=True, plotly=True, return_fig=True)
    fig4 = st.plot_branches(adata,
                            show_text=True, plotly=True, return_fig=True)
    return html.Div([
        dbc.Row([dcc.Graph(figure=fig1), dcc.Graph(figure=fig2)]),
        dbc.Row([dcc.Graph(figure=fig3), dcc.Graph(figure=fig4)])
    ])


###--- Update Stream Plot ---###
@app.callback(
    Output('st-graphic', 'children'),
    Input('color_cells_by_st', 'value'),
    Input('Alpha_sp', 'value'))
def update_dr(colors, alpha):
    fig5 = st.plot_stream_sc(adata, root='S1', color=[colors], alpha=alpha,
                             dist_scale=0.3, show_graph=True, show_text=True, plotly=True, return_fig=True)
    fig6 = st.plot_stream(adata, root='S1', color=[colors], return_svg=True)

    return html.Div([
        dbc.Row([dbc.Col(dcc.Graph(figure=fig5), align="center"),
                 dbc.Col(html.Img(src="data:image/svg+xml;base64,{}".format(fig6)), align="center")])
    ])
