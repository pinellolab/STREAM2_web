# -*- coding: utf-8 -*-
### import web components related packages
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import time
import plotly.graph_objects as go

### import calculation related packages
import numpy as np
import stream as st
import matplotlib

matplotlib.use('Agg')

from app import app

### Load default sample data and do some basic modification to it
adata_computed = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016.pkl', workdir='./stream_result')
adata = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016-raw.h5ad', workdir='./stream_result')
adata.uns[
    'discription'] = 'This scRNA-seq dataset contains 1656 cells and 4768 genes from mouse hematopoietic stem and progenitor cell differentiation. A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation. Blood 128, e20-31 (2016).'

available_normalization = ['Library size correction', 'TF-IDF transformation', 'None']
available_DR = ['Spectral embedding(SE)', 'Modified locally linear embedding(MLLE)', 'UMAP', 'PCA']

### Set optionals
available_samples = ['Nestorowa, S. et al. 2016', 'Harrison, S. et al. 2021', 'Trapnell, C. et al. 2014',
                     ' Tang, Q. et al. 2017']

### Basic layout of this visulization page
layout = html.Div([
    dbc.Container([

        dbc.Col(html.H1("Computation", className="text-center"), className="mb-5 mt-5"),
        dbc.Col(html.H5("Click the following tabs for step-by-step computation.", className="text-center"),
                className="mb-5 mt-5"),

        dbc.Card([
            dbc.CardHeader(
                dbc.Tabs(
                    [
                        dbc.Tab(label='1. Data Selection', tab_id="tab-C-ds", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}, disabled=False),

                        dbc.Tab(label='2. Preprocessing', tab_id="tab-C-pre", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}, disabled=False),

                        dbc.Tab(label='3. Dimension Reduction', tab_id="tab-C-dr", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}, disabled=True),

                        dbc.Tab(label='4. Trajectory Inference', tab_id="tab-C-ti", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}, disabled=True),

                        dbc.Tab(label='5. Features Identification', tab_id="tab-C-fi",
                                label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}, disabled=True)
                    ],
                    id="card-C-tabs", active_tab='tab-C-ds', card=True)
            ),
            dbc.CardHeader(id="card-C-header"),
            dbc.CardBody(id="card-C-content", style={"height": "40rem"}, className="w-90"),
            dbc.CardFooter("If you like STREAM2 please support us by citing it in your work \N{TWO HEARTs}")
        ], color="dark", outline=True),

    ])
])


### tab activate
@app.callback(
    Output("card-C-tabs", "active_tab"),
    [Input("Next1-button", "n_clicks"),
     Input("Next2-button", "n_clicks")],
)
def next_page(n1, n2):
    if n2:
        return 'tab-C-pre'
    elif n1:
        return 'header_C_ds'
    else:
        return 'header_C_ds'


@app.callback(
    Output("card-C-header", "children"), Output("card-C-content", "children"), [Input("card-C-tabs", "active_tab")]
)
def tab_content(active_tab):
    if active_tab == "tab-C-pre":
        return header_C_pre, html.Div(id='card_C_pre')
    elif active_tab == "tab-C-dr":
        return header_C_dr, html.Div(id='card_C_dr')
    elif active_tab == "tab-C-ti":
        return header_C_ti, html.Div(id='card_C_ti')
    elif active_tab == "tab-C-fi":
        # return header_C_fi, card_C_fi
        return header_C_pre, html.Div(id='card_C_fi')
    else:
        return header_C_ds, html.Div(id='card_C_ds')


### Fill in contents for each header
# DS
header_C_ds = dbc.FormGroup([
    dbc.Alert(
        "You can simply choose from our sample datasets or upload your own files (10X Genomics output, Scanpy object and Seurat object are supported)",
        color="success"),

    dcc.RadioItems(
        id='data_C_source',
        options=[
            {'label': 'Sample Datasets', 'value': 'sample'},
            {'label': 'Personal Datasets', 'value': 'personal'}
        ],
        inputStyle={"margin-right": "20px", "margin-left": "20px"},
        value='sample',
        labelStyle={'display': 'inline-block'}
    ),

    html.Div(id='data-C-selection'),

    dbc.Row([
        dbc.Col(dbc.Button(
            dbc.Spinner(id="C1-spinners"),
            id="Compute1-button",
            className="mb-3",
            color="success",
            block=True)),

        dbc.Col(dbc.Button("Go To Next Step -- Preprocessing",
                           id="Next1-button",
                           className="mb-3",
                           color="success",
                           block=True,
                           disabled=True))
    ])
])


###### Data Selection
@app.callback(
    Output('data-C-selection', 'children'),
    Input('data_C_source', 'value'))
def update_ds(Source):
    if Source == "sample":
        return dcc.Dropdown(id='ds_C_sample', options=[{'label': i, 'value': i} for i in available_samples],
                            value='Nestorowa, S. et al. 2016')
    else:
        return dbc.FormGroup([
            dcc.RadioItems(
                id='ds_C_type',
                options=[
                    {'label': '10X Genomics output', 'value': '10x'},
                    {'label': '*.loom', 'value': 'loom'},
                    {'label': 'Scanpy *.h5ad', 'value': 'h5ad'},
                    {'label': 'Seurat *.rds', 'value': 'rds'}
                ],
                inputStyle={"margin-right": "20px", "margin-left": "20px"},
                value='10x',
                labelStyle={'display': 'inline-block'}
            ),

            dcc.Upload([
                'Drag and Drop or',
                html.A('Select Files')
            ], style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
            }, id='ds_C_personal')
        ])


@app.callback(
    [Output("card_C_ds", "children"),
     Output("C1-spinners", "children"),
     Output("Compute1-button", "disabled"),
     Output("Next1-button", "disabled")],
    Input("Compute1-button", "n_clicks")
)
def load_data(n):
    if n is None:
        return "Data is not ready yet. Click the CONFIRM button above and waite it to be uploaded.", "Confirm Data Selection", False, True
    else:
        time.sleep(1)
        return card_C_ds, "Data is successfully uploaded", True, False


# PRE
header_C_pre = dbc.FormGroup([
    dbc.Alert(
        "You can modify parameters for Step preprocessing here!",
        color="success"),

    dbc.Row([
        dbc.Col(
            [dcc.Markdown("""**Expresion Cutoff**"""),
             dcc.Input(
                 id='pre_expr_cutoff',
                 type='number',
                 value=1
             )]),

        dbc.Col(
            [dcc.Markdown("""**min number of features**"""),
             dcc.Input(
                 id='pre_n_features',
                 type='number',
                 value=100
             )]),

        dbc.Col(
            [dcc.Markdown("""**min number of cells**"""),
             dcc.Input(
                 id='pre_n_cells',
                 type='number',
                 value=5
             )]),

        dbc.Col(
            [dcc.Markdown("""**Log2 Transformation**"""),
             dcc.Dropdown(
                 id='pre_log2',
                 options=[{'label': i, 'value': i} for i in ['Yes', 'No']],
                 value='Yes')]),

        dbc.Col(
            [dcc.Markdown("""**Normalization methods**"""),
             dcc.Dropdown(
                 id='pre_nom',
                 options=[{'label': i, 'value': i} for i in available_normalization],
                 value=available_normalization[0])])
    ]),

    dbc.Row([
        dbc.Col(
            [dcc.Markdown("""**Feature Selection methods**"""),
             dcc.Dropdown(
                 id='pre_fs',
                 options=[{'label': i, 'value': i} for i in ['LOESS', 'PCA', 'All']],
                 value='LOESS')]),
        dbc.Col(html.Div(id='pre_av'))
    ]),

    dbc.Row([
        dbc.Col(dbc.Button(
            dbc.Spinner(id="C2-spinners"),
            id="Compute2-button",
            className="mb-3",
            color="success",
            block=True)),

        dbc.Col(dbc.Button("Go To Next Step -- Dimension Reduction",
                           id="Next2-button",
                           className="mb-3",
                           color="success",
                           block=True,
                           disabled=True))
    ])
])


@app.callback(
    [Output("card_C_pre", "children"),
     Output("C2-spinners", "children"),
     Output("Compute2-button", "disabled"),
     Output("Next2-button", "disabled")],
    Input("Compute2-button", "n_clicks")
)
def C2_parameters(n):
    if n is None:
        return "Parameters are not confirmed yet. Click the CONFIRM button above and waite Step 2 to be finished.", "Confirm Parameters", False, True
    else:
        time.sleep(1)
        return card_C_pre_finished, "Step 2 is successfully finished", True, False


C_LOESS = dbc.FormGroup([
    dbc.Row([
        dbc.Col([
            dcc.Markdown("""**Loess Fraction**"""),
            dcc.Input(id='pre_loess_frac',
                      type='number',
                      value=0.01)]),
        dbc.Col([
            dcc.Markdown("""**Percentile**"""),
            dcc.Input(id='pre_percentile',
                      type='number',
                      value=95)])
    ])
])

C_PCA = dbc.FormGroup([
    dbc.Row([
        dbc.Col([
            dcc.Markdown("""**Number of PCs**"""),
            dcc.Input(id='pre_pc',
                      type='number',
                      value=15)]),
        dbc.Col([
            dcc.Markdown("""**First PC**"""),
            dcc.Dropdown(id='pre_first',
                         options=[{'label': i, 'value': i} for i in ['Yes', 'No']],
                         value='Yes')]),
        dbc.Col([
            dcc.Markdown("""**Features**"""),
            dcc.Dropdown(id='pre_pc_f',
                         options=[{'label': i, 'value': i} for i in ['Variable genes', 'All']],
                         value='Variable genes')])
    ])
])


@app.callback(
    Output("pre_av", "children"),
    Input("pre_fs", "value")
)
def pre_dr_p(dr):
    if dr == 'LOESS':
        return C_LOESS
    elif dr == 'PCA':
        return C_PCA
    return C_LOESS


# Dimension Reduction
header_C_dr = dbc.FormGroup([
    dbc.Alert(
        "You can modify parameters for Step Dimension Reduction here!",
        color="success"),

    dbc.Row([

        dbc.Col(
            [dcc.Markdown("""**Number of neighbors**"""),
             dcc.Input(
                 id='pre_n_features',
                 type='number',
                 value=10
             )]),

        dbc.Col(
            [dcc.Markdown("""**Number of components**"""),
             dcc.Input(
                 id='pre_n_cells',
                 type='number',
                 value=5
             )]),

        dbc.Col(
            [dcc.Markdown("""**Feature**"""),
             dcc.Dropdown(
                 id='pre_log2',
                 options=[{'label': i, 'value': i} for i in ['Variable genes', 'Top PCs', 'All']],
                 value='Variable genes')]),

        dbc.Col(
            [dcc.Markdown("""**Dimension Reduction Methods**"""),
             dcc.Dropdown(
                 id='pre_nom',
                 options=[{'label': i, 'value': i} for i in available_DR],
                 value=available_DR[0])])
    ]),

    dbc.Row([
        dbc.Button(
            "Click for advanced parameters",
            id="collapse-button",
            className="mb-3",
            color="success",
            outline=True,
            size="sm",
            style={"text-align": 'center'}
        )
    ]),

    dbc.Row([
        dbc.Collapse(
            dbc.Card(dbc.CardBody("This content is for advanced parameters")),
            id="collapse",
        )
    ]),

    dbc.Row([
        dbc.Button(
            dbc.Spinner(id="C3-spinners"),
            id="Compute3-button",
            className="mb-3",
            color="success",
            block=True,
            disabled=True
        )
    ])
])


@app.callback(
    Output("collapse", "is_open"),
    Input("collapse-button", "n_clicks"),
    State("collapse", "is_open")
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("card_C_dr", "children"),
    Output("C3-spinners", "children"),
    Input("Compute3-button", "n_clicks")
)
def loading_output(n):
    if n:
        time.sleep(7)
        return card_C_dr_finished, 'Finished'
    return "This step is not ready yet. Click the COMPUTE button above and waite it to finish.", "Compute Step 3"


# Trajectory Inference
header_C_ti = dbc.FormGroup([
    dbc.Alert(
        "You can modify parameters for Step Trajectory Inference here!",
        color="success"),

    dbc.Row([

        dbc.Col(
            [dcc.Markdown("""**Number of neighbors**"""),
             dcc.Input(
                 id='pre_n_features',
                 type='number',
                 value=50
             )]),

        dbc.Col(
            [dcc.Markdown("""**lambda**"""),
             dcc.Input(
                 id='epg_lambda',
                 type='number',
                 value=0.02
             )]),

        dbc.Col(
            [dcc.Markdown("""**mu**"""),
             dcc.Input(
                 id='epg_mu',
                 type='number',
                 value=0.1
             )]),

        dbc.Col(
            [dcc.Markdown("""**alpha**"""),
             dcc.Input(
                 id='epg_alpha',
                 type='number',
                 value=0.02
             )]),

        dbc.Col(
            [dcc.Markdown("""**beta**"""),
             dcc.Input(
                 id='epg_beta',
                 type='number',
                 value=0.0
             )]),
    ]),

    dbc.Row([
        dbc.Button(
            "Click for advanced parameters",
            id="C4-collapse-button",
            className="mb-3",
            color="success",
            outline=True,
            size="sm",
            style={"text-align": 'center'}
        )
    ]),

    dbc.Row([
        dbc.Collapse(
            dbc.Card(dbc.CardBody("This content is for advanced parameters")),
            id="C4-collapse",
        )
    ]),

    dbc.Row([
        dbc.Button(
            dbc.Spinner(id="C4-spinners"),
            id="Compute4-button",
            className="mb-3",
            color="success",
            block=True,
            disabled=True
        )
    ])
])


@app.callback(
    Output("C4-collapse", "is_open"),
    Input("C4-collapse-button", "n_clicks"),
    State("C4-collapse", "is_open")
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("card_C_ti", "children"),
    Output("C4-spinners", "children"),
    Input("Compute4-button", "n_clicks")
)
def loading_output(n):
    if n:
        time.sleep(10)
        return card_C_ti_finished, 'Finished'
    return "This step is not ready yet. Click the COMPUTE button above and waite it to finish.", "Compute Step 4"


### Fill in contents for each card body
# DS
fig_ds = st.plot_stream(adata_computed, root='S1', return_svg=True)

card_assay = dbc.Card(
    [
        dbc.CardImg(src="/assets/Data_selection/Assays.png", top=True,
                    style={'width': '110px', 'margin-left': 'auto', 'margin-right': 'auto', 'marginTop': 5}),
        dbc.CardBody(
            html.H5("RNA", className="card-text", style={"text-align": 'center'})
        ),
    ], color="dark", outline=True
)

card_cell = dbc.Card(
    [
        dbc.CardImg(src="/assets/Data_selection/NumCells.png", top=True,
                    style={'width': '110px', 'margin-left': 'auto', 'margin-right': 'auto', 'marginTop': 5}),
        dbc.CardBody(
            html.H5([len(adata.obs), " Cells"], className="card-text", style={"text-align": 'center'})
        ),
    ], color="dark", outline=True
)

card_feature = dbc.Card(
    [
        dbc.CardImg(src="/assets/Data_selection/NumFeatures.png", top=True,
                    style={'width': '110px', 'margin-left': 'auto', 'margin-right': 'auto', 'marginTop': 5}),
        dbc.CardBody(
            html.H5([len(adata.var), " Features"], className="card-text", style={"text-align": 'center'})
        ),
    ], color="dark", outline=True
)

card_C_ds = dbc.FormGroup(
    [
        html.Img(src="data:image/svg+xml;base64,{}".format(fig_ds),
                 style={'width': '100%', 'height': '300px',
                        'lineHeight': '300px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        }),

        dbc.CardBody(
            [
                html.H4("Data Description", className="card-title", style={'textAlign': 'center'}),
                html.P(adata.uns['discription']),
                dbc.Row(
                    [
                        dbc.Col(card_assay), dbc.Col(card_cell), dbc.Col(card_feature)
                    ],
                    className="mb-4",
                )
            ]
        ),
    ],
)

# Pre
fig_qc = st.plot_qc(adata_computed, size=0, jitter=0, fig_size=(2, 2), return_svg=True)

card_C_pre_finished = dbc.FormGroup(
    [
        html.Img(src="data:image/svg+xml;base64,{}".format(fig_qc),
                 style={'width': '100%', 'height': '400px',
                        'lineHeight': '400px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        }),

        dbc.CardBody(
            [
                html.H4("Preprocessing", className="card-title", style={'textAlign': 'center'}),
                html.P("QC Calculation and Feature Selection are successfully finished!")
            ]
        ),
    ],
)

fig_dr = st.plot_dimension_reduction(adata_computed, color=['label'], n_components=3, show_graph=False, show_text=False,
                                     plotly=True, return_fig=True)

card_C_dr_finished = dbc.FormGroup(
    [
        dcc.Graph(figure=fig_dr),

        dbc.CardBody(
            [
                html.H4("Dimension Reduction", className="card-title", style={'textAlign': 'center'}),
                html.P("Yeah! Dimension reduction step is successfully completed!")
            ]
        ),
    ],
)

fig_ti = st.plot_dimension_reduction(adata_computed, color=['label'], n_components=3, show_graph=True,
                                     show_text=False, plotly=True, return_fig=True)
fig_st = st.plot_stream(adata_computed, root='S1', color=['label'], return_svg=True)

card_C_ti_finished = dbc.FormGroup(
    [
        dbc.Row([dcc.Graph(figure=fig_ti, style={'margin-left': 'auto', 'margin-right': 'auto',
                                                 'lineHeight': '60px',
                                                 'borderWidth': '1px',
                                                 'borderStyle': 'dashed',
                                                 'borderRadius': '5px'}),
                 html.Img(src="data:image/svg+xml;base64,{}".format(fig_st),
                          style={'margin-left': 'auto', 'margin-right': 'auto', 'width': '350px%', 'height': '350px',
                                 'lineHeight': '60px',
                                 'borderWidth': '1px',
                                 'borderStyle': 'dashed',
                                 'borderRadius': '5px'})]),

        dbc.CardBody(
            [
                html.H4("Trajectory Inference", className="card-title", style={'textAlign': 'center'}),
                html.P("Yeah! Trajectory inference step is successfully completed!")
            ]
        ),
    ],
)
