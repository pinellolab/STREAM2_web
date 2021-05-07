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

### preset data
adata_computed = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016.pkl', workdir='./stream_result')
adata = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016-raw.h5ad', workdir='./stream_result')
adata.uns[
    'discription'] = 'This scRNA-seq dataset contains 1656 cells and 40594 genes from mouse hematopoietic stem and progenitor cell differentiation. A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation. Blood 128, e20-31 (2016).'
fig_ds = st.plot_stream(adata_computed, root='S1', return_svg=True)

### Set optionals
available_samples = ['Nestorowa, S. et al. 2016', 'Harrison, S. et al. 2021', 'Trapnell, C. et al. 2014',
                     ' Tang, Q. et al. 2017']
available_normalization = ['Library size correction', 'TF-IDF transformation', 'None']
available_DR = ['Spectral embedding(SE)', 'Modified locally linear embedding(MLLE)', 'UMAP', 'PCA']

###----- Header for data selection -----###
header_C_ds = dbc.FormGroup([
    html.Div(
        [
            dbc.Progress(value=20, color="success", striped=True),
            dbc.Alert("STEP 1. DATA SELECTION",color="success"),
            dbc.Alert("You can simply choose from our sample datasets or upload your own files",color="success")
        ]
    ),

    dcc.RadioItems(
        id='C-data-source',
        options=[
            {'label': 'Sample Datasets', 'value': 'sample'},
            {'label': 'Personal Datasets', 'value': 'personal'}
        ],
        inputStyle={"margin-right": "20px", "margin-left": "20px"},
        value='sample',
        labelStyle={'display': 'inline-block'}
    ),

    html.Div(id='C-data-selection')
])


@app.callback(
    Output('C-data-selection', 'children'),
    Input('C-data-source', 'value'))
def update_ds_source(Source):
    if Source == "sample":
        return dcc.Dropdown(id='ds_C_sample', options=[{'label': i, 'value': i} for i in available_samples],
                            value='Nestorowa, S. et al. 2016', style={'color': '#000000'})
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


###----- Header for Quality Control -----###
header_C_qc = dbc.FormGroup([
    html.Div(
        [
            dbc.Progress(value=40, color="success", striped=True),
            dbc.Alert("STEP 2. QUALITY CONTROL",color="success"),
            dbc.Alert(
                "You can modify parameters for Step2 Quality Control here!",
                color="success")
        ]
    ),

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
                 value='Yes', style={'color': '#000000'})]),

        dbc.Col(
            [dcc.Markdown("""**Normalization methods**"""),
             dcc.Dropdown(
                 id='pre_nom',
                 options=[{'label': i, 'value': i} for i in available_normalization],
                 value=available_normalization[0], style={'color': '#000000'})])
    ])
])


###----- Header for Feature Selection -----###
header_C_fs = dbc.FormGroup([
    html.Div(
        [
            dbc.Progress(value=60, color="success", striped=True),
            dbc.Alert("STEP 3. FEATURE SELECTION",color="success"),
            dbc.Alert(
                "You can do feature selection based on variable genes, top components, or all features",
                color="success")
        ]
    ),
     dbc.Col(
            [dcc.Markdown("""**Feature Selection Methods**"""),
             dcc.Dropdown(id='fs_C_type',
                          options=[
                            {'label': 'Variable Genes', 'value': 'vg'},
                            {'label': 'Top Components', 'value': 'tp'}],
                            value='vg', style={'color': '#000000'})]),
    dbc.Row(html.Hr()),
    dbc.Col(html.Div(id='C-feature-selection')),
])

@app.callback(
    Output('C-feature-selection', 'children'),
    Input('fs_C_type', 'value'))
def update_ds_source(Source):
    if Source == "vg":
        return dbc.Row([
            dbc.Col([dcc.Markdown("""**Loess Fraction**"""),
                     dcc.Input(id='fs_loess_frac', type='number',value=0.01)]),
            dbc.Col([dcc.Markdown("""**Percentile**"""),
                     dcc.Input(id='fs_percentile',type='number',value=95)])
            ])
    else:
        return dbc.FormGroup([
            dbc.Row([dbc.Col([dcc.Markdown("""**Feature**"""),
                              dcc.RadioItems(id='fs_feature',options=[{'label': i, 'value': i} for i in ['variable genes','all genes']],
                                             value='variable genes',labelStyle={'display': 'inline-block'})]),
                     dbc.Col(html.Div(id='C-tpc-feature')),
                     dbc.Col([dcc.Markdown("""**Number of components**"""),
                             dcc.Input(id='fs_n_pc', type='number', value=15)]),
                     dbc.Col([dcc.Markdown("""**Include First PC**"""),
                             dcc.RadioItems(
                                 id='fs_fpc',
                                 options=[{'label': 'Yes', 'value': 'True'},{'label': 'No', 'value': 'False'}],
                                 inputStyle={"margin-right": "20px", "margin-left": "20px"},
                                 value='True',
                                 labelStyle={'display': 'inline-block'})])])
            ])

@app.callback(
    Output('C-tpc-feature', 'children'),
    Input('fs_feature', 'value'))
def update_tpc_feature(Feature):
    if Feature == 'variable genes':
        return dbc.Row([
            dbc.Col([dcc.Markdown("""**Loess Fraction**"""),
                     dcc.Input(id='fs_tpc_loess_frac', type='number',value=0.01)]),
            dbc.Col([dcc.Markdown("""**Percentile**"""),
                     dcc.Input(id='fs_tpc_percentile',type='number',value=95)]) ])

###----- Basic layout of this Computation page -----###
Buttons = dbc.Row([
    dbc.Col(dbc.Button('<< previous step', id="C-Previous-button", className="mb-3", color="success", disabled=False,block=True)),
    dbc.Col(dbc.Button(dbc.Spinner(id="C-CONFIRM-spinners"), id="C-Confirm-button", className="mb-3", color='success', disabled=False, outline=True, block=True)),
    dbc.Col(dbc.Button('NEXT STEP >>', id="C-Next-button", className="mb-3", color="success", disabled=False, block=True))]
)

layout = html.Div([
    dbc.Container([

        dbc.Col(html.H1("Computation", className="text-center"), className="mb-5 mt-5"),
        dbc.Col(html.H5("Click the following tabs for step-by-step computation.", className="text-center"),
                className="mb-5 mt-5"),

        dbc.Card([
            dbc.CardHeader(html.Div(header_C_ds, id='header_parameters')),
            dbc.CardHeader(Buttons),
            dbc.CardBody(id="card-C-content", style={"height": "40rem"},
                         className="w-90"),
            dbc.CardFooter("If you like STREAM2 please support us by citing it in your work \N{TWO HEARTs}")
        ], color="dark", inverse=True, outline=False)
    ])
])

##### Card body for Data Selection
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

card_C_qc_finished = dbc.FormGroup(
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

headers = [header_C_ds, header_C_qc, header_C_fs, header_C_qc, header_C_ds, header_C_qc]
cards = [card_C_ds, card_C_qc_finished, card_C_ds, card_C_qc_finished, card_C_ds, card_C_qc_finished]
###----- main update -----###
@app.callback(
    Output("header_parameters", "children"),
    Output("C-Previous-button", "disabled"),
    Output("C-Next-button", "disabled"),
    Input("C-Previous-button", "n_clicks"),
    Input("C-Next-button", "n_clicks"),
)
def Update_header_Content(n_Previous, n_next):
    status_previous = False
    status_next = False

    if n_Previous is None or n_Previous == 'None':
        min = 0
    else:
        min = n_Previous
    if n_next is None or n_next == 'None':
        max = 0
    else:
        max = n_next
    count = max - min

    if count < 1:
        current_header = headers[0]
        status_previous = True
    if count == 5:
        status_next = True
        current_header = headers[count]
    else:
        current_header = headers[count]
    return current_header, status_previous, status_next

@app.callback(
    Output("card-C-content", "children"),
    Output("C-CONFIRM-spinners", "children"),
    Input("C-Confirm-button", "n_clicks"),
    State("C-Previous-button", "n_clicks"),
    State("C-Next-button", "n_clicks"),
)
def Update_Card_Content(n_confirm, n_Previous, n_next):
    if n_Previous is None or n_Previous == 'None':
        min = 0
    else:
        min = n_Previous
    if n_next is None or n_next == 'None':
        max = 0
    else:
        max = n_next
    count = max - min

    if n_confirm:
        time.sleep(5)
        return cards[count], 'CONFIRM'
    else:
        return 'Select data and click CONFIRM button above to start with STREAM2!','CONFIRM'






