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
global adata_computed
adata_computed = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016.pkl', workdir='./stream_result')
global adata
adata = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016-raw.h5ad', workdir='./stream_result')
adata.uns['discription'] = 'This scRNA-seq dataset contains 1656 cells and 40594 genes from mouse hematopoietic stem and progenitor cell differentiation. A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation. Blood 128, e20-31 (2016).'
global fig_ds
fig_ds = st.plot_stream(adata_computed, root='S1', return_svg=True)


### Set optionals
available_samples = ['Nestorowa, S. et al. 2016', 'Harrison, S. et al. 2021', 'Trapnell, C. et al. 2014', ' Tang, Q. et al. 2017']
available_normalization = ['Library size correction', 'TF-IDF transformation', 'None']
available_DR = ['Spectral embedding(SE)', 'Modified locally linear embedding(MLLE)', 'UMAP', 'PCA']



###----- Header for data selection -----###
header_C_ds = dbc.FormGroup([
    dbc.Alert(
        "You can simply choose from our sample datasets or upload your own files (10X Genomics output, Scanpy object and Seurat object are supported)",
        color="success"),

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

    dbc.Col(html.Div(id='C-data-selection')),

    dbc.Row(html.Hr()),

    dbc.Row([
        dbc.Col(dbc.Button('<< previous step',
                           id="C1-Previous-button",
                           className="mb-3",
                           color="success",
                           disabled=True,
                           block=True)),

        dbc.Col(dbc.Button(
            dbc.Spinner(id="C1-spinners"),
            id="C1-Confirm-button",
            className="mb-3",
            color='success',
            outline=True,
            block=True)),

        dbc.Col(dbc.Button('NEXT STEP >>',
                           id="C1-Next-button",
                           className="mb-3",
                           color="success",
                           disabled=True,
                           block=True))
        ])
])

###----- Basic layout of this Computation page -----###
default_header = dbc.FormGroup([html.Div(header_C_ds, id='header_parameters')])

layout = html.Div([
    dbc.Container([

        dbc.Col(html.H1("Computation", className="text-center"), className="mb-5 mt-5"),
        dbc.Col(html.H5("Click the following tabs for step-by-step computation.", className="text-center"),
                className="mb-5 mt-5"),

        dbc.Card([
            dbc.CardHeader(default_header),
            dbc.CardBody(id="card-C-content", style={"height": "40rem"}, className="w-90"),
            dbc.CardFooter("If you like STREAM2 please support us by citing it in your work \N{TWO HEARTs}")
        ], color="dark", inverse=True, outline=False)
    ])
])



###----- update Data Selection Source -----###
@app.callback(
    Output('C-data-selection', 'children'),
    Input('C-data-source', 'value'))
def update_ds_source(Source):
    if Source == "sample":
        return dcc.Dropdown(id='ds_C_sample', options=[{'label': i, 'value': i} for i in available_samples],
                            value='Nestorowa, S. et al. 2016', style={'color':'#000000'})
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
    [Output("card-C-content", "children"),
     Output("C1-spinners", "children"),
     Output("C1-Confirm-button", "disabled"),
     Output("C1-Next-button", "disabled")],
    Input("C1-Confirm-button", "n_clicks"),
)
def Update_Card_Content(n_confirm):

    if n_confirm:
        time.sleep(2)
        return card_C_ds, 'Successfully finished', True, False
    else:
        return 'Data is not confirmed yet! Click the CONFIRM button above to start!', 'Confirm Data', False, True

@app.callback(
    Output("header_parameters", "children"),
    Input("C1-Next-button", "n_clicks"),
)
def Update_Card_header(n_next):
    if n_next:
        return header_C_qc
    else:
        return header_C_ds


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







###----- Update Quality Control -----###
header_C_qc = dbc.FormGroup([
    dbc.Alert(
        "You can modify parameters for Step Quality Control here!",
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

    dbc.Row(html.Hr()),

    dbc.Row([
        dbc.Col(dbc.Button('<< previous step',
                           id="C2-Previous-button",
                           className="mb-3",
                           color="success",
                           disabled=True,
                           block=True)),

        dbc.Col(dbc.Button(
            dbc.Spinner(id="C2-spinners"),
            id="C2-Confirm-button",
            className="mb-3",
            color = 'success',
            outline=True,
            block=True)),

        dbc.Col(dbc.Button('NEXT STEP >>',
                            id="C2-Next-button",
                           className="mb-3",
                           color="success",
                           disabled=True,
                           block=True))
    ])
])

@app.callback(
    [Output("C2-spinners", "children"),
     Output("C2-Confirm-button", "disabled"),
     Output("C2-Next-button", "disabled")],
    Input("C2-Confirm-button", "n_clicks")
)
def C2_confirm(n):
    if n is None:
        return "Confirm parameters", False, True
    else:
        time.sleep(1)
        return  "Step 2 Quality Control is successfully finished", True, False

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
