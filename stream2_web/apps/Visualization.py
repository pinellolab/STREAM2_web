# -*- coding: utf-8 -*-
### import web components related packages
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

### import calculation related packages
import numpy as np
import stream as st
import matplotlib

matplotlib.use('Agg')

from app import app

### Load default sample data and do some basic modification to it
adata = st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016.pkl', workdir='./stream_result')
adata.uns[
    'discription'] = 'This scRNA-seq dataset contains 1656 cells and 4768 genes from mouse hematopoietic stem and progenitor cell differentiation. A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation. Blood 128, e20-31 (2016).'

### Set optionals
available_samples = ['Nestorowa, S. et al. 2016', 'Harrison, S. et al. 2021', 'Trapnell, C. et al. 2014',
                     ' Tang, Q. et al. 2017']
available_projections = ['dimension_reduction', 'visualization_2D', 'flat_tree', 'branches']
available_colors = adata.obs.columns.values.tolist()
# available_colors.insert(0, 'Feature Expression')
avaliable_features = adata.var.index.values.tolist()
available_stream = ['single_cell_stream', 'stream']
avaliable_divergingB = list(adata.uns['de_markers_greater'].keys())
avaliable_transitionB = list(adata.uns['transition_markers'].keys())
avaliable_leafB = list(adata.uns['leaf_markers'].keys())

### Basic layout of this visulization page
layout = html.Div([
    dbc.Container([

        dbc.Col(html.H1("Visualization", className="text-center"), className="mb-5 mt-5"),
        dbc.Col(html.H5("Click the following tabs for step-by-step visualization.", className="text-center"),
                className="mb-5 mt-5"),

        dbc.Card([
            dbc.CardHeader(
                dbc.Tabs(
                    [
                        dbc.Tab(label='1. Data Selection', tab_id="tab-ds", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='2. Quality Control', tab_id="tab-qc", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='3. Dimension Reduction', tab_id="tab-dr", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='4. Stream Plot', tab_id="tab-sp", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"}),

                        dbc.Tab(label='5. Gene visualization', tab_id="tab-gv", label_style={"font-size": "17px"},
                                tab_style={"margin-left": "auto", "margin-right": "auto"})
                    ],
                    id="card-tabs", card=True, active_tab='tab-ds')
            ),
            dbc.CardHeader(id="card-header"),
            dbc.CardBody(id="card-content", style={"height": "48rem"}, className="w-90"),
            dbc.CardFooter("If you like STREAM2 please support us by citing it in your work \N{TWO HEARTs}")
        ], color="dark", outline=True),

    ])
])


### tab activate
@app.callback(
    Output("card-header", "children"), Output("card-content", "children"), [Input("card-tabs", "active_tab")]
)
def tab_content(active_tab):
    if active_tab == "tab-qc":
        return header_qc, card_qc
    elif active_tab == "tab-dr":
        return header_dr, card_dr
    elif active_tab == "tab-sp":
        return header_sp, card_sp
    elif active_tab == "tab-gv":
        return header_gv, card_gv
    else:
        return header_ds, card_ds


### Fill in contents for each header
# DS
header_ds = dbc.FormGroup([
    dbc.Alert(
        "You can simply choose from our sample datasets or upload your own STREAM object (currently only *.pkl and *.h5ad formats are supported)",
        color="success"),

    dcc.RadioItems(
        id='data_source',
        options=[
            {'label': 'Sample Datasets', 'value': 'sample'},
            {'label': 'Personal Datasets', 'value': 'personal'}
        ],
        inputStyle={"margin-right": "20px", "margin-left": "20px"},
        value='sample',
        labelStyle={'display': 'inline-block'}
    ),

    html.Div(id='data-selection')
])


###### Data Selection
@app.callback(
    Output('data-selection', 'children'),
    Input('data_source', 'value'))
def update_ds(Source):
    if Source == "sample":
        return dcc.Dropdown(id='ds_sample', options=[{'label': i, 'value': i} for i in available_samples],
                            value='Nestorowa, S. et al. 2016')
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
            'textAlign': 'center',
        }, id='ds_personal')


# QC
header_qc = dbc.FormGroup([
    dbc.Alert("You can explore QC metrics here!", color="success")
])

# DR
header_dr = dbc.FormGroup([
    dbc.Alert(
        "You can modify detailed dimension reduction figures here!",
        color="success"),

    dbc.Row([
        dbc.Col(
            [dcc.Markdown("""**Alpha**"""),
             dcc.Slider(id='dr_Alpha', min=0, max=1, value=0.8,
                        marks={'%2.1f' % a: '%2.1f' % a for a in np.linspace(0, 1, 11)}, step=None)]),

        dbc.Col(
            [dcc.Markdown("""**Color Cells by**"""),
             dcc.RadioItems(id='dr_Color_By',
                            options=[
                                {'label': 'Feature Expression', 'value': 'feature'},
                                {'label': 'Others', 'value': 'others'}
                            ],
                            inputStyle={"margin-right": "20px", "margin-left": "20px"},
                            value='feature',
                            labelStyle={'display': 'inline-block'}),
             dcc.Dropdown(
                 id='dr_color_cells_by',
                 options=[{'label': i, 'value': i} for i in avaliable_features],
                 value=avaliable_features[0])])
    ])
])


@app.callback(
    Output('dr_color_cells_by', 'options'),
    Input('dr_Color_By', 'value'))
def update_dr(ColorBy):
    if ColorBy == "feature":
        return [{'label': i, 'value': i} for i in avaliable_features]
    else:
        return [{'label': i, 'value': i} for i in available_colors]


# SP
header_sp = header_dr

# GV
header_gv = dbc.FormGroup([
    dbc.Alert(
        "Start with your interested gene sets!",
        color="success"),

    dcc.RadioItems(
        id='gene_source',
        options=[
            {'label': 'Diverging Genes', 'value': 'diverging'},
            {'label': 'Transition Genes', 'value': 'transition'},
            {'label': 'Leaf Genes', 'value': 'leaf'},
        ],
        inputStyle={"margin-right": "20px", "margin-left": "20px"},
        value='diverging',
        labelStyle={'display': 'inline-block'}
    ),

    dbc.Row([
        dbc.Col(html.Div(children=[
            dbc.Row([
                dbc.Col([dcc.Markdown("""**Branch for Diverging Gene Analysis**"""),
                         dcc.Dropdown(id='gv_branch',
                                      options=[
                                          {'label': 'Branch ' + str(i[0][0] + '-' + i[0][1] + ' and ' + 'Branch ' + i[1][0] + '-' + i[1][1]),
                                           'value': str(i)} for i in
                                          avaliable_divergingB],
                                      value=str(avaliable_divergingB[0]))]),

                dbc.Col([dcc.Markdown("""**Relatively Highly Expressed On:**"""),
                         dcc.Dropdown(id='gv_diverging_high')])
            ])
        ], id='gene-visualization'), width=8),

        dbc.Col([dcc.Markdown("""**Gene for visualization on STREAM plots**"""),
                 dcc.Dropdown(id='gv_gene',
                              options=[{'label': i, 'value': i} for i in avaliable_features],
                              value=avaliable_features[0])], width=4)
    ])
])


###### Gene Visualization

@app.callback(
    Output('gene-visualization', 'children'),
    Output('gv_branch', 'value'),
    Input('gene_source', 'value'))
def update_gv(Source):
    if Source == "transition":
        return html.Div([
            dbc.Row([
                dbc.Col([dcc.Markdown("""**Branch for Transition Gene Analysis**"""),
                         dcc.Dropdown(id='gv_branch',
                                      options=[{'label': 'Node ' +str(i[0] + ' to Node ' + i[1]), 'value': str(i)} for i in
                                               avaliable_transitionB],
                                      value=str(avaliable_transitionB[0]))]),
                dbc.Col([html.Div(id='gv_diverging_high')])
            ])
        ]), str(avaliable_transitionB[0])

    elif Source == "leaf":
        return html.Div([
            dbc.Row([
                dbc.Col([dcc.Markdown("""**Branch for Leaf Gene Analysis**"""),
                         dcc.Dropdown(id='gv_branch',
                                      options=[{'label': 'Branch ' + str(i[0] + '-' + i[1]), 'value': str(i)} for i in
                                               avaliable_leafB],
                                      value=str(avaliable_leafB[0]))]),

                dbc.Col([html.Div(id='gv_diverging_high')])
            ]),
        ]), str(avaliable_leafB[0])

    else:
        return html.Div([
            dbc.Row([
                dbc.Col([dcc.Markdown("""**Branch for Diverging Gene Analysis**"""),
                         dcc.Dropdown(id='gv_branch',
                                      options=[
                                          {'label': 'Branch ' + str(i[0][0] + '-' + i[0][1] + ' and ' + 'Branch ' + i[1][0] + '-' + i[1][1]),
                                           'value': str(i)} for i in
                                          avaliable_divergingB],
                                      value=str(avaliable_divergingB[0]))]),

                dbc.Col([dcc.Markdown("""**Relatively Highly Expressed On:**"""),
                         dcc.Dropdown(id='gv_diverging_high')])
            ])
        ]), str(avaliable_divergingB[0])


@app.callback(
    Output('gv_diverging_high', 'options'),
    Output('gv_diverging_high', 'value'),
    Input('gene_source', 'value'),
    Input('gv_branch', 'value'))
def update_gv(Source, Branch):
    if Source == "diverging":
        return [{'label': 'Branch ' + str(i[0] + '-' + i[1]), 'value': str(i)} for i in eval(Branch)], str(eval(Branch)[0])
    else:
        return None, None


### Fill in contents for each card body
# DS
fig_ds = st.plot_stream(adata, root='S1', return_svg=True)

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

card_ds = dbc.FormGroup(
    [
        html.Img(src="data:image/svg+xml;base64,{}".format(fig_ds),
                 style={'width': '100%', 'height': '400px',
                        'lineHeight': '400px',
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

# QC
fig_qc = st.plot_qc(adata, size=0, jitter=0, fig_size=(2, 2), return_svg=True)

card_qc = dbc.FormGroup(
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
                html.H4("Quality Control", className="card-title", style={'textAlign': 'center'}),
                html.P("In quality control, you can see the distribution of number of cells and number of features!")
            ]
        ),
    ],
)

# DR

card_dr = html.Div(id="dr-graphic")


@app.callback(
    Output('dr-graphic', 'children'),
    Input('dr_color_cells_by', 'value'),
    Input('dr_Alpha', 'value'))
def update_dr(color, alpha):
    fig_dr = st.plot_dimension_reduction(adata, color=[color], n_components=3, alpha=alpha, show_graph=True,
                                         show_text=False, plotly=True, return_fig=True)
    fig_dr.update_layout(
        autosize=False,
        width=450,
        height=350,
        margin=dict(
            l=30,
            r=30,
            b=5,
            t=5,
            pad=4
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig_ft = st.plot_visualization_2D(adata, method='umap', n_neighbors=50, alpha=alpha, color=[color],
                                      use_precomputed=True, plotly=True, return_fig=True)
    fig_ft.update_layout(
        autosize=False,
        width=450,
        height=350,
        margin=dict(
            l=30,
            r=30,
            b=5,
            t=5,
            pad=4
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig_2d = st.plot_flat_tree(adata, color=[color], alpha=alpha, dist_scale=0.5, show_graph=True, show_text=True,
                               plotly=True, return_fig=True)
    fig_2d.update_layout(
        autosize=False,
        width=450,
        height=350,
        margin=dict(
            l=30,
            r=30,
            b=5,
            t=5,
            pad=4
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )

    return html.Div([
        dbc.Row([dcc.Graph(figure=fig_dr, style={'margin-left': 'auto', 'margin-right': 'auto',
                                                 'lineHeight': '60px',
                                                 'borderWidth': '1px',
                                                 'borderStyle': 'dashed',
                                                 'borderRadius': '5px'}),
                 dcc.Graph(figure=fig_dr, style={'margin-left': 'auto', 'margin-right': 'auto',
                                                 'lineHeight': '60px',
                                                 'borderWidth': '1px',
                                                 'borderStyle': 'dashed',
                                                 'borderRadius': '5px'})]),
        html.Br(),
        dbc.Row([dcc.Graph(figure=fig_ft, style={'margin-left': 'auto', 'margin-right': 'auto',
                                                 'lineHeight': '60px',
                                                 'borderWidth': '1px',
                                                 'borderStyle': 'dashed',
                                                 'borderRadius': '5px'}),
                 dcc.Graph(figure=fig_2d, style={'margin-left': 'auto', 'margin-right': 'auto',
                                                 'lineHeight': '60px',
                                                 'borderWidth': '1px',
                                                 'borderStyle': 'dashed',
                                                 'borderRadius': '5px'})])
    ])


# SP
card_sp = html.Div(id="sp-graphic")


@app.callback(
    Output('sp-graphic', 'children'),
    Input('dr_color_cells_by', 'value'),
    Input('dr_Alpha', 'value'))
def update_dr(color, alpha):
    fig_sc = st.plot_stream_sc(adata, root='S1', color=[color], alpha=alpha,
                               dist_scale=0.3, show_graph=True, show_text=True, plotly=True, return_fig=True)
    fig_sc.update_layout(
        autosize=False,
        width=450,
        height=350,
        margin=dict(
            l=30,
            r=30,
            b=5,
            t=5,
            pad=4
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig_st = st.plot_stream(adata, root='S1', color=[color], return_svg=True)

    fig_dr = st.plot_dimension_reduction(adata, color=[color], n_components=3, alpha=alpha, show_graph=True,
                                         show_text=False, plotly=True, return_fig=True)
    fig_dr.update_layout(
        autosize=False,
        width=450,
        height=350,
        margin=dict(
            l=30,
            r=30,
            b=5,
            t=5,
            pad=4
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )

    return html.Div([
        dbc.Row([dcc.Graph(figure=fig_sc, style={'margin-left': 'auto', 'margin-right': 'auto',
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
        html.Br(),
        dbc.Row([dcc.Graph(figure=fig_dr, style={'margin-left': 'auto', 'margin-right': 'auto',
                                                 'lineHeight': '60px',
                                                 'borderWidth': '1px',
                                                 'borderStyle': 'dashed',
                                                 'borderRadius': '5px'})
                 ])
    ])


# GV

card_gv = html.Div(id="gv-graphic")


@app.callback(
    Output('gv-graphic', 'children'),
    Input('gene_source', 'value'),
    Input('gv_branch', 'value'),
    Input('gv_diverging_high', 'value'),
    Input('gv_gene', 'value')
)
def update_dr(gene_source, branch, high, gene):
    if high is None:
        high = str(eval(branch)[0])
    if gene is None:
        gene = avaliable_features[0]

    if gene_source == "transition":
        data_df = adata.uns['transition_markers'][eval(branch)].round(2)

    elif gene_source == "leaf":
        data_df = adata.uns['leaf_markers'][eval(branch)].astype(float).round(2)

    elif gene_source == "diverging":
        if high == str(eval(branch)[0]):
            data_df = adata.uns['de_markers_greater'][eval(branch)].round(2)
        else:
            data_df = adata.uns['de_markers_less'][eval(branch)].round(2)

    data_df = data_df.reset_index()
    data_df.rename(columns={'index': 'gene'}, inplace=True)

    fig = go.Figure(data=[go.Table(
        header=dict(
            values=list(data_df.columns),
            line_color='darkslategray', fill_color='lightgrey',
            align='center', font=dict(color='black', size=12)
        ),
        cells=dict(values=data_df.T.values,
                   line_color='darkslategray', align='center', font=dict(color='black', size=11)
                   ))
    ])

    fig_sc = st.plot_stream_sc(adata, root='S1', color=[gene],
                               dist_scale=0.3, show_graph=True, show_text=True, plotly=True, return_fig=True)
    fig_sc.update_layout(
        autosize=False,
        width=450,
        height=350,
        margin=dict(
            l=30,
            r=30,
            b=5,
            t=5,
            pad=4
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig_st = st.plot_stream(adata, root='S1', color=[gene], return_svg=True)

    return html.Div([
        dbc.Row([dcc.Graph(figure=fig,
                           style={"maxHeight": "350px", 'width': '100%', "overflow": "scroll", 'margin-left': 'auto',
                                  'margin-right': 'auto'})]),
        html.Br(),
        dbc.Row([
            dcc.Graph(figure=fig_sc, style={'margin-left': 'auto', 'margin-right': 'auto',
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
    ])
