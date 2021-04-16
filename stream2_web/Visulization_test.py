import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from dash.dependencies import Input, Output

import numpy as np
import stream as st
import matplotlib

matplotlib.use('Agg')

external_stylesheets = [dbc.themes.LUX]

adata=st.read(file_name='./SampleData/Nestorowa_2016/Nestorowa-2016.pkl',workdir='./stream_result')
fig_qc = st.plot_qc(adata,jitter=0.3,fig_size=(2,2),return_svg=True)

available_samples = ['Nestorowa, S. et al. 2016', 'Trapnell, C. et al. 2014', ' Tang, Q. et al. 2017']
available_projections = ['dimension_reduction', 'visualization_2D', 'flat_tree', 'branches']
available_colors = adata.obs.columns
available_stream = ['single_cell_stream', 'stream']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

card = dbc.Card(
    [
        dbc.CardHeader(
            dbc.Tabs(
                [
                    dbc.Tab(label="Data Selection", tab_id="tab-ds", label_style={"font-size":"22px"}, tab_style= {"margin-left":"auto","margin-right":"auto"}),
                    dbc.Tab(label="Quality Control", labelClassName="text-success", tab_id="tab-qc", label_style={"font-size":"22px"}, tab_style= {"margin-left":"auto","margin-right":"auto"}),
                    dbc.Tab(label="Dimension Reduction", labelClassName="text-success", tab_id="tab-dr", label_style={"font-size":"22px"}, tab_style= {"margin-left":"auto","margin-right":"auto"}),
                    dbc.Tab(label="Stream Plot", labelClassName="text-success", tab_id="tab-sp", label_style={"font-size":"22px"}, tab_style= {"margin-left":"auto","margin-right":"auto"}),
                ],
                id="card-tabs",
                card=True,
                active_tab='tab-ds'
            )
        ),
        dbc.CardBody(html.P(id="card-content", className="card-text"),style={"height": "40rem"}),
        dbc.CardFooter("This is the footer")
    ],
    color="dark", outline=True
)

app.layout = html.Div([
    dbc.Container([
        dbc.Col(html.H1("Visualization", className="text-center",style={"color": "#2F6B45"}), className="mb-5 mt-5"),
        dbc.Col(html.H5("Click the following tabs for step-by-step visualization.", className="text-center"),
                className="mb-5 mt-5"),
        dbc.Col(html.Hr(), className="mb-5 mt-5"),
        card,
        html.Div(id='page-content')
    ])
])

card_ds = dbc.FormGroup(
    [
        dbc.Row([
            dbc.Col(html.Div(
                dbc.Alert("You can simply choose from our sample datasets or upload your own STREAM object (currently only *.pkl and *.h5ad formats are supported)",color="success")
            ))
        ],align="center"),

        dbc.Row([
            dbc.Col(html.Hr(), className="mb-5 mt-5")
        ]),

        dbc.Row([
            dbc.Col(
                dcc.RadioItems(
                    id='Data',
                    options=[
                        {'label': 'Sample Datasets', 'value': 'sample'},
                        {'label': 'Personal Datasets', 'value': 'personal'}
                    ],
                    value='sample',
                    labelStyle={'display': 'inline-block',"font-size":"18px"}
                ),align="center"
            )
        ]),

        dbc.Row(dbc.Col(html.Div(id='data-selection')))
    ],
    row=True,
)


### Data Selection
@ app.callback(
    Output('data-selection', 'children'),
    Input('Data', 'value'))

def update_dr(Data):
    if Data == "sample":
        return dcc.Dropdown(options=[{'label': i, 'value': i} for i in available_samples],
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
            'textAlign': 'center'
        })

### tab activate
@app.callback(
    Output("card-content", "children"), [Input("card-tabs", "active_tab")]
)
def tab_content(active_tab):
    if active_tab == "tab-qc":
        return html.Img(src="data:image/svg+xml;base64,{}".format(fig_qc))
    else:
        return card_ds


if __name__ == '__main__':
    app.run_server(debug=True, host="0.0.0.0", port=8080)
