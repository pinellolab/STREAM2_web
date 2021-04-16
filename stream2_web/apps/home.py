import dash_html_components as html
import dash_bootstrap_components as dbc

# needed only if running this as a single page app
# external_stylesheets = [dbc.themes.LUX]

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

card_scRNA = [
    dbc.CardImg(src="/assets/scRNA-seq.png", top=True),
    dbc.CardBody(
        [
            html.H5("scRNA-seq", className="card-title"),
            html.P(
                "Lead to the visualization for scRNA-seq assay",
                className="card-text",
            ),
            dbc.Button("scRNA-seq", color="warning",href="/Visualization"),
        ]
    ),
]

card_scATAC = [
    dbc.CardImg(src="/assets/scATAC-seq.png", top=True),
    dbc.CardBody(
        [
            html.H5("scATAC-seq", className="card-title"),
            html.P(
                "Lead to the visualization for scATAC-seq assay",
                className="card-text",
            ),
            dbc.Button("scATAC-seq", color="warning"),
        ]
    ),
]

card_scProteomics = [
    dbc.CardImg(src="/assets/scProteomics.png", top=True),
    dbc.CardBody(
        [
            html.H5("scProteomics", className="card-title"),
            html.P(
                "Lead to the visualization for scProteomics assay",
                className="card-text",
            ),
            dbc.Button("scProteomics", color="warning"),
        ]
    ),
]

card_scDNA = [
    dbc.CardImg(src="/assets/scDNA.png", top=True),
    dbc.CardBody(
        [
            html.H5("scDNA", className="card-title"),
            html.P(
                "Lead to the visualization for scDNA methylation assay",
                className="card-text",
            ),
            dbc.Button("scDNA", color="warning"),
        ]
    ),
]

# change to app.layout if running as single page app instead
layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col(html.H1("Welcome to the STREAM2 webpage", className="text-center")
                    , className="mb-5 mt-5")
        ]),

        dbc.Row([
            dbc.Col(html.H5(
                children='STREAM2: Fast, scalable, and interactive trajectory analysis of single-cell omics data! ')
                , className="mb-4")
        ]),

        dbc.Row([
            dbc.Col(html.P(
                children='It consists of two main pages: Computation, which contains the end-to-end analysis pipline for multi-omics data, '
                         'Visualization, which gives the interactive multi-omics and assay specific data visualization.')
                , className="mb-4")
        ]),

        dbc.Row([
            dbc.Col(html.Img(src="/assets/HomeCarton.png", height="300px")
                    , className="mb-4")
        ], style={'textAlign': 'center'}),

        dbc.Row([
            dbc.Col(html.Hr(style={'borderWidth': "0.3vh", "width": "25%", "borderColor": "#0070C5"}), className="mb-4")
        ]),

        dbc.Row([
            dbc.Col(html.H4("Computation", className="text-center"), className="mb-4")
        ]),

        dbc.Row([
            dbc.Card(children=[html.Img(src="/assets/Computation.png", height="150px"),
                               dbc.Button("Computation",
                                          href="/Computation",
                                          color="info",
                                          className="mt-3"),
                               ],
                     body=True, color="secondary", outline=True)
        ]),

        dbc.Row([
            dbc.Col(html.Hr(style={'borderWidth': "0.3vh", "width": "25%", "borderColor": "#6DD400"}), className="mb-4")
        ]),

        dbc.Row([
            dbc.Col(html.H4("Visualization", className="text-center"), className="mb-4")
        ]),

        dbc.Row([
            dbc.Card(children=[html.Img(src="/assets/Visualization.png", height="150px"),
                               dbc.Button("Visualization",
                                          href="/Visualization",
                                          color="success",
                                          className="mt-3"),
                               ],
                     body=True, color="secondary", outline=True)
        ]),

        dbc.Row([
            dbc.Col(html.Hr(style={'borderWidth': "0.3vh", "width": "25%", "borderColor": "#FEC700"}), className="mb-4")
        ]),

        dbc.Row([
            dbc.Col(html.H4("Gallery", className="text-center"), className="mb-4")
        ]),

        dbc.Row([
            dbc.Col(html.P(
                children='Samples for applying STREAM2 on multi-omics assay')
                , className="mb-4")
        ]),

        dbc.Row([
            dbc.CardDeck([
                dbc.Card(card_scRNA, color="light"),
                dbc.Card(card_scATAC, color="light"),
                dbc.Card(card_scProteomics, color="light"),
                dbc.Card(card_scDNA, color="light"),
            ], className="mb-4")
        ]),

        html.A("Something else want to be put here.",
               href="#")

    ])

])

# needed only if running this as a single page app
# if __name__ == '__main__':
#     app.run_server(host='127.0.0.1', debug=True)
