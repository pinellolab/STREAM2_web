import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

# must add this line in order for the app to be deployed successfully on Heroku
from app import server
from app import app
# import all pages in the app
from apps import home, Computation, Visualization, Tutorials, Gallery

# building the navigation bar
# https://github.com/facultyai/dash-bootstrap-components/blob/master/examples/advanced-component-usage/Navbars.py

Gallery_dropdown = dbc.DropdownMenu(
    children=[
        dbc.DropdownMenuItem("scRNA-seq"),
        dbc.DropdownMenuItem("scATAC-seq"),
        dbc.DropdownMenuItem("scProteomics"),
        dbc.DropdownMenuItem("scDNA Methylation")
    ],
    nav = True,
    in_navbar = True,
    label = "Gallery",
)

Contact_dropdown = dbc.DropdownMenu(
    children=[
        dbc.DropdownMenuItem("Pinello Lab", href="https://main.pinellolab.partners.org/"),
        dbc.DropdownMenuItem("Zinovyev Lab", href="https://auranic.github.io/"),
        dbc.DropdownMenuItem("Github", href="https://github.com/pinellolab/STREAM2"),
        dbc.DropdownMenuItem("Twitter", href="https://twitter.com/lucapinello"),
        dbc.DropdownMenuItem("EMAIL", href="mailto:LPINELLO@MGH.HARVARD.EDU"),
    ],
    nav = True,
    in_navbar = True,
    label = "Contact Us",
)

navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/logo.png", height="50px")),
                        dbc.Col(dbc.NavbarBrand("STREAM2", className="ml-1")),
                    ],
                    align="center",
                    no_gutters=True,
                ),
                href="/home",
            ),
            dbc.NavbarToggler(id="navbar-toggler2"),
            dbc.Collapse(
                dbc.Nav(
                    # right align dropdown menu with ml-auto className
                    [
                        dbc.NavItem(dbc.NavLink("Home", href="/home"), className="ml-2"),
                        dbc.NavItem(dbc.NavLink("Computation", href="/Computation"), className="ml-2"),
                        dbc.NavItem(dbc.NavLink("Visualization", href="/Visualization"), className="ml-2"),
                        dbc.NavItem(dbc.NavLink("Tutorials", href="https://readthedocs.org/"), className="ml-2"),
                        Gallery_dropdown,
                        Contact_dropdown],
                    className="ml-auto", navbar=True, pills=True
                ),
                id="navbar-collapse2",
                navbar=True,
            ),
        ]
    ),
    color="dark",
    dark=True,
    className="mb-4",
)

def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

for i in [2]:
    app.callback(
        Output(f"navbar-collapse{i}", "is_open"),
        [Input(f"navbar-toggler{i}", "n_clicks")],
        [State(f"navbar-collapse{i}", "is_open")],
    )(toggle_navbar_collapse)

# embedding the navigation bar
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/Computation':
        return Computation.layout
    elif pathname == '/Visualization':
        return Visualization.layout
    elif pathname == '/Tutorials':
        return Tutorials.layout
    else:
        return home.layout

if __name__ == '__main__':
    app.run_server(debug=True, host="0.0.0.0", port=8080)