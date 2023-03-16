from ast import Import
import plotly.graph_objects as go

import json
from textwrap import dedent as d

import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State

import dash_bootstrap_components as dbc

import math

#from app import app

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
# external_stylesheets = [dbc.themes.GRID]
external_stylesheets = [dbc.themes.BOOTSTRAP]

try:
    from jupyter_dash import JupyterDash
except ImportError:
    print('jupyter_dash required: install using "conda install -c plotly conda install -c plotly jupyter-dash"')
    sys.exit()
    
app = JupyterDash(__name__, external_stylesheets=external_stylesheets)
# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

server = app.server
app.config.suppress_callback_exceptions = True

import json

###########################################e

import pandas as pd
import networkx as nx

import arg
from arg import Coalescent, Recombination, Leaf, interval_sum, interval_diff, interval_intersect, get_breakpoints, get_child_lineages, rescale_positions, marginal_arg, traverse_marginal, marginal_trees

import plotly.colors

def get_continuous_color(colorscale, intermed):
    """
    Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
    color for any value in that range.

    Plotly doesn't make the colorscales directly accessible in a common format.
    Some are ready to use:
    
        colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

    Others are just swatches that need to be constructed into a colorscale:

        viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
        colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

    :param colorscale: A plotly continuous colorscale defined with RGB string colors.
    :param intermed: value in the range [0, 1]
    :return: color in rgb string format
    :rtype: str
    """
    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    if intermed <= 0 or len(colorscale) == 1:
        return colorscale[0][1]
    if intermed >= 1:
        return colorscale[-1][1]

    for cutoff, color in colorscale:
        if intermed > cutoff:
            low_cutoff, low_color = cutoff, color
        else:
            high_cutoff, high_color = cutoff, color
            break

    # noinspection PyUnboundLocalVariable
    return plotly.colors.find_intermediate_color(
        lowcolor=low_color, highcolor=high_color,
        intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
        colortype="rgb")

layout = html.Div(
    [
        # Hidden div inside the app that stores the intermediate value
        html.Div(id='intermediate-value', style={'display': 'none'}),

        # row for arg and marginal trees
        dbc.Row(
            [
                # column for arg 
                dbc.Col(
                    [
                        # arg
                        dbc.Container(
                            [
                                dbc.Container(
                                    [

                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        html.B("Simulation:"),
                                                        # "Simulation:",
                                                        dcc.Dropdown(
                                                            id='sim-dropdown',
                                                            options=[
                                                                {'label': "ARG", 'value': 'arg'},
                                                                # {'label': "SMC", 'value': 'smc'},
                                                                {'label': "SMC'", 'value': 'smcprime'},
                                                                {'label': "SMC", 'value': 'smc'}
                                                            ],
                                                            value='arg', searchable=False, clearable=False
                                                        ),
                                                    ], width=2
                                                ),                                                        
                                                dbc.Col(
                                                    [
                                                        html.B("Nr samples:"),
                                                        dcc.Dropdown(
                                                            id='samples-dropdown',
                                                            options=[
                                                                {'label': "3", 'value': 3},
                                                                {'label': "4", 'value': 4},
                                                                {'label': "5", 'value': 5}
                                                            ],
                                                            value=5, searchable=False, clearable=False,
                                                            style={
                                                                    # 'height': '20px', 
                                                                    # 'width': '80px', 
                                                                    'font-size': "0.85rem",
                                                                    # 'min-height': '1px',
                                                                    },
                                                        ),
                                                    ], width=2
                                                ),  
                                                dbc.Col(
                                                    [
                                                        html.B("Length:"),
                                                        dcc.Dropdown(
                                                            id='seqlen-dropdown',
                                                            options=[
                                                                {'label': "1kb", 'value': 1e+3},
                                                                {'label': "2kb", 'value': 2e+3},
                                                                {'label': "4kb", 'value': 4e+3}
                                                            ],
                                                            value=2e+3, searchable=False, clearable=False,
                                                            style={
                                                                    # 'height': '20px', 
                                                                    # 'width': '80px', 
                                                                    'font-size': "0.85rem",
                                                                    # 'min-height': '1px',
                                                                    },
                                                        ),
                                                    ], width=2
                                                ),                                                                                                        
                                                dbc.Col(
                                                    [ 
                                                        # html.Button('New simulation', id='new-arg-button')
                                                        dbc.Button('New', id='new-arg-button', 
                                                            color="primary", #size="sm", #outline=True,
                                                            style={'height': 35, 'font-size': "0.85rem"},
                                                            className="mr-1"
                                                            )
                                                    ], width=3
                                                ),    
                                                dbc.Col(
                                                    [
                                                        html.Div(id='arg-header'),

                                                        # dcc.Markdown(d("""
                                                        # **Ancestral recombination graph:**   
                                                        # Nodes are colored by amount of ancestral sequence.
                                                        # """), ),                    
                                                    ], width=3
                                                ),                                                                                                                                                                                                        
                                            ], justify="between", align="end", style={'padding': 3} 
                                        ),
                                        dcc.Graph(id='arg-figure',
                                                clear_on_unhover=True,
                                                figure={'layout': {
                                                            'height': 570,
                                                            # 'margin': {'l': 0, 'b': 0, 't': 0, 'r': 0},
                                                                }
                                                            },
                                                    ),
                                    ], className='pretty_container', fluid=True,
                                ),
                            ], style={'padding': 20}
                        ),
                    ], width=8, 
                ),

                # column for marginal trees
                dbc.Col(
                    [
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dcc.Markdown(d("""
                                        **Marginal tree(s):** Hover over an ARG node.
                                        """), ),                    
                                        dcc.Graph(id='marginal-tree',
                                                    figure={'layout': {
                                                            # 'title': 'Marginal tree',
                                                            'height': 250,
                                                            # 'margin': {'l': 10, 'b': 0, 't': 0, 'r': 0},
                                                                }
                                                            },),
                                    ], className='pretty_container'
                                ),
                            ], style={'padding': 20, 'padding-left': 0, 'padding-bottom': 0}
                        ),
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dcc.Markdown(d("""
                                        **Ancestral sequences:** Hover over an ARG node.
                                        """), ),                    
                                        dcc.Graph(id='ancestral-sequence',
                                                figure={'layout': {
                                                    'height': 250,
                                                        }
                                                        },),
                                    ], className='pretty_container'
                                ),            
                            ], style={'padding': 20, 'padding-left': 0}
                        ),            
                    ], width=4, 
                ),        
            ], 
            className="g-0"
        ),

        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dbc.Container(
                                            [
                                                dcc.Markdown(d("""
                                                **Coalesce and recombination events:** 
                                                Slide to see progression of events.
                                                """)),
                                            ], 
                                        ),
                                        dbc.Container(
                                            [
                                                dcc.Slider(
                                                    id='event-slider',
                                                    min=0, max=40, value=0, 
                                                    marks={str(i): str(i) for i in range(0, 40)}, 
                                                    step=None,
                                                    ),
                                            ], style={'padding-bottom': 20}
                                        )
                                    ], className='pretty_container'
                                )
                            ], style={'padding': 20, 'padding-top': 0},
                        )
                    ], width=6
                ),  
                dbc.Col(
                    [
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dbc.Container(
                                            [
                                                dcc.Markdown(d("""
                                                    **Recombination points:** 
                                                    Slide to see graph for only part of the sequence.
                                                """)),
                                            ]
                                        ),
                                        dbc.Container(
                                            [
                                                dcc.RangeSlider(
                                                    id='seq-slider',
                                                    min=0,
                                                    max=1000,
                                                    value=[0, 1000],
                                                    # step=None,
                                                    marks={0: '0', 1000: '1'},
                                                    pushable=30,
                                                )
                                            ], style={'padding-bottom': 20}
                                        ),
                                    ], className='pretty_container', 
                                ),
                            ], style={'padding': 20, 'padding-left': 0, 'padding-top': 0}
                        ),
                    ], width=6, align='start',
                ),
            ],
            className="g-0"
        ),
    ], style={'padding': 20}
)


def get_bezier_points(x1, y1, x2, y2, offset=1, absolute=True):
    mid_x = x1 + (x2 - x1) / 2
    mid_y = y1 + (y2 - y1) / 2
    length = math.sqrt((x2-x1)**2 + (y2-y1)**2)
    if absolute:
        hyp = offset
    else:
        hyp = length / offset 
        
    if y2 == y1:
        b11, b12 = mid_x, mid_y + hyp
        b21, b22 = mid_x, mid_y - hyp    
    elif x2 == x1:
        b11, b12 = mid_x + hyp, mid_y
        b21, b22 = mid_x - hyp, mid_y
    else:
        slope = (y2-y1) / (x2-x1)
        if slope < 0:
            angle_in_radians = math.atan(1/slope)
        else:
            angle_in_radians = math.atan(slope)   
        kat = math.sin(angle_in_radians) * hyp
        b11, b12 = mid_x + kat*slope, mid_y + kat/slope
        b21, b22 = mid_x - kat*slope, mid_y - kat/slope
    
    return b11, b12, b21, b22


def arg_figure_data(nodes):

    traces = []

    edge_x = []
    edge_y = []
    
    diamond_shapes = []
    
    # for lineage in get_parent_lineages(nodes, root=False):
    for lineage in get_child_lineages(nodes):
        
        if type(lineage.down) is Recombination and \
            type(lineage.up) is Coalescent and \
            set(lineage.up.children) == set([lineage.down.right_parent, lineage.down.left_parent]):

            # diamond recombination:
            x1 = lineage.down.xpos
            y1 = lineage.down.height

            x2 = lineage.up.xpos
            y2 = lineage.up.height
            
            b11, b12, b21, b22 = get_bezier_points(x1, y1, x2, y2, offset=0.02)

            diamond_shapes.append(
                dict(
                    type="path",
                    path=f"M {x1},{y1} Q {b11},{b12} {x2},{y2}",
                    # line_color="lightgray",
                    layer='below',
                    line= {'width': 2, 'color': 'gray'}
                )
            )
            diamond_shapes.append(
                dict(
                    type="path",
                    path=f"M {x1},{y1} Q {b21},{b22} {x2},{y2}",
                    # line_color="lightgray",
                    layer='below',
                    line= {'width': 2, 'color': 'gray'},
                )
            )            

        else:        
            # start
            edge_x.append(lineage.down.xpos)
            edge_y.append(lineage.down.height)
            # end
            edge_x.append(lineage.up.xpos)
            edge_y.append(lineage.up.height)
            # gap
            edge_x.append(None)
            edge_y.append(None)

    traces.append(dict(
        x=edge_x,
        y=edge_y,
        mode='lines',
        opacity=1,
        hoverinfo = 'skip',
        line={
            'color': 'grey',
        },
        name=''
    ))

    node_x = []
    node_y = []    
    node_text = []
    node_color = []
    for node in nodes:
        node_x.append(node.xpos)
        node_y.append(node.height)
        prop_ancestral = 1
        if type(node) is Coalescent:
            prop_ancestral = interval_sum(node.parent.intervals)
        elif type(node) is Recombination:
            prop_ancestral = interval_sum(node.child.intervals)
        node_text.append(f"Fraction ancestral: {round(prop_ancestral, 2)}<br>Event: {type(node).__name__}")

        node_color.append(prop_ancestral)

    traces.append(dict(
        x=node_x,
        y=node_y,
        text=node_text,
        # range_color=[0, 1],
        # cmin=0,
        # cmax=1,
        mode='markers',
        opacity=1,
        hoverinfo ='text',
        marker={
            'size': 10,
            'color': node_color,
            'cmin': 0,
            'cmax': 1,
            'line': {'width': 0.7, 'color': 'white'},
            # 'colorscale': 'Viridis',
            'colorscale': 'Rainbow',
            'colorbar': {'title': 'Fraction<br>ancestral<br>sequence',
                        'titleside': 'top',
                        'thickness': 15,
                        'len': 0.5,
                        # 'tickmode': 'array',
                        'tickvals': [0, 0.5, 1],
                        # 'ticktext': ['0', '1'],
                        'ticks': 'outside',
                        },
        },
        name=''
    ))

    return dict(data=traces,
                layout=dict(xaxis=dict(fixedrange=True, 
                                       range=[-0.1, 1.1], #title='Samples',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            yaxis=dict(fixedrange=True, 
                                       range=[-0.1, 1.1], #title='Time',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            hovermode='closest',
                            range_color=[0,1],
                            margin= {'l': 50, 'b': 20, 't': 20, 'r': 20},
                            transition = {'duration': 0},
                            showlegend=False,
                            shapes=diamond_shapes,
                            )
                )


def tree_figure_data(node_lists):

    traces = []

    edge_x = []
    edge_y = [] 
    node_x = []
    node_y = []    
    node_color = []

    for i, nodes in enumerate(node_lists):

        # for lineage in get_parent_lineages(nodes, root=False):
        for lineage in get_child_lineages(nodes):
            # start
            edge_x.append(lineage.down.xpos)
            edge_y.append(lineage.down.height)
            # end
            edge_x.append(lineage.up.xpos)
            edge_y.append(lineage.up.height)
            # gap
            edge_x.append(None)
            edge_y.append(None)

        for node in nodes:
            node_x.append(node.xpos)
            node_y.append(node.height)

            node_color.append(i/len(node_lists))

    traces.append(dict(
        x=edge_x,
        y=edge_y,
        mode='lines',
        opacity=1,
        hoverinfo = 'skip',
        line={
            'color': 'grey',
        },
        name=''
    ))

    traces.append(dict(
        x=node_x,
        y=node_y,
        mode='markers',
        opacity=1,
        hoverinfo ='text',
        marker={
            'size': 7,
            'color': node_color,
            'cmin': 0,
            'cmax': 1,
            'colorscale': 'Rainbow',
            'line': {'width': 0.3, 'color': 'white'},
        },
        name=''
    ))

    return dict(data=traces,
                layout=dict(xaxis=dict(fixedrange=True, 
                                       range=[-0.05, 1.05], #title='Samples',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            yaxis=dict(fixedrange=True, 
                                       range=[-0.1, 1.1], #title='Time',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            hovermode='closest',
                            range_color=[0,1],
                            margin= {'l': 50, 'b': 20, 't': 20, 'r': 20},
                            transition = {'duration': 0},
                            showlegend=False,
                            )
                )


@app.callback(
    Output('arg-header', 'children'),
    [Input('new-arg-button', 'n_clicks')])
def update_header(n_clicks):

    if n_clicks is None:
        n_sim = 1
    else:
        n_sim = n_clicks + 1

    return dcc.Markdown(d("""
                **Simulation #{}:**   
                """.format(n_sim)))

@app.callback(Output('intermediate-value', 'children'), 
    [Input('new-arg-button', 'n_clicks'),
     Input('sim-dropdown', 'value'),
     Input('samples-dropdown', 'value'),
     Input('seqlen-dropdown', 'value')])
def new_data(n_clicks, sim, samples, length):

    nodes = arg.get_arg_nodes(L=length, n=samples, simulation=sim)
#    rescale_positions(nodes)
    json_str = arg.arg2json(nodes)
    return json_str

@app.callback(
    [Output(component_id='event-slider', component_property='min'),
     Output(component_id='event-slider', component_property='max'),
     Output(component_id='event-slider', component_property='step'),
     Output(component_id='event-slider', component_property='value')],
    [Input('intermediate-value', 'children')])    
def update_event_slider(jsonified_data):
    if jsonified_data:
        nodes = arg.json2arg(jsonified_data)
    else:
        nodes = []

    nr_leaves = len([n for n in nodes if type(n) is arg.Leaf])
    nr_events = len(nodes)-nr_leaves
    return 0, nr_events, 1, nr_events

@app.callback(
    [Output(component_id='seq-slider', component_property='min'),
     Output(component_id='seq-slider', component_property='max'),
     Output(component_id='seq-slider', component_property='value'),
     Output(component_id='seq-slider', component_property='marks')],
    [Input('intermediate-value', 'children')])    
def update_seq_slider(jsonified_data):
    if jsonified_data:
        nodes = arg.json2arg(jsonified_data)
    else:
        nodes = []
    breakpoints = get_breakpoints(nodes)
    marks = dict((b*1000, str(i+1)) for i, b in enumerate(breakpoints))
    if not marks:
        marks = None
    return 0, 1000, [0, 1000], marks


@app.callback(
    Output('arg-figure', 'figure'),
    [Input('intermediate-value', 'children'),
     Input('event-slider', 'value'),
     Input('seq-slider', 'value')])
def update_arg_figure(jsonified_data, event, interval):

    if jsonified_data:
        nodes = arg.json2arg(jsonified_data)

        interval = [i/1000 for i in interval]

        # Get marginal arg for interval
        marg_arg_nodes = marginal_arg(nodes, interval)
        # print(interval)
        # get only subset of events
        nr_leaves = len([n for n in nodes if type(n) is arg.Leaf])
        new_nodes = marg_arg_nodes[:nr_leaves+event]
    else:
        new_nodes = []

    return arg_figure_data(new_nodes)


@app.callback(
    Output('marginal-tree', 'figure'),
    [Input('intermediate-value', 'children'),
     Input('arg-figure', 'hoverData'),
     Input('seq-slider', 'value')])
def update_marg_tree_figure(jsonified_data, hover, slider_interval):

    marg_tree_list = []
    if hover and jsonified_data:
        nodes = arg.json2arg(jsonified_data)
        focus_node_idx = hover['points'][0]['pointIndex']
        focus_node = nodes[focus_node_idx]

        if type(focus_node) is Recombination:
            intervals = focus_node.child.intervals
        else:
            intervals = focus_node.parent.intervals

        # slider interval is 0-1000 not 0-1:
        slider_interval = [x/1000 for x in slider_interval]
        # get part of intervals that intersect slider interval:
        intervals = interval_intersect([slider_interval], intervals)

        for interval in intervals:
            # get marginal arg under focus node            
            new_nodes = traverse_marginal(focus_node, interval)
            new_nodes = list(new_nodes)
            new_nodes.sort(key=lambda x: x.height)

            marg_trees, _ = marginal_trees(new_nodes, interval)
            marg_tree_list.extend(marg_trees)
            
        nr_cols = len(marg_tree_list)

        space = 0.5
    
        for i in range(nr_cols):
            tree = marg_tree_list[i]
            for node in tree:
                node.xpos = node.xpos/(nr_cols+(nr_cols-1)*space) + i/nr_cols
            marg_tree_list[i] = tree

    # TODO: Maybe keep "dangling root" branch here

    if marg_tree_list:
        return(tree_figure_data(marg_tree_list))
    else:
        return(tree_figure_data([]))


@app.callback(
    Output('ancestral-sequence', 'figure'),
    [Input('intermediate-value', 'children'),
     Input('arg-figure', 'hoverData'),
     Input('seq-slider', 'value')])    
def update_ancestral_seq_figure(jsonified_data, hover, slider_interval):

    traces = []
    shape_list = []
    
    if hover and jsonified_data:
        nodes = arg.json2arg(jsonified_data)
        focus_node_idx = hover['points'][0]['pointIndex']
        focus_node = nodes[focus_node_idx]

        # slider interval is 0-1000 not 0-1:
        slider_interval = [x/1000 for x in slider_interval]

        gray_segments = list(map(tuple, interval_diff([[0, 1]], [slider_interval])))

        def get_segments(focus_node, intervals):
            segments = list()
            marg_tree_list = list()
            for interval in intervals:  
                new_nodes = traverse_marginal(focus_node, interval)
                new_nodes = list(new_nodes)
                new_nodes.sort(key=lambda x: x.height)
                marg_trees, marg_segm = marginal_trees(new_nodes, interval)
#                 print(marg_trees, marg_segm)
                marg_tree_list.extend(marg_trees)
                segments.extend(marg_segm)
            return segments

        def get_shapes(segments, gray_segments, x, y, color_map):
            shape_list = list()
            shape = dict(type='rect', xref='x', yref='y', fillcolor=None, line= {'width': 1},
                        x0=x, y0=y, x1=x+2/5, y1=y+0.1)
            shape_list.append(shape)                
            for i, segment in enumerate(segments):
                color=color_map[segment]
                shape = dict(type='rect', xref='x', yref='y', fillcolor=color, line= {'width': 1},
                            x0=x+segment[0]*2/5, y0=y, x1=x+segment[1]*2/5, y1=y+0.1)
                shape_list.append(shape)
            for i, segment in enumerate(gray_segments):
                shape = dict(type='rect', xref='x', yref='y', fillcolor='lightgray', line= {'width': 1},
                            x0=x+segment[0]*2/5, y0=y, x1=x+segment[1]*2/5, y1=y+0.1)
                shape_list.append(shape)
            return shape_list

        if type(focus_node) is Leaf:

            colors, _ = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Rainbow)
            colorscale = plotly.colors.make_colorscale(colors)
            color = get_continuous_color(colorscale, intermed=0)

            shape_list = [dict(type='rect', xref='x', yref='y', fillcolor=color, line= {'width': 1},
                x0=1.5/5, y0=0.25, x1=3.5/5, y1=0.35)]
            for segment in gray_segments:
                shape = dict(type='rect', xref='x', yref='y', fillcolor='lightgray', line= {'width': 1},
                            x0=1.5/5+segment[0]*2/5, y0=0.25, x1=1.5/5+segment[1]*2/5, y1=0.35)
                shape_list.append(shape)                

        elif type(focus_node) is Recombination:
            # print("###", focus_node.left_parent.intervals, focus_node.right_parent.intervals)
            segments1 = get_segments(focus_node, focus_node.left_parent.intervals)
            segments2 = get_segments(focus_node, focus_node.right_parent.intervals)
            segments3 = get_segments(focus_node, focus_node.child.intervals)

            # get part of intervals that intersect slider interval:
            segments1 = list(map(tuple, interval_intersect([slider_interval], segments1)))
            segments2 = list(map(tuple, interval_intersect([slider_interval], segments2)))
            segments3 = list(map(tuple, interval_intersect([slider_interval], segments3)))

            unique_segments = sorted(set().union(segments1, segments2, segments3))
            color_map = dict()
            colors, _ = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Rainbow)
            colorscale = plotly.colors.make_colorscale(colors)
            for i, s in enumerate(unique_segments):
                color_map[s] = get_continuous_color(colorscale, intermed=i/len(unique_segments))

            shape_list = \
                get_shapes(segments1, gray_segments, x=0, y=0.75, color_map=color_map) + \
                get_shapes(segments2, gray_segments, x=3/5, y=0.75, color_map=color_map) + \
                get_shapes(segments3, gray_segments, x=1.5/5, y=0.25, color_map=color_map) + \
                   [dict(type='line', xref='x', yref='y', line= {'width': 2, 'color': 'gray'},
                            x0=0.5, y0=0.55, x1=0.5, y1=0.35),
                    dict(type='line', xref='x', yref='y', line= {'width': 2, 'color': 'gray'},
                            x0=0.5, y0=0.55, x1=1/5, y1=0.75),                            
                    dict(type='line', xref='x', yref='y', line= {'width': 2, 'color': 'gray'},
                            x0=0.5, y0=0.55, x1=4/5, y1=0.75)]

            # print("slider", slider_interval)
            # shape_list.append(dict(type='rect', xref='x', yref='y', fillcolor='grey', line= {'width': 1},
            #             x0=slider_interval[0], y0=0.5, x1=slider_interval[1], y1=0.5+0.1))
                
        else:
            segments1 = get_segments(focus_node, focus_node.children[0].intervals)
            segments2 = get_segments(focus_node, focus_node.children[1].intervals)
            segments3 = get_segments(focus_node, focus_node.parent.intervals)

            # get part of intervals that intersect slider interval:
            segments1 = list(map(tuple, interval_intersect([slider_interval], segments1)))
            segments2 = list(map(tuple, interval_intersect([slider_interval], segments2)))
            segments3 = list(map(tuple, interval_intersect([slider_interval], segments3)))

            unique_segments = sorted(set().union(segments1, segments2, segments3))
            color_map = dict()
            colors, _ = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Rainbow)
            colorscale = plotly.colors.make_colorscale(colors)
            for i, s in enumerate(unique_segments):
                color_map[s] = get_continuous_color(colorscale, intermed=i/len(unique_segments))

            shape_list = \
                get_shapes(segments1, gray_segments, x=0, y=0.25, color_map=color_map) + \
                get_shapes(segments2, gray_segments, x=3/5, y=0.25, color_map=color_map) + \
                get_shapes(segments3, gray_segments, x=1.5/5, y=0.75, color_map=color_map) + \
                    [dict(type='line', xref='x', yref='y', line= {'width': 2, 'color': 'gray'},
                            x0=0.5, y0=0.55, x1=0.5, y1=0.75),
                    dict(type='line', xref='x', yref='y', line= {'width': 2, 'color': 'gray'},
                            x0=0.5, y0=0.55, x1=1/5, y1=0.35),                            
                    dict(type='line', xref='x', yref='y', line= {'width': 2, 'color': 'gray'},
                            x0=0.5, y0=0.55, x1=4/5, y1=0.35)]

    figure_data = dict(
                    data=traces,
                    layout=dict(xaxis=dict(fixedrange=True,
                                           range=[-0.01, 1.01], #title='Samples',
                                           showgrid=False, showline=False, 
                                           zeroline=False, showticklabels=False
                                           ),
                               yaxis=dict(fixedrange=True,
                                           range=[0, 1], #title='Time',
                                           showgrid=False, showline=False, 
                                           zeroline=False, showticklabels=False
                                           ),
                               margin= {'l': 0, 'b': 0, 't': 20, 'r': 0},
                               transition = {'duration': 0},
                               showlegend=False,
                               shapes=shape_list,
                               )
                    )

#        figure_data['layout']['shapes'].extend(shape_list)

    return figure_data
    


app.layout = layout

if __name__ == '__main__':
    app.run_server(debug=True)