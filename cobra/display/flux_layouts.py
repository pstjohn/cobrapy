"""
A D3-based network layout engine for plotting model solutions in the jupyter window
"""

import os
import json
import itertools

from jinja2 import Environment, FileSystemLoader
from IPython.display import HTML

from ..io.json import _to_dict

env = Environment(loader=FileSystemLoader(
    os.path.join(os.path.dirname(__file__), 'templates')))


def flux_map(cobra_model, included_metabolites=None,
             excluded_metabolites=None, excluded_reactions=None,
             excluded_compartments=None):
    """Create a flux map representation of the cobra.Model, including or
    excluding the given metabolites, reactions, and/or compartments 
    
    """
    # build cofactor metabolites from strings
    cobra_metabolites = []
    if excluded_metabolites:
        compartments = cobra_model.get_compartments()
        metabolite_list = [
            cf + '_' + co for cf, co in itertools.product(
                excluded_metabolites, compartments)] + excluded_metabolites
        for cofactor in metabolite_list:
            try: 
                cobra_metabolites += [
                    cobra_model.metabolites.get_by_id(cofactor)]
            except KeyError: pass
            # what if its already a cobra metabolite?

    cobra_rxns = []
    if excluded_reactions:
        for rxnid in excluded_reactions:
            try: 
                cobra_rxns += [
                    cobra_model.reactions.get_by_id(rxnid)]
            except KeyError: pass


    # Exclude metabolites and reactions in the given comparment
    excluded_metabolites = set(cobra_metabolites)
    excluded_reactions = set(cobra_rxns)

    if excluded_compartments:
        met_compartments = set((m for m in cobra_model.metabolites.iquery(
            lambda x: set(x.compartment).intersection(
                set(excluded_compartments)))))
        excluded_metabolites |= met_compartments

        # Do I want to redo this not to include excluded metabolites? 
        rxn_compartments = set((r for r in cobra_model.reactions.iquery(
            lambda x: set(x.compartments).intersection(
                set(excluded_compartments)))))
        excluded_reactions |= rxn_compartments

    for reaction in excluded_reactions:
        reaction.notes['map_info'] = {'hidden' : True}

    for metabolite in excluded_metabolites:
        metabolite.notes['map_info'] = {'hidden' : True}

    return render_model(cobra_model)







    



def create_model_json(cobra_model):
    """ Convert a cobra.Model object to a json string for d3. Adds flux
    information if the model has been solved
    
    """
    # Add flux info
    for reaction in cobra_model.reactions:
        try:
            reaction.notes['map_info']['flux'] = reaction.x

        except KeyError:
            # Create a new "map_info" object
            reaction.notes['map_info'] = {'flux' : reaction.x}

        except AttributeError:
            # Model likely hasn't been solved, get out now
            return json.dumps(_to_dict(cobra_model), allow_nan=False)

    for metabolite in cobra_model.metabolites:
        carried_flux = sum([abs(r.x * r.metabolites[metabolite]) for r in
                            metabolite.reactions])/2

        try:
            metabolite.notes['map_info']['flux'] = carried_flux

        except KeyError:
            # Create a new "map_info" object
            metabolite.notes['map_info'] = {'flux' : carried_flux}

    return json.dumps(_to_dict(cobra_model), allow_nan=False)



def render_model(cobra_model):
    """ Render a cobra.Model object in the current window """

    # Increment figure counter
    render_model._fignum += 1

    # Get figure name and JSON string for the cobra model
    figure_id = 'd3flux{:0>3d}'.format(render_model._fignum)
    modeljson = create_model_json(cobra_model)

    # Initialize the jinja templates
    template_css = env.get_template('network_style.css')
    template_html = env.get_template('output_template.html')
    template_js = env.get_template('d3flux.js')

    # Render the jinja templates with the given variables
    css = template_css.render()
    html = template_html.render(figure_id=figure_id)
    js = template_js.render(figure_id=figure_id, modeljson=modeljson)
    
    # Define wrappers to compile to single HTML string
    css_head = '<style type="text/css">\n'
    css_foot = '\n</style>'
    js_head = '<script type="text/Javascript">\n'
    js_foot = '\n</script>'

    # wrap css and js elements
    css = css_head + css + css_foot
    js = js_head + js + js_foot

    # compile and return HTML
    return HTML(css + html + js)

# Initialize figure counter
render_model._fignum = 0


common_cofactors = cofactors = ['coa', 'nadh', 'nad', 'nadph', 'nadp', 'atp',
                                'adp', 'amp', 'q8', 'q8h2', 'pi', 'co2', 'h2o',
                                'h', 'o2', 'h2', 'nh4']
