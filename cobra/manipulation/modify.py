from copy import deepcopy
from warnings import warn
from itertools import chain
from ast import NodeTransformer

from contextlib import contextmanager

from six import iteritems

from .. import Reaction, Metabolite
from .delete import get_compiled_gene_reaction_rules
from ..core.Gene import ast2str
from ..io.sbml3 import _renames


def _escape_str_id(id_str):
    """make a single string id SBML compliant"""
    for c in ("'", '"'):
        if id_str.startswith(c) and id_str.endswith(c) \
                and id_str.count(c) == 2:
            id_str = id_str.strip(c)
    for char, escaped_char in _renames:
        id_str = id_str.replace(char, escaped_char)
    return id_str


class _GeneRenamer(NodeTransformer):

    def visit_Name(self, node):
        node.id = _escape_str_id(node.id)
        return node


def escape_ID(cobra_model):
    """makes all ids SBML compliant"""
    for x in chain([cobra_model],
                   cobra_model.metabolites,
                   cobra_model.reactions,
                   cobra_model.genes):
        x.id = _escape_str_id(x.id)
    cobra_model.repair()
    gene_renamer = _GeneRenamer()
    for rxn, rule in iteritems(get_compiled_gene_reaction_rules(cobra_model)):
        if rule is not None:
            rxn._gene_reaction_rule = ast2str(gene_renamer.visit(rule))


def get_growth_medium(cobra_model):
    """Searches a cobra model for the currently active exchange reactions which
    allow metabolites to be created.

    Parameters
    ----------
    cobra_model : cobra.Model


    Returns
    -------
    medium : dict
        A dictionary containing pairs of reaction_id's : lower_bounds for the
        exchange reactions present in the model

    """
    exchange_reactions = cobra_model.reactions.query("system_boundary", 'boundary')
    
    def is_active(reaction):
        """ Test whether the reaction's bounds allow feasible flux to create
        metabolite
        """
        if reaction.reactants:
            if reaction.lower_bound < 0: return True
            else: return False

        elif reaction.products:
            raise RuntimeWarning(
                'Reaction {} has unconventional direction'.format(
                    reaction.id))
            if reaction.upper_bound > 0: return True
            else: return False

        else: raise RuntimeError('Empty Reaction: {}'.format(reaction.id))

    medium = {r.id : r.lower_bound for r in exchange_reactions if
              is_active(r)}

    return medium



def initialize_growth_medium(cobra_model, the_medium='MgM',
                             external_boundary_compartment='e',
                             external_boundary_reactions=None,
                             reaction_lower_bound=0.,
                             reaction_upper_bound=1000.,
                             irreversible=False,
                             reactions_to_disable=None):
    """Sets all of the input fluxes to the model to zero and then will
    initialize the input fluxes to the values specified in the_medium if
    it is a dict or will see if the model has a composition dict and use
    that to do the initialization.

    cobra_model: A cobra.Model object.


    the_medium: A string, or a dictionary.
    If a string then the initialize_growth_medium function expects that
    the_model has an attribute dictionary called media_compositions, which is a
    dictionary of dictionaries for various medium compositions.  Where a medium
    composition is a dictionary of external boundary reaction ids for the
    medium components and the external boundary fluxes for each medium
    component.


    external_boundary_compartment: None or a string.
    If not None then it specifies the compartment in which to disable all of
    the external systems boundaries.

    external_boundary_reactions: None or a list of external_boundaries that are
    to have their bounds reset.  This acts in conjunction with
    external_boundary_compartment.


    reaction_lower_bound: Float.  The default value to use for the lower
    bound for the boundary reactions.

    reaction_upper_bound: Float.  The default value to use for the upper
    bound for the boundary.

    irreversible: Boolean.  If the model is irreversible then the medium
    composition is taken as the upper bound

    reactions_to_disable: List of reactions for which the upper and lower
    bounds are disabled.  This is superceded by the contents of
    media_composition

    """
    # Zero all of the inputs to the model
    if hasattr(the_medium, 'keys'):
        medium_composition = the_medium
    else:
        if hasattr(cobra_model, 'media_compositions'):
            if the_medium in cobra_model.media_compositions:
                medium_composition = cobra_model.media_compositions[the_medium]
            else:
                raise Exception("%s is not in the model's media list" %
                                the_medium)
        else:
            raise Exception("the model doesn't have attribute "
                            "media_compositions and the medium is not a dict")
    if external_boundary_reactions is not None:
        if isinstance(external_boundary_reactions[0], str):
            external_boundary_reactions = map(cobra_model.reactions.get_by_id,
                                              external_boundary_reactions)
    elif external_boundary_compartment is None:
            warn("We are initializing the medium without first adjusting all"
                 "external boundary reactions")

    # Select the system_boundary reactions to reset
    if external_boundary_compartment is not None:
        _system_boundaries = dict([(x, x.get_compartments())
                                   for x in cobra_model.reactions
                                   if x.boundary == 'system_boundary'])
        [_system_boundaries.pop(k) for k, v in list(_system_boundaries.items())
         if len(v) == 1 and external_boundary_compartment not in v]
        if external_boundary_reactions is None:
            external_boundary_reactions = _system_boundaries.keys()
        else:
            external_boundary_reactions += _system_boundaries.keys()

    for the_reaction in external_boundary_reactions:
        the_reaction.lower_bound = reaction_lower_bound
        if the_reaction.upper_bound == 0:
            the_reaction.upper_bound = reaction_upper_bound
    # Disable specified reactions
    if reactions_to_disable is not None:
        if isinstance(reactions_to_disable[0], str):
            reactions_to_disable = map(cobra_model.reactions.get_by_id,
                                       reactions_to_disable)
        for the_reaction in reactions_to_disable:
            the_reaction.lower_bound = the_reaction.upper_bound = 0.

    # Update the model inputs based on the_medium
    for the_component in medium_composition.keys():
        the_reaction = cobra_model.reactions.get_by_id(the_component)
        if irreversible:
            the_reaction.upper_bound = medium_composition[the_component]
        else:
            the_reaction.lower_bound = medium_composition[the_component]


def convert_to_irreversible(cobra_model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    cobra_model: A Model object which will be modified in place.

    """
    reactions_to_add = []
    for reaction in cobra_model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            reverse_reaction.objective_coefficient = \
                reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in iteritems(reaction._metabolites)}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    cobra_model.add_reactions(reactions_to_add)


def revert_to_reversible(cobra_model, update_solution=True, allow_missing=False):
    """This function will convert a reversible model made by
    convert_to_irreversible into a reversible model.

    cobra_model: A cobra.Model which will be modified in place.

    """
    reverse_reactions = [x for x in cobra_model.reactions
                         if "reflection" in x.notes and
                         x.id.endswith('_reverse')]

    # If there are no reverse reactions, then there is nothing to do
    if len(reverse_reactions) == 0:
        return

    update_solution = update_solution and cobra_model.solution is not None \
        and cobra_model.solution.status not in ["NA", "infeasible"]

    if update_solution:
        x_dict = cobra_model.solution.x_dict

    missing_reactions = []
    for reverse in reverse_reactions:
        forward_id = reverse.notes.pop("reflection")

        try:
            forward = cobra_model.reactions.get_by_id(forward_id)
            forward.lower_bound = -reverse.upper_bound

            # update the solution dict
            if update_solution:
                if reverse.id in x_dict:
                    x_dict[forward_id] -= x_dict.pop(reverse.id)

            if "reflection" in forward.notes:
                forward.notes.pop("reflection")

        except KeyError:
            if not allow_missing: 
                raise KeyError(
                    "Forward reaction {} not found in model".format(
                        forward_id))
            else:
                # Remove 'reflection' tag, as forward reaction is no longer in
                # the model
                reverse_reaction.notes.pop("reflection")
                
                # Add the reaction to a list of reactions not to remove.
                missing_reactions += [reverse]


    # Since the metabolites and genes are all still in
    # use we can do this faster removal step.  We can
    # probably speed things up here.
    cobra_model.remove_reactions(
        set(reverse_reactions).difference(set(missing_reactions)))

    # update the solution vector
    if update_solution:
        cobra_model.solution.x_dict = x_dict
        cobra_model.solution.x = [x_dict[r.id] for r in cobra_model.reactions]

@contextmanager
def knocked_out(model, ko_list):
    """Convenience function to temporarily knockout reactions.

    model: a cobra.Model model
        The original, wild-type model

    ko_list: a list of cobra.Reactions
        Reactions to knock out

    This function will temporary knockout reactions, and reset them to their
    original bounds using context management. For example:

    >>> with knocked_out(wt_model, [rxn1, rxn2]) as ko_model:
    >>>     my_knockout_test(ko_model)
    >>>
    >>> my_wt_test(wt_model)

    will knockout reactions rxn1 and rxn2 only within the scope of the with
    statement.

    """
    # Correct non-list input
    if type(ko_list) is not list: ko_list = [ko_list]

    # Store original bounds for later use
    try:
        wt_bounds = {rxn : model.reactions.get_by_id(rxn.id).bounds 
                     for rxn in ko_list}
    except AttributeError:
        wt_bounds = {rxn : model.reactions.get_by_id(rxn).bounds 
                     for rxn in ko_list}


    # Knockout reactions (set bounds to zero)
    for rxn in ko_list: model.reactions.get_by_id(rxn).knock_out()

    # Yield statement for context manager, with-code gets executed here
    yield model

    # Reset reaction bounds to their original state
    for rxn in ko_list: model.reactions.get_by_id(rxn).bounds = wt_bounds[rxn]

def canonical_form(model, objective_sense='maximize',
                   already_irreversible=False, copy=True):
    """Return a model (problem in canonical_form).

    Converts a minimization problem to a maximization, makes all variables
    positive by making reactions irreversible, and converts all constraints to
    <= constraints.


    model: class:`~cobra.core.Model`. The model/problem to convert.

    objective_sense: str. The objective sense of the starting problem, either
    'maximize' or 'minimize'. A minimization problems will be converted to a
    maximization.

    already_irreversible: bool. If the model is already irreversible, then pass
    True.

    copy: bool. Copy the model before making any modifications.

    """
    if copy:
        model = model.copy()

    if not already_irreversible:
        convert_to_irreversible(model)

    if objective_sense == "minimize":
        # if converting min to max, reverse all the objective coefficients
        for reaction in model.reactions:
            reaction.objective_coefficient = - reaction.objective_coefficient
    elif objective_sense != "maximize":
        raise Exception("Invalid objective sense '%s'. "
                        "Must be 'minimize' or 'maximize'." % objective_sense)

    # convert G and E constraints to L constraints
    for metabolite in model.metabolites:
        if metabolite._constraint_sense == "G":
            metabolite._constraint_sense = "L"
            metabolite._bound = - metabolite._bound
            for reaction in metabolite.reactions:
                coeff = reaction.get_coefficient(metabolite)
                # reverse the coefficient
                reaction.add_metabolites({metabolite: -2 * coeff})
        elif metabolite._constraint_sense == "E":
            # change existing constraint to L
            metabolite._constraint_sense = "L"
            # add new constraint
            new_constr = Metabolite("%s__GE_constraint" % metabolite.id)
            new_constr._constraint_sense = "L"
            new_constr._bound = - metabolite._bound
            for reaction in metabolite.reactions:
                coeff = reaction.get_coefficient(metabolite)
                reaction.add_metabolites({new_constr: -coeff})

    # convert lower bounds to LE constraints
    for reaction in model.reactions:
        if reaction.lower_bound < 0:
            raise Exception("Bounds of irreversible reactions should be >= 0,"
                            " for %s" % reaction.id)
        elif reaction.lower_bound == 0:
            continue
        # new constraint for lower bound
        lb_constr = Metabolite("%s__LB_constraint" % reaction.id)
        lb_constr._constraint_sense = "L"
        lb_constr._bound = - reaction.lower_bound
        reaction.add_metabolites({lb_constr: -1})
        reaction.lower_bound = 0

    return model
