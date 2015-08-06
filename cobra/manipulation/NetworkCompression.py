"""
Preforms network simplification as outlined in Quek, L.-E., & Nielsen, L. K.
(2014). A depth-first search algorithm to compute elementary flux modes by
linear programming. BMC Systems Biology, 8(1), 94.
http://doi.org/10.1186/s12918-014-0094-2.

Uses a pandas.DataFrame to describe to describe and manipulate model stochiometry
"""
import numpy as np
import pandas as pd

from ..core import DataframeBasedModel, Model
from .modify import convert_to_irreversible, revert_to_reversible

def condense_model(cobra_model, deepcopy_model=True):
    """Condense a cobra.Model object into a minimum set of reactants and
    metabolites that retains the fundimental complexity of the full model
    
    deepcopy_model: bool
        Whether or not to create a deep copy of the input model
    """
    
    input_model = cobra_model.copy()

    # Known bug in compression: don't pass reactions with 0 bounds.
    for reaction in input_model.reactions.iquery(lambda x: x.bounds == (0,0)):
        reaction.remove_from_model()

    convert_to_irreversible(input_model)
    df_model = DataframeBasedModel(input_model)
    
    # Create a set of metabolites not to compress for model comprehension
    # purposes
    ignored_metabolites = set((r.metabolites.keys()[0].id for r in
                               df_model.reactions.iquery('system_boundary',
                                                      'boundary')))

    def objective_precurors():
        for reaction in df_model.objective.keys():
            for metabolite in reaction.metabolites.keys():
                yield metabolite.id

    ignored_metabolites.union(set(objective_precurors()))

    # Run the compression algorithm
    compressed_df = condense(df_model.S, ignored_metabolites)
    df_model.S = compressed_df

    # Revert to a reversible model
    revert_to_reversible(df_model, update_solution=False, allow_missing=True)

    # Revert to a standard cobrapy model
    df_model.__class__ = Model

    return df_model



def condense(df, ignored_metabolites=None, verbose=True):
    """Iterate over condensing algorithms until the dataframe size remains constant

    ignored_metabolites: set
        A set of metabolites to ignore, typically those exchanged with the
        external environment or used in the objective function.

    verbose: bool
        whether or not network size (n_met, r_rxns) and called function should
        be called on each iteration. (These functions aren't optimized --
        recommended)

    TODO:
        * keep a list of merged/deleted reactions. Use this list to
            * propagate gene reaction rules 
            * update solution dictionaries
            * decompress model
            * check and update reaction bounds - errors pop up if knocked-out
            reactions are passed
    """

    func_order = [trim, eliminate, lump, match]

    while True:
        df_size = np.array(df.shape)
        for func in func_order:
            if func == lump:
                out = func(df, ignored_metabolites)
            else:
                out = func(df)
            if type(out) is not bool:
                df = out
                break
        new_df_size = np.array(df.shape)
        diff_df_size = new_df_size - df_size
        if not np.any(diff_df_size):
            return df
        elif verbose: 
            print "{0:>12}: {1:>12}".format(func.func_name, new_df_size)





def trim(df):
    """Remove all metabolites which are only involved in less than two
    reactions, as they cannot support flux. Also remove all reactions
    associated with these metabolites, as they are necessarily flux=0.

    """
    unused_metabolites = df.loc[
        df.apply(lambda x: len(x[np.abs(x) > 0]), axis=1) < 2].index

    unused_reactions = df.loc[:, df.loc[unused_metabolites].any(0)].columns

    if len(unused_metabolites | unused_reactions) > 0:
        # print unused_metabolites
        return df.drop(unused_metabolites, axis=0).drop(unused_reactions,
                                                        axis=1)
    else: return False


def eliminate(df):
    """Remove all reactions which do not have any metabolites

    """
    empty_reactions = df.loc[:, (df == 0).all(0)].columns
    if len(empty_reactions) > 0: 
        # print empty_reactions
        return df.drop(empty_reactions, axis=1)
    else: return False
    


def match(df):
    """Given a stochiometric dataframe, find reactions (columns) which are
    isoenzymes. Scales each reaction to ensure matches occur regardless of scaling.

    """

    duplicate_reactions = list(df.loc[:,df.apply(lambda x: x/np.abs(x).sum(),
                                                 axis=0).T.duplicated()].columns)

    if len(duplicate_reactions) > 0: 
        return df.drop(duplicate_reactions, axis=1)
    
    else: return False


def lump(df, ignored_metabolites=None):
    """Finds metabolites which are produced or consumed by a single reaction,
    and updates the dataframe so that reactions bypass the given metabolite.
    
    ignored_metabolites: set
        A set of metabolites to ignore, typically those exchanged with the
        external environment.
    """

    if ignored_metabolites is None:
        ignored_metabolites = set()
    
    # A list of metabolites which are consumed by only one reaction
    consumed_metabolites = df.index[
        df.apply(lambda x: len(x[x < 0]), axis=1) == 1].difference(
            ignored_metabolites)
    if len(consumed_metabolites) > 0:
        return skip_metabolite(df, consumed_metabolites[0],
                               direction='consuming')

    else:
        produced_metabolites = df.index[
            df.apply(lambda x: len(x[x > 0]), axis=1) == 1].difference(
                ignored_metabolites)
        if len(produced_metabolites) > 0:
            return skip_metabolite(df, produced_metabolites[0],
                                   direction='producing')

        else: return False


def skip_metabolite(df, metabolite, direction='consuming'):
    """ Update df to direct all fluxes into metabolite to the next node in the
    network.  metabolite should be consumed by only one reaction. For use with
    the lump routine.

    direction: 'consuming' or 'producing'
        whether or not the metabolite is consumed or produced by a single
        reaction

    """

    if direction == 'consuming':
        # get reactions which produce the metabolite
        multiple_reactions = {rxn : stoch for rxn, stoch in
                               df.loc[metabolite, df.loc[metabolite] >
                                      0].iteritems()}

        # Get information on the single reaction which consumes the metabolite
        sole_reaction = df.loc[
            metabolite, df.loc[metabolite] < 0]
        sole_stoch = -sole_reaction.values[0]

    elif direction == 'producing':
        # get reactions which produce the metabolite
        multiple_reactions = {rxn : stoch for rxn, stoch in
                               df.loc[metabolite, df.loc[metabolite] <
                                      0].iteritems()}

        # Get information on the single reaction which consumes the metabolite
        sole_reaction = df.loc[
            metabolite, df.loc[metabolite] > 0]
        sole_stoch = -sole_reaction.values[0]

    else: raise TypeError(
        "{} is not in ('producing', 'consuming')".format(direction))

    assert len(sole_reaction) == 1, "{0} produced by {1}".format(
        metabolite, ', '.join(sole_reaction.index))

    assert len(multiple_reactions) > 0, "Dead-end reaction, aborting"

    if len(multiple_reactions) > 1:
        # Add sole reaction to the multiple reaction, scaled such that the
        # metabolite is removed.
        for reaction, prod_stoch in multiple_reactions.iteritems():
            df[reaction] += (prod_stoch/sole_stoch) * df[sole_reaction.index[0]]

        # Single sole reaction is no longer needed and can be removed.
        # print sole_reaction, metabolite
        return df.drop(sole_reaction.index[0], axis=1).drop(metabolite, axis=0)

    else:
        # Here either the consuming reaction can be merged into the producing,
        # or vice versa. We typically want to keep the reaction which is
        # already larger, as this preserves larger Biomass reactions
        rxn1 = multiple_reactions.keys()[0]
        stoch1 = multiple_reactions.values()[0]
        rxn2 = sole_reaction.index[0]
        stoch2 = sole_reaction.values[0]

        rxn1_size = (df[rxn1] != 0).sum()
        rxn2_size = (df[rxn2] != 0).sum()

        if rxn1_size > rxn2_size:
            df[rxn1] += -(stoch1/stoch2) * df[rxn2]
            return df.drop(rxn2, axis=1).drop(metabolite, axis=0)
        else:
            df[rxn2] += -(stoch2/stoch1) * df[rxn1]
            return df.drop(rxn1, axis=1).drop(metabolite, axis=0)




