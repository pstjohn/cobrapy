"""
Functions to automate common analysis operations on elementary flux mode
calculations. Will allow the ElementaryFluxModes.py file to contain only the
necessary functions for the wrapping of efmtool, and likewise for MinimalCutSets.py.
"""

import numpy as np
import pandas as pd
from ..core import Reaction

def boundary_efms(cobra_model, efms, tol=1E-4):
    """Shortcut method to condense elementary flux modes to a unique set of
    input-output modes.

    cobra_model: a cobra.Model object

    efms: a pandas.DataFrame 
        represents the elementary flux modes of the system (output of
        calculate_elementary_modes)

    tol: float
        Tolerance for combining similar elementary flux modes

    """

    def scale_df(df):
        df_ss = np.sqrt((df**2).sum(1))
        to_drop = df_ss[df_ss == 0].index
        df_t = df.drop(to_drop)
        df_ss = df_ss.drop(to_drop)
        return df_t.divide(df_ss.values, axis='rows'), to_drop

    def remove_duplicates(df, tol=1E-4):
        df_round = df.fillna(0)
        df_round, to_drop = scale_df(df_round)
        df = df.drop(to_drop)
        df_round = np.round(df_round, decimals=int(-np.log10(1E-4)))
        return df.ix[~(df_round*(1/tol)).astype('int64').duplicated()]
    
    # Get the boundary reactions from the cobra.Model
    boundary_rxns = [r.id for r in cobra_model.reactions.query(
        'system_boundary', 'boundary') if r.id in efms.columns]

    ioefms = efms[boundary_rxns]

    return remove_duplicates(ioefms)

def check_mfa_results(model, cutsets, targets):
    """Run metabolic flux analysis for each of the calculated cutsets,
    returning optimal growth rates and fluxes through target reactions.

    model: a cobra.Model object
    cutsets: a pandas.DataFrame containing the results of `calculate_minimum_cut_sets`
    targets: a list of reaction strings indicating the desired results of MCS

    """
    sorted_cutsets = cutsets.iloc[cutsets.sum(1).argsort()]
    rxn_kos = sorted_cutsets.T.apply(lambda x: list(x.loc[x].index))

    target_dict = {key : lambda model, key=key: model.reactions.get_by_id(key).x 
                   for key in targets}
    target_dict.update({'growth' : lambda model: model.objective.keys()[0].x})

    def find_growth_rates():
        for cutset, knockouts in rxn_kos.iteritems():
            ko_model = model.copy()
            for r in knockouts:
                ko_model.reactions.get_by_id(r).knock_out()
            try: 
                s = ko_model.optimize()
                assert np.isfinite(s.f)
            except Exception:
                yield (cutset, {key : np.NaN for key, val in
                                target_dict.iteritems()})
                continue

            yield (cutset, {key : val(ko_model) for key, val in
                            target_dict.iteritems()})

    out = pd.DataFrame(rxn_kos, columns=['Knockouts']).join(
        pd.DataFrame(dict(find_growth_rates())).T)

    return out[out['growth'] > 1E-6]

def format_reaction(cobra_model, efm, tol=1E-8):
    """Build a cobra.Reaction representation of an elementary mode.

    cobra_model: a cobra.Model object from which the efm was calculated
    
    efm: a pandas.Series object representing the desired elementary mode
    
    """
    
    efm = np.round(efm, 2)
    
    efm_rxn = Reaction(efm.name)
        
    for rxn_id, flux in efm[efm != 0].iteritems():
        model_rxn = cobra_model.reactions.get_by_id(rxn_id)
        efm_rxn.add_metabolites((model_rxn * flux).metabolites)
            
    efm_rxn._metabolites = {
        metabolite : stoich for metabolite, stoich in
        efm_rxn.metabolites.iteritems() if abs(stoich) >= tol}
    
    return efm_rxn

def get_maximum_reactions(first, data, cobra_model, second=None, exclude=None, tol=1E-8):
    """Print reactions representing the maximum production of a given metabolite.

    first: string
        A string representing the desired metabolite, with the '_' compartment
        suffix omitted.

    data: a pandas.DataFrame
        Dataframe representing the EFMs to search.

    cobra_model: the cobra.Model object from which to build the reaction

    second: string or None
        If present, select from the best modes of 'first' that are also the
        maximum for 'second'

    exclude: list or None
        A list of metabolite strings. If present, EFMs containing this
        metabolite will be omitted from the maximum calculation

    tol: float
        Value below which to consider an EFM flux inactive.
        
    """
    
    first_str = data.columns[data.columns.str.contains(first + '_')][0]
    if exclude: 
        exclude_list = [
            data.columns[data.columns.str.contains(e + '_')][0]
            for e in exclude]
    else: exclude_list = []
    
    efmlist = data[data[first_str] <= (data[first_str].min() * (1 - tol))]
    efmlist = efmlist.iloc[(efmlist == 0).sum(1).argsort()[::-1]]
    efmlist = efmlist.loc[(efmlist[exclude_list] == 0).all(1)]
    
    if second:
        second_str = data.columns[data.columns.str.contains(second + '_')][0]
        efmlist = efmlist[efmlist[second_str] <= 
                          (efmlist[second_str].min() * (1 - tol))]
    
    for id_, efm in efmlist.iterrows():
        print format_reaction(cobra_model, efm).reaction



