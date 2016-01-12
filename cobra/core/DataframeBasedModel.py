import pandas as pd

from .Model import Model

class DataframeBasedModel(Model):
    """DataframeBasedModel adds pandas dataframes to a cobra.Model in order to
    enable metabolic network compression.

    """

    def __init__(self, description=None, deepcopy_model=False):
        """
        description: None | String | cobra.Model

        deepcopy_model: Boolean.  If True and description is
        a cobra.Model then make a deepcopy of the Model before
        creating the ArrayBasedModel.

        """

        if deepcopy_model and isinstance(description, Model):
            description = description.copy()
        Model.__init__(self, description)
        self._S = None
        self.update()

    @property
    def S(self):
        """Stoichiometric dataframe of the model

        """
        return self._S

    @S.setter
    def S(self, df):
        """Update the model's reactions and metabolites according to the new
        dataframe. Remove unused metabolites and reactions not found in the new
        dataframe. Raise an error if new metabolite or reaction IDs are seen.

        df: pandas.DataFrame
            New stochiometric dataframe
        """

        metabolite_set = frozenset((m.id for m in self.metabolites))
        reaction_set = frozenset((r.id for r in self.reactions))

        # Assert no new reactions or metabolites are being added.
        new_reactions = set(df.columns).difference(reaction_set)
        new_metabolites = set(df.index).difference(metabolite_set)

        assert not new_reactions, "{} not found in model".format(
            ', '.join(new_reactions))
        assert not new_metabolites, "{} not found in model".format(
            ', '.join(new_metabolites))

        # Get a list of metabolites and reactions to remove, remove them
        missing_metabolites = metabolite_set.difference(set(df.index))
        missing_reactions = reaction_set.difference(set(df.columns))
        self.remove_metabolites(missing_metabolites)
        self.remove_reactions(missing_reactions)
        
        # Iterate over each reaction in the network, updating the metabolite
        # dictionaries
        for rxn_id, col in df.iteritems():
            column_stochiometry = col[col != 0]
            reaction = self.reactions.get_by_id(rxn_id)

            # Clear the existing metabolite dictionary
            reaction.subtract_metabolites(reaction.metabolites)

            # Create a new dictionary from the dataframe
            new_metabolite_dict = {
                self.metabolites.get_by_id(met_id) : stochiometry for met_id,
                stochiometry in column_stochiometry.iteritems()}

            # Add new dictionary to the reaction
            reaction.add_metabolites(new_metabolite_dict)

        self._S = df
        self.repair()

    def remove_metabolites(self, metabolites, method='subtractive'):
        """ Remove given metabolites from the model
        
        metabolites: list
            which metabolites to remove
        
        method: 'subtractive' or 'destructive'
            see cobra.Metabolite.remove_from_model()

        """
        for metabolite in metabolites:
            try:
                self.metabolites.get_by_id(
                    metabolite).remove_from_model(method=method)
            except KeyError:
                self.metabolites.get_by_id(
                    metabolite.id).remove_from_model(method=method)

    def update(self):
        """ Regenerate the stochiometric dataframe. Not efficient! Shouldn't
        get placed inside a loop
        
        """

        mlabels = (m.id for m in self.metabolites)
        rlabels = (r.id for r in self.reactions)
        self._S = pd.DataFrame(index=mlabels, columns=rlabels)

        def metabolite_dictionaries():
            """ Iterator to yeild reaction, metabolite pairs with string IDs
            """
            for reaction in self.reactions:
                metabolites = {
                    m.id : stoch for m, stoch in
                    reaction.metabolites.iteritems()}
                yield reaction.id, metabolites

        self._S.update(pd.DataFrame(
            {key : val for key, val in metabolite_dictionaries()}
        ))

        self._S.replace(pd.np.NaN, 0, inplace=True)


    @property
    def lower_bounds(self):
        return pd.Series([r.lower_bound for r in self.reactions],
                         index=[r.id for r in self.reactions])

    @property
    def upper_bounds(self):
        return pd.Series([r.upper_bound for r in self.reactions],
                         index=[r.id for r in self.reactions])

    @property
    def fluxes(self):
        return pd.Series([r.x for r in self.reactions],
                         index=[r.id for r in self.reactions])

    @property
    def objectives(self):
        return pd.Series([r.objective_coefficient for r in self.reactions],
                         index=[r.id for r in self.reactions])


        
