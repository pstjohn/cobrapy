from warnings import warn
from copy import deepcopy, copy
import pandas as pd
import itertools

from six import iteritems, string_types

from ..solvers import optimize
from .Object import Object
from .Solution import Solution
from .Reaction import Reaction
from .DictList import DictList

# Note, when a reaction is added to the Model it will no longer keep personal
# instances of its Metabolites, it will reference Model.metabolites to improve
# performance.  When doing this, take care to monitor metabolite coefficients.
# Do the same for Model.reactions[:].genes and Model.genes

class Model(Object):
    """Metabolic Model

    Refers to Metabolite, Reaction, and Gene Objects.
    """

    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model"""
        self.__dict__.update(state)
        for y in ['reactions', 'genes', 'metabolites']:
            for x in getattr(self, y):
                x._model = self
        if not hasattr(self, "name"):
            self.name = None

    def __init__(self, id_or_model=None, name=None):
        if isinstance(id_or_model, Model):
            Object.__init__(self, name=name)
            self.__setstate__(id_or_model.__dict__)
            if not hasattr(self, "name"):
                self.name = None
        else:
            Object.__init__(self, id_or_model, name=name)
            self._trimmed = False
            self._trimmed_genes = []
            self._trimmed_reactions = {}
            self.genes = DictList()
            self.reactions = DictList()  # A list of cobra.Reactions
            self.metabolites = DictList()  # A list of cobra.Metabolites
            self.pathways = DictList()
            self.compartments = {} # For SBML input-output
            # genes based on their ids {Gene.id: Gene}
            self.solution = Solution(None)
            self.solution.model = self
            # self.media_compositions = {}

    @property
    def description(self):
        warn("description deprecated")
        return self.name if self.name is not None else ""

    @description.setter
    def description(self, value):
        self.name = value
        warn("description deprecated")

    def __add__(self, other_model):
        """Adds two models. +

        The issue of reactions being able to exists in multiple Models now
        arises, the same for metabolites and such.  This might be a little
        difficult as a reaction with the same name / id in two models might
        have different coefficients for their metabolites due to pH and whatnot
        making them different reactions.

        """
        new_model = self.copy()
        new_reactions = deepcopy(other_model.reactions)
        new_model.add_reactions(new_reactions)
        new_model.id = self.id + '_' + other_model.id
        return new_model

    def __iadd__(self, other_model):
        """Adds a Model to this model +=

        The issue of reactions being able to exists in multiple Models now
        arises, the same for metabolites and such.  This might be a little
        difficult as a reaction with the same name / id in two models might
        have different coefficients for their metabolites due to pH and whatnot
        making them different reactions.

        """
        new_reactions = deepcopy(other_model.reactions)
        self.add_reactions(new_reactions)
        self.id = self.id + '_' + other_model.id
        return self

    def copy(self):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite,
        Gene, and Reaction objects are created anew but in a faster fashion
        than deepcopy
        """
        new = self.__class__()
        do_not_copy = {"metabolites", "reactions", "genes"}
        for attr in self.__dict__:
            if attr not in do_not_copy:
                new.__dict__[attr] = self.__dict__[attr]

        new.metabolites = DictList()
        do_not_copy = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in iteritems(metabolite.__dict__):
                if attr not in do_not_copy:
                    new_met.__dict__[attr] = copy(
                        value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr not in do_not_copy:
                    new_gene.__dict__[attr] = copy(
                        value) if attr == "formula" else value
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in iteritems(reaction.__dict__):
                if attr not in do_not_copy:
                    new_reaction.__dict__[attr] = value
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in iteritems(reaction._metabolites):
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            for gene in reaction._genes:
                new_gene = new.genes.get_by_id(gene.id)
                new_reaction._genes.add(new_gene)
                new_gene._reaction.add(new_reaction)
        return new

    def add_metabolites(self, metabolite_list):
        """Will add a list of metabolites to the the object, if they do not
        exist and then expand the stochiometric matrix

        metabolite_list: A list of :class:`~cobra.core.Metabolite` objects

        """
        if not hasattr(metabolite_list, '__iter__'):
            metabolite_list = [metabolite_list]
        # First check whether the metabolites exist in the model
        metabolite_list = [x for x in metabolite_list
                           if x.id not in self.metabolites]
        for x in metabolite_list:
            x._model = self
        self.metabolites += metabolite_list

    def add_reaction(self, reaction):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction: A :class:`~cobra.core.Reaction` object

        """
        self.add_reactions([reaction])

    def add_reactions(self, reaction_list):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction_list: A list of :class:`~cobra.core.Reaction` objects

        """
        # Only add the reaction if one with the same ID is not already
        # present in the model.

        # This function really should not used for single reactions
        if not hasattr(reaction_list, "__len__"):
            reaction_list = [reaction_list]
            warn("Use add_reaction for single reactions")

        reaction_list = DictList(reaction_list)
        reactions_in_model = [
            i.id for i in reaction_list if self.reactions.has_id(
                i.id)]

        if len(reactions_in_model) > 0:
            raise Exception("Reactions already in the model: " +
                            ", ".join(reactions_in_model))

        # Add reactions. Also take care of genes and metabolites in the loop
        for reaction in reaction_list:
            reaction._model = self  # the reaction now points to the model
            # keys() is necessary because the dict will be modified during
            # the loop
            for metabolite in list(reaction._metabolites.keys()):
                # if the metabolite is not in the model, add it
                # should we be adding a copy instead.
                if not self.metabolites.has_id(metabolite.id):
                    self.metabolites.append(metabolite)
                    metabolite._model = self
                    # this should already be the case. Is it necessary?
                    metabolite._reaction = set([reaction])
                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(
                        metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)

            for gene in list(reaction._genes):
                # If the gene is not in the model, add it
                if not self.genes.has_id(gene.id):
                    self.genes.append(gene)
                    gene._model = self
                    # this should already be the case. Is it necessary?
                    gene._reaction = set([reaction])
                # Otherwise, make the gene point to the one in the model
                else:
                    model_gene = self.genes.get_by_id(gene.id)
                    if model_gene is not gene:
                        reaction._dissociate_gene(gene)
                        reaction._associate_gene(model_gene)

        self.reactions += reaction_list

    def to_array_based_model(self, deepcopy_model=False, **kwargs):
        """Makes a :class:`~cobra.core.ArrayBasedModel` from a cobra.Model which
        may be used to perform linear algebra operations with the
        stoichiomatric matrix.

        deepcopy_model: Boolean.  If False then the ArrayBasedModel points
        to the Model

        """
        from .ArrayBasedModel import ArrayBasedModel
        return ArrayBasedModel(self, deepcopy_model=deepcopy_model, **kwargs)

    def optimize(self, objective_sense='maximize',
                 minimize_absolute_flux=False,
                 **kwargs):

        r"""Optimize model using flux balance analysis

        objective_sense: 'maximize' or 'minimize'

        minimize_absolute_flux: False, or float
            Whether or not to invoke optimize_maf in order to generate
            parsimonious flux balance solutions. Float indicates the percent of
            the maximum objective function to maintain while minimizing overall
            flux.

        solver: 'glpk', 'cglpk', 'gurobi', 'cplex' or None

        quadratic_component: None or :class:`scipy.sparse.dok_matrix`
            The dimensions should be (n, n) where n is the number of reactions.

            This sets the quadratic component (Q) of the objective coefficient,
            adding :math:`\\frac{1}{2} v^T \cdot Q \cdot v` to the objective.

        tolerance_feasibility: Solver tolerance for feasibility.

        tolerance_markowitz: Solver threshold during pivot

        time_limit: Maximum solver time (in seconds)

        .. NOTE :: Only the most commonly used parameters are presented here.
                   Additional parameters for cobra.solvers may be available and
                   specified with the appropriate keyword argument.

        """
        if minimize_absolute_flux:
            from ..flux_analysis import optimize_maf
            solution = optimize_maf(self, float(minimize_absolute_flux),
                                        **kwargs)
        else: 
            solution = optimize(self, objective_sense=objective_sense, **kwargs)
        self.solution = solution
        return solution

    def remove_reactions(self, reactions, delete=True,
                         remove_orphans=False):
        """remove reactions from the model

        reactions: [:class:`~cobra.core.Reaction.Reaction`] or [str]
            The reactions (or their id's) to remove

        delete: Boolean
            Whether or not the reactions should be deleted after removal.
            If the reactions are not deleted, those objects will be
            recreated with new metabolite and gene objects.

        remove_orphans: Boolean
            Remove orphaned genes and metabolites from the model as well

        """
        if isinstance(reactions, string_types) or hasattr(reactions, "id"):
            warn("need to pass in a list")
            reactions = [reactions]
        for reaction in reactions:
            try:
                reaction = self.reactions[self.reactions.index(reaction)]
            except ValueError:
                warn('%s not in %s' % (reaction, self))
            else:
                if delete:
                    reaction.delete(remove_orphans=remove_orphans)
                else:
                    reaction.remove_from_model(remove_orphans=remove_orphans)

    def repair(self, rebuild_index=True, rebuild_relationships=True):
        """Update all indexes and pointers in a model"""

        # Assert all id's are unique
        metabolite_ids = set((m.id for m in self.metabolites))
        assert len(metabolite_ids) == len(self.metabolites), \
            "Duplicate found in metabolite IDs"

        if rebuild_index:  # DictList indexes
            self.reactions._generate_index()
            self.metabolites._generate_index()
            self.genes._generate_index()
        if rebuild_relationships:
            for met in self.metabolites:
                met._reaction.clear()
            for gene in self.genes:
                gene._reaction.clear()
            for rxn in self.reactions:
                for met in rxn._metabolites:
                    met._reaction.add(rxn)
                for gene in rxn._genes:
                    gene._reaction.add(rxn)
        # point _model to self
        for l in (self.reactions, self.genes, self.metabolites):
            for e in l:
                e._model = self
        if self.solution is None:
            self.solution = Solution(None)
            self.solution.model = self
        return

    def change_objective(self, objectives):
        """Change the model objective"""
        self.objective = objectives

    def get_compartments(self):
        return frozenset((x.compartment for x in self.metabolites))

    @property
    def objective(self):
        return {reaction: reaction.objective_coefficient
                for reaction in self.reactions
                if reaction.objective_coefficient != 0}

    @objective.setter
    def objective(self, objectives):
        # set all objective coefficients to 0 initially
        for x in self.reactions:
            x.objective_coefficient = 0.
        # case of a single reaction
        if isinstance(objectives, string_types) or \
                isinstance(objectives, Reaction):
            self.reactions.get_by_id(str(objectives)).objective_coefficient = 1
        elif isinstance(objectives, int):
            self.reactions[objectives].objective_coefficient = 1

        # case of an iterable
        else:
            for reaction_id in objectives:
                if isinstance(reaction_id, int):  # index in a list
                    reaction = self.reactions[reaction_id]
                else:
                    reaction = self.reactions.get_by_id(str(reaction_id))
                # objective coefficient obtained from a dict, and is 1. if
                # from a list.
                reaction.objective_coefficient = objectives[reaction_id] \
                    if hasattr(objectives, "items") else 1.


    def summary(self, tol=1E-8, round=2, fva=None):
        """Print a summary of the input and output fluxes of the model.
        
        tol: float
            tolerance for determining if a flux is zero (not printed)

        round: int
            number of digits after the decimal place to print

        fva: int or None
            Whether or not to calculate and report flux variability in the
            output summary


        perhaps I should remove the pandas requirement? ive used it pretty much
        everywhere else though, so once more shouldnt hurt.

        """

        obj_fluxes = pd.Series({'{:<15}'.format(r.id): '{:.3f}'.format(r.x)
                                for r in self.objective.iterkeys()})

        if not fva:

            out_rxns = self.reactions.query(
                lambda rxn: rxn.x > tol).query('system_boundary', 'boundary')
            in_rxns = self.reactions.query(
                lambda rxn: rxn.x < -tol).query('system_boundary', 'boundary')
            
            out_fluxes = pd.Series({r.reactants[0] : r.x for r in out_rxns})
            in_fluxes = pd.Series({r.reactants[0] : r.x for r in in_rxns})

            out_fluxes = pd.np.round(out_fluxes.sort_values(ascending=False), round)
            in_fluxes  = pd.np.round(in_fluxes.sort_values(), round)

            table = pd.np.array(
                [((a if a else ''), (b if b else ''), (c if c else ''))
                 for a, b, c in itertools.izip_longest(
                         ['IN FLUXES'] + in_fluxes.to_string().split('\n'), 
                         ['OUT FLUXES'] + out_fluxes.to_string().split('\n'),
                         ['OBJECTIVES'] + obj_fluxes.to_string().split('\n'))])


        else:
            from ..flux_analysis.variability import flux_variability_analysis
            fva_results = pd.DataFrame(
                flux_variability_analysis(self, fraction_of_optimum=fva)).T

            half_span = (fva_results.maximum - fva_results.minimum)/2
            median = fva_results.minimum + half_span

            out_rxns = self.reactions.query(
                lambda rxn: median.loc[rxn.id] > tol
            ).query('system_boundary', 'boundary')

            in_rxns = self.reactions.query(
                lambda rxn: median.loc[rxn.id] < -tol
            ).query('system_boundary', 'boundary')

            out_fluxes = pd.DataFrame(
                {r.reactants[0] : {'x'   : median.loc[r.id], 
                                   'err' : half_span.loc[r.id]}
                 for r in out_rxns}).T

            in_fluxes = pd.DataFrame(
                {r.reactants[0] : {'x'   : median.loc[r.id], 
                                   'err' : half_span.loc[r.id]}
                 for r in in_rxns}).T

            out_fluxes = pd.np.round(
                out_fluxes.sort_values(by='x', ascending=False), round)
            in_fluxes  = pd.np.round(
                in_fluxes.sort_values(by='x'), round)

            in_fluxes_s = in_fluxes.apply(
                lambda x: u'{0:0.2f} \u00B1 {1:0.2f}'.format(x.x, x.err),
                axis=1)
            out_fluxes_s = out_fluxes.apply(
                lambda x: u'{0:0.2f} \u00B1 {1:0.2f}'.format(x.x, x.err),
                axis=1)
            out_fluxes_s = out_fluxes.apply(lambda x: unicode(x.x) + u" \u00B1 "
                                          + unicode(x.err), axis=1)
            # out_fluxes_s = out_fluxes.apply(lambda x: unicode(x.x) + u" \u00B1 "
            #                               + unicode(x.err), axis=1)


            table = pd.np.array(
                [((a if a else ''), (b if b else ''), (c if c else ''))
                 for a, b, c in itertools.izip_longest(
                         ['IN FLUXES'] + in_fluxes_s.to_string().split('\n'), 
                         ['OUT FLUXES'] + out_fluxes_s.to_string().split('\n'),
                         ['OBJECTIVES'] + obj_fluxes.to_string().split('\n'))])



        print u'\n'.join([u"{a:<30}{b:<30}{c:<20}".format(a=a, b=b, c=c) for a,b,c in table])


    def to_json(self, filename, pretty=False):
        """ Save the model to a json file.
        
        """
        from ..io import save_json_model
        save_json_model(self, filename, pretty=pretty)


    def add_pathway(self, pathway):
        """Add a cobra.Pathway to the model. Requires that all metabolites,
        reactions are already present. Pathways will not change
        any model behavior, but are merely a convienience feature for
        understanding and grouping reactions 

        pathway: a cobra.Pathway object

        """

        # Assert that pathway is not already in the model, and that all
        # reactions and metabolites in the pathway are already in the model.
        if self.pathways.has_id(pathway.id):
            raise Exception("Pathway already in model: " +
                            pathway.id)

        reactions_not_in_model = [
            i.id for i in pathway.reactions if not self.reactions.has_id(
                i.id)]

        if len(reactions_not_in_model) > 0:
            raise Exception("Reactions not in the model: " +
                            ", ".join(reactions_not_in_model))


        metabolites_not_in_model = [
            i.id for i in pathway.metabolites if not self.metabolites.has_id(
                i.id)]

        if len(metabolites_not_in_model) > 0:
            raise Exception("Metabolites not in the model: " +
                            ", ".join(metabolites_not_in_model))

        pathway._model = self

        # Superpathways will need to be added later, assuming that connections
        # between pathways have been broken.
        pathway.superpathways = DictList()

        # Associate subpathways with those already in the model
        model_subpathways = []
        for subpathway in pathway.subpathways:
            try:
                model_subpathway = self.pathways.get_by_id(subpathway.id)
            except KeyError:
                self.add_pathway(subpathway)
                model_subpathway = self.pathways.get_by_id(subpathway.id)

            model_subpathway.superpathways.append(pathway)

            model_subpathways += [model_subpathway]

        pathway.subpathways = model_subpathways


        # Associate reactions with those already in the model
        model_reactions = {self.reactions.get_by_id(reaction.id) : stoich for
                           reaction, stoich in pathway.reactions.iteritems()}
        pathway.reactions = model_reactions

        # Associate metabolites with those already in the model
        for metabolite in list(pathway._metabolites.keys()):
            stoichiometry = pathway._metabolites.pop(metabolite)
            model_metabolite = self.metabolites.get_by_id(
                metabolite.id)
            pathway._metabolites[model_metabolite] = stoichiometry

        self.pathways.append(pathway)

        
    def add_pathways(self, pathways):
        """Add multiple pathways to the model.

        """
        for pathway in pathways:
            self.add_pathway(pathway)


