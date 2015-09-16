from warnings import warn

try:
    import graphviz as gv
except ImportError:
    gv = False

from .Reaction import Reaction
from .DictList import DictList

class Pathway(Reaction):
    """Pathway is a class for holding meta-information on metabolic pathways in
    a cobra.Model object.  Pathways can be constructed from both reactions and
    other pathways (superpathways).

    """

    def __init__(self, name=None):
        """An object for housing reaction pathways, sets of related, linked
        reactions which accomplish metabolic functions 
        
        """

        Reaction.__init__(self, name)

        self.reactions = {}
        self.subpathways = DictList()
        self.superpathways = DictList()

    @property
    def metabolites(self):

        metabolites = {}
        for reaction, stoich in self.reactions.iteritems():
            for metabolite, stoich in reaction.metabolites.iteritems():
                try: 
                    metabolites[metabolite] += stoich
                except KeyError:
                    metabolites[metabolite] = stoich

                if metabolites[metabolite] == 0:
                    metabolites.pop(metabolite)
                        
        return metabolites
        

    def add_reaction(self, reaction, stoich=1.0):
        """Add a cobra.Reaction to the pathway

        """
        
        try: 
            self.reactions[reaction] += stoich
        except KeyError:
            self.reactions[reaction] = stoich

    def add_reactions(self, reactions, stoichiometries=None):
        """Add an iterable of reactions for the pathway

        """

        if stoichiometries == None: 
            stoichiometries = len(reactions) * [1.0]

        for reaction, stoich in zip(reactions, stoichiometries):
            self.add_reaction(reaction, stoich)

    def add_pathway(self, other_pathway, stoich=1.0):
        """Add a cobra.Pathway to the pathway to create a super pathway.

        """

        # We don't need to recursively descend into the sub-pathway tree, as
        # self.reactions should contain the complete list of reaction leaves

        self.add_reactions(other_pathway.reactions)

        # Subpathways only contains a list of the direct pathway children.
        self.subpathways.append(other_pathway)

        # Make sub-pathway aware of this super-pathway parent
        other_pathway.superpathways.append(self)


    def add_pathways(self, pathways, stoichiometries=None):
        """Add an iterable of subpathways to the pathway

        """

        if stoichiometries == None: 
            stoichiometries = len(pathways) * [1.0]

        for pathway, stoich in zip(pathways, stoichiometries):
            self.add_pathway(pathway, stoich)


    def visualize(self, view='ipython', filename=None):
        """Visualize the pathway using graphviz. Works best for smaller pathways.

        view: 'ipython' or 'default'
            How to show the resulting SVG graph:
                * 'ipython' displays the SVG inline using IPython.display.SVG.
                * 'default' uses the default graphviz view option.

        filename: string
            Directory in which to save SVG. defaults to a temporary file which
            is deleted.

        """
        if not gv: raise RuntimeError('Graphviz not available')

        # Create digraph
        dot = gv.Digraph(comment=self.id)
        
        visited_metabolites = set()

        nodes = []

        # Iterate over reactions in the pathway
        for reaction in self.reactions.iterkeys():

            reactants = []
            products = []

            # Get metabolite primaries from metacyc, if available
            reactants_p, products_p = self.notes['PRIMARIES'][
                reaction.notes['BIOCYC'][0]]

            # Create reaction node
            dot.node(reaction.id, _attributes={'shape' : 'diamond'})

            for metabolite in reaction.reactants:
                if metabolite.notes['BIOCYC'][0] in reactants_p:
                    # Metabolite is a primary reactant, add it to nodes if not present
                    if metabolite not in visited_metabolites:
                        dot.node(metabolite.id, _attributes={'shape' : 'oval'})
                        visited_metabolites.add(metabolite)

                    # Add edge between metabolite and reaction
                    dot.edge(metabolite.id, reaction.id, _attributes={
                        'dir' : 'both' if reaction.reversibility else None,
                        'weight' : '1',
                    })

                    reactants += [metabolite.id]

                else:
                    # Add a unique metabolite node to avoid complicated overlaps
                    dot.node(metabolite.id + '_' + reaction.id, _attributes={
                        'shape' : 'plaintext',
                        'label' : metabolite.id,
                    })

                    dot.edge(
                        metabolite.id + '_' + reaction.id, reaction.id,
                        _attributes={
                            'dir' : 'both' if reaction.reversibility else None,
                            'weight' : '0',
                        })

                    reactants += [metabolite.id + '_' + reaction.id]



            for metabolite in reaction.products:
                if metabolite.notes['BIOCYC'][0] in products_p:
                    # Metabolite is a primary reactant, add it to nodes if not present
                    if metabolite not in visited_metabolites:
                        dot.node(metabolite.id, _attributes={'shape' : 'oval'})
                        visited_metabolites.add(metabolite)

                    # Add edge between metabolite and reaction
                    dot.edge(reaction.id, metabolite.id, _attributes={
                        'dir' : 'both' if reaction.reversibility else None,
                        'weight' : '1',
                    })

                    products += [metabolite.id]

                else:
                    # Add a unique metabolite node to avoid complicated overlaps
                    dot.node(metabolite.id + '_' + reaction.id, _attributes={
                        'shape' : 'plaintext',
                        'label' : metabolite.id,
                    })

                    dot.edge(
                        reaction.id, metabolite.id + '_' + reaction.id,
                        _attributes={
                            'dir' : 'both' if reaction.reversibility else None,
                            'weight' : '0',
                        })

                    products += [metabolite.id + '_' + reaction.id]

            nodes += [reactants, products]

        # Ensure ranks are the same for common reactants and metabolites
        for nodeset in nodes:
            if len(nodeset) > 1:
                dot.body.append(u'\t{{rank=same; "{}"}}'.format('" "'.join(nodeset)))

        if view == 'ipython': return dot
        else: raise Exception('other views not yet implemented')

