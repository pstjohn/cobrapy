from warnings import warn

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

    def add_reactions(self, reactions, stoichometries=None):
        """Add an iterable of reactions for the pathway

        """

        if stoichometries == None: 
            stoichometries = len(reactions) * [1.0]

        for reaction, stoich in zip(reactions, stoichometries):
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


    def add_pathways(self, pathways, stoichometries=None):
        """Add an iterable of subpathways to the pathway

        """

        if stoichometries == None: 
            stoichometries = len(pathways) * [1.0]

        for pathway, stoich in zip(pathways, stoichometries):
            self.add_pathway(pathway, stoich)



