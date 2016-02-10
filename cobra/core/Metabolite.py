from warnings import warn
import re
from copy import deepcopy

import pandas as pd

from .Species import Species

# Numbers are not required because of the |(?=[A-Z])? block. See the
# discussion in https://github.com/opencobra/cobrapy/issues/128 for
# more details.
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


class Metabolite(Species):
    """Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

    """

    def __init__(self, id=None, formula=None,
                 name="", compartment=None):
        """
        id: str

        formula: str
            Chemical formula (i.e. H2O)

        name: str
            A human readable name.

        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object

        """
        Species.__init__(self, id, name)
        self.formula = formula
        # because in a Model a metabolite may participate in multiple Reactions
        self.compartment = compartment
        self.charge = None

        self._constraint_sense = 'E'
        self._bound = 0.

    @property
    def formula(self):
        """Describes a Chemical Formula
        A legal formula string contains only letters and numbers.
        """
        try: return self._formula
        except AttributeError:
            # Handle loading old pickled classes
            self._formula = self.__dict__['formula']
            return self._formula

    @formula.setter
    def formula(self, formula):
        formula = str(formula)
        if "*" in formula:
            warn("invalid character '*' found in formula '%s'" % formula)
            formula = formula.replace("*", "")
        if "(" in formula or ")" in formula:
            warn("invalid formula (has parenthesis) in '%s'" % formula)
        self._formula = formula

    @property
    def elements(self):
        """A dictionary breaking the chemical formula down by element."""

        def int_or_float(s):
            try:
                return int(s)
            except ValueError:
                # Do we want float warnings?
                # I don't, but they would go here.
                return float(s)

        return {element : int_or_float(stoich) if stoich else 1 
                for element, stoich in element_re.findall(self.formula)}


    @elements.setter
    def elements(self, elements_dict):

        def formula_items():
            for element, stoich in elements_dict.items():
                # Ignore stoichiometric entries close to machine precision.
                if abs(stoich) < 1E-10: continue
                yield ''.join((element, str(stoich) if stoich != 1 else ''))

        self.formula = ''.join(formula_items())

    @property
    def formula_weight(self):
        """Calculate the formula weight"""
        try:
            return sum([count * elements_and_molecular_weights[element]
                        for element, count in self.elements.items()])
        except KeyError as e:
            warn("The element %s does not appear in the peridic table" % e)

    @property
    def y(self):
        """The shadow price for the metabolite in the most recent solution

        Shadow prices are computed from the dual values of the bounds in
        the solution.

        """
        try:
            return self._model.solution.y_dict[self.id]
        except Exception as e:
            if self._model is None:
                raise Exception("not part of a model")
            if not hasattr(self._model, "solution") or \
                    self._model.solution is None or \
                    self._model.solution.status == "NA":
                raise Exception("model has not been solved")
            if self._model.solution.status != "optimal":
                raise Exception("model solution was not optimal")
            raise e  # Not sure what the exact problem was

    def remove_from_model(self, method='subtractive', **kwargs):
        """Removes the association from self.model

        method: 'subtractive' or 'destructive'.
            If 'subtractive' then the metabolite is removed from all
            associated reactions.  If 'destructive' then all associated
            reactions are removed from the Model.

        """
        # why is model being taken in as a parameter? This plays
        # back to the question of allowing a Metabolite to be associated
        # with multiple Models
        if "model" in kwargs:
            warn("model argument deprecated")

        self._model.metabolites.remove(self)
        self._model = None
        if method.lower() == 'subtractive':
            for the_reaction in list(self._reaction):
                the_coefficient = the_reaction._metabolites[self]
                the_reaction.subtract_metabolites({self: the_coefficient})
        elif method.lower() == 'destructive':
            rxns_to_remove = list(self._reaction)
            for x in rxns_to_remove:
                x.remove_from_model()
        else:
            raise Exception(method + " is not 'subtractive' or 'destructive'")


    def summary(self, ignore_inactive=0.01, element=None, ret_df=False,
                fva=None):
        """ Print a summary of the reactions which produce and consume this
        metabolite 
        
        ignore_inactive: float
            A [0, 1] value which specifies the fraction beneath which fluxes
            are ignored

        element: None or str
            If not None, multiplies fluxes the by the composition of the given
            element.
            
        ret_df: bool
            If true, return the flux_summary dataframe

        """

        def rxn_generator():
            for rxn in self.reactions:
                return_dict = {
                    'id' : rxn.id,
                }

                # Correct the direction of the reaction. A positive flux
                # produces the metabolite
                return_dict['flux'] = rxn.x * rxn.metabolites[self]

                # Correct to moles of element if one was provided
                if element: return_dict['flux'] *= self.elements[element]
                

                return_dict['reaction'] = rxn.reaction
                # # Correct the reaction direction for reactions running in
                # # reverse
                # if rxn.x >= 0:
                #     return_dict['reaction'] = rxn.reaction

                # elif rxn.x < 0:
                #     # Invert reaction direction
                #     return_dict['reaction'] = (-rxn).reaction

                yield return_dict

        flux_summary = pd.DataFrame(rxn_generator())

        assert flux_summary.flux.sum() < 1E-6, "Error in flux balance"

        if ret_df:
            # Sort and return the flux dataframe
            flux_summary.sort_values('flux', ascending=False, inplace=True)
            return flux_summary

        producing = flux_summary[flux_summary.flux > 0].copy()
        consuming = flux_summary[flux_summary.flux < 0].copy()

        for df in [producing, consuming]:
            df['percent'] = df.flux / df.flux.sum()
            df.drop(df[df['percent'] < ignore_inactive].index, axis=0, inplace=True)

            # df['%'] = df['%'].map('{:.1f}%'.format)

        producing.sort_values('percent', ascending=False, inplace=True)
        consuming.sort_values('percent', ascending=False, inplace=True)

        if not fva: 

            producing.flux = producing.flux.apply(
                lambda x: '{:6.2g}'.format(x))
            consuming.flux = consuming.flux.apply(
                lambda x: '{:6.2g}'.format(x))

            flux_len = 6

        else:
            from ..flux_analysis.variability import flux_variability_analysis
            fva_results = pd.DataFrame(
                flux_variability_analysis(self.model, self.reactions,
                                          fraction_of_optimum=fva)).T
            half_span = (fva_results.maximum - fva_results.minimum)/2
            median = fva_results.minimum + half_span

            producing.flux = producing.id.apply(
                lambda x, median=median, err=half_span: 
                u'{0:0.2f} \u00B1 {1:0.2f}'.format(median[x], err[x]))
            consuming.flux = consuming.id.apply(
                lambda x, median=median, err=half_span: 
                u'{0:0.2f} \u00B1 {1:0.2f}'.format(median[x], err[x]))

            flux_len = max(producing.flux.apply(len).max(), 
                           consuming.flux.apply(len).max()) + 1


        for df in [producing, consuming]:

           df['reaction'] = df['reaction'].map(lambda x: x[:54])
           df['id'] = df['id'].map(lambda x: x[:8])



        head = "PRODUCING REACTIONS -- " + self.name[:55]
        print head
        print "-"*len(head)
        print ("{0:^6} {1:>" + str(flux_len) + "} {2:>8} {3:^54}").format(
                   '%', 'FLUX', 'RXN ID', 'REACTION')

        for row in producing.iterrows():
            print (u"{0.percent:6.1%} {0.flux:>" + str(flux_len) + 
                   "} {0.id:>8} {0.reaction:>54}").format(row[1])


        print
        print "CONSUMING REACTIONS -- " + self.name[:55]
        print "-"*len(head)
        print ("{0:^6} {1:>" + str(flux_len) + "} {2:>8} {3:^54}").format(
                   '%', 'FLUX', 'RXN ID', 'REACTION')

        for row in consuming.iterrows():
            print (u"{0.percent:6.1%} {0.flux:>" + str(flux_len) + 
                   "} {0.id:>8} {0.reaction:>54}").format(row[1])


    def __add__(self, other_metabolite):
        """ Create a metabolite pool by adding together two metabolites. Useful
        for seeing flux into and out of groups of metabolites.
        
        Using a bunch of internal class methods here. Honestly not sure why
        some of them are protected, but this should likely only get used to
        call Metabolite.summary() on a pooled metabolite.

        """
        pool = deepcopy(self)
        other = deepcopy(other_metabolite)

        # Pop off the compartment subscript, if applicable.
        pool.id = re.sub('_.$', '', self.id) + '_' + re.sub('_.$', '', other.id)

        # (A + B + C) should yeild a name of "A and B and C Pool". Need to pull
        # off the pool each time a new one is added.
        old_name = re.sub(' Pool', '', self.name)
        if ' and ' in old_name:
            old_name = re.sub(' and ', ', ', old_name)
            old_name += ',' # Oxford comma
        pool.name = (' and '.join([old_name, other.name])) + ' Pool'

        model = self.model.copy()

        for reaction in set(self.reactions).union(other.reactions):
            pooled_reaction = deepcopy(reaction)
            pooled_reaction._model = model

            # Find the combined stochiometry for the pooled metabolite in each
            # reaction
            stoich_self = 0 # Initialize stoich counters
            stoich_other = 0
            for metabolite in pooled_reaction.metabolites.keys():
                if metabolite.id == self.id:
                    stoich_self = pooled_reaction._metabolites.pop(metabolite)
                if metabolite.id == other.id:
                    stoich_other = pooled_reaction._metabolites.pop(metabolite)
            pooled_reaction.add_metabolites({pool : stoich_self + stoich_other})
    
        return pool
        

        







elements_and_molecular_weights = {
    'H':   1.007940,
    'He':  4.002602,
    'Li':  6.941000,
    'Be':  9.012182,
    'B':   10.811000,
    'C':   12.010700,
    'N':   14.006700,
    'O':   15.999400,
    'F':   18.998403,
    'Ne':  20.179700,
    'Na':  22.989770,
    'Mg':  24.305000,
    'Al':  26.981538,
    'Si':  28.085500,
    'P':   30.973761,
    'S':   32.065000,
    'Cl':  35.453000,
    'Ar':  39.948000,
    'K':   39.098300,
    'Ca':  40.078000,
    'Sc':  44.955910,
    'Ti':  47.867000,
    'V':   50.941500,
    'Cr':  51.996100,
    'Mn':  54.938049,
    'Fe':  55.845000,
    'Co':  58.933200,
    'Ni':  58.693400,
    'Cu':  63.546000,
    'Zn':  65.409000,
    'Ga':  69.723000,
    'Ge':  72.640000,
    'As':  74.921600,
    'Se':  78.960000,
    'Br':  79.904000,
    'Kr':  83.798000,
    'Rb':  85.467800,
    'Sr':  87.620000,
    'Y':   88.905850,
    'Zr':  91.224000,
    'Nb':  92.906380,
    'Mo':  95.940000,
    'Tc':  98.000000,
    'Ru':  101.070000,
    'Rh':  102.905500,
    'Pd':  106.420000,
    'Ag':  107.868200,
    'Cd':  112.411000,
    'In':  114.818000,
    'Sn':  118.710000,
    'Sb':  121.760000,
    'Te':  127.600000,
    'I':   126.904470,
    'Xe':  131.293000,
    'Cs':  132.905450,
    'Ba':  137.327000,
    'La':  138.905500,
    'Ce':  140.116000,
    'Pr':  140.907650,
    'Nd':  144.240000,
    'Pm':  145.000000,
    'Sm':  150.360000,
    'Eu':  151.964000,
    'Gd':  157.250000,
    'Tb':  158.925340,
    'Dy':  162.500000,
    'Ho':  164.930320,
    'Er':  167.259000,
    'Tm':  168.934210,
    'Yb':  173.040000,
    'Lu':  174.967000,
    'Hf':  178.490000,
    'Ta':  180.947900,
    'W':   183.840000,
    'Re':  186.207000,
    'Os':  190.230000,
    'Ir':  192.217000,
    'Pt':  195.078000,
    'Au':  196.966550,
    'Hg':  200.590000,
    'Tl':  204.383300,
    'Pb':  207.200000,
    'Bi':  208.980380,
    'Po':  209.000000,
    'At':  210.000000,
    'Rn':  222.000000,
    'Fr':  223.000000,
    'Ra':  226.000000,
    'Ac':  227.000000,
    'Th':  232.038100,
    'Pa':  231.035880,
    'U':   238.028910,
    'Np':  237.000000,
    'Pu':  244.000000,
    'Am':  243.000000,
    'Cm':  247.000000,
    'Bk':  247.000000,
    'Cf':  251.000000,
    'Es':  252.000000,
    'Fm':  257.000000,
    'Md':  258.000000,
    'No':  259.000000,
    'Lr':  262.000000,
    'Rf':  261.000000,
    'Db':  262.000000,
    'Sg':  266.000000,
    'Bh':  264.000000,
    'Hs':  277.000000,
    'Mt':  268.000000,
    'Ds':  281.000000,
    'Rg':  272.000000,
    'Cn':  285.000000,
    'Uuq': 289.000000,
    'Uuh': 292.000000
}
