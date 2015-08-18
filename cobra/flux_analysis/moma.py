"""
As a non-academic, I'm sadly not able to get free licenses for CPLEX or Gurobi,
making the current MOMA implementation inaccessible.  In general, I'm somewhat
surprised COBRApy doesn't just rely on a third-party python optimization repo
(like CVXPY, PuLP, Pyomo) to keep a consistent solver interface. As it stands,
its fairly hard to do any optimizations outside of the cobra.Model structure,
for instance including customized objective functions.  The
cobra.flux_analysis.moma script or the parsimonious fba seem to suffer from
this - this will be my attempt to fix these issues by introducing additional
dependencies.


# TODO: make `x_wildtype` a dict, so that reactions can be completely removed
# from the model without destroying the moma evalutions.
"""

import itertools
import numpy as np
import cvxpy as cvx

from ..core import Model, Reaction, Solution


class MOMAModel(Model):

    def __init__(self, wt_model, deepcopy_model=False):
        """ Calculates the results of reaction knockouts using the technique of
        Minimization of Metabolic Adjustment. This class is intended to act as
        a drop-in replacement for cobra.Model, where knockout simulations
        performed by 

        >>> moma_model = MOMAModel(wt_model)
        >>> moma_model.reactions.XX.knockout()
        >>> moma_model.optimize()

        Are instead simulated via MOMA to minimize distance from the flux state
        specified in wt_model's objective.
    
        Parameters:
        ===========

        wt_model: cobra.Model
            A model describing the wild-type organism. Should likely be
            optimized, but otherwise will get optimized with parisomonius FBA.

        deepcopy_model: bool
            Whether or not to make a copy of the model

        For details and motiviation of the algorithm, see:
        Segre, D., Vitkup, D., & Church, G. M. (2002). Analysis of optimality
        in natural and perturbed metabolic networks. Proceedings of the
        National Academy of Sciences of the United States of America, 99(23),
        15112-15117.  http://doi.org/10.1073/pnas.232349399

        """

        if deepcopy_model: wt_model = wt_model.copy()

        # Optimize the wt_model to get a list of desired fluxes
        if not wt_model.solution.x:
            wt_model.optimize(minimize_absolute_flux=1.0)

        self.x_wildtype = np.array(wt_model.solution.x) # numpy frozenarray?

        Model.__init__(self, wt_model)

        # create cvxpy problem
        self.create_cvx_problem()

        # create + update MOMAreactions
        for reaction in self.reactions:
            lb = float(reaction.lower_bound)
            ub = float(reaction.upper_bound)
            reaction.__class__ = MOMAReaction
            reaction._lower_bound = lb
            reaction._upper_bound = ub



    def create_cvx_problem(self, method='quadratic'):
        """ Create CVXPY problem from the model. Using 'v' to denote the
        symbolic flux values.

        method: 'quadratic' or 'linear'
            The method used to compare flux values before and after knockout.

        """

        A_red, self._c, ub, lb, nr, nm = generate_model_matricies(self)

        # Unknown fluxes
        self._v = cvx.Variable(nr, name='flux')

        # Initialize bound parameters (allows these to easily change)
        self._lower_bound = cvx.Parameter(nr)
        self._upper_bound = cvx.Parameter(nr)
        self._lower_bound.value = lb
        self._upper_bound.value = ub
        
        constraints = [
            A_red*self._v == 0,     # Stoichiometric constraint
            self._v >= self._lower_bound, # Flux contraints
            self._v <= self._upper_bound,
        ] 
        
        if method == 'quadratic':
            objective = cvx.Minimize(cvx.sum_squares(self._v - self.x_wildtype))

        elif method == 'linear':
            objective = cvx.minimize(cvx.sum_entries(
                cvx.abs(self._v - self.x_wildtype)))

        else: raise RuntimeError('Method {} not supported'.format(method))

        self._cvx_prob = cvx.Problem(objective, constraints)


        
    def _update_lower_bound(self, reaction):
        """ Utility function for propogating bound changes to cvx_prob """
        reaction_index = self.reactions.index(reaction.id)
        self._lower_bound.value[reaction_index] = reaction.lower_bound

    def _update_upper_bound(self, reaction):
        """ Utility function for propogating bound changes to cvx_prob """
        reaction_index = self.reactions.index(reaction.id)
        self._upper_bound.value[reaction_index] = reaction.upper_bound

    
    def optimize(self, MOMA=True, method=None, **kwargs):
        """ Run the minimization of metabolic adjustment. If `method`, update
        the corresponding objective funciton.

        method: None, 'quadratic' or 'linear'
            The method used to compare flux values before and after knockout.
            If None, use pre-specified objective

        """

        # Only update the objective function if desired
        if method:
            if method == 'quadratic':
                self._cvx_prob.objective = cvx.Minimize(
                    cvx.sum_squares(self._v - self.x_wildtype))

            elif method == 'linear':
                self._cvx_prob.objective = cvx.Minimize(
                    cvx.sum_entries(cvx.abs(self._v - self.x_wildtype)))

        # Solve the MOMA problem
        optimum_dist = self._cvx_prob.solve()

        # Parse results for a cobrapy Solution object
        x_opt = np.array(self._v.value).flatten()
        f = self._c.dot(self._v.value)
        x_dict = {reaction.id : xi for reaction, xi in
                  itertools.izip(self.reactions, x_opt)}

        # Create the solution object. I use the original model's objective in
        # calculating f, as this is more representative of the effect of the
        # knockout.
        cobra_solution = Solution(f=f, x=x_opt, x_dict=x_dict,
                                  status=self._cvx_prob.status)

        # Regardless, we may want access to the achieved MOMA objective.
        # Attaching this to the Solution.
        cobra_solution.optimum_dist = optimum_dist

        self.solution = cobra_solution

        return cobra_solution

        



         



class MOMAReaction(Reaction):

    def __init__(self, reaction):
        """ Wrapper class for cobra.Reaction to allow the updating of the
        MOMAModel problem. Changes to bounds will automagically trigger updates
        in the owning model's cvx_prob.

        NOTE: Because MOMAModel lazily edits Reaction.__class__, this method
        will likely never be called.
        
        """
        Reaction.__init__(reaction)
        self._lower_bound = reaction.lower_bound
        self._upper_bound = reaction.upper_bound
        self.model = reaction.model


    @property
    def lower_bound(self):
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, val):
        self._lower_bound = val
        self.model._update_lower_bound(self)

    @property
    def upper_bound(self):
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, val):
        self._upper_bound = val
        self.model._update_upper_bound(self)




def independent_rows(A, tol = 1e-05):
    """
    Return an array composed of independent rows of A.

    Note the answer may not be unique; this function returns one of many
    possible answers.

    http://stackoverflow.com/q/13312498/190597 (user1812712)
    http://math.stackexchange.com/a/199132/1140 (Gerry Myerson)
    http://mail.scipy.org/pipermail/numpy-discussion/2008-November/038705.html
        (Anne Archibald)

    """
    Q, R = np.linalg.qr(A.T)
    independent = np.where(np.abs(R.diagonal()) > tol)[0]
    return A.T[:, independent].T


def generate_model_matricies(cobra_model):
    """ Generate stochiometric and bound matricies without creating the
    ArrayBasedModel for simplicity """

    nr = len(cobra_model.reactions)
    nm = len(cobra_model.metabolites)

    A = np.zeros((nm, nr))
    c = np.zeros(nr) # Objective function
    ub = np.zeros(nr) # reaction upper bounds
    lb = np.zeros(nr) # reaction upper bounds

    for i, reaction in enumerate(cobra_model.reactions):
        c[i] = reaction.objective_coefficient
        ub[i] = reaction.upper_bound
        lb[i] = reaction.lower_bound
        for metabolite, stoich in reaction.metabolites.iteritems():
            j = cobra_model.metabolites.index(metabolite.id)
            A[j,i] = stoich

    # We can (and need to for CVXOPT) remove metabolite balances which do not
    # add additional constraints to the system.
    A_red = independent_rows(A)

    return A_red, c, ub, lb, nr, nm




# TODO move these to cobra.tests
if __name__ == "__main__":

    import cobra.test
    model = cobra.test.create_test_model('textbook')
    test_moma_model = MOMAModel(model, deepcopy_model=True)

    # knocks out a reaction and runs the moma evaluation
    test_moma_model.reactions.PDH.knock_out()
    test_moma_model.optimize()
