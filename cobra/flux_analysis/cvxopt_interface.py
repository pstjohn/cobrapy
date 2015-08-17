"""
As a non-academic, I'm sadly not able to get free licenses for CPLEX or Gurobi.
This file will contain methods for optimizing lp and qp problems (perhaps MILP)
via CVXOPT.

In general, I'm somewhat surprised COBRApy doesn't just rely on a third-party
python optimization repo (like CVXOPT, PuLP, Pyomo) to keep a consistent solver
interface. As it stands, its fairly hard to do any optimizations outside of the
cobra.Model structure, for instance including customized objective functions. 
The cobra.flux_analysis.moma script or the parsimonious fba script could
probably both be simplified.

These functions likely won't be optimized for speed, but rather quickly provide
access to common FBA methods without the currently wrapped solvers.
"""

import itertools

import cvxopt as cvx
import numpy as np

cvx.solvers.options['show_progress'] = False
from ..core import Solution

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

    


def minimize_absolute_value(cobra_model, xg, solver='glpk'):
    """ Solves the linear program

    minimize: sum(| x - x_goal |)
    subject to:
        A x = b
        lb < x < ub

    where x is the vector of reaction fluxes, xg (x goal) are the desired values
    (i.e., 0 for parsimonious FBA, x_wt for MOMA), and A, ub, and lb are taken
    directly from the cobra_model.

    """

    A, c, ub, lb, nr, nm = generate_model_matricies(cobra_model)

    # Create matricies for CVXOPT
    # x_star = [ x z ]
    c_abs = np.concatenate((np.zeros(nr), np.ones(nr))) # Minimize all z variables
    A_abs = np.hstack((A, np.zeros(A.shape)))

    G_abs = np.vstack((
        np.hstack((-np.eye(nr),  np.zeros((nr,nr)))),
        np.hstack(( np.eye(nr),  np.zeros((nr, nr)))),
        np.hstack(( np.eye(nr), -np.eye(nr))),
        np.hstack((-np.eye(nr), -np.eye(nr)))
        ))

    h_abs = np.concatenate((
        -lb,  ub,
         xg, -xg,
    ))

    cvx_c_abs = cvx.matrix(c_abs)
    cvx_A_abs = cvx.matrix(A_abs)
    cvx_b_abs = cvx.matrix(np.zeros(A_abs.shape[0]))
    cvx_G_abs = cvx.matrix(G_abs)
    cvx_h_abs = cvx.matrix(h_abs)

    sol = cvx.solvers.lp(cvx_c_abs, cvx_G_abs, cvx_h_abs, cvx_A_abs, cvx_b_abs,
                         solver='glpk')
    x = np.array(sol['x'][:nr])

    x_dict = {reaction.id : x for reaction, x in
              itertools.izip(cobra_model.reactions, x)}
    
    # Sadly optimum dual values are not possible with this method since the
    # stochiometric matrix's rows are compressed to full rank.

    return Solution(f=c.dot(x), x=x, x_dict=x_dict, status=sol['status'])



def minimize_quadratic(cobra_model, xg):
    """ Solves the quadratic program

    minimize: sum( (x - x_goal)**2 )
    subject to:
        A x = b
        lb < x < ub

    where x is the vector of reaction fluxes, xg (x goal) are the desired values
    (i.e., 0 for parsimonious FBA, x_wt for MOMA), and A, ub, and lb are taken
    directly from the cobra_model.

    """

    A, c, ub, lb, nr, nm = generate_model_matricies(cobra_model)

    # Create matricies for CVXOPT
    cvx_P_quad = cvx.matrix(np.eye(nr))
    cvx_q_quad = cvx.matrix(xg)
    cvx_G_quad = cvx.matrix(np.vstack([-np.eye(nr), np.eye(nr)]))
    cvx_h_quad = cvx.matrix(np.concatenate((-lb, ub)))
    cvx_A = cvx.matrix(A)
    cvx_b = cvx.matrix(np.zeros(A.shape[0]))

    qpsol = cvx.solvers.qp(cvx_P_quad, cvx_q_quad, cvx_G_quad, cvx_h_quad,
                           cvx_A, cvx_b)
    x = np.array(qpsol['x'])

    x_dict = {reaction.id : x for reaction, x in
              itertools.izip(cobra_model.reactions, x)}
    
    # Sadly optimum dual values are not possible with this method since the
    # stochiometric matrix's rows are compressed to full rank.

    return Solution(f=c.dot(x), x=x, x_dict=x_dict, status=qpsol['status'])

