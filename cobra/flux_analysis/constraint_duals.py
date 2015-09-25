import cvxpy as cvx
import pandas as pd
import numpy as np


from .moma import generate_model_matricies


def calculate_constraint_duals(cobra_model, tol=1E-6, solver='GLPK'):
    """Calculate the dual values for the upper and lower bounds of a cobra.Model.

    cobra_model: a cobra.Model object

    tol: float
        Numerical tolerance for setting a dual value to zero.

    solver: string
        The solver to use (from cvxpy)
    
    Dual variables for linear programming constraints indicate how much the
    objective function could be increased by a corresponding relaxation in the
    constraint. This function uses cvxpy to calculate dual variables (shadow
    prices) for the upper and lower bounds of the model, respectively. 
    
    This function returns a pandas DataFrame with upper and lower dual values,
    respectively.

    """

    A_red, c, ub, lb, nr, nm = generate_model_matricies(cobra_model)

    # Unknown fluxes
    v = cvx.Variable(nr, name='flux')

    # Initialize bound parameters (allows these to easily change)
    lower_bound = cvx.Parameter(nr)
    upper_bound = cvx.Parameter(nr)
    lower_bound.value = lb
    upper_bound.value = ub

    constraints = [
        A_red*v == 0,     # Stoichiometric constraint
        v >= lower_bound, # Flux contraints
        v <= upper_bound,
    ] 

    objective_coeffs = np.array([rxn.objective_coefficient for rxn in
                                cobra_model.reactions])

    otol = 1E-8
    objective = cvx.Maximize(objective_coeffs * v - otol*cvx.sum_entries(cvx.abs(v)))

    cvx_prob = cvx.Problem(objective, constraints)
    cvx_prob.solve(solver=solver)

    dual_vals = pd.DataFrame(
        np.array([constraints[1].dual_value.tolist(),
                  constraints[2].dual_value.tolist()]).squeeze().T,
        index=[r.id for r in cobra_model.reactions],
        columns=['lower_bound', 'upper_bound'])

    return dual_vals.applymap(lambda x: x if x > tol else 0.)



    
