from ..manipulation.modify import convert_to_irreversible, revert_to_reversible


def optimize_maf(cobra_model, fraction_of_optimum=1.0, **kwargs):
    """optimize while minimizing absolute flux
    
    This function implements the strategy of maximizing an objective (ie
    Biomass) while minimizing the summed absolute fluxes in the model. In
    addition to the biological motivations (Scheutz et al., Science 2012), this
    function returns a more understandable flux vector (similar to loopless
    COBRA) while staying in the realm of LP.

    Just realized this is the same functionality as optimize_minimum_flux in
    .parsimonious, but I like my implementation better.

    Updates everything in-place, returns model to original state at end.
    """

    # Store objective for later use
    original_objective = dict(cobra_model.objective)

    # Convert to irreversible, so all reactions will have a positive flux
    convert_to_irreversible(cobra_model)

    solution = None
    try: 
        # Maximize the existing objective function.
        solution = cobra_model.optimize(minimize_absolute_flux=False, **kwargs)

        # Minimizing absolute flux not relevant for infeasible solution
        if cobra_model.solution.f is not None:

            # Set the lower bound of the objective function to its 
            # maximum value x fraction_of_optimum
            original_lower_bounds = {}
            for reaction in original_objective.iterkeys():
                original_lower_bounds[reaction] = float(reaction.lower_bound)
                reaction.lower_bound = reaction.x * fraction_of_optimum

            # Set the new objective as the minimization of all fluxes
            cobra_model.objective = {r : -1 for r in cobra_model.reactions}
            solution = cobra_model.optimize(minimize_absolute_flux=False,
                                            **kwargs)

            # Reset lower bounds
            for reaction, lb in original_lower_bounds.iteritems():
                reaction.lower_bound = lb

    finally:
        # Return the model to its original state
        revert_to_reversible(cobra_model)
        cobra_model.objective = original_objective
        
        try:
            solution.f = sum([coeff * reaction.x for reaction, coeff in
                              cobra_model.objective.iteritems()])
        except Exception:
            # Sadly cobra will raise the ambiguous "Exception" when trying to
            # access a reaction.x if none exists, requiring this incredibly
            # vague exception block.
            pass

    return solution
