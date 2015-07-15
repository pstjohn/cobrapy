from ..manipulation.modify import convert_to_irreversible, revert_to_reversible


def optimize_maf(cobra_model, fraction_of_optimum=1.0):
    """optimize while minimizing absolute flux
    
    This function implements the strategy of maximizing an objective (ie
    Biomass) while minimizing the summed absolute fluxes in the model. In
    addition to the biological motivations (Scheutz et al., Science 2012), this
    function returns a more understandable flux vector (similar to loopless
    COBRA) while staying in the realm of LP.

    Updates everything in-place, returns model to original state at end.
    """

    # Convert to irreversible, so all reactions will have a positive flux
    convert_to_irreversible(cobra_model)

    # Maximize the existing objective function. (Store it for later use)
    original_objective = dict(cobra_model.objective)
    cobra_model.optimize()

    # Set the lower bound of the objective function to its 
    # maximum value x fraction_of_optimum
    for reaction in original_objective.iterkeys():
        reaction.lower_bound = reaction.x * fraction_of_optimum

    # Set the new objective as the minimization of all fluxes
    cobra_model.objective = {r : -1 for r in cobra_model.reactions}
    cobra_model.optimize()

    # Return the model to its original state
    revert_to_reversible(cobra_model)
    cobra_model.objective = original_objective
