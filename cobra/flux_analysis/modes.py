from contextlib import contextmanager
from optlang.symbolics import add

@contextmanager
def flux_mode_constraints(model, c=1):
    """ Add indicator constraints to allow the enumeration of k-shortest
    elementary flux modes, as described in [1]_.

    Parameters
    ----------
    model: cobra.core.model
        Model for which to compute elementary flux modes
    c : float
        The minimum allowable flux through a given elementary flux mode.
        Defaults to 1.

    Yields
    ------
    model
        The input model with additional constraints added.
    indicator_variables : dict
        A dictionary of optlang.Variable objects corresponding to the added
        indicator variables.

    """
    with model:

        indicator_variables = dict()
        indicator_constraints = list()

        Variable = model.problem.Variable
        Constraint = model.problem.Constraint

        for reaction in model.reactions:
        
            if reaction.upper_bound > 0:
                
                fwd_id = 'y_fwd_' + reaction.id
                y_fwd = Variable(fwd_id, type='binary')
                indicator_variables[fwd_id] = y_fwd

                indicator_constraints += [Constraint(
                    reaction.forward_variable,
                    indicator_variable=y_fwd, active_when=0, lb=0,
                    ub=0,
                    name='indicator_constraint_fwd_1_{}'.format(reaction.id)
                )]

                indicator_constraints += [Constraint(
                    reaction.forward_variable,
                    indicator_variable=y_fwd, active_when=1, lb=c,
                    name='indicator_constraint_fwd_2_{}'.format(reaction.id)
                )]
        
            # Only if y is reversible
            if reaction.lower_bound < 0:
                
                rev_id = 'y_rev_' + reaction.id
                y_rev = Variable(rev_id, type='binary')
                indicator_variables[rev_id] = y_rev

                indicator_constraints += [Constraint(
                    reaction.reverse_variable,
                    indicator_variable=y_rev, active_when=0, lb=0,
                    ub=0,
                    name='indicator_constraint_rev_1_{}'.format(reaction.id)
                )]

                indicator_constraints += [Constraint(
                    reaction.reverse_variable,
                    indicator_variable=y_rev, active_when=1, lb=c,
                    name='indicator_constraint_rev_2_{}'.format(reaction.id)
                )]
            
            if reaction.reversibility:
                
                indicator_constraints += [Constraint(
                    y_fwd + y_rev, lb=0, ub=1,
                    name='one_direction_constraint_{}'.format(reaction.id)
                )]

        indicator_constraints += [Constraint(
            add(*indicator_variables.values()), lb=1,
            name='an_EM_must_constain_at_least_one_active_reaction'
        )]

        model.add_cons_vars(
            list(indicator_variables.values()) + indicator_constraints)

        model.objective = add(*indicator_variables.values())
        model.objective_direction = 'min'

        model._indicator_variables = indicator_variables
        
        yield model

        # Clean up additional model attribute
        del model._indicator_variables


def elementary_flux_modes(model, c=1):
    """ Calculate elementary flux modes, from shortest to longest, using the
    K-shortest approach [1]_. 

    Parameters
    ----------
    model: cobra.core.model
        Model for which to compute elementary flux modes
    c : float
        The minimum allowable flux through a given elementary flux mode.
        Defaults to 1.

    Yields
    ------
    efm : pandas.Series
        A pandas series corresponding to the current elementary flux mode

    References
    ----------
    .. [1] de Figueiredo LF, Podhorski A, Rubio A, et al. Computing the
       shortest elementary flux modes in genome-scale metabolic networks.
       Bioinformatics.  2009;25(23):3158-3165.
       doi:10.1093/bioinformatics/btp564.
    
    """

    with flux_mode_constraints(model, c) as model:
        pass



