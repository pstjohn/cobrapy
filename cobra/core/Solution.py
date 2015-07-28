import json

class Solution(object):
    """Stores the solution from optimizing a cobra.Model. This is
    used to provide a single interface to results from different
    solvers that store their values in different ways.

    f: The objective value

    solver: A string indicating which solver package was used.

    x: List or Array of the values from the primal.

    x_dict: A dictionary of reaction ids that maps to the primal values.

    y: List or Array of the values from the dual.

    y_dict: A dictionary of reaction ids that maps to the dual values.

    """

    def __init__(self, f, x=None,
                 x_dict=None, y=None, y_dict=None,
                 solver=None, the_time=0, status='NA'):
        self.solver = solver
        self.f = f
        self.x = x
        self.x_dict = x_dict
        self.status = status
        self.y = y
        self.y_dict = y_dict

    def save(self, filename=None):
        if filename == None:
            # Will raise an attribute error if there's no attached solution
            filename = self.model.id + '_solution.json'
        with open(filename, 'wb') as f:
            json.dump(self.x_dict, f)
        return

    def __repr__(self):
        if self.f is None:
            return "<Solution '%s' at 0x%x>" % (self.status, id(self))
        return "<Solution %.2f at 0x%x>" % (self.f, id(self))
