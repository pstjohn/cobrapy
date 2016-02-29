class Object(object):
    """Defines common behavior of object in cobra.core"""

    def __init__(self, id=None, name=""):
        """
        id: None or a string

        """
        self.id = id
        self.name = name

        self.notes = {'map_info' : {}}
        self.annotation = {}

    def __getstate__(self):
        """To prevent excessive replication during deepcopy."""
        state = self.__dict__.copy()
        if '_model' in state:
            state['_model'] = None
        return state

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.__str__(), id(self))

    _str_method = 'id'
    def __str__(self):
        return str(self.__getattribute__(self._str_method))
