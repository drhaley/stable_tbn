from math import inf as infinity

class Constraints:
    """
    User-defined constraints are placed into this container class
    """
    def __init__(self):
        pass

    def optimize(self):
        return True

    def sort(self):
        return True

    def max_polymers(self):
        return infinity
