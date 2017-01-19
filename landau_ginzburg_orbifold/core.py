class LandauGinzburgOrbifold:
    def __init__(self, group, polynomial=None, charges=None):
        self._group_ = group
        self._polynomial_ = polynomial
        self._charges_ = charges
