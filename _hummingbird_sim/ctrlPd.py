import numpy as np
import hummingbirdParam as P


class ctrlPD:
    def __init__(self):
        pass

    def update(self):
        
        return np.array([[P.Fe],[0]])


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








