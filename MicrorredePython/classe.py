class Microgrid:
    def __init__(self, A, B, C, D, E):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E


class Controller:
    def __init__(self, Np, Nu, delta, lambda_, alpha, Ymax, Ymin, 
                        Umax, Umin, DeltaUmax, DeltaUmin, 
                        Pbatmax, Pbatmin, DeltaPbatmax, DeltaPbatmin, ref):
        self.Np = Np
        self.Nu = Nu
        self.delta = delta
        self.lambda_ = lambda_
        self.alpha = alpha
        self.Ymax = Ymax
        self.Ymin = Ymin
        self.Umax = Umax
        self.Umin = Umin
        self.DeltaUmax = DeltaUmax
        self.DeltaUmin = DeltaUmin
        self.Pbatmax = Pbatmax
        self.Pbatmin = Pbatmin
        self.DeltaPbatmax = DeltaPbatmax
        self.DeltaPbatmin = DeltaPbatmin
        self.ref = ref
