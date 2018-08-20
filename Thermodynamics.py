__author__ = 'brianmendoza'


class NearestNeighbors:

    def __init__(self, sequence, offsequence):
        self.Ttable = {}
        self.Ttable.fromkeys(["AA", "TT"], (-7.9, -22.2))
        self.Ttable.fromkeys(["AT", "TA"])
        self.Ttable.fromkeys(["TA", ""])
        self.Ttable.fromkeys(["TG","CA"])
        self.Ttable.fromkeys(["AC", "GT"])
        self.Ttable.fromkeys(["",""])
        self.Ttable.fromkeys()
        self.Ttable = {"AATT": (-7.9, -22.2), "ATTA": (-7.2, -20.4), "TAAT": (-7.2, -21.3), "CAGT": (-8.5, -22.7), "GTCA": (-8.4, -22.4),
                  "CTGA": (-7.8, -21.0), "GACT": (-8.2, -22.2), "CGGC": (-10.6, -27.2), "GCCG": (-9.8, -24.4), "GGCC": (-8.0, -19.9),
                  "AT": (2.3, 4.1), "GC": (0.1, -2.8)}
        self.MatchTable = {}
        #these keys are the RNA sequence values from Sugimoto et al. 1995
        self.MatchTable["AA"] = (-7.8, -21.9, -1.0)  #goes enthalpy, entropy, and free energy in kcal-cal-kcal units
        self.MatchTable["AC"] = (-5.9, -12.3, -2.1)
        self.MatchTable["AG"] = (-9.1, -23.5, -1.8)
        self.MatchTable["AU"] = (-8.3, -23.9, -0.9)
        self.MatchTable["CA"] = (-9.0, -26.1, -0.9)
        self.MatchTable["CC"] = (-9.3, -23.2, -2.1)
        self.MatchTable["CG"] = (-16.3, -47.1, -1.7)
        self.MatchTable["CU"] = (-7.0, -19.7, -0.9)
        self.MatchTable["GA"] = (-5.5, -13.5, -1.3)
        self.MatchTable["GC"] = (-8.0, -17.1, -2.7)
        self.MatchTable["GG"] = (-12.8, -31.9, -2.9)
        self.MatchTable["GU"] = (-7.8, -21.6, -1.1)
        self.MatchTable["UA"] = (-7.8, -23.2, -0.6)
        self.MatchTable["UC"] = (-8.6, -22.9, -1.5)
        self.MatchTable["UG"] = (-10.4, -28.4, -1.6)
        self.MatchTable["UU"] = (-11.5, -36.4, -0.2)
        self.initiation = (1.9, -3.9, 3.1)

    def enthalpy_calc(self, seq):
        enthalpy = 0
        start = seq[0]
        if start == 'G' or 'C':
            enthalpy += self.Ttable["GC"][0]
        else:
            enthalpy += self.Ttable["AT"][0]
        for i in range(0,len(seq)-1):
            first = seq[i]
            second = seq[i+1]
        end = seq[-1]

    def free_energy_calc(self):



    def mismatch_calc(self):




n = NearestNeighbors("ATTCGAG", "None")

