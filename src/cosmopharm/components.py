from typing import Optional
from numbers import Number

class Component:
    def __init__(self,
                 name: Optional[str] = None,
                 Mw: Optional[Number] = None,  # Positive number expected
                 T_fus: Optional[Number] = None,  # Positive number expected
                 H_fus: Number = 0,
                 Cp_fus_a_fit: Number = 0,
                 Cp_fus_bT_fit: Number = 0,
                 v_298: Optional[Number] = None,
                 v_hc: Optional[Number] = None,
                 ):

        self.name = name
        self.Mw = Mw  # molar weight in g/mol

        # for SLE calculations
        self.T_fus = T_fus
        self.H_fus = H_fus
        self.Cp_fus_a_fit = Cp_fus_a_fit
        self.Cp_fus_bT_fit = Cp_fus_bT_fit
        self.v_298 = v_298
        self.v_hc = v_hc

    def __repr__(self):
        return f"<Component('{self.name}')>"

    def __str__(self):
        return f"{self.name}"
