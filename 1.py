import math

class Material:
    def __init__(self, name, density, young_modulus, yield_strength):
        self.name = name #nazwa materiału
        self.density = density #gestosc w g/cm^3
        self.young_modulus = young_modulus * 1e9 #moduł youga w Pa (z GPa)
        self.yield_strength = yield_strength

class Section:
    def __init__(self, shape, dimensions):
        self.shape = shape #kształt przekroju - prostokąt albo koło
        self.dimension = dimensions #wymiary; przy prostokątnym {"width": x, "height": y}, przy okrągłym {"radius": r}

    def pole_przekroju (self):
        if self.shape == 'rectangle':
            ppr = self.dimension['width'] * self.dimension['height']
            return ppr
        elif self.shape == 'circle':
            ppc = math.pi * self.dimension['radius']**2
            return ppc
        else:
            print('Nieobsługiwane pole przekroju')

    def moment_bezwladnosci (self):
        if self.shape == 'rectangle':
            mbr = self.dimension['width'] * self.dimension['height']**3/12
            return mbr
        elif self.shape == 'circle':
            mbc = math.pi * self.dimension ['radius']**4/4
        else:
            print('Nieobsługiwane pole przekroju')
