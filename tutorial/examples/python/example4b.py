import py_rootbox as rb


# User tropism 1: print input arguments to command line
class My_Info_Tropism(rb.Tropism):

    def tropismObjective(self, pos, old, a, b, dx, root):
        print("Postion \t", pos)
        print("Heading \t", old.column(0))
        print("Test for angle alpha = \t", a)
        print("Test for angle beta = \t", b)
        print("Eesolution of next segment \t", dx)
        print("Root id", root.getId())
        print()
        return 0.


# User tropism 2: depending on root age use plagio- or gravitropism
class My_Age_Tropism(rb.Tropism):

    def __init__(self, n, sigma, age):
        super(My_Age_Tropism, self).__init__()
        self.plagio = rb.Plagiotropism(0., 0.)
        self.gravi = rb.Gravitropism(0., 0.)
        self.setTropismParameter(n, sigma)
        self.age = age

    def tropismObjective(self, pos, old, a, b, dx, root):
        age = root.getAge()
        if age < self.age:
            d = self.plagio.tropismObjective(pos, old, a, b, dx, root)
            return d
        else:
            return self.gravi.tropismObjective(pos, old, a, b, dx, root)


# set up the root system
rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010"
rs.openFile(name)
rs.initialize()

# Set useer defined after initialize
mytropism1 = My_Info_Tropism()
mytropism1.setTropismParameter(2., 0.2)
mytropism2 = My_Age_Tropism(2., 0.5, 5.)  # after 5 days switch from plagio- to gravitropism
rs.setTropism(mytropism2, 2)  # 2 for laterals, -1 for all root types

# Simulate
simtime = 100  # e.g. 30 or 60 days
dt = 1
N = round(simtime / dt)
for _ in range(0, N):
    rs.simulate(dt)

# Export results (as vtp)
rs.write("../results/example_4b.vtp")

print("done.")
