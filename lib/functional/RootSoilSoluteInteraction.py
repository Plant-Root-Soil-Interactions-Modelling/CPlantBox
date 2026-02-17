

class RootSoilSoluteInteraction:

    def __init__(self, rs, vmax, km, cmin, exu):

        self.rs = rs

        self.vmax = vmax  # kg/(m2 day) -> York et. al (2016) (Lynch group)
        self.km = km  # kg/m3
        self.cmin = cmin  # kg/m3
        self.exu = exu  # kg /(m2 day)

    def set_defaults(self):

            self.vmax = 45.25 * 62 * 1.e-11 * (24.*3600.)  # kg/(m2 day) -> York et. al (2016) (Lynch group)
            self.km = 10.67 * 62 * 1.e-6  # kg/m3
            self.cmin = 4.4 * 62 * 1.e-6  # kg/m3
            self.exu = 1.e-2  # data Eva Oburger kg /(m2 day)

    def solute_fluxes(self, c):
        """ concentrations @param c [kg/m3], returns [g/day]
        self.Vmax [kg/(m2 day)]
        self.Km [kg/m3]
        self.Cmin [kg/m3]
        """
        segs = self.rs.segments
        a = self.rs.radii
        l = self.rs.segLength()
        assert(len(segs) == len(c))
        sf = np.zeros(len(segs),)
        c = np.maximum(c, self.CMin)  ###############################################
        for i, s in enumerate(segs):
            sf[i] = -2 * np.pi * a[i] * l[i] * 1.e-4 * ((c[i] - self.CMin) * self.Vmax / (self.Km + (c[i] - self.CMin)))  # kg/day
        sf = np.minimum(sf, 0.)
        return sf * 1.e3  # kg/day -> g/day

    def exudate_fluxes(self, rs_age, kex):
        """ returns [g/day]
        self.Exu [kg/(m2 day)]
        """
        segs = self.rs.segments
        tipI = np.array(self.rs.getRootTips()) - 1
        polylengths = self.rs.getParameter("length")
        types = self.rs.getParameter("type")

        a = self.rs.radii
        l = self.rs.segLength()
        sf = np.zeros((len(segs)))
        kex_all = np.zeros((len(segs)))

        for i, s in enumerate(tipI):
            if types[i] != 0:
                l_ = 0
                s_ = s
                while l_ <= kex[0][1]:
                    l_ = l_ + l[s_]
                    kexu = (kex[1][1] - kex[1][0]) / (kex[0][1] - kex[0][0]) * l_ + kex[1][0]  # linear decrease
                    kex_all[s_] = max(0, kexu)
                    sf[s_] = 2 * np.pi * a[s_] * l[s_] * 1.e-4 * kex_all[s_] * 1.e3  # g/day/root
                    s_ = s_ - 1
                    if l_ > polylengths[i]:
                        break

        return sf, kex_all  #  g/day
