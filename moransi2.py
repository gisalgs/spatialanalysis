def moransi2(z, wlist):
    """
    Input
      z: list of values
      w: weight list
    Output
      I: Moran's I measure
    """
    n = len(z)
    d = 0.0
    var = 0.0
    mean = float(sum(z))/n
    for i in range(n):
        var += (z[i]-mean)**2
    S0 = len(wlist)
    for e in wlist:
        d += (z[e[0]]-mean)*(z[e[1]]-mean)
    I = n*d/S0/var
    return I
