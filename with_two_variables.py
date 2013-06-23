#!/usr/bin/python2
# -*- coding: utf-8 -*-

from scipy import integrate, misc
import numpy as np
from numpy import linalg
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


def prawa(x, y):
    """prawa strona równania"""
    return g[0]*x + g[1]*y + g[2]*x*y + g[3]*x*x + g[4]*y*y


def lewa(x, y):
    return h[0]*x + h[1]*y + h[2]*x*y + h[3]*x*x + h[4]*y*y


def fi(nr, x, y, x1, x2, x3, y1, y2, y3):
    """funkcja kształtu - równanie płaszczyzny przechodzącej przez trzy punkty,
    gdzie punktowi nr przyporządkowuje się wartość 1, pozostałym 0"""
    if nr == 1:
        nominator = (x-x2) * (y3-y2) - (y-y2) * (x3-x2)
        denominator = (x1-x2) * (y3-y2) - (y1-y2) * (x3-x2)
    elif nr == 2:
        nominator = (x-x1) * (y3-y1) - (y-y1) * (x3-x1)
        denominator = (x2-x1) * (y3-y1) - (y2-y1) * (x3-x1)
    else:
        nominator = (x-x1) * (y2-y1) - (y-y1) * (x2-x1)
        denominator = (x3-x1) * (y2-y1) - (y3-y1) * (x2-x1)
    return nominator / denominator


def straight_line(y, x1, x2, y1, y2):
    """zwraca prostą przechodzącą przez dwa określone pukty
    o współrzędnych (x1, y1) i (x2, y2)"""
    return x1 + (y-y1)*(x2-x1)/(y2-y1)


def help_local(inner_function, xp, xq, xr, yp, yq, yr):
    xmin = min(xp, xq, xr)
    ymin = min(yp, yq, yr)

    def upper(y):
        return straight_line(y, xmin+ddx, xmin, ymin, ymin+ddy)

    def lower(y):
        return xmin

    return integrate.dblquad(inner_function, ymin, ymin+ddy, lower, upper)[0]


def local_right(xp, xq, xr, yp, yq, yr):
    """element wektoru prawych stron"""
    def inner_function(x, y):
        return prawa(x, y)*fi(1, x, y, xp, xq, xr, yp, yq, yr)

    return help_local(inner_function, xp, xq, xr, yp, yq, yr)


def local_main(nr, xp, xq, xr, yp, yq, yr):
    """macierz lokalna"""
    def pox(nr, x, y):
        """pochodna funkcji kształtu po x"""
        if nr == 1:
            return (yr-yq) / ((xp-xq)*(yr-yq) - (yp-yq)*(xr-xq))
        elif nr == 2:
            return (yr-yp) / ((xq-xp)*(yr-yp) - (yq-yp)*(xr-xp))
        return (yq-yp) / ((xr-xp)*(yq-yp) - (yr-yp)*(xq-xp))

    def poy(nr, x, y):
        """pochodna funkcji kształtu po y"""
        if nr == 1:
            return -(xr-xq) / ((xp-xq)*(yr-yq) - (yp-yq)*(xr-xq))
        elif nr == 2:
            return -(xr-xp) / ((xq-xp)*(yr-yp) - (yq-yp)*(xr-xp))
        return -(xq-xp) / ((xr-xp)*(yq-yp) - (yr-yp)*(xq-xp))

    def inner_function(x, y):
        result = - a[0]*pox(1, x, y)*pox(nr, x, y)
        result -= a[1]*pox(1, x, y)*poy(nr, x, y)
        result -= a[2]*poy(1, x, y)*poy(nr, x, y)
        result -= a[3]*poy(1, x, y)*pox(nr, x, y)
        rest = a[4] * pox(nr, x, y) + a[5] * poy(nr, x, y)
        rest += a[6] * lewa(x, y) * fi(nr, x, y, xp, xq, xr, yp, yq, yr)
        result += rest*fi(1, x, y, xp, xq, xr, yp, yq, yr)
        return result
    return help_local(inner_function, xp, xq, xr, yp, yq, yr)


def divide_into_triangles():
    """podzielenie powierzchni na trójkąty,
    zwraca krotkę z numerami węzłów odpowiadającymi wierzchołkom trójkąta"""
    siatka = []
    for x in range(nnx-1):
        for y in range(nny-1):
            siatka.append([nny*x + y+1, nny*(x+1) + y, nny*x + y])
            siatka.append([nny*x + y+1, nny*(x+1) + y, nny*(x+1) + y+1])
    return siatka


def set_nodes():
    """ustanowienie węzłów podziału,
    zwraca krotkę z ich współrzędnymi"""
    return [(xa+ddx*x, ya+ddy*y) for x in range(nnx) for y in range(nny)]


def agregate(wezel, siatka):
    """agregacja elementów z macierzy lokalnych do macierzy globalnych"""
    main_matrix = [[0 for i in range(nodes_number)] for j in range(nodes_number)]
    right_matrix = [0 for i in range(nodes_number)]
    for ii in range(triangles_number):
        for jj in range(3):
            i = siatka[ii][jj]
            j = siatka[ii][(jj+1)%3]
            k = siatka[ii][(jj+2)%3]
            xp = wezel[i][0]
            yp = wezel[i][1]
            xq = wezel[j][0]
            yq = wezel[j][1]
            xr = wezel[k][0]
            yr = wezel[k][1]
            right_matrix[i] += local_right(xp, xq, xr, yp, yq, yr)
            main_matrix[i][i] += local_main(1, xp, xq, xr, yp, yq, yr)
            main_matrix[i][j] += local_main(2, xp, xq, xr, yp, yq, yr)
            main_matrix[i][k] += local_main(3, xp, xq, xr, yp, yq, yr)
    # uwzgędnienie warunków brzegowych
    for ii in range(nodes_number):
        if wezel[ii][0] == xa and xaa != 'N':
            right_matrix[ii] = xaa
        elif wezel[ii][0] == xb and xbb != 'N':
            right_matrix[ii] = xbb
        elif wezel[ii][1] == ya and yaa != 'N':
            right_matrix[ii] = yaa
        elif wezel[ii][1] == yb and ybb != 'N':
            right_matrix[ii] = ybb
        else:
            continue
        for i in range(nodes_number):
            main_matrix[ii][i] = 0
        main_matrix[ii][ii] = 1
    return (main_matrix, right_matrix)


def shape_line(x, y, xc, yc):
    """funkcja przybliżająca rozwiązanie w otoczeniu węzła (xc, yc)"""
    if xc+ddx >= x >= xc and yc+ddy >= y >= yc and \
            ddx*(y-yc) + ddy*(x-xc) <= ddy*ddx:
        return fi(1, x, y, xc, xc+ddx, xc, yc, yc, yc+ddy)
    if xc+ddx >= x >= xc and yc-ddy <= y <= yc:
        if ddx*(y+ddy-yc) + ddy*(x-xc) > ddy*ddx:
            return fi(1, x, y, xc, xc+ddx, xc+ddx, yc, yc, yc-ddy)
        else:
            return fi(1, x, y, xc, xc+ddx, xc, yc, yc-ddy, yc-ddy)
    if xc-ddx <= x <= xc and yc-ddy <= y <= yc and \
            ddx*(y+ddy-yc) + ddy*(x+ddx-xc) >= ddy*ddx:
        return fi(1, x, y, xc, xc-ddx, xc, yc, yc, yc-ddx)
    if xc-ddx <= x <= xc and yc+ddy >= y >= yc:
        if ddx*(y-yc) + ddy*(x+ddx-xc) > ddy*ddx:
            return fi(1, x, y, xc, xc-ddx, xc, yc, yc+ddy, yc+ddy)
        else:
            return fi(1, x, y, xc, xc-ddx, xc-ddx, yc, yc, yc+ddy)
    return 0


def approximate(x, y, nodes, wsp):
    """funkcja przybliżająca rozwiązanie w całym zakresie"""
    result = 0
    for n in range(len(nodes)):
        result += wsp[n]*shape_line(x, y, nodes[n][0], nodes[n][1])
    return result


def plot(X, Y, Z, z_min=0):
    """przedstawienie wyniku na wykresie"""
    fig = pyplot.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z)
    cset = ax.contourf(X, Y, Z, zdir='z', offset=z_min, cmap=cm.coolwarm)
    pyplot.show()


def get_val(xy, ab):
    print """podaj warunek brzegowy Dirichleta dla {xy} = {xy}{ab}
            (podaj wartość funkcji) lub napisz N, jeśli to warunek Neumana
            (UWAGA: pochodna musi być tam równa 0)""".format(xy=xy, ab=ab)
    val = raw_input()
    try:
        val = int(val)
    except ValueError:
        if val != 'N':
            raise ValueError('''invalid literal - it should be 'N' or a number: '{}' provided'''.format(val))
    return val


if __name__ == "__main__":
    print "podaj granice poszukiwanego obszaru pierwszej zmiennej - 'xa xb'"
    xa, xb = [int(i) for i in raw_input().split()]
    print "podaj granice poszukiwanego obszaru dla drugiej zmiennej - 'ya yb'"
    ya, yb = [int(i) for i in raw_input().split()]
    print "podaj współczynniki równania - 'a1 a2 a3 a4 a5 a6 a7', gdzie\
            równanie wyraża się wzorem a1*d^2f(x,y)/dx^2 + a2*d^f(x,y)/dxdy + \
            a3*d^f(x,y)/dy^2 + a4*d^f(x,y)/dydx + a5*df(x,y)/dx + a6*df(x,y)/dy+\
            a7f(x,y) + g(x) = h(x)"
    a = [int(i) for i in raw_input().split()]
    print "podaj współczynniki funkcji g(x) - 'a b c d e', gdzie\
            g(x) = ax+by+cxy+dx^2+ey^2"
    g = [int(i) for i in raw_input().split()]
    print "podaj współczynniki funkcji h(x) - 'a b c d e', gdzie\
            h(x) = ax+by+cxy+dx^2+ey^2"
    h = [int(i) for i in raw_input().split()]
    xaa = get_val('x', 'a')
    xbb = get_val('x', 'b')
    yaa = get_val('y', 'a')
    ybb = get_val('y', 'b')
    print "podaj jak bardzo dokładny ma być wynik (1-10)"
    nnx = nny = 10*int(raw_input())
    ddx = float(xb - xa)/float(nnx - 1)
    ddy = float(yb - ya)/float(nny - 1)
    nodes_number = nnx*nny
    triangles_number = 2*(nnx-1)*(nny-1)
    wezel = set_nodes()
    siatka = divide_into_triangles()
    main_matrix, right_matrix = agregate(wezel, siatka)
    wsp = linalg.solve(main_matrix, right_matrix)
    X = np.arange(xa, xb+ddx, ddx)
    Y = np.arange(ya, yb+ddy, ddy)
    X, Y = np.meshgrid(X, Y)
    Z = np.array([approximate(x, y, wezel, wsp)
        for x, y in zip(np.ravel(X), np.ravel(Y))])
    Z = Z.reshape(X.shape)
    z_min = min(wsp)
    plot(X, Y, Z, z_min)
