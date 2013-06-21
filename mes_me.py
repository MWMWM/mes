#!/usr/bin/python2
# -*- coding: utf-8 -*-

import math
import time
from scipy import integrate, misc
import numpy as np
from numpy import linalg
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages


def fi(nr, x, y, x1, x2, x3, y1, y2, y3):
    """funkcja kształtu - równanie płaszczyzny przechodzącej przez trzy punkty,
    gdzie punktowi nr przyporządkowuje się wartość 1, pozostałym 0"""
    if nr == 1:
        return ((x-x2)*(y3-y2)-(y-y2)*(x3-x2)) / ((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2))
    elif nr == 2:
        return ((x-x1)*(y3-y1)-(y-y1)*(x3-x1)) / ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
    return ((x-x1)*(y2-y1)-(y-y1)*(x2-x1)) / ((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))


def f(x, y):
    """prawa strona równania"""
    return x*y


def straight_line(y, x1, x2, y1, y2):
    """zwraca prostą przechodzącą przez dwa określone pukty
    o współrzędnych (x1, y1) i (x2, y2)"""
    return x1 + (y-y1)*(x2-x1)/(y2-y1)


def help_local(inner_function, xp, xq, xr, yp, yq, yr, ddx, ddy):
    xmin = min(xp, xq, xr)
    ymin = min(yp, yq, yr)

    def upper(y):
        return straight_line(y, xmin+ddx, xmin, ymin, ymin+ddy)

    def lower(y):
        return xmin
    return integrate.dblquad(inner_function, ymin, ymin+ddy, lower, upper)[0]


def local_right(xp, xq, xr, yp, yq, yr, ddx, ddy):
    """element wektoru prawych stron"""
    def inner_function(x, y):
        return f(x, y)*fi(1, x, y, xp, xq, xr, yp, yq, yr)
    return help_local(inner_function, xp, xq, xr, yp, yq, yr, ddx, ddy)


def local_main_old(nr, xp, xq, xr, yp, yq, yr, ddx, ddy):
    """macierz lokalna korzystająca z wbudowanej funkcji liczącej pochodne"""
    def deriv1(x, y, nr):
        return fi(nr, x, y, xp, xq, xr, yp, yq, yr)

    def deriv2(y, x, nr):
        return fi(nr, x, y, xp, xq, xr, yp, yq, yr)

    def inner_function(x, y):
        result = -misc.derivative(deriv1, x, args=(y, 1))*misc.derivative(
                deriv1, x, args=(y, nr))
        result -= misc.derivative(deriv2, y, args=(x, 1))*misc.derivative(
                deriv2, y, args=(x, nr))
        result -= 3*x*y*fi(1, x, y, xp, xq, xr, yp, yq, yr)*fi(
                nr, x, y, xp, xq, xr, yp, yq, yr)
        return result
    return help_local(inner_function, xp, xq, xr, yp, yq, yr, ddx, ddy)


def local_main(nr, xp, xq, xr, yp, yq, yr, ddx, ddy):
    """macierz lokalna"""
    def pox(nr, x, y, x1, x2, x3, y1, y2, y3):
        """pochodna funkcji kształtu po x"""
        if nr == 1:
            return (y3-y2) / ((x1-x2)*(y3-y2) - (y1-y2)*(x3-x2))
        elif nr == 2:
            return (y3-y1) / ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
        return (y2-y1) / ((x3-x1)*(y2-y1) - (y3-y1)*(x2-x1))

    def poy(nr, x, y, x1, x2, x3, y1, y2, y3):
        """pochodna funkcji kształtu po y"""
        if nr == 1:
            return -(x3-x2) / ((x1-x2)*(y3-y2) - (y1-y2)*(x3-x2))
        elif nr == 2:
            return -(x3-x1) / ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
        return -(x2-x1) / ((x3-x1)*(y2-y1) - (y3-y1)*(x2-x1))

    def inner_function(x, y):
        result = -pox(1, x, y, xp, xq, xr, yp, yq, yr)*pox(nr, x, y, xp, xq, xr, yp, yq, yr)
        result -= poy(1, x, y, xp, xq, xr, yp, yq, yr)*poy(nr, x, y, xp, xq, xr, yp, yq, yr)
        result -= 3*x*y*fi(1, x, y, xp, xq, xr, yp, yq, yr)*fi(nr, x, y, xp, xq, xr, yp, yq, yr)
        return result
    return help_local(inner_function, xp, xq, xr, yp, yq, yr, ddx, ddy)


def divide_into_triangles(nodes_number, nnx, nny):
    """podzielenie powierzchni na trójkąty,
    zwraca krotkę z numerami węzłów odpowiadającymi wierzchołkom trójkąta"""
    siatka = []
    for x in range(nnx-1):
        for y in range(nny-1):
            siatka.append([nny*x + y+1, nny*(x+1) + y, nny*x + y])
            siatka.append([nny*x + y+1, nny*(x+1) + y, nny*(x+1) + y+1])
    return siatka


def set_nodes(xa, ya, nnx, nny, ddx, ddy):
    """ustanowienie węzłów podziału,
    zwraca krotkę z ich współrzędnymi"""
    return [(xa+ddx*x, ya+ddy*y) for x in range(nnx) for y in range(nny)]


def agregacja(wezel, siatka, nodes_number, triangles_number, ddx, ddy):
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
            right_matrix[i] += local_right(xp, xq, xr, yp, yq, yr, ddx, ddy)
            main_matrix[i][i] += local_main(1, xp, xq, xr, yp, yq, yr, ddx, ddy)
            main_matrix[i][j] += local_main(2, xp, xq, xr, yp, yq, yr, ddx, ddy)
            main_matrix[i][k] += local_main(3, xp, xq, xr, yp, yq, yr, ddx, ddy)
    # uwzgędnienie warunków brzegowych
    for ii in range(nodes_number):
        if wezel[ii][0] == 0 or wezel[ii][1] == 0 or wezel[ii][1] == 1:
            for i in range(nodes_number):
                main_matrix[ii][i] = 0
            main_matrix[ii][ii] = 1
            right_matrix[ii] = 0
    return (main_matrix, right_matrix)


def shape_line(x, y, xc, yc, ddx, ddy):
    """funkcja przybliżająca rozwiązanie w otoczeniu węzła (xc, yc)"""
    if xc+ddx >= x >= xc and yc+ddy >= y >= yc and ddx*(y-yc) + ddy*(x-xc) <= ddy*ddx:
        return fi(1, x, y, xc, xc+ddx, xc, yc, yc, yc+ddy)
    if xc+ddx >= x >= xc and yc-ddy <= y <= yc:
        if ddx*(y+ddy-yc) + ddy*(x-xc) > ddy*ddx:
            return fi(1, x, y, xc, xc+ddx, xc+ddx, yc, yc, yc-ddy)
        else:
            return fi(1, x, y, xc, xc+ddx, xc, yc, yc-ddy, yc-ddy)
    if xc-ddx <= x <= xc and yc-ddy <= y <= yc and ddx*(y+ddy-yc) + ddy*(x+ddx-xc) >= ddy*ddx:
        return fi(1, x, y, xc, xc-ddx, xc, yc, yc, yc-ddx)
    if xc-ddx <= x <= xc and yc+ddy >= y >= yc:
        if ddx*(y-yc) + ddy*(x+ddx-xc) > ddy*ddx:
            return fi(1, x, y, xc, xc-ddx, xc, yc, yc+ddy, yc+ddy)
        else:
            return fi(1, x, y, xc, xc-ddx, xc-ddx, yc, yc, yc+ddy)
    return 0


def approximate(x, y, ddx, ddy, nodes, wsp):
    """funkcja przybliżająca rozwiązanie w całym zakresie"""
    result = 0
    for n in range(len(nodes)):
        result += wsp[n]*shape_line(x, y, nodes[n][0], nodes[n][1], ddx, ddy)
    return result


def plot(X, Y, Z, nnx, nny):
    """przedstawienie wyniku na wykresie"""
    fig = pyplot.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z)
    cset = ax.contourf(X, Y, Z, zdir='z', offset=-0.08, cmap=cm.coolwarm)
    p = PdfPages('mes{}na{}.pdf'.format(nnx, nny))
    p.savefig(fig)
    p.close()


def count_difference(xa, xb, ya, yb, ddx1, ddy1, ddx2, ddy2, nodes1, wsp1, nodes2, wsp2, nn):
    suma = 0
    norma = 0
    for i in np.arange(xa, xb+ddx1, ddx1):
        for j in np.arange(ya, yb+ddy1, ddy1):
            suma += (approximate(i, j, ddx1, ddy1, wezel1, wsp1) - approximate(i, j, ddx2, ddy2, wezel2, wsp2))**2
            norma += approximate(i, j, ddx1, ddy1, wezel1, wsp1)**2
    blad_bezwzgledny = math.sqrt(suma)/(nn)
    norma = math.sqrt(norma)
    blad_wzgledny = blad_bezwzgledny / norma
    return blad_bezwzgledny, blad_wzgledny


def plot_2D(x, y, nazwa, xlabel, ylabel):
    fig = pyplot.figure()
    pyplot.plot(x, y)
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    p = PdfPages('{}.pdf'.format(nazwa))
    p.savefig(fig)
    p.close()


if __name__ == "__main__":
    #określenie granic
    xa = 0
    xb = 2
    ya = 0
    yb = 1
    #określenie początkowej liczby podziałów
    nnx1 = 5
    nny1 = 4
    blad = []
    czas = []
    N = 26
    for test_number in range(N):
        ddx1 = float(xb - xa)/float(nnx1 - 1)
        ddy1 = float(yb - ya)/float(nny1 - 1)
        liczbawezlow1 = nnx1*nny1
        liczbatrojkatow1 = 2*(nnx1-1)*(nny1-1)
        started = time.time()
        wezel1 = set_nodes(xa, ya, nnx1, nny1, ddx1, ddy1)
        siatka1 = divide_into_triangles(liczbawezlow1, nnx1, nny1)
        main_matrix1, right_matrix1 = agregacja(wezel1, siatka1, liczbawezlow1, liczbatrojkatow1, ddx1, ddy1)
        wsp1 = linalg.solve(main_matrix1, right_matrix1)
        X = np.arange(xa, xb+ddx1, ddx1)
        Y = np.arange(ya, yb+ddy1, ddy1)
        X, Y = np.meshgrid(X, Y)
        Z = np.array([approximate(x, y, ddx1, ddy1, wezel1, wsp1)
            for x, y in zip(np.ravel(X), np.ravel(Y))])
        Z = Z.reshape(X.shape)
        plot(X, Y, Z, nnx1, nny1)
        if test_number > 0:
            Z_differ = np.array([(approximate(x, y, ddx1, ddy1, wezel1, wsp1) - approximate(x, y, ddx2, ddy2, wezel2, wsp2))
                for x, y in zip(np.ravel(X), np.ravel(Y))])
            b1, b2 = count_difference(xa, xb, ya, yb, ddx1, ddy1, ddx2, ddy2, wezel1, wsp1, wezel2, wsp2, nnx1*nny1)
            blad.append(b2)
        ended = time.time()
        czas.append(ended - started)
        nnx1 += 1
        nny1 += 1
        ddx2 = ddx1
        ddy2 = ddy1
        wezel2 = wezel1
        wsp2 = wsp1
    iksy = [(5+x)*(4+x) for x in range(N)]
    plot_2D(iksy[1:], blad, "blad", u"liczba podziałów", u"błąd")
    plot_2D(iksy, czas, "czas", u"liczba podziałów", "czas")
    # zobrazowanie funkcji kształtu
    x = [[0, 1], [0, 1]]
    y = [[0, 0], [1, 1]]
    z = [[0, 0], [0, 1]]
    plot(x, y, z, "shape", "line")
