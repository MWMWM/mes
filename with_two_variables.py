#!/usr/bin/python2
# -*- coding: utf-8 -*-

from scipy import integrate, misc
import numpy as np
from numpy import linalg
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

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

class DifferentialEquation(object):
    def __init__(self, values):
        self.nnx = int(values['nnx']) * 10
        self.nny = int(values['nny']) * 10
        self.xa, self.xb = [int(x) for x in values['x'].split()]
        self.ya, self.yb = [int(y) for y in values['y'].split()]
        self.a = [int(aa) for aa in values['a'].split()]
        self.g = [int(gg) for gg in values['g'].split()]
        self.h = [int(hh) for hh in values['h'].split()]
        boundary = values['boundary'].split()
        self.xaa = self.check_validity(boundary[0])
        self.xbb = self.check_validity(boundary[1])
        self.yaa = self.check_validity(boundary[2])
        self.ybb = self.check_validity(boundary[3])
        self.ddx = float(self.xb - self.xa)/float(self.nnx - 1)
        self.ddy = float(self.yb - self.ya)/float(self.nny - 1)
        self.nodes_number = self.nnx * self.nny
        self.triangles_number = 2*(self.nnx-1)*(self.nny-1)

    def check_validity(self, val):
        try:
            val = int(val)
        except ValueError:
            if val != 'N':
                raise ValueError('''invalid literal - it should be 'N' or a number: '{}' provided'''.format(val))
        return val

    def set_nodes(self):
        """ustanowienie węzłów podziału, zwraca krotkę z ich współrzędnymi"""
        return [(self.xa + self.ddx*x, self.ya + self.ddy*y) 
                for x in range(self.nnx) for y in range(self.nny)]

    def divide_into_triangles(self):
        """podzielenie powierzchni na trójkąty,
        zwraca krotkę z numerami węzłów odpowiadającymi wierzchołkom trójkąta"""
        s = []
        for x in range(self.nnx-1):
            for y in range(self.nny-1):
                s.append([self.nny*x + y+1, self.nny*(x+1) + y, self.nny*x + y])
                s.append([self.nny*x + y+1, self.nny*(x+1) + y, self.nny*(x+1) + y+1])
        return s
    
    def agregate(self):
        """agregacja elementów z macierzy lokalnych do macierzy globalnych"""
        self.main_matrix = [[0 for i in range(self.nodes_number)] 
                for j in range(self.nodes_number)]
        self.right_matrix = [0 for i in range(self.nodes_number)]
        for ii in range(self.triangles_number):
            for jj in range(3):
                i = self.siatka[ii][jj]
                j = self.siatka[ii][(jj+1)%3]
                k = self.siatka[ii][(jj+2)%3]
                xp = self.nodes[i][0]
                yp = self.nodes[i][1]
                xq = self.nodes[j][0]
                yq = self.nodes[j][1]
                xr = self.nodes[k][0]
                yr = self.nodes[k][1]
                self.right_matrix[i] += self.local_right(xp, xq, xr, yp, yq, yr)
                self.main_matrix[i][i] += self.local_main(1, xp, xq, xr, yp, yq, yr)
                self.main_matrix[i][j] += self.local_main(2, xp, xq, xr, yp, yq, yr)
                self.main_matrix[i][k] += self.local_main(3, xp, xq, xr, yp, yq, yr)

    def consider_boundary_conditions(self):
        """uwzgędnienie warunków brzegowych"""
        for ii in range(self.nodes_number):
            if self.nodes[ii][0] == self.xa and self.xaa != 'N':
                self.right_matrix[ii] = self.xaa
            elif self.nodes[ii][0] == self.xb and self.xbb != 'N':
                self.right_matrix[ii] = self.xbb
            elif self.nodes[ii][1] == self.ya and self.yaa != 'N':
                self.right_matrix[ii] = self.yaa
            elif self.nodes[ii][1] == self.yb and self.ybb != 'N':
                self.right_matrix[ii] = self.ybb
            else:
                continue
            for i in range(self.nodes_number):
                self.main_matrix[ii][i] = 0
            self.main_matrix[ii][ii] = 1
    
    def local_right(self, xp, xq, xr, yp, yq, yr):
        """element wektoru prawych stron"""
        def inner_function(x, y):
            return self.prawa(x, y)*fi(1, x, y, xp, xq, xr, yp, yq, yr)
        return self.help_local(inner_function, xp, xq, xr, yp, yq, yr)

    def local_main(self, nr, xp, xq, xr, yp, yq, yr):
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
            result = - self.a[0]*pox(1, x, y)*pox(nr, x, y)
            result -= self.a[1]*pox(1, x, y)*poy(nr, x, y)
            result -= self.a[2]*poy(1, x, y)*poy(nr, x, y)
            result -= self.a[3]*poy(1, x, y)*pox(nr, x, y)
            rest = self.a[4] * pox(nr, x, y) + self.a[5] * poy(nr, x, y)
            rest += self.a[6] * self.lewa(x, y) * fi(nr, x, y, xp, xq, xr, yp, yq, yr)
            result += rest*fi(1, x, y, xp, xq, xr, yp, yq, yr)
            return result
        return self.help_local(inner_function, xp, xq, xr, yp, yq, yr)

    def help_local(self, inner_function, xp, xq, xr, yp, yq, yr):
        xmin = min(xp, xq, xr)
        ymin = min(yp, yq, yr)
        xmax = max(xp, xq, xr)
        ymax = max(yp, yq, yr)
        def upper(y):
            return straight_line(y, xmax, xmin, ymin, ymax)

        def lower(y):
            return xmin
        return integrate.dblquad(inner_function, ymin, ymax, lower, upper)[0]

    def prawa(self, x, y):
        """prawa strona równania"""
        return self.g[0]*x + self.g[1]*y + self.g[2]*x*y + self.g[3]*x*x + self.g[4]*y*y

    def lewa(self, x, y):
        return self.h[0]*x + self.h[1]*y + self.h[2]*x*y + self.h[3]*x*x + self.h[4]*y*y

    def shape_line(self, x, y, xc, yc):
        """funkcja przybliżająca rozwiązanie w otoczeniu węzła (xc, yc)"""
        if xc+self.ddx >= x >= xc and yc+self.ddy >= y >= yc and \
                self.ddx*(y-yc) + self.ddy*(x-xc) <= self.ddy*self.ddx:
            return fi(1, x, y, xc, xc+self.ddx, xc, yc, yc, yc+self.ddy)
        if xc+self.ddx >= x >= xc and yc-self.ddy <= y <= yc:
            if self.ddx*(y+self.ddy-yc) + self.ddy*(x-xc) > self.ddy*self.ddx:
                return fi(1, x, y, xc, xc+self.ddx, xc+self.ddx, yc, yc, yc-self.ddy)
            else:
                return fi(1, x, y, xc, xc+self.ddx, xc, yc, yc-self.ddy, yc-self.ddy)
        if xc-self.ddx <= x <= xc and yc-self.ddy <= y <= yc and \
                self.ddx*(y+self.ddy-yc) + self.ddy*(x+self.ddx-xc) >= self.ddy*self.ddx:
            return fi(1, x, y, xc, xc-self.ddx, xc, yc, yc, yc-self.ddx)
        if xc-self.ddx <= x <= xc and yc+self.ddy >= y >= yc:
            if self.ddx*(y-yc) + self.ddy*(x+self.ddx-xc) > self.ddy*self.ddx:
                return fi(1, x, y, xc, xc-self.ddx, xc, yc, yc+self.ddy, yc+self.ddy)
            else:
                return fi(1, x, y, xc, xc-self.ddx, xc-self.ddx, yc, yc, yc+self.ddy)
        return 0

    def approximate(self, x, y):
        """funkcja przybliżająca rozwiązanie w całym zakresie"""
        result = 0
        for n in range(len(self.nodes)):
            result += self.wsp[n]*self.shape_line(x, y, self.nodes[n][0], self.nodes[n][1])
        return result

    def solve_equation(self):
        self.nodes = self.set_nodes()
        self.siatka = self.divide_into_triangles()
        self.agregate()
        self.consider_boundary_conditions()
        self.wsp = linalg.solve(self.main_matrix, self.right_matrix)
        X = np.arange(self.xa, self.xb + self.ddx, self.ddx)
        Y = np.arange(self.ya, self.yb + self.ddy, self.ddy)
        X, Y = np.meshgrid(X, Y)
        Z = np.array([self.approximate(x, y)
            for x, y in zip(np.ravel(X), np.ravel(Y))])
        Z = Z.reshape(X.shape)
        z_min = min(self.wsp)
        plot(X, Y, Z, z_min)


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
    return val


def get_values_from_command_line():
    print "podaj granice poszukiwanego obszaru pierwszej zmiennej - 'xa xb'"
    x = raw_input()
    print "podaj granice poszukiwanego obszaru dla drugiej zmiennej - 'ya yb'"
    y = raw_input()
    print "podaj współczynniki równania - 'a1 a2 a3 a4 a5 a6 a7', gdzie\
            równanie wyraża się wzorem a1*d^2f(x,y)/dx^2 + a2*d^f(x,y)/dxdy + \
            a3*d^f(x,y)/dy^2 + a4*d^f(x,y)/dydx + a5*df(x,y)/dx + a6*df(x,y)/dy+\
            a7f(x,y) + g(x) = h(x)"
    a = raw_input()
    print "podaj współczynniki funkcji g(x) - 'a b c d e', gdzie\
            g(x) = ax+by+cxy+dx^2+ey^2"
    g = raw_input()
    print "podaj współczynniki funkcji h(x) - 'a b c d e', gdzie\
            h(x) = ax+by+cxy+dx^2+ey^2"
    h = raw_input()
    xaa = get_val('x', 'a')
    xbb = get_val('x', 'b')
    yaa = get_val('y', 'a')
    ybb = get_val('y', 'b')
    print "podaj jak bardzo dokładny ma być wynik (1-10)"
    nnx = nny = raw_input()
    return {'x': x, 'y': y, 'a': a, 'g': g, 'h': h, 'nnx': nnx, 'nny': nny,
            'boundary': xaa + ' ' + xbb + ' ' + yaa + ' ' + ybb}

if __name__ == "__main__":
    val = get_values_from_command_line()
    d = DifferentialEquation(val)
    d.solve_equation()
