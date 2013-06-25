#!/usr/bin/python2
# -*- coding: utf-8 -*-


import sys
from PyQt4 import QtGui
from with_two_variables import DifferentialEquation


def Pack(data):
    data = data.split()
    data = {'x': [data[0], data[1]], 'y': [data[2], data[3]], 
            'boundary': [data[i] for i in range(4, 8)], 
            'a': [data[i] for i in range(8, 15)],
            'g': [data[i] for i in range(15, 20)],
            'h': [data[i] for i in range(20, 25)],
           'nnx': data[25], 'nny': data[25]}
    return data


class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.widgets = MyWindow()
        self.setCentralWidget(self.widgets)
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        openFile = QtGui.QAction('Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.triggered.connect(self.PutDataFromFile)
        fileMenu.addAction(openFile)
        self.widgets.open_button.clicked.connect(self.PutDataFromFile)

        getResult = QtGui.QAction('Execute', self)
        getResult.setShortcut('Ctrl+E')
        getResult.triggered.connect(self.GetAnswer)
        fileMenu.addAction(getResult)
        self.widgets.answer_button.clicked.connect(self.GetAnswer)

        self.show()

    def PutDataFromFile(self):
        dialog_window = QtGui.QFileDialog()
        fname = dialog_window.getOpenFileName(self, directory='/home/')
        with open(fname, 'r') as f:
            data = f.read().split()
        self.widgets.PasteData(data)
        self.show()

    def GetAnswer(self):
        data = self.widgets.CopyData()
        data = Pack(data)
        diff = DifferentialEquation(data)
        diff.solve_equation()


class MyWindow(QtGui.QWidget):
    def __init__(self):
        super(MyWindow, self).__init__()
        self.setWindowTitle('Differential Equation Solver')
        self.grid = QtGui.QGridLayout()
        self.grid.setSpacing(2)
        self.AddLabel(u'Podaj granice poszukiwanego obszaru i odpowiadające im warunki brzegowe', 0, 0, 5)
        self.xa = self.AddLabelAndEditLine(u' - dolna granica I arg', 1)
        self.xb = self.AddLabelAndEditLine(u' - górna granica I arg', 2)
        self.ya = self.AddLabelAndEditLine(u' - dolna granica II arg', 3)
        self.yb = self.AddLabelAndEditLine(u' - górna granica II arg', 4)
        self.xaa = self.AddLabelAndEditLine(u' - warunek brzegowy', 1, 2)
        self.xbb = self.AddLabelAndEditLine(u' - warunek brzegowy', 2, 2)
        self.yaa = self.AddLabelAndEditLine(u' - warunek brzegowy', 3, 2)
        self.ybb = self.AddLabelAndEditLine(u' - warunek brzegowy', 4, 2)
        self.AddLabel(u'Podaj postać poszukiwanego równania', 5, 0, 5)
        self.a = self.AddEquationLine(6)
        self.g = self.AddCoefficientsLine('g', 7)
        self.h = self.AddCoefficientsLine('h', 8)
        self.AddLabel(u'Podaj jak dokładne ma być rozwiązanie', 9, 0, 5)
        self.AddLabelAndEditLine('- liczba w granicach [1-10]', 10)
        self.answer_button = self.AddAnswerButton()
        self.open_button = self.AddOpenButton()
        self.setLayout(self.grid)
        self.setGeometry(300, 300, 350, 300)

    def AddLabelAndEditLine(self, text, line_nb, column_nb=0):
        label = QtGui.QLabel(text)
        line_edit = QtGui.QLineEdit()
        self.grid.addWidget(label, line_nb, column_nb+1)
        self.grid.addWidget(line_edit, line_nb, column_nb)
        return line_edit

    def AddLabel(self, text, line_nb, column_nb=0, how_long=0):
        label = QtGui.QLabel(text)
        if how_long:
            self.grid.addWidget(label, line_nb, column_nb, 1, how_long)
        else:
            self.grid.addWidget(label, line_nb, column_nb)

    def AddEquationLine(self, line_nb):
        a1 = self.AddLabelAndEditLine(
                u'\u2202 ^2 f(x,y) / \u2202 x^2 + ', line_nb, 0)
        a2 = self.AddLabelAndEditLine(
                u'\u2202 ^2 f(x,y) / \u2202 x \u2202 y + ', line_nb, 2)
        a3 = self.AddLabelAndEditLine(
                u'\u2202 ^2 f(x,y) / \u2202 y^2 + ', line_nb, 4)
        a4 = self.AddLabelAndEditLine(
                u'\u2202 ^2 f(x,y) / \u2202 y \u2202 x + ', line_nb, 6)
        a5 = self.AddLabelAndEditLine(
                u'\u2202 f(x,y) / \u2202 x + ', line_nb, 8)
        a6 = self.AddLabelAndEditLine(
                u'\u2202 f(x,y) / \u2202 y + ', line_nb, 10)
        a7 = self.AddLabelAndEditLine(u' f(x,y) *g(x) ', line_nb, 12)
        a8 = self.AddLabel(u' = h(x),', line_nb, 14)
        return [a1, a2, a3, a4, a5, a6]

    def AddCoefficientsLine(self, func, line_nb):
        self.AddLabel('gdzie {}(x) = '.format(func), line_nb)
        a1 = self.AddLabelAndEditLine(u'x + ', line_nb, 2)
        a2 = self.AddLabelAndEditLine(u'y + ', line_nb, 4)
        a3 = self.AddLabelAndEditLine(u'xy + ', line_nb, 6)
        a4 = self.AddLabelAndEditLine(u'x^2 + ', line_nb, 8)
        a5 = self.AddLabelAndEditLine(u'y^2', line_nb, 10)
        return [a1, a2, a3, a4, a5]

    def AddOpenButton(self):
        b = QtGui.QPushButton(u'Wczytaj dane', self)
        b.setToolTip('Click to read data from .txt file')
        self.grid.addWidget(b, 11, 0, 1, 3)
        return b    

    def AddAnswerButton(self):
        b = QtGui.QPushButton(u'Narysuj rozwiązanie', self)
        b.setToolTip('Click to get the answer')
        self.grid.addWidget(b, 11, 4, 1, 3)
        return b

    def PasteData(self, data):
        i = 0
        for child in self.children():
            if isinstance(child, QtGui.QLineEdit):
                child.setText(str(data[i]))
                i +=1

    def CopyData(self):
        data = ''
        for child in self.children():
            if isinstance(child, QtGui.QLineEdit):
                data += str(child.text()) + ' '
        return data


def main():
    app = QtGui.QApplication(sys.argv)
    w = MainWindow()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
