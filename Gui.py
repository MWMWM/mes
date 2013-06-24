#!/usr/bin/python2
# -*- coding: utf-8 -*-

import sys
from PyQt4 import QtGui
from with_two_variables import DifferentialEquation


class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.form_widget = MyWindow()
        self.setCentralWidget(self.form_widget)
        openFile = QtGui.QAction(QtGui.QIcon('open.png'), 'Solve using .txt file', self)
        openFile.setShortcut('Ctrl+S')
        openFile.triggered.connect(self.GetAnswer)
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)
        self.show()

    def ChangeToDictFormat(self, data):
        return {'x': data[0], 'y': data[1], 'a': data[2], 'g': data[3], 'h': data[4],
                'boundary': data[5], 'nnx': data[6], 'nny': data[6]}

    def GetAnswer(self):
        dialog_window = QtGui.QFileDialog()
        fname = dialog_window.getOpenFileName(self, 'Open file', '/home/m/Documents/studyjne/mes')
        with open(fname, 'r') as f:        
            data = f.read().split('\n')
            data = self.ChangeToDictFormat(data)
            diff = DifferentialEquation(data)
        dialog_window.close()
        diff.solve_equation()


class MyWindow(QtGui.QWidget):
    def __init__(self):
        super(MyWindow, self).__init__()
        self.setWindowTitle('Differential Equation Solver')
        self.line_nb = 0
        self.grid = QtGui.QGridLayout()
        self.grid.setSpacing(10)
        self.xa, self.xxa = self.AddLabelAndEditLine('xa')
        self.xb, self.xxb = self.AddLabelAndEditLine('xb')
        self.ya, self.yya = self.AddLabelAndEditLine('ya')
        self.yb, self.yyb = self.AddLabelAndEditLine('yb')
        self.a1, self.aa1 = self.AddLabelAndEditLine(u'współczynnik a1')
        self.setLayout(self.grid)
        self.setGeometry(300, 300, 350, 300)
        self.AddAnswerButton()


    def AddLabelAndEditLine(self, text, where=0):
        label = QtGui.QLabel(text)
        line_edit = QtGui.QLineEdit()
        self.grid.addWidget(label, self.line_nb, where)
        self.grid.addWidget(line_edit, self.line_nb, where+1)
        self.line_nb += 1
        return label, line_edit

    def AddAnswerButton(self):
        b = QtGui.QPushButton('Answer', self)
        b.setToolTip('Click to get the answer')
        b.clicked.connect(self.GetAnswer)
        self.grid.addWidget(b, self.line_nb, 0)

    def GetAnswer(self):
        pass
        #DifferentialEquation.solve_equation()

def main():
    app = QtGui.QApplication(sys.argv)
    w = MainWindow()
    sys.exit(app.exec_())
        
        
if __name__ == '__main__':
    main()

