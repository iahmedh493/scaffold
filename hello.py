
import sys

from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtWidgets import *

from MainWindow_ui import Ui_MainWindow
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import uic
import subprocess
import os




class MyMainWindow(QMainWindow):

    def __init__(self):
        super(MyMainWindow, self).__init__()

        # uic.loadUi('ui/home.ui', self)
        # uic.loadUi('toyota.ui', self)
        self.main_ui = Ui_MainWindow()
        self.main_ui.setupUi(self)
        # ******************************************
        # initialize elements

        self.start_btn = self.main_ui.pushButton
        self.start_btn.clicked.connect(self.run_deseq2)

    
    def run_deseq2(self):#(list1, path_):
        # os.chdir(path_)
        command = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
        script = "data/deseq2_run.R"
    
        cmd = [command, script] #+ args1
        subprocess.check_output(cmd, universal_newlines=True)
        print("complete")
       


if __name__ == "__main__":

    app = QApplication(sys.argv)

    window = MyMainWindow()
    window.show()

    sys.exit(app.exec_())







