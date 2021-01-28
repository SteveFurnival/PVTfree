#!/usr/bin/python


import sys

from PySide6.QtWidgets import QApplication, QMessageBox


# Create the application object

app = QApplication(sys.argv)


# Create a simple dialog box

msg_box = QMessageBox()

msg_box.setText("Hello World!")

msg_box.show()


sys.exit(msg_box.exec_())
