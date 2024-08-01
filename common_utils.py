from PyQt5 import QtWidgets

def show_message(fontSize, icon, title, message, button=QtWidgets.QMessageBox.StandardButton.Close):
    msgBox = QtWidgets.QMessageBox()
    msgBox.setStyleSheet(f"font: {fontSize}pt 'Arial'")
    msgBox.setIcon(icon)
    msgBox.setWindowTitle(title)
    msgBox.setText(message)
    msgBox.addButton(button)
    msgBox.exec()