from PyQt5 import QtWidgets, QtGui, QtCore
import traceback
import math
import models.GlobalSettings as GlobalSettings

logger = GlobalSettings.logger

def show_message(fontSize, icon, title, message, button=QtWidgets.QMessageBox.StandardButton.Close):
    try:
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet(f"font: {fontSize}pt 'Arial'")
        msgBox.setIcon(icon)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.addButton(button)
        msgBox.exec()
    except Exception as e:
        print(e)

def show_error(message, e):
    try:
        logger.critical(message)
        logger.critical(e)
        logger.critical(traceback.format_exc())

        show_message(
            fontSize=12,
            icon=QtWidgets.QMessageBox.Icon.Critical,
            title="Fatal Error",
            message=f"Fatal Error:\n{str(e)}\n\nFor more information on this error, look at CASPER.log in the application folder."
        )
    except Exception as e:
        print(e)

    exit(-1)

def scale_ui(window, base_width=1920, base_height=1080, font_size=12, header_font_size=30, custom_scale_width=None, custom_scale_height=None):
    try:
        window.repaint()
        QtWidgets.QApplication.processEvents()

        screen = window.screen()
        width = screen.geometry().width()
        height = screen.geometry().height()

        # Font scaling
        window.centralWidget().setStyleSheet(f"font: {font_size}pt 'Arial';")
        if hasattr(window, 'menuBar'):
            window.menuBar().setStyleSheet(f"font: {font_size}pt 'Arial';")

        # Header scaling
        # if hasattr(window, 'label_8'):
            # window.label_8.setStyleSheet(f"font: bold {header_font_size}pt 'Arial';")

        if hasattr(window, 'title'):
            scaled_title_font_size = int(header_font_size * (width / base_width))
            window.title.setStyleSheet(f"font: bold {scaled_title_font_size}pt 'Arial';")

        window.adjustSize()

        currentWidth = window.size().width()
        currentHeight = window.size().height()

        # Set logo if it's a StartupWindow
        if isinstance(window, QtWidgets.QMainWindow) and hasattr(window, 'logo'):
            window.logo.setPixmap(QtGui.QPixmap(GlobalSettings.appdir + "assets/CASPER-logo.jpg"))

        # Window resize and center
        scaledWidth = int((width * (custom_scale_width if custom_scale_width else 1150)) / base_width)
        scaledHeight = int((height * (custom_scale_height if custom_scale_height else 650)) / base_height)

        if scaledHeight < currentHeight:
            scaledHeight = currentHeight
        if scaledWidth < currentWidth:
            scaledWidth = currentWidth

        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        x = centerPoint.x() - (math.ceil(scaledWidth / 2))
        y = centerPoint.y() - (math.ceil(scaledHeight / 2))
        window.setGeometry(x, y, scaledWidth, scaledHeight)

        window.repaint()
        QtWidgets.QApplication.processEvents()

        return

    except Exception as e:
        show_error(f"Error in scaleUI() in {type(window).__name__}.", e)

def center_ui(window):
    try:
        window.repaint()
        QtWidgets.QApplication.processEvents()

        # Get the dimensions of the window
        width = window.width()
        height = window.height()

        # Get the screen number based on the cursor's current position
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()

        # Calculate new x and y coordinates
        x = centerPoint.x() - (width // 2)
        y = centerPoint.y() - (height // 2)

        # Set the new geometry for the window
        window.setGeometry(x, y, width, height)
        window.repaint()
    except Exception as e:
        print(f"Error centering window: {e}")
    