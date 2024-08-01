from PyQt5 import QtWidgets, QtGui, QtCore
import math
from error_handling import show_error
import GlobalSettings

def scale_ui(window, base_width=1920, base_height=1080, font_size=12, header_font_size=30):
    try:
        
        window.repaint()
        QtWidgets.QApplication.processEvents()

        screen = window.screen()
        width = screen.geometry().width()
        height = screen.geometry().height()

        # Font scaling
        window.centralWidget().setStyleSheet(f"font: {font_size}pt 'Arial';")
        window.menuBar().setStyleSheet(f"font: {font_size}pt 'Arial';")

        # Header scaling
        if hasattr(window, 'label_8'):  # Check if the window has 'label_8' attribute
            window.label_8.setStyleSheet(f"font: bold {header_font_size}pt 'Arial';")

        window.adjustSize()

        currentWidth = window.size().width()
        currentHeight = window.size().height()

        # Set logo if it's a StartupWindow
        if isinstance(window, QtWidgets.QMainWindow) and hasattr(window, 'logo'):
            window.logo.setPixmap(QtGui.QPixmap(GlobalSettings.appdir + "assets/CASPER-logo.jpg"))

        # Window resize and center
        scaledWidth = int((width * 1150) / base_width)
        scaledHeight = int((height * 650) / base_height)

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
        # Log the error or handle it as needed
        print(f"Error centering window: {e}")
    