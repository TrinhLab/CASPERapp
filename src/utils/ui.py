from PyQt6 import QtWidgets, QtGui, QtCore
import traceback
import models.GlobalSettings as GlobalSettings


def show_message(title, message, fontSize=12, icon=QtWidgets.QMessageBox.Icon.Information, button=QtWidgets.QMessageBox.StandardButton.Close):
    try:
        msgBox = QtWidgets.QMessageBox()
        msgBox.setStyleSheet(f"font: {fontSize}pt 'Arial'")
        msgBox.setIcon(icon)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.addButton(button)
        msgBox.exec()
    except Exception as e:
        print(f"Error showing message: {e}") 

def show_error(settings, message, e):
    try:
        logger = settings.get_logger()
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
        print(f"Error showing error message: {e}") 

    exit(-1)

def scale_ui(window, base_width=1920, base_height=1080, font_size=12, header_font_size=30, custom_scale_width=None, custom_scale_height=None):
    try:
        window.repaint()
        QtWidgets.QApplication.processEvents()

        # Get the primary screen
        screen = QtGui.QGuiApplication.primaryScreen()
        
        # Get the geometry of the screen
        screen_geometry = screen.geometry()
        width = screen_geometry.width()
        height = screen_geometry.height()

        # Font scaling
        window.centralWidget().setStyleSheet(f"font: {font_size}pt 'Arial';")

        if hasattr(window, 'title'):
            scaled_title_font_size = int(header_font_size * (width / base_width))
            window.title.setStyleSheet(f"font: bold {scaled_title_font_size}pt 'Arial';")

        window.adjustSize()

        currentWidth = window.size().width()
        currentHeight = window.size().height()

        # Window resize and center
        scaledWidth = int((width * (custom_scale_width if custom_scale_width else 1150)) / base_width)
        scaledHeight = int((height * (custom_scale_height if custom_scale_height else 650)) / base_height)

        if scaledHeight < currentHeight:
            scaledHeight = currentHeight
        if scaledWidth < currentWidth:
            scaledWidth = currentWidth

        # Center the window on the screen
        centerPoint = screen_geometry.center()
        x = centerPoint.x() - (scaledWidth // 2)
        y = centerPoint.y() - (scaledHeight // 2)
        window.setGeometry(x, y, scaledWidth, scaledHeight)

        window.repaint()
        QtWidgets.QApplication.processEvents()

    except Exception as e:
        print(f"Error in scale_ui: {e}") 

def center_ui(window):
    try:
        window.repaint()
        QtWidgets.QApplication.processEvents()

        # Get the dimensions of the window
        width = window.width()
        height = window.height()

        # Get the primary screen
        screen = QtGui.QGuiApplication.primaryScreen()
        
        # Get the geometry of the screen
        screen_geometry = screen.geometry()
        centerPoint = screen_geometry.center()

        # Calculate new x and y coordinates
        x = centerPoint.x() - (width // 2)
        y = centerPoint.y() - (height // 2)

        # Set the new geometry for the window
        window.setGeometry(x, y, width, height)
        window.repaint()
    except Exception as e:
        print(f"Error centering window: {e}") 

def position_window(new_window, parent_window=None):
    # Check if the window is already visible and active
    if new_window.view.isVisible() and new_window.view.isActiveWindow():
        # If the window is already visible and active, just ensure it's in the foreground
        new_window.view.raise_()
        new_window.view.activateWindow()
        QtWidgets.QApplication.setActiveWindow(new_window.view)
        return

    if parent_window is None:
        parent_window = QtWidgets.QApplication.activeWindow()
    
    if parent_window:
        if hasattr(parent_window, 'last_position') and parent_window.last_position:
            new_window.view.move(parent_window.last_position)
        else:
            parent_geo = parent_window.geometry()
            new_window.view.move(parent_geo.x() + 50, parent_geo.y() + 50)
    else:
        center_ui(new_window)
    
    new_window.view.show()
    new_window.view.raise_()
    new_window.view.activateWindow()
    QtWidgets.QApplication.setActiveWindow(new_window.view)
    
    # Force the window to be active and in the foreground
    new_window.view.setWindowState(new_window.view.windowState() & ~QtCore.Qt.WindowState.WindowMinimized | QtCore.Qt.WindowState.WindowActive)
    new_window.view.raise_()
    new_window.view.activateWindow()
