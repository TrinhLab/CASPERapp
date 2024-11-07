from PyQt6.QtWidgets import (
    QMainWindow, QPushButton, QWidget, QVBoxLayout,
    QHBoxLayout, QLabel, QFrame, 
)
from PyQt6.QtGui import QIcon, QAction
from PyQt6.QtCore import Qt
from PyQt6 import uic, QtWidgets, QtCore, QtGui
from utils.ui import show_error
from utils.LoggingMixin import LoggingMixin
import os
from typing import Optional
import qdarktheme
from views.CloseableTabWidget import CloseableTabWidget

class MainWindowView(QMainWindow, LoggingMixin):
    def __init__(self, global_settings):
        QMainWindow.__init__(self)
        LoggingMixin.__init__(self)
        self.settings = global_settings
        self._init_ui()
        self.oldPos = None

    def _init_ui(self) -> None:
        self.log_method_call("_init_ui")
        
        screen = QtGui.QGuiApplication.primaryScreen()
        screen_geometry = screen.geometry()
        centerPoint = screen_geometry.center()

        # Load UI file
        uic.loadUi(self.settings.get_ui_dir_path() + "/main_window.ui", self)
        self._init_window_properties()
        self._init_ui_elements()
        self.apply_theme()
        
        # Calculate and set position
        final_size = self.size()
        x = centerPoint.x() - (final_size.width() // 2)
        y = centerPoint.y() - (final_size.height() // 2)
        
        self.setGeometry(x, y, final_size.width(), final_size.height())
        self.setUpdatesEnabled(True)
        self.show()
        self.repaint()
        
        self.log_debug(f"Window initialized at position ({x}, {y}) with size {final_size}")

    def _init_window_properties(self) -> None:
        self.setWindowFlags(Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)

        toolbars = self.findChildren(QtWidgets.QToolBar)
        for toolbar in toolbars:
            toolbar.hide()

    def _init_ui_elements(self) -> None:
        self._init_menuBar()
        self._init_custom_title_bar()

        main_widget = QWidget()
        main_layout = QVBoxLayout(main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        main_layout.addWidget(self.title_bar, 0)
        main_layout.addWidget(self._init_divider(), 0)

        # Create and set up tab container
        tab_container = QWidget()
        tab_container_layout = QVBoxLayout(tab_container)
        tab_container_layout.setContentsMargins(0, 0, 0, 0)
        tab_container_layout.setSpacing(0)

        # Initialize and add CloseableTabWidget
        self.tab_widget = CloseableTabWidget(self)
        self.tab_widget.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, 
                                    QtWidgets.QSizePolicy.Policy.Expanding)
        self.tab_widget.setStyleSheet("""
            QTabWidget::pane {
                border: 1px solid #444444;
                padding: 10px;
            }
        """)
        tab_container_layout.addWidget(self.tab_widget)
        main_layout.addWidget(tab_container, 1)

        # Configure tab widget properties
        self.tab_widget.setMovable(True)

        # Set main widget as central widget
        self.setCentralWidget(main_widget)

    def _init_menuBar(self) -> None:
        self.action_change_database_directory = self._find_widget("actChangeDatabaseDirectory", QAction)
        self.action_open_genome_browser = self._find_widget("actOpenGenomeBrowser", QAction)
        self.action_open_repository = self._find_widget("actionGoToCASPERRepository", QAction)
        self.action_open_NCBI_BLAST = self._find_widget("actionGoToNCBIBLAST", QAction)
        self.action_open_NCBI = self._find_widget("actGoToNCBI", QAction)

    def _find_widget(self, name: str, widget_type: type) -> Optional[QtWidgets.QWidget]:
        """Find a widget by name and type"""
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.log_warning(f"Widget '{name}' not found in UI file")
        return widget

    def _init_custom_title_bar(self) -> None:
        self.title_bar = QWidget(self)
        self.title_bar.setObjectName("custom_title_bar")
        self.title_bar.setFixedHeight(32) 

        # Create the main horizontal layout for the title bar
        layout = QHBoxLayout(self.title_bar)
        layout.setContentsMargins(10, 0, 10, 0) 
        layout.setSpacing(5) 

        # ----- Window Control Buttons -----
        self.minimize_window_button = QPushButton("-", self.title_bar)
        self.minimize_window_button.setObjectName("minimize_window_button")
        self.minimize_window_button.setFixedSize(20, 20) 

        self.maximize_window_button = QPushButton("⛶", self.title_bar)
        self.maximize_window_button.setObjectName("maximize_window_button")
        self.maximize_window_button.setFixedSize(20, 20) 

        self.close_window_button = QPushButton("✕", self.title_bar)
        self.close_window_button.setObjectName("close_window_button")
        self.close_window_button.setFixedSize(20, 20) 

        button_style = """
        QPushButton {
            padding: 0px;
            margin: 0px;
            line-height: 20px;
            text-align: center;
        }
        """
        self.minimize_window_button.setStyleSheet(button_style)
        self.maximize_window_button.setStyleSheet(button_style)
        self.close_window_button.setStyleSheet(button_style)

        # ----- Left Widget (Minimize, Maximize, Close) -----
        left_widget = QWidget()
        left_layout = QHBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(5)
        left_layout.addWidget(self.close_window_button)
        left_layout.addWidget(self.minimize_window_button)
        left_layout.addWidget(self.maximize_window_button)

        # ----- Theme Toggle Button -----
        self.theme_toggle_button = QPushButton(self.title_bar)
        self.theme_toggle_button.setObjectName("theme_toggle_button")
        self.theme_toggle_button.setFixedSize(20, 20) 
        self.theme_toggle_button.setStyleSheet("border: none;")
        self.update_theme_icon()

        # ----- Right Widget (Theme Toggle + Stretch) -----
        right_widget = QWidget()
        right_layout = QHBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(5)

        # Add a stretch to push the toggle button to the far right within right_widget
        right_layout.addStretch()
        right_layout.addWidget(self.theme_toggle_button)

        # Adjust left_widget to calculate its required width
        left_widget.adjustSize()
        left_width = left_widget.sizeHint().width()

        # Set right_widget's fixed width to match left_widget's width
        right_widget.setFixedWidth(left_width)

        self.title_label = QLabel("CASPER", self.title_bar)
        self.title_label.setObjectName("title_label")
        self.title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center the text in the label

        # Add Widgets to the Main Title Bar Layout
        layout.addWidget(left_widget)       
        layout.addStretch(1)                
        layout.addWidget(self.title_label) 
        layout.addStretch(1)                
        layout.addWidget(right_widget)     

        # Add mouse tracking to the title bar
        self.title_bar.mousePressEvent = self.mousePressEvent
        self.title_bar.mouseMoveEvent = self.mouseMoveEvent
        self.title_bar.mouseReleaseEvent = self.mouseReleaseEvent
        self.title_bar.setMouseTracking(True)

    def _init_divider(self):
        divider = QFrame()
        divider.setObjectName("custom_divider")
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)
        return divider

    def _handle_init_error(self, e: Exception) -> None:
        error_msg = f"Error initializing MainWindowView: {str(e)}"
        self.log_error("_init_ui", e)
        show_error(self.settings, "Initialization Error", error_msg)
        raise

    def update_theme_icon(self) -> None:
        try:
            icon_name = "dark_mode.png" if self.settings.get_theme() == "dark" else "light_mode.png"
            icon_path = os.path.join(self.settings.get_assets_dir_path(), icon_name)
            icon = QIcon(icon_path)
            self.theme_toggle_button.setIcon(icon)
            self.theme_toggle_button.setIconSize(QtCore.QSize(16, 16))
        except Exception as e:
            self.log_error("update_theme_icon", e)
            show_error(self.settings, "Theme Error", "Failed to update theme icon")

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.log_debug(f"Window resized. New size: {self.size()}")

    def apply_theme(self):
        themes = {
            "dark": {
                "bg_color": "#2b2b2b",
                "fg_color": "#ffffff",
                "button_bg_color": "#3a3a3a",
                "button_border_color": "#5a5a5a",
                "button_hover_bg_color": "#4a4a4a",
                "input_bg_color": "#3a3a3a",
                "input_border_color": "#5a5a5a",
                "menu_bg_color": "#2b2b2b",
                "menu_item_hover_bg_color": "#3a3a3a",
                "tab_bg_color": "#2b2b2b",
                "tab_border_color": "#444444",
                "tab_selected_border_color": "#51b85e",
                "tab_hover_bg_color": "#3b3b3b",
                "divider_color": "#444444"
            },
            "light": {
                "bg_color": "#f0f0f0",
                "fg_color": "#000000",
                "button_bg_color": "#e0e0e0",
                "button_border_color": "#c0c0c0",
                "button_hover_bg_color": "#d0d0d0",
                "input_bg_color": "#ffffff",
                "input_border_color": "#c0c0c0",
                "menu_bg_color": "#f0f0f0",
                "menu_item_hover_bg_color": "#e0e0e0",
                "tab_bg_color": "#f0f0f0",
                "tab_border_color": "#c0c0c0",
                "tab_selected_border_color": "#51b85e",
                "tab_hover_bg_color": "#e0e0e0",
                "divider_color": "#c0c0c0"
            }
        }

        current_theme = self.settings.get_theme()
        theme = themes["dark"] if current_theme == "dark" else themes["light"]
        qdarktheme.setup_theme(current_theme)

        # Set the stylesheet
        self.setStyleSheet(f"""
            QWidget {{ background-color: {theme['bg_color']}; color: {theme['fg_color']}; }}
            QPushButton {{ background-color: {theme['button_bg_color']}; border: 1px solid {theme['button_border_color']}; }}
            QPushButton:hover {{ background-color: {theme['button_hover_bg_color']}; }}
            QLineEdit, QTextEdit, QPlainTextEdit {{ background-color: {theme['input_bg_color']}; border: 1px solid {theme['input_border_color']}; }}
            QComboBox {{ background-color: {theme['input_bg_color']}; border: 1px solid {theme['input_border_color']}; }}
            QMenuBar {{ background-color: {theme['menu_bg_color']}; }}
            QMenuBar::item:selected {{ background-color: {theme['menu_item_hover_bg_color']}; }}
            QMenu {{ background-color: {theme['menu_bg_color']}; }}
            QMenu::item:selected {{ background-color: {theme['menu_item_hover_bg_color']}; }}
            QFrame#custom_divider {{ border-bottom: 1px solid {theme['divider_color']}; }}
        """)

        # Set the tab widget stylesheet
        self.tab_widget.setStyleSheet(f"""
            QTabWidget::pane {{
                border-top: 1px solid {theme['divider_color']};
                top: -1px;
            }}
            QTabWidget::tab-bar {{
                alignment: left;
            }}
            QTabBar::tab {{
                background: {theme['tab_bg_color']};
                color: {theme['fg_color']};
                padding: 5px 5px 3px 5px;
                border: none;
                border-bottom: 1px solid {theme['tab_border_color']};
                border-right: 1px solid {theme['tab_border_color']};
                margin: 0px;
            }}
            QTabBar::tab:first {{
                border-left: none;
            }}
            QTabBar::tab:selected {{
                background: {theme['tab_bg_color']};
                border-top: 2px solid {theme['tab_selected_border_color']};
                border-bottom: 1px solid {theme['tab_bg_color']};
                color: {theme['fg_color']};
                padding-top: 3px;
            }}
            QTabBar::tab:hover {{
                background: {theme['tab_hover_bg_color']};
            }}
        """)

    def mousePressEvent(self, event):
        """Handle mouse press events for window dragging"""
        if event.button() == Qt.MouseButton.LeftButton:
            self.oldPos = event.globalPosition().toPoint()

    def mouseMoveEvent(self, event):
        """Handle mouse move events for window dragging"""
        if self.oldPos is not None:
            delta = event.globalPosition().toPoint() - self.oldPos
            self.move(self.x() + delta.x(), self.y() + delta.y())
            self.oldPos = event.globalPosition().toPoint()

    def mouseReleaseEvent(self, event):
        """Handle mouse release events for window dragging"""
        if event.button() == Qt.MouseButton.LeftButton:
            self.oldPos = None