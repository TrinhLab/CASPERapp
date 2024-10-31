from PyQt6.QtWidgets import (
    QMainWindow, QPushButton, QRadioButton, QComboBox, QPlainTextEdit,
    QProgressBar, QMenuBar, QMenu, QStackedWidget, QWidget, QVBoxLayout,
    QHBoxLayout, QLabel, QFrame, QTabWidget, QToolButton, QTabBar
)
from PyQt6.QtGui import QIcon, QAction, QFont, QCursor 
from PyQt6.QtCore import Qt, QPoint, pyqtSignal, QSize
from PyQt6 import uic, QtWidgets, QtCore, QtGui
from utils.ui import scale_ui, show_error
import os
from typing import Optional
from functools import partial
import qdarktheme

class CloseableTabWidget(QTabWidget):
    tab_closed = pyqtSignal(QWidget)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTabsClosable(False)
        self.tabCloseRequested.connect(self.closeTab)
        self._tabs = {}  # Dictionary to keep track of tab widgets
        self.tabBar().tabMoved.connect(self._handle_tab_moved)

    def closeTab(self, index):
        try:
            if self.count() > 1 and index != 0:
                widget = self.widget(index)
                if widget:
                    # Get tab text before removal
                    tab_text = self.tabText(index)
                    
                    # Clean up the controller if it exists
                    controller = getattr(widget, 'controller', None)
                    if controller and hasattr(controller, 'model') and hasattr(controller.model, 'cleanup'):
                        controller.model.cleanup()
                    
                    # Remove from tracking dictionary
                    if tab_text in self._tabs:
                        del self._tabs[tab_text]
                    
                    # Remove the tab
                    self.removeTab(index)
                    
                    # Emit signal before deletion
                    self.tab_closed.emit(widget)
                    
                    # Schedule widget for deletion
                    widget.deleteLater()
                    
                    # Update all remaining tabs
                    self._update_all_tabs()
        except Exception as e:
            print(f"Error closing tab: {e}")

    def addTab(self, widget, label):
        try:
            if widget and label:
                # Store widget reference with unique identifier
                tab_id = f"{label}_{id(widget)}"
                self._tabs[tab_id] = {
                    'widget': widget,
                    'label': label,
                    'close_button': None
                }
                
                # Add the tab
                index = super().addTab(widget, label)
                
                if index != 0:
                    # Create and setup close button
                    close_button = self._create_close_button(index, label)
                    self._tabs[tab_id]['close_button'] = close_button
                    self.tabBar().setTabButton(index, QTabBar.ButtonPosition.RightSide, close_button)
                
                return index
        except Exception as e:
            print(f"Error adding tab: {e}")
            return -1

    def _create_close_button(self, index, label):
        """Create a new close button for a tab"""
        close_button = QToolButton(self.tabBar())
        close_button.setObjectName(f"close_button_{label}")
        close_icon = self.style().standardIcon(QtWidgets.QStyle.StandardPixmap.SP_TitleBarCloseButton)
        close_button.setIcon(close_icon)
        close_button.setIconSize(QSize(16, 16))
        close_button.setAutoRaise(True)
        close_button.setStyleSheet("""
            QToolButton {
                border: none;
                padding: 0px;
            }
            QToolButton:hover {
                background: #c42b1c;
            }
        """)
        close_button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        close_button.setFixedSize(18, 18)
        close_button.clicked.connect(lambda checked, idx=index: self.safely_close_tab(idx))
        return close_button

    def safely_close_tab(self, index):
        """Safely handle tab closure with error checking"""
        try:
            if 0 <= index < self.count():
                current_widget = self.widget(index)
                if current_widget and index != 0:
                    self.closeTab(index)
        except Exception as e:
            print(f"Error in safely_close_tab: {e}")

    def _handle_tab_moved(self, from_index: int, to_index: int):
        """Handle tab movement and update close buttons"""
        try:
            self._update_all_tabs()
        except Exception as e:
            print(f"Error handling tab movement: {e}")

    def _update_all_tabs(self):
        """Update all tabs and their close buttons"""
        try:
            for i in range(1, self.count()):  # Skip index 0 (home tab)
                widget = self.widget(i)
                if widget:
                    label = self.tabText(i)
                    tab_id = f"{label}_{id(widget)}"
                    
                    # Create new close button if needed
                    if tab_id not in self._tabs or not self._tabs[tab_id].get('close_button'):
                        close_button = self._create_close_button(i, label)
                        self._tabs[tab_id] = {
                            'widget': widget,
                            'label': label,
                            'close_button': close_button
                        }
                        self.tabBar().setTabButton(i, QTabBar.ButtonPosition.RightSide, close_button)
                    else:
                        # Update existing close button's click connection
                        close_button = self._tabs[tab_id]['close_button']
                        close_button.clicked.disconnect()
                        close_button.clicked.connect(lambda checked, idx=i: self.safely_close_tab(idx))
        except Exception as e:
            print(f"Error updating tabs: {e}")

    def moveTab(self, from_index, to_index):
        """Override moveTab to safely handle tab movement"""
        try:
            if (0 <= from_index < self.count() and 
                0 <= to_index < self.count() and 
                from_index != 0 and 
                to_index != 0):
                
                super().moveTab(from_index, to_index)
                self._update_all_tabs()
                
        except Exception as e:
            print(f"Error moving tab: {e}")


class MainWindowView(QMainWindow):
    def __init__(self, global_settings):
        super().__init__()
        self.global_settings = global_settings
        self.logger = global_settings.get_logger()
        self._init_ui()
        self.oldPos = None

    def _init_ui(self) -> None:
        # Hide the window and disable updates during initialization
        self.hide()
        self.setUpdatesEnabled(False)
        try:
            # Calculate center position first
            screen = QtGui.QGuiApplication.primaryScreen()
            screen_geometry = screen.geometry()
            centerPoint = screen_geometry.center()
            
            # Load and initialize UI
            self._load_ui_file()
            self._init_window_properties()
            self._init_ui_elements()
            self.apply_theme()
            self._scale_ui()
            
            # Get final size
            final_size = self.size()
            
            # Calculate position only once
            x = centerPoint.x() - (final_size.width() // 2)
            y = centerPoint.y() - (final_size.height() // 2)
            
            # Set position and size in a single operation
            self.setGeometry(x, y, final_size.width(), final_size.height())
            
            # Re-enable updates and show window
            self.setUpdatesEnabled(True)
            self.show()
            self.repaint()  # Force immediate repaint
            
            self.logger.debug(f"Window initialized at position ({x}, {y}) with size {final_size}")
        except Exception as e:
            self._handle_init_error(e)

    def _load_ui_file(self) -> None:
        ui_file = os.path.join(self.global_settings.get_ui_dir_path(), "main_window.ui")
        uic.loadUi(ui_file, self)

    def _init_window_properties(self) -> None:
        """
        Creates a frameless, translucent window without a toolbar.
        """
        # Set window flags before other properties
        self.setWindowFlags(Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        # Hide toolbars
        toolbars = self.findChildren(QtWidgets.QToolBar)
        for toolbar in toolbars:
            toolbar.hide()
        # Ensure window starts hidden
        self.setVisible(False)

    def _init_ui_elements(self) -> None:
        # Initialize menu bar and custom title bar
        self._init_menuBar()
        self._init_custom_title_bar()

        # Create main widget and layout
        main_widget = QWidget()
        main_layout = QVBoxLayout(main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Add title bar and divider
        main_layout.addWidget(self.title_bar, 0)
        main_layout.addWidget(self._init_divider(), 0)

        # Create and set up tab container
        tab_container = QWidget()
        tab_container_layout = QVBoxLayout(tab_container)
        tab_container_layout.setContentsMargins(0, 0, 0, 0)
        tab_container_layout.setSpacing(0)

        # Use the _add_new_tab method to add a new tab
        # self._add_new_tab()

        # Initialize and add CloseableTabWidget
        self.tab_widget = CloseableTabWidget(self)
        self.tab_widget.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)
        self.tab_widget.setStyleSheet("""
            QTabWidget::pane {
                border: 1px solid #444444;
                padding: 10px;  /* Add padding here */
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
        widget = self.findChild(widget_type, name)
        if widget is None:
            self.global_settings.logger.warning(f"Widget '{name}' not found in UI file.")
        return widget
    
    def _init_custom_title_bar(self) -> None:
        self.title_bar = QWidget(self)
        self.title_bar.setObjectName("custom_title_bar")
        self.title_bar.setFixedHeight(32)  # Reduced height

        # Create the main horizontal layout for the title bar
        layout = QHBoxLayout(self.title_bar)
        layout.setContentsMargins(10, 0, 10, 0)  # Equal margins on left and right
        layout.setSpacing(5)  # Reduced spacing between items

        # ----- Window Control Buttons -----
        button_font = QFont("Arial", 8) 
        
        self.minimize_window_button = QPushButton("-", self.title_bar)
        self.minimize_window_button.setObjectName("minimize_window_button")
        self.minimize_window_button.setFixedSize(20, 20) 
        self.minimize_window_button.setFont(button_font)

        self.maximize_window_button = QPushButton("⛶", self.title_bar)
        self.maximize_window_button.setObjectName("maximize_window_button")
        self.maximize_window_button.setFixedSize(20, 20) 
        self.maximize_window_button.setFont(button_font)

        self.close_window_button = QPushButton("✕", self.title_bar)
        self.close_window_button.setObjectName("close_window_button")
        self.close_window_button.setFixedSize(20, 20) 
        self.close_window_button.setFont(button_font)

        # Apply a style to center the text vertically and horizontally
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

        # ----- Synchronize Widths of Left and Right Widgets -----
        # Adjust left_widget to calculate its required width
        left_widget.adjustSize()
        left_width = left_widget.sizeHint().width()

        # Set right_widget's fixed width to match left_widget's width
        right_widget.setFixedWidth(left_width)

        # ----- Title Label -----
        self.title_label = QLabel("CASPER", self.title_bar)
        self.title_label.setObjectName("title_label")
        self.title_label.setFont(QFont("Arial", 10, QFont.Weight.Bold)) 
        self.title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center the text in the label

        # Add Widgets to the Main Title Bar Layout
        layout.addWidget(left_widget)       
        layout.addStretch(1)                
        layout.addWidget(self.title_label) 
        layout.addStretch(1)                
        layout.addWidget(right_widget)     


    def _init_divider(self):
        divider = QFrame()
        divider.setObjectName("custom_divider")
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)
        return divider

    # def _add_new_tab(self):
    #     new_tab = QWidget()
    #     layout = QVBoxLayout(new_tab)
        
    #     label = QLabel("This is a new tab", new_tab)
    #     layout.addWidget(label)
        
    #     new_tab_button = QPushButton("Open New Tab", new_tab)
    #     new_tab_button.clicked.connect(self._add_new_tab)
    #     layout.addWidget(new_tab_button)
    #     tab_index = self.tab_widget.addTab(new_tab, f"Tab {self.tab_widget.count() + 1}")
    #     self.tab_widget.setCurrentIndex(tab_index)


    def _scale_ui(self):
        """Modified scale_ui to only handle sizing, not positioning"""
        try:
            screen = QtGui.QGuiApplication.primaryScreen()
            screen_geometry = screen.geometry()
            width = screen_geometry.width()
            height = screen_geometry.height()

            # Font scaling
            self.centralWidget().setStyleSheet(f"font: 12pt 'Arial';")

            if hasattr(self, 'title'):
                scaled_title_font_size = int(30 * (width / 1920))
                self.title.setStyleSheet(f"font: bold {scaled_title_font_size}pt 'Arial';")

            # Calculate size only
            scaledWidth = int((width * 575) / 1920)
            scaledHeight = int((height * 400) / 1080)

            # Ensure minimum size
            self.adjustSize()
            currentWidth = self.size().width()
            currentHeight = self.size().height()

            if scaledHeight < currentHeight:
                scaledHeight = currentHeight
            if scaledWidth < currentWidth:
                scaledWidth = currentWidth

            # Only resize, don't reposition
            self.resize(scaledWidth, scaledHeight)

        except Exception as e:
            self.logger.error(f"Error in _scale_ui: {str(e)}")

    def _handle_init_error(self, e: Exception) -> None:
        error_msg = f"Error initializing MainWindowView: {str(e)}"
        self.global_settings.logger.error(error_msg, exc_info=True)
        show_error(self.global_settings, "Initialization Error", error_msg)
        raise

    def update_theme_icon(self) -> None:
        icon_name = "dark_mode.png" if self.global_settings.get_theme() == "dark" else "light_mode.png"
        icon_path = os.path.join(self.global_settings.get_assets_dir_path(), icon_name)
        icon = QIcon(icon_path)
        self.theme_toggle_button.setIcon(icon)
        self.theme_toggle_button.setIconSize(QtCore.QSize(16, 16)) 

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.oldPos = event.globalPosition().toPoint()

    def mouseMoveEvent(self, event):
        if self.oldPos:
            delta = event.globalPosition().toPoint() - self.oldPos
            self.move(self.x() + delta.x(), self.y() + delta.y())
            self.oldPos = event.globalPosition().toPoint()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.oldPos = None

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.logger.debug(f"Window resized. New size: {self.size()}")

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

        # Get the current theme
        current_theme = self.global_settings.get_theme()
        theme = themes["dark"] if current_theme == "dark" else themes["light"]

        # Apply the selected theme using qdarktheme
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

    def show_window(self) -> None:
        """Shows the window without repositioning"""
        self.show()
        self.repaint()






