from PyQt6.QtWidgets import QTabWidget, QTabBar, QToolButton, QWidget
from PyQt6.QtCore import pyqtSignal, QSize, Qt
from PyQt6.QtGui import QCursor
from PyQt6 import QtWidgets
import logging

class CloseableTabWidget(QTabWidget):
    tab_closed = pyqtSignal(QWidget)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTabsClosable(False)
        self.tabCloseRequested.connect(self.closeTab)
        self._tabs = {}  # Dictionary to keep track of tab widgets
        self.tabBar().tabMoved.connect(self._handle_tab_moved)
        self.logger = logging.getLogger(__name__)

    def closeTab(self, index):
        """Close a tab at the given index"""
        self.logger.debug(f"Attempting to close tab at index {index}")
        
        if not (self.count() > 1 and index != 0):
            self.logger.debug("Tab closure conditions not met")
            return
            
        widget = self.widget(index)
        if not widget:
            self.logger.warning(f"No widget found at index {index}")
            return
            
        # Critical operations need try-catch
        try:
            tab_text = self.tabText(index)
            
            # Cleanup controller if exists
            controller = getattr(widget, 'controller', None)
            if controller and hasattr(controller, 'model') and hasattr(controller.model, 'cleanup'):
                controller.model.cleanup()
            
            # Remove from tracking and emit signal
            if tab_text in self._tabs:
                del self._tabs[tab_text]
            
            self.removeTab(index)
            self.tab_closed.emit(widget)
            widget.deleteLater()
            self._update_all_tabs()
            
            self.logger.debug(f"Successfully closed tab '{tab_text}'")
        except Exception as e:
            self.logger.error(f"Failed to close tab: {e}", exc_info=True)
            raise

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
                
                # Create and setup close button for all tabs except Home and Startup
                if label not in ["Home", "Startup"]:
                    close_button = self._create_close_button(index, label)
                    self._tabs[tab_id]['close_button'] = close_button
                    self.tabBar().setTabButton(index, QTabBar.ButtonPosition.RightSide, close_button)
                
                return index
        except Exception as e:
            self.logger.error(f"Error adding tab: {e}")
            return -1

    def _create_close_button(self, index, label):
        """Create a new close button for a tab"""
        try:
            close_button = QToolButton(self.tabBar())
            close_button.setObjectName(f"close_button_{label}")
            
            # Always use a custom close icon style
            close_button.setText("×")  # Using multiplication symbol as close icon
            close_button.setStyleSheet("""
                QToolButton {
                    border: none;
                    padding: 0px;
                    color: #666666;
                    background: transparent;
                    font-size: 16px;
                    font-weight: bold;
                }
                QToolButton:hover {
                    color: #ffffff;
                    background: #c42b1c;
                }
            """)
            
            close_button.setAutoRaise(True)
            close_button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
            close_button.setFixedSize(20, 20)
            close_button.clicked.connect(lambda checked, idx=index: self.safely_close_tab(idx))
            
            return close_button
            
        except Exception as e:
            self.logger.error(f"Error creating close button: {e}")
            raise

    def safely_close_tab(self, index):
        """Safely handle tab closure with error checking"""
        try:
            if 0 <= index < self.count():
                current_widget = self.widget(index)
                if current_widget and index != 0:
                    self.closeTab(index)
        except Exception as e:
            self.logger.error(f"Error in safely_close_tab: {e}")

    def _handle_tab_moved(self, from_index: int, to_index: int):
        """Handle tab movement and update close buttons"""
        try:
            self._update_all_tabs()
        except Exception as e:
            self.logger.error(f"Error handling tab movement: {e}")

    def _update_all_tabs(self):
        """Update all tabs and their close buttons"""
        try:
            for i in range(self.count()):
                widget = self.widget(i)
                if widget:
                    label = self.tabText(i)
                    tab_id = f"{label}_{id(widget)}"
                    
                    # Create new close button if needed and if not Home or Startup tab
                    if (label not in ["Home", "Startup"] and 
                        (tab_id not in self._tabs or not self._tabs[tab_id].get('close_button'))):
                        close_button = self._create_close_button(i, label)
                        self._tabs[tab_id] = {
                            'widget': widget,
                            'label': label,
                            'close_button': close_button
                        }
                        self.tabBar().setTabButton(i, QTabBar.ButtonPosition.RightSide, close_button)
                    elif label not in ["Home", "Startup"] and tab_id in self._tabs:
                        # Update existing close button's click connection
                        close_button = self._tabs[tab_id]['close_button']
                        close_button.clicked.disconnect()
                        close_button.clicked.connect(lambda checked, idx=i: self.safely_close_tab(idx))
        except Exception as e:
            self.logger.error(f"Error updating tabs: {e}")

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
            self.logger.error(f"Error moving tab: {e}") 