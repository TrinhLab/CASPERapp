import logging
from PyQt6.QtWidgets import QGraphicsScene, QGraphicsLineItem, QGraphicsSimpleTextItem
from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QPen, QColor, QFont, QBrush

class Ruler(QGraphicsScene):
    """Component for displaying the ruler with position markers"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.logger = logging.getLogger(__name__)
        self.base_width = 15
        self.strand_margin = 40  # Width for 5' and 3' indicators
        self.ruler_height = 30
        self.ruler_color = QColor(0, 120, 215)  # Blue 
        
        # Set background color to match interface
        background_color = QColor(240, 240, 240)
        self.setBackgroundBrush(QBrush(background_color))
        
        # Ensure the background fills the entire scene
        self.setSceneRect(QRectF(0, 0, 1000, self.ruler_height - 15))
        
        # Tick mark settings
        self.major_tick_height = 10  # Height for major ticks (multiples of 10)
        self.medium_tick_height = 7  # Height for medium ticks (multiples of 5)
        self.minor_tick_height = 4  # Height for minor ticks
        
        self.ruler_y = 15  # Y position of main ruler line

        # Initialize tracker line as None - we'll create it when needed
        self.tracker_line = None

    def create_ruler(self, bases_per_line):
        """Create ruler with position markers
        
        Args:
            bases_per_line (int): Number of bases per line in sequence viewer
        """
        try:
            self.clear()  # Clear everything

            # Create horizontal blue line aligned with sequence
            ruler_line = QGraphicsLineItem(
                self.strand_margin, self.ruler_y,
                bases_per_line * self.base_width + self.strand_margin, self.ruler_y
            )
            ruler_line.setPen(QPen(self.ruler_color, 1))
            self.addItem(ruler_line)
            
            # Add tick marks and numbers
            for i in range(0, bases_per_line):
                x_pos = i * self.base_width + self.strand_margin + self.base_width/2
                pos_1_based = i + 1
                
                # Determine tick properties based on position
                if pos_1_based % 10 == 0:  # Major ticks (every 10)
                    tick_height = self.major_tick_height
                    tick_start = 10
                    # Add number
                    text = QGraphicsSimpleTextItem(str(pos_1_based))
                    text.setFont(QFont("Arial", 8))
                    text_width = text.boundingRect().width()
                    text.setPos(x_pos - text_width/2, 0)
                    self.addItem(text)
                elif pos_1_based % 5 == 0:  # Medium ticks (every 5)
                    tick_height = self.medium_tick_height
                    tick_start = 11
                else:  # Minor ticks
                    tick_height = self.minor_tick_height
                    tick_start = 13
                
                # Create tick mark
                tick = QGraphicsLineItem(
                    x_pos, tick_start,
                    x_pos, tick_start + tick_height
                )
                tick.setPen(QPen(self.ruler_color, 1))
                self.addItem(tick)
            
            # Create new tracker line
            self.tracker_line = QGraphicsLineItem()
            pen = QPen(self.ruler_color)
            pen.setWidth(2)
            pen.setStyle(Qt.PenStyle.SolidLine)
            self.tracker_line.setPen(pen)
            self.tracker_line.setZValue(100)
            self.tracker_line.hide()
            self.addItem(self.tracker_line)
            
        except Exception as e:
            self.logger.error(f"Error creating ruler: {str(e)}")

    def boundingRect(self):
        """Return the bounding rectangle of the ruler"""
        return self.sceneRect()

    def update_tracker_position(self, x_pos):
        """Update the position of the tracker line"""
        try:
            if self.tracker_line is None:
                # Create tracker line if it doesn't exist
                self.tracker_line = QGraphicsLineItem()
                pen = QPen(self.ruler_color)
                pen.setWidth(2)
                pen.setStyle(Qt.PenStyle.SolidLine)
                self.tracker_line.setPen(pen)
                self.tracker_line.setZValue(100)
                self.addItem(self.tracker_line)

            # Ensure x_pos is within bounds
            scene_rect = self.sceneRect()
            x_pos = max(self.strand_margin, min(x_pos, scene_rect.width() - self.strand_margin))
            
            # Set tracker line position
            self.tracker_line.setLine(x_pos, 0, x_pos, self.ruler_height)
            self.tracker_line.show()
            self.update()
            
        except Exception as e:
            self.logger.error(f"Error updating tracker position: {str(e)}")
