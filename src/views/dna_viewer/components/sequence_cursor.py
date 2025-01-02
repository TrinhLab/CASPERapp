from PyQt6.QtWidgets import QGraphicsLineItem
from PyQt6.QtCore import Qt, QRectF, QPointF
from PyQt6.QtGui import QPen, QColor, QPainter

class SequenceCursor(QGraphicsLineItem):
    """A vertical I-beam cursor for DNA sequence that spans between strands"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Use a solid blue line 
        self.pen = QPen(QColor(0, 100, 255)) 
        self.pen.setWidth(2)
        self.pen.setStyle(Qt.PenStyle.SolidLine)
        self.setPen(self.pen)
        
        # Keep cursor on top
        self.setZValue(9999)
        
        # Don't make cursor selectable/focusable
        self.setAcceptHoverEvents(False)
        self.setAcceptedMouseButtons(Qt.MouseButton.NoButton)
        
        # Disable caching for immediate updates
        self.setCacheMode(QGraphicsLineItem.CacheMode.NoCache)

    def set_position(self, x, y, height):
        """Position the cursor at given coordinates"""
        self.setLine(x, y, x, y + height)
        self.update()

    def paint(self, painter, option, widget):
        """Draw the I-beam cursor"""
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # Draw vertical line
        painter.setPen(self.pen)
        line = self.line()
        painter.drawLine(line)
        
        # Draw small horizontal lines at top and bottom using QPointF
        bar_width = 6
        
        # Top bar
        top_start = QPointF(line.x1() - bar_width/2, line.y1())
        top_end = QPointF(line.x1() + bar_width/2, line.y1())
        painter.drawLine(top_start, top_end)
        
        # Bottom bar
        bottom_start = QPointF(line.x1() - bar_width/2, line.y2())
        bottom_end = QPointF(line.x1() + bar_width/2, line.y2())
        painter.drawLine(bottom_start, bottom_end)

    def boundingRect(self):
        """Bounding rectangle for the cursor"""
        line = self.line()
        padding = 4 
        return QRectF(line.x1() - padding, 
                     line.y1() - padding,
                     padding * 2, 
                     line.y2() - line.y1() + padding * 2)