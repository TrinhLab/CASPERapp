from PyQt6.QtWidgets import QGraphicsObject
from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QFont, QColor

class NucleotideItem(QGraphicsObject):
    """Component for displaying individual nucleotides in the gene sequence viewer"""
    clicked = pyqtSignal(object)
    
    def __init__(self, nucleotide, x, y, width, is_uppercase=True, is_complement=False, is_padding=False, parent=None):
        """Initialize nucleotide item
        
        Args:
            nucleotide (str): The nucleotide base (A, T, G, C)
            x (float): X position
            y (float): Y position 
            width (float): Width of nucleotide
            is_uppercase (bool): Whether nucleotide should be uppercase
            is_complement (bool): Whether this is a complement strand nucleotide
            is_padding (bool): Whether this nucleotide is part of padding sequence
            parent (QGraphicsObject): Parent item
        """
        super().__init__(parent)
        self.nucleotide = nucleotide
        self.setPos(x, y)
        self.width = width
        self.base_width = width
        self.height = 20
        self.is_uppercase = is_uppercase
        self.is_complement = is_complement
        self.is_padding = is_padding
        self.is_highlighted = False
        self.highlight_color = None
        
        # Enable hover events
        self.setAcceptHoverEvents(True)

    def paint(self, painter, option, widget):
        """Paint the nucleotide"""
        try:
            # Draw highlight background if highlighted
            if self.is_highlighted and self.highlight_color:
                painter.save()
                painter.fillRect(self.boundingRect(), self.highlight_color)
                painter.restore()
            
            # Draw nucleotide centered
            painter.setFont(QFont("Courier", 12))
            
            # Use grey for padding sequence, black for all other nucleotides
            if self.is_padding:
                painter.setPen(QColor(100, 100, 100))
            else:
                painter.setPen(Qt.GlobalColor.black)
            
            # Get complement nucleotide if needed
            display_nucleotide = self._get_complement() if self.is_complement else self.nucleotide
            
            # Draw text centered
            painter.drawText(self.boundingRect(), Qt.AlignmentFlag.AlignCenter, display_nucleotide)
                
        except Exception as e:
            print(f"Error in paint: {str(e)}")

    def boundingRect(self):
        """Return the bounding rectangle for the nucleotide"""
        return QRectF(0, 0, self.width, self.height)

    def mousePressEvent(self, event):
        """Handle mouse press events"""
        if event.button() == Qt.MouseButton.LeftButton:
            sequence_viewer = self.parent()
            if sequence_viewer:
                pos = sequence_viewer.get_nucleotide_position(self)
                local_x = event.pos().x()
                
                # Calculate position relative to letter boundaries
                text_x = (self.width - self.base_width) / 2
                relative_x = local_x - text_x
                
                # Determine cursor position
                if relative_x <= 0:  # Before letter
                    cursor_pos = pos
                elif relative_x >= self.base_width:  # After letter
                    cursor_pos = pos + 1
                else:  # On letter
                    cursor_pos = pos + (1 if relative_x > self.base_width / 2 else 0)
                
                # Update selection and cursor
                sequence_viewer.selection_start = cursor_pos
                sequence_viewer.selection_end = cursor_pos
                sequence_viewer.selection_active = True
                
                sequence_viewer.cursor_position_changed.emit(cursor_pos)
                sequence_viewer._update_selection()
                self.update()
        
        event.accept()

    def _get_complement(self):
        """Get complement nucleotide"""
        complement_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
            'K': 'M', 'M': 'K', 'R': 'Y', 'Y': 'R',
            'k': 'm', 'm': 'k', 'r': 'y', 'y': 'r',
            'S': 'S', 's': 's'
        }
        return complement_map.get(self.nucleotide, self.nucleotide)
