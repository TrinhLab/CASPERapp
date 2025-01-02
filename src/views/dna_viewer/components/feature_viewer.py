from PyQt6.QtWidgets import QGraphicsObject
from PyQt6.QtCore import Qt, QRectF, pyqtSignal, QPointF
from PyQt6.QtGui import QBrush, QColor, QFont

class FeatureViewer(QGraphicsObject):
    """Component for displaying DNA features like genes and exons"""
    cursor_position_changed = pyqtSignal(int)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence = ""
        self.features = []
        self.start_pos = 0
        self.base_width = 15
        self.bases_per_line = 70
        self.feature_height = 20
        self.line_height = 25
        self.feature_spacing = 2  # Spacing between features and strands
        self.setAcceptHoverEvents(True)
        
        # Add logger
        import logging
        self.logger = logging.getLogger(__name__)

    def set_data(self, sequence, features, start_pos):
        """Set sequence and features data"""
        self.sequence = sequence
        self.features = sorted(features, key=lambda x: x['start'])
        self.start_pos = start_pos
        self.update()

    def paint(self, painter, option, widget):
        """Paint the features"""
        if not self.features or not self.sequence:
            return

        # Process each line of sequence
        current_pos = 0
        while current_pos < len(self.sequence):
            line_text = self.sequence[current_pos:current_pos + self.bases_per_line]
            line_num = current_pos // self.bases_per_line
            
            # Calculate y position to be directly below negative strand
            y_pos = line_num * self.line_height * 2
            feature_y = y_pos + self.line_height * 2  # Position below negative strand
            
            # Calculate sequence width for this line
            sequence_width = len(line_text) * self.base_width
            
            # Draw features
            for feature in self.features:
                try:
                    # Calculate relative positions within current line
                    feature_start = feature['start'] - current_pos
                    feature_end = feature['end'] - current_pos
                    
                    # Skip if feature is not in current line
                    if feature_end < 0 or feature_start >= self.bases_per_line:
                        continue
                    
                    # Clip to line boundaries
                    feature_start = max(0, feature_start)
                    feature_end = min(self.bases_per_line, feature_end)
                    
                    # Calculate pixel positions
                    x_start = feature_start * self.base_width
                    x_end = feature_end * self.base_width
                    
                    # Create rectangle for feature
                    feature_rect = QRectF(
                        x_start,
                        feature_y,
                        x_end - x_start,
                        self.feature_height
                    )
                    
                    # Set color based on feature type
                    if feature.get('type') == 'exon':
                        color = QColor(100, 180, 255)  # Light blue for exons
                    else:
                        color = QColor(255, 140, 0)  # Orange for genes
                    
                    # Draw feature rectangle
                    painter.setBrush(QBrush(color))
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.drawRect(feature_rect)
                    
                    # Draw label if enough space
                    label = feature.get('name', '')
                    text_width = painter.fontMetrics().horizontalAdvance(label)
                    if (x_end - x_start) > text_width:
                        text_x = x_start + ((x_end - x_start) - text_width) / 2
                        text_y = feature_y + self.feature_height/2 + 4
                        painter.setPen(Qt.GlobalColor.white)
                        painter.setFont(QFont("Arial", 8))
                        painter.drawText(QPointF(text_x, text_y), label)
                        
                except Exception as e:
                    self.logger.error(f"Error drawing feature: {str(e)}")
                    continue

            current_pos += self.bases_per_line

    def boundingRect(self):
        """Return the bounding rectangle of the component"""
        if not self.sequence:
            return QRectF()
            
        # Calculate exact width based on sequence length
        last_line_length = len(self.sequence) % self.bases_per_line
        if last_line_length == 0:
            last_line_length = self.bases_per_line
        width = max(self.base_width * self.bases_per_line, 
                   self.base_width * last_line_length) + 100
        
        # Calculate height for actual sequence lines
        total_lines = (len(self.sequence) + self.bases_per_line - 1) // self.bases_per_line
        height = total_lines * self.line_height * 2
        
        return QRectF(0, 0, width, height)

    def mousePressEvent(self, event):
        """Handle mouse press to show insertion point"""
        if event.button() == Qt.MouseButton.LeftButton:
            # Calculate position based on click location
            local_pos = event.pos()
            line_number = int(local_pos.y() // (self.line_height * 2))
            base_position = int(local_pos.x() // self.base_width)
            
            # Calculate absolute position
            position = self.start_pos + line_number * self.bases_per_line + base_position
            
            # Emit cursor position
            self.cursor_position_changed.emit(position)
            
        event.accept()

    def mouseMoveEvent(self, event):
        """Handle mouse move to update insertion point"""
        if event.buttons() & Qt.MouseButton.LeftButton:
            # Calculate position based on mouse location
            local_pos = event.pos()
            line_number = int(local_pos.y() // (self.line_height * 2))
            base_position = int(local_pos.x() // self.base_width)
            
            # Calculate absolute position
            position = self.start_pos + line_number * self.bases_per_line + base_position
            
            # Emit cursor position
            self.cursor_position_changed.emit(position)
            
        event.accept()
