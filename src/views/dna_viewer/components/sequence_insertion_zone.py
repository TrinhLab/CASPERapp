from PyQt6.QtWidgets import QGraphicsObject
from PyQt6.QtCore import Qt, QRectF, pyqtSignal, QPointF, QEvent
from PyQt6.QtGui import QColor, QPainterPath
from .sequence_cursor import SequenceCursor
import logging

class SequenceInsertionZone(QGraphicsObject):
    """Handles the interactive zone between base pairs for insertion/deletion"""
    
    insertion_point_selected = pyqtSignal(int)  # Emits position when zone is clicked
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.logger = logging.getLogger(__name__)
        
        self.zone_width = 8 
        
        # State tracking
        self.hover_position = None
        self.active_zones = []  # List of (x, y, width, height) for each zone
        
        self.sequence_length = 0
        self.bases_per_line = 70 # Default value
        self.base_width = 15     
        self.strand_margin = 40 
        self.line_height = 25    
        
        self.sequence_cursor = SequenceCursor(self)
        self.sequence_cursor.hide()
        self.current_cursor_pos = None
        
        self.setAcceptHoverEvents(True)
        self.setAcceptedMouseButtons(Qt.MouseButton.LeftButton)
        self.setFlag(QGraphicsObject.GraphicsItemFlag.ItemIsSelectable, True)
        self.setFlag(QGraphicsObject.GraphicsItemFlag.ItemIsFocusable, True)
        self.setFlag(QGraphicsObject.GraphicsItemFlag.ItemClipsToShape, False)
        self.setZValue(100)  # Keep zones above sequence
        
        # Add tracking for visible area
        self.visible_rect = QRectF()
        self.scale_factor = 1.0
        
        # Add tracking for hover state
        self.last_valid_hover_pos = None
        self.hover_active = False
        
    def create_zones(self, sequence_length, base_width, strand_margin, line_height, bases_per_line, line_spacing=None):
        """Create insertion zones between bases"""
        self.active_zones.clear()
        self.line_height = line_height
        self.line_spacing = line_spacing if line_spacing is not None else line_height * 3
        self.bases_per_line = bases_per_line
        self.base_width = base_width
        self.strand_margin = strand_margin
        self.sequence_length = sequence_length
        self.update()

    def contains(self, point):
        """Override contains to better handle hover detection"""
        # Convert point to local coordinates if needed
        local_point = point
        if isinstance(point, QPointF):
            local_point = self.mapFromScene(point)
        
        # Calculate position
        line_number = int(local_point.y() / self.line_spacing)
        x_relative = local_point.x() - self.strand_margin
        position_in_line = int(x_relative / self.base_width)
        absolute_position = (line_number * self.bases_per_line) + position_in_line
        
        # Calculate total lines
        total_lines = (self.sequence_length + self.bases_per_line - 1) // self.bases_per_line
        
        # Check if position is valid with more precise bounds
        in_x_range = -self.zone_width <= x_relative <= (self.bases_per_line * self.base_width + self.zone_width)
        in_y_range = 0 <= line_number < total_lines
        position_valid = 0 <= absolute_position <= self.sequence_length
        
        # Calculate vertical position within line with more generous zones
        y_in_line = local_point.y() % self.line_spacing
        
        # Define zones with appropriate overlap and coverage
        zone_height = self.line_height * 1.8 
        middle_point = self.line_spacing / 2
        
        # Upper zone covers from start to middle + overlap
        upper_zone_start = 0
        upper_zone_end = middle_point + (zone_height / 2) 
        
        # Lower zone covers from middle - overlap to end
        lower_zone_start = middle_point - (zone_height / 2)  
        lower_zone_end = self.line_spacing  
        
        # Allow interaction in both upper and lower strand regions with overlap
        upper_strand_zone = upper_zone_start <= y_in_line <= upper_zone_end
        lower_strand_zone = lower_zone_start <= y_in_line <= lower_zone_end
        in_vertical_zone = upper_strand_zone or lower_strand_zone
        
        return in_x_range and in_y_range and position_valid and in_vertical_zone

    def hoverMoveEvent(self, event):
        """Handle hover using direct coordinate calculation"""
        pos = event.pos()
        
        if self.contains(pos):
            # Calculate snapped position for cursor only
            line_number = int(pos.y() / self.line_spacing)
            x_relative = pos.x() - self.strand_margin
            position_in_line = int(x_relative / self.base_width)
            absolute_position = (line_number * self.bases_per_line) + position_in_line
            
            # Update cursor behavior
            if self.hover_position != absolute_position:
                self.hover_position = absolute_position
                self.setCursor(Qt.CursorShape.IBeamCursor)
            
            # Update ruler tracker with exact mouse position - no snapping
            if self.scene() and self.scene().parent():
                dna_viewer = self.scene().parent()
                if hasattr(dna_viewer, 'ruler_scene'):
                    # Get exact mouse position in scene coordinates
                    scene_pos = self.mapToScene(pos)
                    
                    # Adjust for any offset between sequence view and ruler
                    ruler_x = scene_pos.x()
                    
                    # Account for margin differences between sequence and ruler
                    margin_diff = self.strand_margin - dna_viewer.ruler_scene.strand_margin
                    if margin_diff != 0:
                        ruler_x -= margin_diff
                    
                    # Update tracker immediately
                    dna_viewer.ruler_scene.update_tracker_position(ruler_x)
                    
                    # Force immediate updates
                    dna_viewer.ruler_scene.update()
                    dna_viewer.ruler_view.viewport().update()
                    
                    # Update scene for smooth rendering
                    if self.scene():
                        self.scene().update()
        else:
            if self.hover_position is not None:
                self.hover_position = None
                self.setCursor(Qt.CursorShape.ArrowCursor)
                
                # Hide ruler tracker when not hovering
                if self.scene() and self.scene().parent():
                    dna_viewer = self.scene().parent()
                    if hasattr(dna_viewer, 'ruler_scene'):
                        dna_viewer.ruler_scene.tracker_line.hide()
                        dna_viewer.ruler_scene.update()
                        dna_viewer.ruler_view.viewport().update()
        
        return super().hoverMoveEvent(event)


    def mousePressEvent(self, event):
        """Handle mouse press to only place cursor without highlighting"""
        pos = event.pos()
        
        # Calculate position relative to sequence
        x_relative = pos.x() - self.strand_margin
        line_number = int(pos.y() / self.line_spacing)
        
        # Calculate nearest space between bases for cursor
        raw_position = x_relative / self.base_width
        cursor_position = round(raw_position)  # Round to nearest space
        absolute_position = (line_number * self.bases_per_line) + cursor_position
        
        # Confine cursor position within sequence boundaries
        absolute_position = max(0, min(absolute_position, self.sequence_length))
        cursor_position = absolute_position % self.bases_per_line
        line_number = absolute_position // self.bases_per_line
        
        # Store selection start but don't highlight yet
        self.selection_start = absolute_position
        
        # Reset selection state
        self.selection_active = False
        
        # Position cursor at nearest space between bases
        cursor_x = self.strand_margin + (cursor_position * self.base_width)
        cursor_y = (line_number * self.line_spacing) + (self.line_height * 0.1)
        cursor_height = self.line_height * 2 + 5
        
        # Show cursor
        self.sequence_cursor.set_position(cursor_x, cursor_y, cursor_height)
        self.sequence_cursor.show()
        
        # Store and emit current cursor position
        self.current_cursor_pos = absolute_position
        
        # Get DNA feature viewer and emit cursor position change
        if self.scene() and self.scene().parent():
            dna_viewer = self.scene().parent()
            sequence_viewer = dna_viewer.sequence_viewer
            
            # Clear selection state and highlights
            sequence_viewer.selection_start = None
            sequence_viewer.selection_end = None
            sequence_viewer.selection_active = False
            
            # Clear selection highlights
            selection_color = QColor(100, 150, 255, 100)
            for nuc in sequence_viewer.nucleotides:
                if nuc.highlight_color == selection_color:
                    nuc.is_highlighted = False
                    nuc.highlight_color = None
                    nuc.update()
            
            # Emit cursor position change immediately
            sequence_viewer.cursor_position_changed.emit(absolute_position)
        
        event.accept()

    def mouseMoveEvent(self, event):
        """Handle mouse drag for selection"""
        # Only start selection if mouse has moved
        if not self.selection_active:
            # Check if mouse has moved enough to start selection
            initial_pos = event.pos()
            x_relative = initial_pos.x() - self.strand_margin
            raw_position = x_relative / self.base_width
            current_position = round(raw_position)
            
            # Only activate selection if mouse has moved to a different position
            if (current_position != self.selection_start // self.bases_per_line):
                self.selection_active = True
        
        if self.selection_active:
            pos = event.pos()
            
            # Calculate current position
            x_relative = pos.x() - self.strand_margin
            line_number = int(pos.y() / self.line_spacing)
            raw_position = x_relative / self.base_width
            cursor_position = round(raw_position)  # For cursor placement
            
            # For highlighting, use the same rounding logic as cursor
            highlight_position = cursor_position
            
            # Convert to sequence positions
            current_pos = (line_number * self.bases_per_line) + highlight_position
            
            # Confine position within sequence boundaries
            current_pos = max(0, min(current_pos, self.sequence_length))
            cursor_position = current_pos % self.bases_per_line
            line_number = current_pos // self.bases_per_line
            
            if 0 <= current_pos <= self.sequence_length:
                scene = self.scene()
                if scene and scene.parent():
                    dna_viewer = scene.parent()
                    sequence_viewer = dna_viewer.sequence_viewer
                    
                    # Update ruler tracker with exact mouse position during selection
                    if hasattr(dna_viewer, 'ruler_scene'):
                        # Get exact mouse position in scene coordinates
                        scene_pos = self.mapToScene(pos)
                        
                        # Adjust for any offset between sequence view and ruler
                        ruler_x = scene_pos.x()
                        
                        # Account for margin differences between sequence and ruler
                        margin_diff = self.strand_margin - dna_viewer.ruler_scene.strand_margin
                        if margin_diff != 0:
                            ruler_x -= margin_diff
                        
                        # Update tracker immediately
                        dna_viewer.ruler_scene.update_tracker_position(ruler_x)
                        dna_viewer.ruler_scene.update()
                        dna_viewer.ruler_view.viewport().update()
                    
                    # Clear previous selection highlights
                    selection_color = QColor(100, 150, 255, 100)
                    for nuc in sequence_viewer.nucleotides:
                        if nuc.highlight_color == selection_color:
                            nuc.is_highlighted = False
                            nuc.highlight_color = None
                            nuc.update()
                    
                    # Determine selection range
                    if current_pos >= self.selection_start:
                        start = self.selection_start
                        end = current_pos - 1
                    else:
                        start = current_pos
                        end = self.selection_start - 1
                    
                    try:
                        # Only highlight if there's a valid range
                        if start <= end:
                            for base_idx in range(start, end + 1):
                                pos_strand_idx = base_idx * 2
                                neg_strand_idx = base_idx * 2 + 1
                                
                                for idx in [pos_strand_idx, neg_strand_idx]:
                                    if idx < len(sequence_viewer.nucleotides):
                                        nuc = sequence_viewer.nucleotides[idx]
                                        if not nuc.is_highlighted or nuc.highlight_color == selection_color:
                                            nuc.is_highlighted = True
                                            nuc.highlight_color = selection_color
                                            nuc.update()
                        
                            # Emit selection signal with adjusted end position for display
                            absolute_start = dna_viewer._original_start_pos + start + 1  # Add 1 only for display
                            absolute_end = dna_viewer._original_start_pos + end + 1  # Add 1 only for display
                            sequence_viewer.sequence_selected.emit(absolute_start, absolute_end)
                        
                        # Update scene and cursor
                        if scene:
                            scene.update()
                            
                        # Position cursor at nearest space between bases
                        cursor_x = self.strand_margin + (cursor_position * self.base_width)
                        cursor_y = (line_number * self.line_spacing) + (self.line_height * 0.1)
                        cursor_height = self.line_height * 2 + 5
                        
                        self.sequence_cursor.set_position(cursor_x, cursor_y, cursor_height)
                        self.sequence_cursor.show()
                        
                        self.current_cursor_pos = current_pos
                            
                    except Exception as e:
                        print(f"Error highlighting: {str(e)}")
                
            event.accept()

    def mouseReleaseEvent(self, event):
        """Handle mouse release to end selection"""
        if self.selection_active:
            self.selection_active = False
            
            # Emit the final selection position
            if hasattr(self, 'selection_start') and self.selection_start is not None:
                self.insertion_point_selected.emit(self.selection_start)
            
            event.accept()

    def boundingRect(self):
        """Return bounding rectangle for entire sequence area"""
        if not hasattr(self, 'sequence_length') or self.sequence_length == 0:
            return QRectF(0, 0, 1, 1)  # Return minimal rect if no sequence
        
        total_lines = (self.sequence_length + self.bases_per_line - 1) // self.bases_per_line
        width = self.strand_margin * 2 + (self.bases_per_line * self.base_width)
        height = total_lines * self.line_spacing
        
        # Add padding
        padding = 20
        return QRectF(-padding, -padding, width + 2*padding, height + 2*padding)

    def paint(self, painter, option, widget):
        """Paint method required by QGraphicsObject"""
        pass