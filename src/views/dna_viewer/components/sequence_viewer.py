import traceback
from PyQt6.QtWidgets import QGraphicsObject, QApplication, QGraphicsSimpleTextItem, QGraphicsLineItem
from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QPen, QFont, QColor
from .nucleotide_item import NucleotideItem
import logging

class SequenceViewer(QGraphicsObject):
    """Component for displaying and interacting with DNA sequence"""
    sequence_selected = pyqtSignal(int, int)  # Emit start and end positions
    cursor_position_changed = pyqtSignal(int)  # Emit cursor position
    
    def __init__(self, parent=None, logger=None):
        super().__init__(parent)
        
        self.logger = logging.getLogger(__name__)
        
        self._init_properties()
        self._init_graphics_storage()

    def _init_properties(self):
        """Initialize properties"""
        # Layout properties
        self.strand_margin = 40
        self.base_width = 15
        self.bases_per_line = 70
        self.line_height = 25
        self.line_spacing = 80
        
        # Sequence properties
        self.sequence = ""
        self.start_pos = 0
        
        # Selection properties
        self.selection_start = None
        self.selection_end = None
        self.drag_start_pos = None
        self.selection_active = False
        
        # Cursor properties
        self.cursor_pos = None
        
        # Clipboard
        self.clipboard = QApplication.clipboard()
        
        # Add tracking for cleared selection highlights
        self.cleared_selection_positions = set()  # Track positions where selection was cleared

    def _init_graphics_storage(self):
        """Initialize storage for graphics items"""
        self.nucleotides = []
        self.highlighted_regions = []
        self.nucleotide_map = {'+': [], '-': []}
        self.plot_lines = []
        self.tick_lines = []

    def set_data(self, sequence, start_pos=None):
        """Set sequence data and create nucleotide items"""
        try:
            if not sequence:
                self.logger.warning("Empty sequence provided")
                return
            
            self.sequence = sequence
            self.start_pos = start_pos if start_pos is not None else 0
            
            # Store current highlights, excluding cleared selection highlights
            current_highlights = []
            selection_blue = QColor(100, 150, 255, 100)
            
            for nuc in self.nucleotides:
                if nuc.is_highlighted:
                    idx = self.get_nucleotide_position(nuc)
                    pos = idx // 2  # Convert to sequence position
                    
                    # Only store if it's not a selection highlight that was cleared
                    if nuc.highlight_color != selection_blue or pos not in self.cleared_selection_positions:
                        current_highlights.append({
                            'position': pos,
                            'color': nuc.highlight_color
                        })
            
            # Batch update
            if self.scene():
                views = self.scene().views()
                for view in views:
                    view.setUpdatesEnabled(False)
            
            try:
                # Create display
                self._create_display()
                
                # Reapply highlights, excluding cleared selections
                for highlight in current_highlights:
                    pos = highlight['position']
                    color = highlight['color']
                    
                    # Skip if this was a cleared selection
                    if color == selection_blue and pos in self.cleared_selection_positions:
                        continue
                        
                    pos_strand_idx = pos * 2
                    neg_strand_idx = pos * 2 + 1
                    
                    for idx in [pos_strand_idx, neg_strand_idx]:
                        if idx < len(self.nucleotides):
                            nuc = self.nucleotides[idx]
                            nuc.is_highlighted = True
                            nuc.highlight_color = color
                            nuc.update()
                
            finally:
                # Re-enable updates
                if self.scene():
                    for view in views:
                        view.setUpdatesEnabled(True)
                        view.viewport().update()

        except Exception as e:
            self.logger.error(f"Error in set_data: {str(e)}")

    def _create_display(self):
        try:
            # Clear existing items
            self.cleanup_graphics()
            
            total_lines = (len(self.sequence) + self.bases_per_line - 1) // self.bases_per_line
            
            nucleotides_batch = []
            plot_lines_batch = []
            tick_lines_batch = []
            position_numbers = [] 
            
            current_pos = 0
            while current_pos < len(self.sequence):
                remaining_bases = len(self.sequence) - current_pos
                bases_this_line = min(self.bases_per_line, remaining_bases)
                line_text = self.sequence[current_pos:current_pos + bases_this_line]
                
                line_num = current_pos // self.bases_per_line
                y_pos = line_num * self.line_spacing
                
                # Create nucleotides for this line
                for i, nuc in enumerate(line_text):
                    x_pos = i * self.base_width + self.strand_margin
                    
                    # Calculate absolute position in sequence
                    abs_pos = current_pos + i
                    
                    # Determine if this nucleotide is part of padding
                    # Only consider it padding if it's at the start or end of the sequence
                    is_padding = (nuc.islower() and nuc in 'atgc' and 
                                (abs_pos < 30 or abs_pos >= len(self.sequence) - 30)) 
                    
                    # Create positive strand nucleotide
                    nuc_item = NucleotideItem(
                        nucleotide=nuc,
                        x=x_pos,
                        y=y_pos + self.line_height * 0.1,
                        width=self.base_width,
                        is_uppercase=nuc.isupper(),
                        is_padding=is_padding,
                        parent=self
                    )
                    nucleotides_batch.append(nuc_item)
                    
                    # Create complement strand nucleotide
                    complement_item = NucleotideItem(
                        nucleotide=nuc,
                        x=x_pos,
                        y=y_pos + self.line_height * 1.45,
                        width=self.base_width,
                        is_uppercase=nuc.isupper(),
                        is_complement=True,
                        is_padding=is_padding,
                        parent=self
                    )
                    nucleotides_batch.append(complement_item)
                
                # Create plot line
                plot_y = y_pos + self.line_height * 1.15
                plot_line = QGraphicsLineItem(
                    self.strand_margin,
                    plot_y,
                    self.strand_margin + bases_this_line * self.base_width,
                    plot_y,
                    self
                )
                plot_line.setPen(QPen(Qt.GlobalColor.black, 1))
                plot_lines_batch.append(plot_line)
                
                # Add position number at end of line
                end_pos = str(self.start_pos + current_pos + bases_this_line)
                pos_item = QGraphicsSimpleTextItem(end_pos, self)
                pos_item.setFont(QFont("Courier", 12))
                
                # Position text with fixed offset
                extra_spacing = 25 if line_num == total_lines - 1 else 0
                pos_x = self.strand_margin + (bases_this_line * self.base_width) + extra_spacing + 10
                pos_y = plot_y - pos_item.boundingRect().height()/2
                pos_item.setPos(pos_x, pos_y)
                position_numbers.append(pos_item)
                
                # Create tick marks
                for i in range(bases_this_line):
                    x_pos = i * self.base_width + self.strand_margin
                    pos_1_based = current_pos + i + 1
                    
                    tick_height = 12 if (i == 0 and current_pos == 0) or \
                                 (i == bases_this_line - 1 and current_pos + bases_this_line == len(self.sequence)) else \
                                 10 if pos_1_based % 10 == 0 else \
                                 8 if pos_1_based % 5 == 0 else 5
                    
                    tick = QGraphicsLineItem(
                        x_pos + self.base_width/2,
                        plot_y - tick_height/2,
                        x_pos + self.base_width/2,
                        plot_y + tick_height/2,
                        self
                    )
                    tick.setPen(QPen(Qt.GlobalColor.black, 1))
                    tick_lines_batch.append(tick)
                
                current_pos += bases_this_line
            
            # Add all items to scene in batches
            self.nucleotides = nucleotides_batch
            self.plot_lines = plot_lines_batch
            self.tick_lines = tick_lines_batch
            
            # Add strand indicators and position numbers
            self._add_strand_indicators(total_lines)
            
        except Exception as e:
            self.logger.error(f"Error in create optimized display: {str(e)}")

    def _add_strand_indicators(self, total_lines):
        # Add first line indicators
        five_prime_pos = QGraphicsSimpleTextItem("5'", self)
        five_prime_pos.setFont(QFont("Arial", 10))
        five_prime_pos.setPos(0, self.line_height * 0.26)
        
        three_prime_neg = QGraphicsSimpleTextItem("3'", self)
        three_prime_neg.setFont(QFont("Arial", 10))
        three_prime_neg.setPos(0, self.line_height * 1.58)
        
        # Add last line indicators
        last_line_width = (len(self.sequence) % self.bases_per_line) * self.base_width
        if last_line_width == 0:
            last_line_width = self.bases_per_line * self.base_width
        
        three_prime_pos = QGraphicsSimpleTextItem("3'", self)
        three_prime_pos.setFont(QFont("Arial", 10))
        three_prime_pos.setPos(
            last_line_width + self.strand_margin + 20,
            (total_lines - 1) * self.line_spacing + self.line_height * 0.26
        )
        
        five_prime_neg = QGraphicsSimpleTextItem("5'", self)
        five_prime_neg.setFont(QFont("Arial", 10))
        five_prime_neg.setPos(
            last_line_width + self.strand_margin + 20,
            (total_lines - 1) * self.line_spacing + self.line_height * 1.58
        )

    def cleanup_graphics(self):
        """Clean up all graphics items"""
        try:
            # Remove plot lines
            for line in self.plot_lines:
                if line.scene():
                    line.scene().removeItem(line)
            self.plot_lines.clear()
            
            # Remove tick lines
            for line in self.tick_lines:
                if line.scene():
                    line.scene().removeItem(line)
            self.tick_lines.clear()

            # Remove nucleotides
            for nuc in self.nucleotides:
                if nuc.scene():
                    nuc.scene().removeItem(nuc)
            self.nucleotides.clear()

            # Remove all text items (including position numbers)
            for item in self.childItems():
                if isinstance(item, QGraphicsSimpleTextItem):
                    if item.scene():
                        item.scene().removeItem(item)

            # Clear nucleotide maps but preserve highlights
            self.nucleotide_map['+'].clear()
            self.nucleotide_map['-'].clear()
            
            self.logger.debug("Cleaned up all graphics items")
            
        except Exception as e:
            self.logger.error(f"Error in cleanup_graphics: {str(e)}")

    def highlight_sequence(self, start_pos, end_pos, color, strand='+'):
        """Highlight sequence with proper strand handling"""
        try:
            # Store highlight information
            self.highlighted_regions.append((start_pos, end_pos, color, strand))
            
            # Calculate which lines contain the sequence
            start_line = start_pos // self.bases_per_line
            end_line = end_pos // self.bases_per_line
            
            # Calculate positions within lines
            start_pos_in_line = start_pos % self.bases_per_line
            end_pos_in_line = end_pos % self.bases_per_line
            
            self.logger.debug(f"Highlighting from line {start_line} to {end_line}")
            self.logger.debug(f"Start pos in line: {start_pos_in_line}, End pos in line: {end_pos_in_line}")
            
            # For each line that contains part of the sequence
            for line_num in range(start_line, end_line + 1):
                # Calculate start and end positions for this line
                line_start = start_pos_in_line if line_num == start_line else 0
                line_end = end_pos_in_line if line_num == end_line else self.bases_per_line - 1
                
                # Calculate base indices for this line
                base_start = line_num * self.bases_per_line + line_start
                base_end = line_num * self.bases_per_line + line_end
                
                # Highlight nucleotides
                for i in range(base_start, base_end + 1):
                    if i >= len(self.nucleotides):
                        break
                    
                    # For negative strand, highlight both strands' nucleotides
                    if strand == '-':
                        # Highlight negative strand nucleotide
                        neg_idx = i * 2 + 1  # Odd indices for negative strand
                        if neg_idx < len(self.nucleotides):
                            nuc = self.nucleotides[neg_idx]
                            nuc.is_highlighted = True
                            nuc.highlight_color = color
                            nuc.update()
                    else:
                        # Highlight positive strand nucleotide
                        pos_idx = i * 2  # Even indices for positive strand
                        if pos_idx < len(self.nucleotides):
                            nuc = self.nucleotides[pos_idx]
                            nuc.is_highlighted = True
                            nuc.highlight_color = color
                            nuc.update()
            
            # Force scene update
            if self.scene():
                self.scene().update()
            
            self.logger.debug(f"Highlighted sequence on strand {strand} from {start_pos} to {end_pos}")
            
        except Exception as e:
            self.logger.error(f"Error in highlight_sequence: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def clear_highlights(self):
        """Clear all highlights with optimized rendering"""
        try:
            # Get view and disable updates
            view = None
            if self.scene():
                views = self.scene().views()
                if views:
                    view = views[0]
                    view.setUpdatesEnabled(False)
            
            try:
                self.highlighted_regions.clear()
                
                # Batch collect nucleotides that need updating
                nucleotides_to_update = []
                
                for nuc in self.nucleotides:
                    if nuc.is_highlighted:
                        nuc.is_highlighted = False
                        nuc.highlight_color = None
                        nucleotides_to_update.append(nuc)
                
                # Update all nucleotides in a single batch
                if nucleotides_to_update:
                    # Use prepareGeometryChange for better performance
                    for nuc in nucleotides_to_update:
                        nuc.prepareGeometryChange()
                    
                    # Force a single scene update
                    if self.scene():
                        self.scene().update()
                    
            finally:
                # Re-enable view updates
                if view:
                    view.setUpdatesEnabled(True)
                    view.viewport().update()
            
        except Exception as e:
            self.logger.error(f"Error in clear_highlights: {str(e)}")
            # Make sure view updates are re-enabled
            if self.scene():
                views = self.scene().views()
                if views:
                    views[0].setUpdatesEnabled(True)
                    views[0].viewport().update()

    def get_nucleotide_position(self, nucleotide):
        """Get the position of a nucleotide in the sequence"""
        try:
            idx = self.nucleotides.index(nucleotide)
            return idx
        except ValueError:
            return -1

    def boundingRect(self):
        """Return the bounding rectangle"""
        if not self.sequence:
            return QRectF(0, 0, 100, 100)  # Return a default size when empty
        
        # Calculate total width including margins
        width = (self.base_width * self.bases_per_line) + (self.strand_margin * 2)
        
        # Calculate height using line spacing
        total_lines = (len(self.sequence) + self.bases_per_line - 1) // self.bases_per_line
        height = total_lines * self.line_spacing
        
        # Add some padding
        height += 50
        width += 50
        
        return QRectF(0, 0, width, height)

    def paint(self, painter, option, widget):
        """Paint method required by QGraphicsObject"""
        pass