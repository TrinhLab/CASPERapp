from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QGraphicsView, QGraphicsScene, 
                            QGraphicsObject, QGraphicsSimpleTextItem, QApplication, 
                            QLabel, QFrame, QGraphicsLineItem)
from PyQt6.QtCore import Qt, QRectF, pyqtSignal, QPointF, QLineF, QSizeF, QTimer
from PyQt6.QtGui import QPainter, QPen, QBrush, QColor, QPainterPath, QFont, QPolygonF, QTransform, QKeySequence 

class DNAFeatureViewer(QWidget):  # Change to QWidget
    """Custom widget for displaying DNA features with sequence"""
    sequence_selected = pyqtSignal(int, int)  # Emit start and end positions when sequence is selected
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Add logger
        import logging
        self.logger = logging.getLogger(__name__)
        
        # Create layout for this widget
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)
        
        # Create graphics view with left alignment
        self.view = QGraphicsView()
        self.scene = QGraphicsScene(self)
        self.view.setScene(self.scene)
        self.view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        
        # Create components
        self.sequence_viewer = SequenceViewer()
        self.feature_viewer = FeatureViewer()
        
        # Add components to scene
        self.scene.addItem(self.sequence_viewer)
        self.scene.addItem(self.feature_viewer)
        
        # Connect signals from both viewers
        self.sequence_viewer.sequence_selected.connect(self._on_sequence_selected)
        self.sequence_viewer.cursor_position_changed.connect(self._on_cursor_position_changed)
        self.feature_viewer.cursor_position_changed.connect(self._on_cursor_position_changed)
        
        # Create status panel
        self.status_panel = QLabel()
        self.status_panel.setFrameStyle(QFrame.Shape.Panel | QFrame.Shadow.Sunken)
        self.status_panel.setStyleSheet("""
            QLabel {
                background-color: #f0f0f0;
                padding: 5px;
                border-top: 1px solid #ccc;
                min-height: 20px;
            }
        """)
        
        # Add widgets to layout
        self.layout.addWidget(self.view)
        self.layout.addWidget(self.status_panel)
        
        # Setup view
        self.view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        self.view.setMinimumHeight(200)

        # Add resize event handling
        self.view.viewport().installEventFilter(self)

        # Add resize timer for debouncing
        self.resize_timer = QTimer()
        self.resize_timer.setSingleShot(True)
        self.resize_timer.timeout.connect(self._delayed_resize)
        self.cached_size = None

        # Create ruler view with left alignment
        self.ruler_view = QGraphicsView()
        self.ruler_scene = QGraphicsScene(self)
        self.ruler_view.setScene(self.ruler_scene)
        self.ruler_view.setFixedHeight(25)
        self.ruler_view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.ruler_view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.ruler_view.setViewportMargins(0, 0, 0, 0)
        self.ruler_view.setFrameStyle(0)
        self.ruler_view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        
        # Set opaque white background to hide sequence behind ruler
        self.ruler_view.setBackgroundBrush(QBrush(Qt.GlobalColor.white))
        self.ruler_view.setAutoFillBackground(True)
        
        # Add ruler to layout before main view
        self.layout.insertWidget(0, self.ruler_view)
        
        # Connect horizontal scroll bars
        self.view.horizontalScrollBar().valueChanged.connect(
            self.ruler_view.horizontalScrollBar().setValue
        )
        
        # Create initial ruler
        self._create_ruler()

    def eventFilter(self, obj, event):
        """Handle viewport resize events with debouncing"""
        if obj == self.view.viewport() and event.type() == event.Type.Resize:
            # Cache the new size
            self.cached_size = event.size()
            # Reset and start the timer
            self.resize_timer.stop()
            self.resize_timer.start(150)  # 150ms delay
        return super().eventFilter(obj, event)

    def _delayed_resize(self):
        """Handle resize after debouncing period"""
        if self.cached_size:
            try:
                self.view.setUpdatesEnabled(False)
                self.sequence_viewer.setVisible(False)
                self.feature_viewer.setVisible(False)
                
                self._handle_resize(self.cached_size)
                
                # Call reapply_highlights on sequence_viewer instead of self
                self.sequence_viewer._reapply_highlights()
                
                self.sequence_viewer.setVisible(True)
                self.feature_viewer.setVisible(True)
                self.view.setUpdatesEnabled(True)
                
                # Force immediate update
                self.view.viewport().update()
                
            except Exception as e:
                self.logger.error(f"Error in delayed resize: {str(e)}")
            finally:
                self.view.setUpdatesEnabled(True)
                self.sequence_viewer.setVisible(True)
                self.feature_viewer.setVisible(True)

    def _handle_resize(self, new_size):
        """Handle viewport resize by adjusting sequence display"""
        try:
            viewport_width = max(1, new_size.width())
            margin = 100  # Space for position numbers
            
            # Calculate available width for bases
            available_width = viewport_width - margin
            base_width = self.sequence_viewer.base_width
            
            # Calculate maximum number of bases that can fit
            max_bases = (available_width // base_width)
            
            # Round down to nearest multiple of 10
            new_bases = (max_bases // 10) * 10
            
            # Ensure minimum of 10 bases
            new_bases = max(10, new_bases)
            
            # Calculate total width needed
            total_width = (new_bases * base_width) + margin
            
            # Only update sequence if bases per line needs to change
            if new_bases != self.sequence_viewer.bases_per_line:
                self.logger.debug(
                    f"Resizing from {self.sequence_viewer.bases_per_line} to {new_bases} bases per line"
                )
                
                # Batch update both viewers
                self.view.setUpdatesEnabled(False)
                self.ruler_view.setUpdatesEnabled(False)
                
                # Update sequence viewer
                self.sequence_viewer.bases_per_line = new_bases
                self.sequence_viewer._create_nucleotide_items()
                
                # Update feature viewer to match
                if hasattr(self, 'feature_viewer'):
                    self.feature_viewer.bases_per_line = new_bases
                    self.feature_viewer.update()
                
                # Calculate total height needed
                sequence_length = len(self.sequence_viewer.sequence)
                total_lines = (sequence_length + new_bases - 1) // new_bases
                total_height = total_lines * self.sequence_viewer.line_spacing
                
                # Set fixed scene rect with left alignment and full height
                scene_rect = QRectF(0, 0, total_width, total_height)
                self.scene.setSceneRect(scene_rect)
                
                # Update ruler with exact viewport width
                self._create_ruler()
                self.ruler_view.setFixedWidth(viewport_width + 10)
                # Set ruler scene rect to exactly match viewport width
                self.ruler_scene.setSceneRect(0, 0, viewport_width + 10, 25)
                
                # Keep scroll positions in sync
                scroll_value = self.view.horizontalScrollBar().value()
                self.ruler_view.horizontalScrollBar().setValue(scroll_value)
                
                # Ensure view shows full content
                self.view.setSceneRect(scene_rect)
                
                self.view.setUpdatesEnabled(True)
                self.ruler_view.setUpdatesEnabled(True)
                
                # Force immediate update
                self.view.viewport().update()
                self.ruler_view.viewport().update()
            else:
                # Even if we don't resize, ensure proper alignment
                self.ruler_view.setFixedWidth(viewport_width + 10)
                # Update ruler scene rect to match viewport exactly
                self.ruler_scene.setSceneRect(0, 0, viewport_width + 10, 25)
                scroll_value = self.view.horizontalScrollBar().value()
                self.ruler_view.horizontalScrollBar().setValue(scroll_value)
                self.ruler_view.viewport().update()

        except Exception as e:
            self.logger.error(f"Error in _handle_resize: {str(e)}")

    def set_data(self, sequence, features, start_pos=None):
        """Set data for both viewers"""
        if start_pos is None:
            start_pos = 0
            
        # Update both components
        self.sequence_viewer.set_data(sequence, start_pos)
        self.feature_viewer.set_data(sequence, features, start_pos)
        
        # Position feature viewer to overlap with sequence viewer
        self.feature_viewer.setY(0)
        
        # Update scene rect to encompass both viewers
        combined_rect = self.sequence_viewer.boundingRect().united(self.feature_viewer.boundingRect())
        self.scene.setSceneRect(combined_rect)
        
        # Update status panel with initial sequence info
        sequence_length = len(sequence)
        self.status_panel.setText(f"Showing: {start_pos}...{start_pos + sequence_length} = {sequence_length} bp")
        
        self.update()

    def _on_sequence_selected(self, start_pos, end_pos):
        """Handle sequence selection"""
        # Calculate selected length
        selected_length = end_pos - start_pos + 1
        
        # Update status panel with selection info
        self.status_panel.setText(f"Selected: {start_pos}...{end_pos} = {selected_length} bp")
        
        # Emit signal for other components
        self.sequence_selected.emit(start_pos, end_pos)

    def clear_selection(self):
        """Clear selection and reset status panel"""
        if hasattr(self, 'sequence_viewer'):
            sequence = self.sequence_viewer.sequence
            start_pos = self.sequence_viewer.start_pos
            sequence_length = len(sequence)
            self.status_panel.setText(f"Showing: {start_pos}...{start_pos + sequence_length} = {sequence_length} bp")

    def _on_cursor_position_changed(self, position):
        """Handle cursor position changes"""
        if position >= 0:
            self.status_panel.setText(f"Insertion Point: {position}")
        else:
            # Reset to showing current sequence range if cursor position is invalid
            if hasattr(self, 'sequence_viewer'):
                sequence = self.sequence_viewer.sequence
                start_pos = self.sequence_viewer.start_pos
                sequence_length = len(sequence)
                self.status_panel.setText(f"Showing: {start_pos}...{start_pos + sequence_length} = {sequence_length} bp")

    def _create_ruler(self):
        """Create ruler with position markers"""
        try:
            self.ruler_scene.clear()
            
            # Get current bases per line
            bases_per_line = self.sequence_viewer.bases_per_line
            base_width = self.sequence_viewer.base_width
            
            # Calculate total width including margin
            total_width = bases_per_line * base_width + 100  # Match sequence viewer width
            
            # Create horizontal blue line aligned with sequence
            ruler_line = QGraphicsLineItem(0, 15, bases_per_line * base_width, 15)
            ruler_line.setPen(QPen(QColor(0, 120, 215), 1))
            self.ruler_scene.addItem(ruler_line)
            
            # Add tick marks and numbers for every base
            for i in range(0, bases_per_line):
                x_pos = i * base_width + base_width/2  # Center tick marks between bases
                
                # Use 1-based indexing for position calculation
                pos_1_based = i + 1
                
                # Determine tick height based on position
                if pos_1_based % 10 == 0:  # Major ticks (every 10)
                    tick_height = 10
                    tick_start = 10
                    # Add number
                    text = QGraphicsSimpleTextItem(str(pos_1_based))
                    text.setFont(QFont("Arial", 8))
                    text_width = text.boundingRect().width()
                    text.setPos(x_pos - text_width/2, 0)  # Position above line
                    self.ruler_scene.addItem(text)
                elif pos_1_based % 5 == 0:  # Medium ticks (every 5)
                    tick_height = 7
                    tick_start = 11
                else:  # Small ticks (every 1)
                    tick_height = 4
                    tick_start = 13
                
                # Create tick mark
                tick = QGraphicsLineItem(x_pos, tick_start, x_pos, tick_start + tick_height)
                tick.setPen(QPen(QColor(0, 120, 215), 1))
                self.ruler_scene.addItem(tick)
            
            # Set scene rect to exactly match sequence viewer width
            self.ruler_scene.setSceneRect(0, 0, total_width, 25)
            
        except Exception as e:
            self.logger.error(f"Error creating ruler: {str(e)}")

class NucleotideItem(QGraphicsObject):
    clicked = pyqtSignal(object)
    
    def __init__(self, nucleotide, x, y, width, is_uppercase=False, is_complement=False, parent=None):
        super().__init__(parent)
        self.nucleotide = nucleotide
        self.spacing = 0
        self.base_width = width
        self.rect = QRectF(0, 0, width, width * 2)
        self.setPos(x, y)
        self.is_uppercase = is_uppercase
        self.is_complement = is_complement
        self.is_highlighted = False
        self.highlight_color = None
        self.show_cursor = False
        self.cursor_side = 'right'
        self.setAcceptHoverEvents(True)
        
        # Get logger from parent
        sequence_viewer = self.parent()
        if sequence_viewer and hasattr(sequence_viewer, 'logger'):
            self.logger = sequence_viewer.logger
        else:
            import logging
            self.logger = logging.getLogger(__name__)

    def boundingRect(self):
        return self.rect

    def paint(self, painter, option, widget):
        try:
            # Draw highlight background if highlighted
            if self.is_highlighted and self.highlight_color:
                painter.fillRect(self.rect, self.highlight_color)
            
            # Draw nucleotide
            painter.setFont(QFont("Courier", 12))
            if self.is_uppercase:
                painter.setPen(Qt.GlobalColor.black)
            else:
                painter.setPen(QColor(100, 100, 100))
            
            # Get complement nucleotide if needed
            display_nucleotide = self.nucleotide
            if self.is_complement:
                complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                                'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
                display_nucleotide = complement_map.get(self.nucleotide, self.nucleotide)
            
            # Draw text centered
            painter.drawText(self.rect, Qt.AlignmentFlag.AlignCenter, display_nucleotide)
            
            # Draw cursor if active
            if self.show_cursor:
                painter.setPen(QPen(Qt.GlobalColor.black, 1))
                if self.cursor_side == 'right':
                    x = self.rect.right() - 1
                else:
                    x = self.rect.left() + 1
                cursor_height = self.rect.height()
                painter.drawLine(QPointF(x, 0), QPointF(x, cursor_height))
                
        except Exception as e:
            self.logger.error(f"Error in paint: {str(e)}")

    def mousePressEvent(self, event):
        try:
            if event.button() == Qt.MouseButton.LeftButton:
                sequence_viewer = self.parent()
                if sequence_viewer:
                    pos = sequence_viewer.get_nucleotide_position(self)
                    local_x = event.pos().x()
                    
                    # Determine cursor position based on click location
                    is_right_side = local_x > self.rect.width() / 2
                    
                    # Always show cursor
                    self.show_cursor = True
                    self.cursor_side = 'right' if is_right_side else 'left'
                    
                    # Calculate cursor position
                    cursor_pos = pos + 1 if is_right_side else pos
                    
                    self.logger.debug(f"""
                    Mouse click details:
                    - Nucleotide: {self.nucleotide}
                    - Click X: {local_x}
                    - Is right side: {is_right_side}
                    - Cursor position: {cursor_pos}
                    """)
                    
                    # Start selection
                    sequence_viewer.selection_active = True
                    sequence_viewer.selection_start = cursor_pos
                    sequence_viewer.selection_end = cursor_pos
                    
                    # Clear other cursors
                    for nuc in sequence_viewer.nucleotides:
                        if nuc != self:
                            nuc.show_cursor = False
                            nuc.update()
                    
                    sequence_viewer.cursor_position_changed.emit(cursor_pos)
                    sequence_viewer._update_selection()
                    self.update()
                    
        except Exception as e:
            self.logger.error(f"Error in mousePressEvent: {str(e)}")
        
        event.accept()

    def mouseMoveEvent(self, event):
        sequence_viewer = self.parent()
        if sequence_viewer and event.buttons() & Qt.MouseButton.LeftButton:
            try:
                pos = sequence_viewer.get_nucleotide_position(self)
                local_x = event.pos().x()
                
                # Calculate position relative to letter boundaries
                text_x = (self.rect.width() - self.base_width) / 2
                relative_x = local_x - text_x
                
                # Force cursor to show
                self.show_cursor = True
                
                # Determine cursor position
                if relative_x <= 0:  # Before letter
                    self.cursor_side = 'left'
                    cursor_pos = pos
                elif relative_x >= self.base_width:  # After letter
                    self.cursor_side = 'right'
                    cursor_pos = pos + 1
                else:  # On letter
                    is_after = relative_x > self.base_width / 2
                    self.cursor_side = 'right' if is_after else 'left'
                    cursor_pos = pos + 1 if is_after else pos
                
                # Update selection
                sequence_viewer.selection_end = pos
                
                # Clear other cursors
                for nuc in sequence_viewer.nucleotides:
                    if nuc != self:
                        nuc.show_cursor = False
                        nuc.update()
                
                sequence_viewer.cursor_position_changed.emit(cursor_pos)
                sequence_viewer._update_selection()
                self.update()  # Force redraw
                
            except Exception as e:
                self.logger.error(f"Error in mouseMoveEvent: {str(e)}")
        
        event.accept()

    def hoverMoveEvent(self, event):
        local_x = event.pos().x()
        mid_point = self.rect.width() / 2
        
        sequence_viewer = self.parent()
        if sequence_viewer:
            # Only show cursor if not selecting
            if not sequence_viewer.selection_active:
                self.show_cursor = True
                self.cursor_side = 'right' if local_x >= mid_point else 'left'
                
                # Clear cursor from other nucleotides
                for nuc in sequence_viewer.nucleotides:
                    if nuc != self:
                        nuc.show_cursor = False
                        nuc.update()
                
                pos = sequence_viewer.get_nucleotide_position(self)
                cursor_pos = pos + 1 if self.cursor_side == 'right' else pos
                sequence_viewer.cursor_position_changed.emit(cursor_pos)
                self.update()

    def hoverLeaveEvent(self, event):
        sequence_viewer = self.parent()
        if sequence_viewer and not sequence_viewer.selection_active:
            self.show_cursor = False
            self.update()
        super().hoverLeaveEvent(event)

class SequenceViewer(QGraphicsObject):
    sequence_selected = pyqtSignal(int, int)
    cursor_position_changed = pyqtSignal(int)  # New signal for cursor position
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence = ""
        self.start_pos = 0
        self.base_width = 15
        self.bases_per_line = 70
        self.line_height = 25
        self.line_spacing = 80  # Increased spacing between gene lines
        self.nucleotides = []
        self.selection_start = None
        self.selection_end = None
        self.setAcceptHoverEvents(True)
        
        # Add clipboard support
        self.clipboard = QApplication.clipboard()
        
        # Add tracking for drag start position
        self.drag_start_pos = None
        self.selection_active = False
        
        # Enable mouse tracking
        self.setAcceptedMouseButtons(Qt.MouseButton.LeftButton)
        self.setFlag(QGraphicsObject.GraphicsItemFlag.ItemIsFocusable)
        
        # Add logger
        import logging
        self.logger = logging.getLogger(__name__)
        
        # Configure logging
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.DEBUG)
        
        # Add highlight state tracking
        self.highlighted_regions = []  # List of (start, end, color, strand) tuples

    def _get_text_width(self):
        """Calculate text width based on current bases per line"""
        return self.base_width * self.bases_per_line

    def set_data(self, sequence, start_pos):
        self.sequence = sequence
        self.start_pos = start_pos
        self._create_nucleotide_items()
        self.update()

    def _create_nucleotide_items(self):
        try:
            # Get the view from the scene
            view = None
            if self.scene():
                views = self.scene().views()
                if views:
                    view = views[0]
                    view.setUpdatesEnabled(False)  # Use view instead of scene
            
            # Clear existing items
            self.cleanup_graphics()
            
            current_pos = 0
            max_width = 0  # Track maximum line width
            
            # Pre-calculate total lines
            total_lines = (len(self.sequence) + self.bases_per_line - 1) // self.bases_per_line
            
            # Pre-allocate lists for better performance
            self.plot_lines = []
            self.tick_lines = []
            self.nucleotides = []
            
            # Track nucleotides by strand and line
            self.nucleotide_map = {
                '+': [],  # List of lists for positive strand nucleotides by line
                '-': []   # List of lists for negative strand nucleotides by line
            }
            current_line = 0
            self.nucleotide_map['+'].append([])  # Initialize first line
            self.nucleotide_map['-'].append([])

            # Create nucleotides and setup display
            self._create_display()

            # Reapply highlights after creating nucleotides
            self._reapply_highlights()
            
        except Exception as e:
            self.logger.error(f"Error in _create_nucleotide_items: {str(e)}")
            # Make sure we re-enable updates even if there's an error
            if view:
                view.setUpdatesEnabled(True)

    def _create_display(self):
        """Create the nucleotide display"""
        try:
            current_pos = 0
            max_width = 0  # Add max_width definition here
            
            while current_pos < len(self.sequence):
                # Calculate exact number of bases for this line
                remaining_bases = len(self.sequence) - current_pos
                bases_this_line = min(self.bases_per_line, remaining_bases)
                line_text = self.sequence[current_pos:current_pos + bases_this_line]
                
                line_num = current_pos // self.bases_per_line
                y_pos = line_num * self.line_spacing

                # Calculate width for position numbers
                max_width = max(max_width, self.bases_per_line * self.base_width)

                # Create positive strand nucleotides
                line_nucleotides_pos = []
                for i, nucleotide in enumerate(line_text):
                    x_pos = i * self.base_width
                    nuc_item = NucleotideItem(
                        nucleotide=nucleotide,
                        x=x_pos,
                        y=y_pos + self.line_height * -0.3,
                        width=self.base_width,
                        is_uppercase=nucleotide.isupper(),
                        parent=self
                    )
                    self.nucleotides.append(nuc_item)
                    line_nucleotides_pos.append(nuc_item)

                # Create complement strand nucleotides
                line_nucleotides_neg = []
                for i, nucleotide in enumerate(line_text):
                    x_pos = i * self.base_width
                    nuc_item = NucleotideItem(
                        nucleotide=nucleotide,
                        x=x_pos,
                        y=y_pos + self.line_height * 1.1,
                        width=self.base_width,
                        is_uppercase=nucleotide.isupper(),
                        is_complement=True,
                        parent=self
                    )
                    self.nucleotides.append(nuc_item)
                    line_nucleotides_neg.append(nuc_item)

                # Store nucleotides by strand and line
                if line_num >= len(self.nucleotide_map['+']):
                    self.nucleotide_map['+'].append([])
                    self.nucleotide_map['-'].append([])
                self.nucleotide_map['+'][line_num].extend(line_nucleotides_pos)
                self.nucleotide_map['-'][line_num].extend(line_nucleotides_neg)

                # Draw plot line matching exactly the sequence width for this line
                plot_y = y_pos + self.line_height
                plot_line = QGraphicsLineItem(0, plot_y,
                                            bases_this_line * self.base_width, plot_y, self)
                plot_line.setPen(QPen(Qt.GlobalColor.black, 1))
                self.plot_lines.append(plot_line)

                # Draw tick marks only for actual bases in this line
                for i in range(bases_this_line):
                    x_pos = i * self.base_width
                    
                    # Convert to 1-based index for position calculation
                    pos_1_based = current_pos + i + 1  # Add 1 for 1-based indexing
                    
                    # Determine tick height based on position
                    if i == 0 and current_pos == 0:  # First base
                        tick_height = 12  # Longest tick for start
                    elif i == bases_this_line - 1 and current_pos + bases_this_line == len(self.sequence):  # Last base
                        tick_height = 12  # Longest tick for end
                    elif pos_1_based % 10 == 0:  # Major ticks (every 10)
                        tick_height = 10
                    elif pos_1_based % 5 == 0:  # Medium ticks (every 5)
                        tick_height = 8
                    else:  # Regular ticks
                        tick_height = 5
                    
                    tick_line = QGraphicsLineItem(
                        x_pos + self.base_width/2, 
                        plot_y - tick_height/2,
                        x_pos + self.base_width/2, 
                        plot_y + tick_height/2,
                        self
                    )
                    self.tick_lines.append(tick_line)

                # Add position number aligned with the plot line - remove the +1
                end_pos = str(self.start_pos + current_pos + bases_this_line)  # Removed +1
                pos_item = QGraphicsSimpleTextItem(end_pos, self)
                pos_item.setFont(QFont("Courier", 12))
                
                # Calculate position for consistent alignment
                text_width = pos_item.boundingRect().width()
                pos_x = max_width + 10  # Fixed position based on maximum width
                pos_y = plot_y - pos_item.boundingRect().height()/2
                pos_item.setPos(pos_x, pos_y)

                current_pos += bases_this_line
            
            # Get parent view if it exists
            view = self.scene().views()[0] if self.scene() and self.scene().views() else None
            if view:
                view.setUpdatesEnabled(True)  # Re-enable updates if we have a view
            self.update()
            
        except Exception as e:
            self.logger.error(f"Error in _create_display: {str(e)}")

    def _reapply_highlights(self):
        """Reapply stored highlights after recreating nucleotides"""
        try:
            # Clear existing highlights from nucleotides
            for nuc in self.nucleotides:
                nuc.is_highlighted = False
                nuc.highlight_color = None

            # Reapply each stored highlight
            for start_pos, end_pos, color, strand in self.highlighted_regions:
                # Calculate which lines contain the sequence
                start_line = start_pos // self.bases_per_line
                end_line = end_pos // self.bases_per_line
                
                # Calculate positions within lines
                start_pos_in_line = start_pos % self.bases_per_line
                end_pos_in_line = end_pos % self.bases_per_line

                # Get the correct strand's nucleotide map
                strand_map = self.nucleotide_map[strand]

                # Handle multi-line sequences
                for line_num in range(start_line, end_line + 1):
                    if line_num >= len(strand_map):
                        continue

                    # Calculate start and end positions for this line
                    if line_num == start_line:
                        line_start = start_pos_in_line
                    else:
                        line_start = 0

                    if line_num == end_line:
                        line_end = end_pos_in_line
                    else:
                        line_end = self.bases_per_line - 1

                    # Get nucleotides for this line segment
                    line_nucleotides = strand_map[line_num]
                    
                    # Calculate the range of nucleotides to highlight
                    start_idx = min(line_start, len(line_nucleotides))
                    end_idx = min(line_end + 1, len(line_nucleotides))
                    
                    # Highlight the nucleotides
                    for i in range(start_idx, end_idx):
                        nuc = line_nucleotides[i]
                        nuc.is_highlighted = True
                        nuc.highlight_color = color
                        nuc.update()

        except Exception as e:
            self.logger.error(f"Error in _reapply_highlights: {str(e)}")

    def _update_selection(self):
        """Update the visual selection"""
        try:
            if self.selection_start is None:
                self.logger.debug("No selection start point")
                return

            start_idx = min(self.selection_start, self.selection_end or self.selection_start)
            end_idx = max(self.selection_start, self.selection_end or self.selection_start)
            
            self.logger.debug(f"Updating selection: start={start_idx}, end={end_idx}")

            # Update highlighting for all nucleotides
            for i, nuc in enumerate(self.nucleotides):
                if start_idx <= i <= end_idx:
                    nuc.is_highlighted = True
                    nuc.highlight_color = QColor(200, 200, 255, 100)
                    self.logger.debug(f"Highlighting nucleotide at position {i}: {nuc.nucleotide}")
                else:
                    nuc.is_highlighted = False
                    nuc.highlight_color = None
                nuc.update()

            # Emit selection signal
            if self.selection_active:
                self.logger.debug(f"Emitting selection signal: {self.start_pos + start_idx} to {self.start_pos + end_idx}")
                self.sequence_selected.emit(
                    self.start_pos + start_idx, 
                    self.start_pos + end_idx
                )
                
        except Exception as e:
            self.logger.error(f"Error in _update_selection: {str(e)}")

    def mousePressEvent(self, event):
        """Handle mouse press for selection start"""
        try:
            # Convert scene position to local coordinates
            local_pos = self.mapFromScene(event.scenePos())
            
            # Calculate base position more precisely
            x_pos = local_pos.x()
            line_number = int(local_pos.y() // (self.line_height * 2))
            
            # Calculate which letter space was clicked
            exact_position = x_pos / self.base_width
            base_position = int(exact_position)
            
            # Calculate if click was in left or right half of the letter space
            fraction = exact_position - base_position
            is_right_side = fraction > 0.5
            
            # Adjust base_position based on where exactly the click occurred
            if is_right_side:
                cursor_index = base_position
                cursor_side = 'right'
            else:
                cursor_index = max(0, base_position - 1)
                cursor_side = 'right' if base_position == 0 else 'left'
            
            # Calculate final index
            index = line_number * self.bases_per_line + cursor_index
            index = max(0, min(index, len(self.nucleotides) - 1))
            
            if 0 <= index < len(self.nucleotides):
                # Start selection
                self.selection_active = True
                self.selection_start = index
                self.selection_end = index
                
                # Update cursor position
                nuc = self.nucleotides[index]
                nuc.show_cursor = True
                nuc.cursor_side = cursor_side
                
                # Clear other cursors
                for i, other_nuc in enumerate(self.nucleotides):
                    if i != index:
                        other_nuc.show_cursor = False
                        other_nuc.update()
                
                # Update selection and cursor
                cursor_pos = index + 1 if cursor_side == 'right' else index
                self.cursor_position_changed.emit(self.start_pos + cursor_pos)
                self._update_selection()
                nuc.update()
                
        except Exception as e:
            print(f"Error in mousePressEvent: {str(e)}")
        
        event.accept()

    def mouseMoveEvent(self, event):
        """Handle mouse move for selection update"""
        if not self.selection_active:
            return
        
        try:
            # Convert scene position to local coordinates
            local_pos = self.mapFromScene(event.scenePos())
            
            # Calculate base position more precisely
            x_pos = local_pos.x()
            line_number = int(local_pos.y() // (self.line_height * 2))
            
            # Calculate which letter space was clicked
            exact_position = x_pos / self.base_width
            base_position = int(exact_position)
            
            # Calculate if mouse is in left or right half of the letter space
            fraction = exact_position - base_position
            is_right_side = fraction > 0.5
            
            # Adjust base_position based on mouse position
            if is_right_side:
                cursor_index = base_position
                cursor_side = 'right'
            else:
                cursor_index = max(0, base_position - 1)
                cursor_side = 'right' if base_position == 0 else 'left'
            
            # Calculate final index
            index = line_number * self.bases_per_line + cursor_index
            index = max(0, min(index, len(self.nucleotides) - 1))
            
            if 0 <= index < len(self.nucleotides):
                # Update selection end point
                self.selection_end = index
                
                # Update cursor position
                nuc = self.nucleotides[index]
                nuc.show_cursor = True
                nuc.cursor_side = cursor_side
                
                # Clear other cursors
                for i, other_nuc in enumerate(self.nucleotides):
                    if i != index:
                        other_nuc.show_cursor = False
                        other_nuc.update()
                
                # Update selection and cursor
                cursor_pos = index + 1 if cursor_side == 'right' else index
                self.cursor_position_changed.emit(self.start_pos + cursor_pos)
                self._update_selection()
                nuc.update()
                
        except Exception as e:
            print(f"Error in mouseMoveEvent: {str(e)}")
        
        event.accept()

    def mouseReleaseEvent(self, event):
        """Handle mouse release"""
        if self.selection_active:
            if self.selection_start is not None and self.selection_end is not None:
                start_pos = min(self.selection_start, self.selection_end)
                end_pos = max(self.selection_start, self.selection_end)
                
                # Get selected sequence
                selected_sequence = self._get_selected_sequence(start_pos, end_pos)
                
                # Copy to clipboard
                self.clipboard.setText(selected_sequence)
                
                # Emit selection signal
                self.sequence_selected.emit(
                    self.start_pos + start_pos, 
                    self.start_pos + end_pos
                )
        else:
            # If no selection was made, keep showing the insertion point
            if hasattr(self, 'drag_start_pos') and self.drag_start_pos is not None:
                self.cursor_position_changed.emit(self.start_pos + self.drag_start_pos)
        
        event.accept()

    def _find_closest_nucleotide(self, pos):
        """Find the closest nucleotide to the given position"""
        closest_item = None
        min_distance = float('inf')
        
        for nuc in self.nucleotides:
            nuc_pos = nuc.scenePos()
            distance = (pos.x() - nuc_pos.x()) ** 2 + (pos.y() - nuc_pos.y()) ** 2
            
            if distance < min_distance:
                min_distance = distance
                closest_item = nuc
                
        return closest_item

    def _get_selected_sequence(self, start_pos, end_pos):
        """Get the sequence of selected nucleotides"""
        selected_nucs = []
        for i in range(start_pos, end_pos + 1):
            if i < len(self.nucleotides):
                selected_nucs.append(self.nucleotides[i].nucleotide)
        return ''.join(selected_nucs)

    def keyPressEvent(self, event):
        """Handle keyboard shortcuts"""
        if event.matches(QKeySequence.StandardKey.Copy):
            if self.selection_start is not None and self.selection_end is not None:
                start_pos = min(self.selection_start, self.selection_end)
                end_pos = max(self.selection_start, self.selection_end)
                selected_sequence = self._get_selected_sequence(start_pos, end_pos)
                self.clipboard.setText(selected_sequence)
        super().keyPressEvent(event)

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

            # Get the correct strand's nucleotide map
            strand_map = self.nucleotide_map[strand]

            # Handle multi-line sequences
            for line_num in range(start_line, end_line + 1):
                if line_num >= len(strand_map):
                    continue

                # Calculate start and end positions for this line
                if line_num == start_line:
                    line_start = start_pos_in_line
                else:
                    line_start = 0

                if line_num == end_line:
                    line_end = end_pos_in_line
                else:
                    line_end = self.bases_per_line - 1

                # Get nucleotides for this line segment
                line_nucleotides = strand_map[line_num]
                
                # Calculate the range of nucleotides to highlight
                start_idx = min(line_start, len(line_nucleotides))
                end_idx = min(line_end + 1, len(line_nucleotides))
                
                # Highlight the nucleotides
                for i in range(start_idx, end_idx):
                    nuc = line_nucleotides[i]
                    nuc.is_highlighted = True
                    nuc.highlight_color = color
                    nuc.update()

                self.logger.debug(
                    f"Highlighted nucleotides on line {line_num} "
                    f"for strand {strand} from {start_idx} to {end_idx}"
                )

        except Exception as e:
            self.logger.error(f"Error in highlight_sequence: {str(e)}")
            self.logger.error(f"Start pos: {start_pos}, End pos: {end_pos}, Strand: {strand}")

    def clear_highlights(self):
        """Clear all highlights"""
        self.highlighted_regions.clear()
        for nuc in self.nucleotides:
            nuc.is_highlighted = False
            nuc.highlight_color = None
            nuc.update()

    def boundingRect(self):
        if not self.sequence:
            return QRectF()
            
        # Calculate exact width based on actual sequence in last line
        last_line_length = len(self.sequence) % self.bases_per_line
        if last_line_length == 0 and len(self.sequence) > 0:
            last_line_length = self.bases_per_line
            
        width = self.base_width * last_line_length + 100  # Add space for position numbers
        
        # Calculate height using line spacing
        total_lines = (len(self.sequence) + self.bases_per_line - 1) // self.bases_per_line
        height = total_lines * self.line_spacing
        
        return QRectF(0, 0, width, height)

    def get_nucleotide_position(self, nucleotide):
        """Get the position of a nucleotide in the sequence"""
        try:
            idx = self.nucleotides.index(nucleotide)
            return self.start_pos + idx
        except ValueError:
            return -1

    def cleanup_graphics(self):
        """Clean up all graphics items"""
        if hasattr(self, 'plot_lines'):
            for line in self.plot_lines:
                if line in self.scene().items():
                    self.scene().removeItem(line)
            self.plot_lines.clear()
        
        if hasattr(self, 'tick_lines'):
            for line in self.tick_lines:
                if line in self.scene().items():
                    self.scene().removeItem(line)
            self.tick_lines.clear()

        for nuc in self.nucleotides:
            if nuc in self.scene().items():
                self.scene().removeItem(nuc)
        self.nucleotides.clear()

        for item in self.scene().items():
            if isinstance(item, QGraphicsSimpleTextItem):
                self.scene().removeItem(item)

class FeatureViewer(QGraphicsObject):
    cursor_position_changed = pyqtSignal(int)  # Add signal for cursor position
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence = ""
        self.features = []
        self.start_pos = 0
        self.base_width = 15
        self.bases_per_line = 70
        self.feature_height = 20
        self.line_height = 25
        self.feature_spacing = 2  # Reduce spacing between features and strands
        self.setAcceptHoverEvents(True)  # Enable hover events

    def set_data(self, sequence, features, start_pos):
        """Updated to accept sequence parameter"""
        self.sequence = sequence
        self.features = sorted(features, key=lambda x: x['start'])
        self.start_pos = start_pos
        self.update()

    def paint(self, painter, option, widget):
        if not self.features or not self.sequence:
            return

        # Process each line of sequence
        current_pos = 0
        while current_pos < len(self.sequence):
            line_text = self.sequence[current_pos:current_pos + self.bases_per_line]
            line_num = current_pos // self.bases_per_line
            
            # Calculate y position to be directly below negative strand
            y_pos = line_num * self.line_height * 2
            feature_y = y_pos + self.line_height * 2  # Position directly below negative strand
            
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
                    
                    # Draw orange rectangle
                    painter.setBrush(QColor(255, 140, 0))
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.drawRect(feature_rect)
                    
                    # Draw label if enough space
                    label = feature.get('name', 'HFL1')
                    text_width = painter.fontMetrics().horizontalAdvance(label)
                    if (x_end - x_start) > text_width:
                        text_x = x_start + ((x_end - x_start) - text_width) / 2
                        text_y = feature_y + self.feature_height/2 + 4
                        painter.setPen(Qt.GlobalColor.white)
                        painter.setFont(QFont("Arial", 8))
                        painter.drawText(QPointF(text_x, text_y), label)
                        
                except Exception as e:
                    if hasattr(self, 'logger'):
                        self.logger.error(f"Error drawing feature: {str(e)}")
                    continue

            current_pos += self.bases_per_line

    def boundingRect(self):
        if not self.sequence:
            return QRectF()
            
        # Calculate exact width based on sequence length
        last_line_length = len(self.sequence) % self.bases_per_line
        if last_line_length == 0:
            last_line_length = self.bases_per_line
        width = max(self.base_width * self.bases_per_line, 
                   self.base_width * last_line_length) + 100
        
        # Calculate height for actual sequence lines only
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