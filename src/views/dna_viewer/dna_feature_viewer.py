from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QGraphicsView, QGraphicsScene, 
                            QLabel, QFrame)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor 
from .components.sequence_viewer import SequenceViewer
from .components.feature_viewer import FeatureViewer
from .components.ruler import Ruler
from .components.sequence_insertion_zone import SequenceInsertionZone
import logging
import traceback

class DNAFeatureViewer(QWidget):
    """Main widget for displaying DNA sequences"""
    sequence_selected = pyqtSignal(int, int)  # Emit start and end positions
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Get logger from parent or global settings
        if parent and hasattr(parent, 'logger'):
            self.logger = parent.logger
        else:
            self.logger = logging.getLogger(__name__)
        
        # Create layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)
        
        # Initialize components
        self._init_views()
        self._init_components()
        self._init_status_panel()
        self._init_connections()
        
        # Set focus policy for the widget itself
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFocus()

    def _init_views(self):
        """Initialize graphics views"""
        # Create ruler view
        self.ruler_view = QGraphicsView()
        self.ruler_scene = Ruler()
        self.ruler_view.setScene(self.ruler_scene)
        self.ruler_view.setFixedHeight(25)
        self.ruler_view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.ruler_view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.ruler_view.setViewportMargins(0, 0, 0, 0)
        self.ruler_view.setFrameStyle(0)
        self.ruler_view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.ruler_view.setEnabled(False)  # Disable all user interaction with ruler
        
        # Set background color for both view and viewport
        background_color = QColor(240, 240, 240)
        self.ruler_view.setBackgroundBrush(QBrush(background_color))
        self.ruler_view.viewport().setStyleSheet(f"background-color: rgb({background_color.red()}, {background_color.green()}, {background_color.blue()})")
        self.ruler_view.setAutoFillBackground(True)
        
        # Create main view with proper event handling
        self.view = QGraphicsView()
        self.scene = QGraphicsScene(self)
        self.view.setScene(self.scene)
        
        # Add these lines to remove the frame from main view
        self.view.setFrameStyle(0)  # Remove frame
        self.view.setViewportMargins(0, 0, 0, 0)  # Remove margins
        
        # Set alignment to force left anchoring
        self.view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        
        # Make sure scene events are handled
        self.scene.setItemIndexMethod(QGraphicsScene.ItemIndexMethod.NoIndex)

    def _init_components(self):
        try:
            # Create sequence and feature viewers 
            if not hasattr(self, 'sequence_viewer'):
                self.sequence_viewer = SequenceViewer(logger=self.logger)
                # Reduce left margin/padding
                self.sequence_viewer.strand_margin = 20 
                self.scene.addItem(self.sequence_viewer)
            
            if not hasattr(self, 'feature_viewer'):
                # Create feature viewer and add to scene
                self.feature_viewer = FeatureViewer()
                # Match margin with sequence viewer
                self.feature_viewer.strand_margin = 20 
                self.scene.addItem(self.feature_viewer)
            
            # Add sequence insertion zone
            if not hasattr(self, 'insertion_zone'):
                self.insertion_zone = SequenceInsertionZone()
                # Match margin with sequence viewer
                self.insertion_zone.strand_margin = 20 
            self.scene.addItem(self.insertion_zone)
            
            # Add views to layout
            self.layout.insertWidget(0, self.ruler_view)
            self.layout.addWidget(self.view)

        except Exception as e:
            self.logger.error(f"Error in _init_components: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def _init_status_panel(self):
        """Initialize status panel"""
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
        self.layout.addWidget(self.status_panel)

    def _init_connections(self):
        """Initialize signal connections"""
        # Connect sequence viewer signals
        self.sequence_viewer.sequence_selected.connect(self._on_sequence_selected)
        self.sequence_viewer.cursor_position_changed.connect(self._on_cursor_position_changed)
        
        # Connect feature viewer signals
        self.feature_viewer.cursor_position_changed.connect(self._on_cursor_position_changed)
        
        # Add viewport resize handling
        self.view.viewport().installEventFilter(self)

    def set_data(self, sequence, features=None, start_pos=None):
        """Set sequence and feature data"""
        try:
            if start_pos is None:
                start_pos = 0
            
            # self.logger.debug(f"Features: {features[:2] if features else None}")
            
            # Store original start position for status panel
            self._original_start_pos = start_pos

            # Update components
            self.sequence_viewer.set_data(sequence, start_pos)
            if features is not None:
                self.feature_viewer.set_data(sequence, features, start_pos)
            
            # Update insertion zone
            self.insertion_zone.create_zones(
                sequence_length=len(sequence),
                base_width=self.sequence_viewer.base_width,
                strand_margin=self.sequence_viewer.strand_margin,
                line_height=self.sequence_viewer.line_height,
                bases_per_line=self.sequence_viewer.bases_per_line,
                line_spacing=self.sequence_viewer.line_spacing
            )
            
            # Position components
            self.feature_viewer.setY(0)
            self.insertion_zone.setY(0)
            
            # Update scene rect
            combined_rect = self.sequence_viewer.boundingRect().united(
                self.feature_viewer.boundingRect()
            ).united(
                self.insertion_zone.boundingRect()
            )
            self.scene.setSceneRect(combined_rect)
            
            # Update ruler
            self.ruler_scene.create_ruler(self.sequence_viewer.bases_per_line)
            
            # Update status panel with original gene positions
            sequence_length = len(sequence) if sequence else 0
            self.status_panel.setText(
                f"Showing: {self._original_start_pos}...{self._original_start_pos + sequence_length} = {sequence_length} bp"
            )
            
            self.update()
            self.logger.debug("Finished setting data")
            
        except Exception as e:
            self.logger.error(f"Error setting data: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")

    def _on_sequence_selected(self, start_pos, end_pos):
        """Handle sequence selection"""
        try:
            # Calculate selected length
            selected_length = end_pos - start_pos + 1
            
            # Update only the status panel
            self.status_panel.setText(f"Selected: {start_pos}...{end_pos} = {selected_length} bp")
            
            # Emit signal without updating line edits
            self.sequence_selected.emit(start_pos, end_pos)
            
        except Exception as e:
            self.logger.error(f"Error handling sequence selection: {str(e)}")

    def _on_cursor_position_changed(self, position):
        """Handle cursor position changes"""
        try:
            if position >= 0:
                # Convert position to be relative to gene's actual start position
                absolute_position = self._original_start_pos + position
                self.status_panel.setText(f"Insertion Point: {absolute_position}")
                self.logger.debug(f"Cursor position changed to absolute position: {absolute_position}")
            else:
                # Reset to showing current sequence range
                if hasattr(self, 'sequence_viewer'):
                    sequence = self.sequence_viewer.sequence
                    sequence_length = len(sequence)
                    self.status_panel.setText(
                        f"Showing: {self._original_start_pos}...{self._original_start_pos + sequence_length} = {sequence_length} bp"
                    )
        except Exception as e:
            self.logger.error(f"Error handling cursor position change: {str(e)}")

    def eventFilter(self, obj, event):
        """Handle viewport resize events"""
        if obj == self.view.viewport() and event.type() == event.Type.Resize:
            try:
                viewport_width = event.size().width()
                margin = 100 
                
                available_width = viewport_width - margin
                base_width = self.sequence_viewer.base_width
                
                max_bases = (available_width // base_width)
                new_bases = (max_bases // 10) * 10
                new_bases = max(10, min(new_bases, 200))
                
                if new_bases != self.sequence_viewer.bases_per_line:
                    # Store current sequence and features
                    current_sequence = self.sequence_viewer.sequence
                    current_start = self.sequence_viewer.start_pos
                    
                    # Store only guide highlights (red/green), not selection highlights (blue)
                    current_highlights = []
                    selection_blue = QColor(100, 150, 255, 100)
                    
                    for nuc in self.sequence_viewer.nucleotides:
                        if nuc.is_highlighted and nuc.highlight_color != selection_blue:
                            idx = self.sequence_viewer.get_nucleotide_position(nuc)
                            current_highlights.append({
                                'position': idx // 2,
                                'color': nuc.highlight_color
                            })
                    
                    # Store cursor position
                    cursor_sequence_pos = None
                    if hasattr(self.insertion_zone, 'current_cursor_pos'):
                        cursor_sequence_pos = self.insertion_zone.current_cursor_pos

                    try:
                        # Update bases per line
                        self.sequence_viewer.bases_per_line = new_bases
                        self.feature_viewer.bases_per_line = new_bases
                        
                        # Update components
                        self.set_data(current_sequence, None, current_start)
                        
                        # Reapply only guide highlights
                        for highlight in current_highlights:
                            pos = highlight['position']
                            pos_strand_idx = pos * 2
                            neg_strand_idx = pos * 2 + 1
                            
                            for idx in [pos_strand_idx, neg_strand_idx]:
                                if idx < len(self.sequence_viewer.nucleotides):
                                    nuc = self.sequence_viewer.nucleotides[idx]
                                    nuc.is_highlighted = True
                                    nuc.highlight_color = highlight['color']
                                    nuc.update()
                        
                        # Restore cursor position
                        if cursor_sequence_pos is not None:
                            # Calculate new visual position based on sequence position
                            line_number = cursor_sequence_pos // new_bases
                            pos_in_line = cursor_sequence_pos % new_bases
                            
                            # Calculate exact pixel coordinates for cursor
                            cursor_x = self.sequence_viewer.strand_margin + (pos_in_line * self.sequence_viewer.base_width)
                            cursor_y = (line_number * self.sequence_viewer.line_spacing) + (self.sequence_viewer.line_height * 0.1)
                            cursor_height = self.sequence_viewer.line_height * 2 + 5
                            
                            # Update cursor position
                            if hasattr(self.insertion_zone, 'sequence_cursor'):
                                self.insertion_zone.sequence_cursor.set_position(cursor_x, cursor_y, cursor_height)
                                self.insertion_zone.sequence_cursor.show()
                                self.insertion_zone.current_cursor_pos = cursor_sequence_pos
                        
                        # Update ruler
                        self.ruler_scene.create_ruler(new_bases)
                        
                        # Update scroll positions
                        scroll_value = self.view.horizontalScrollBar().value()
                        self.ruler_view.horizontalScrollBar().setValue(scroll_value)
                        
                    finally:
                        # Force immediate update
                        self.view.viewport().update()
                        self.ruler_view.viewport().update()
                        
            except Exception as e:
                self.logger.error(f"Error handling resize event: {str(e)}")
                
        return super().eventFilter(obj, event)