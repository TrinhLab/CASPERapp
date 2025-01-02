from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton

class BaseInsertionDialog(QDialog):
    """Dialog for inserting base pairs"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Insert Base Pairs")
        self.setModal(True)
        
        self.valid_bases = set('ATGCRYMKSWBDHVNatgcrymkswbdhvn')
        
        layout = QVBoxLayout(self)

        label = QLabel(
            "Enter base pairs to insert:\n"
            "A (Adenine), T (Thymine), G (Guanine), C (Cytosine)\n"
            "R (A/G), Y (C/T), M (A/C), K (G/T), S (G/C), W (A/T)\n"
            "H (A/C/T), B (G/C/T), V (G/C/A), D (G/A/T), N (Any)"
        )
        layout.addWidget(label)
        
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText("e.g., ATGC, RYKMSWBDHVN")
        self.input_field.textChanged.connect(self._validate_input)
        layout.addWidget(self.input_field)
        
        button_layout = QHBoxLayout()
        
        self.insert_button = QPushButton("Insert")
        self.insert_button.clicked.connect(self.accept)
        self.insert_button.setEnabled(False)  # Disabled until valid input
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        
        button_layout.addWidget(self.insert_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)
        
        self.setMinimumWidth(400)
    
    def _validate_input(self, text):
        """Validate input and filter invalid characters"""
        # Remove any invalid characters
        valid_text = ''.join(c for c in text if c in self.valid_bases)
        
        # If text changed, update the field
        if valid_text != text:
            self.input_field.setText(valid_text)
        
        # Enable/disable insert button based on input
        self.insert_button.setEnabled(bool(valid_text))
        
    def get_bases(self):
        """Return the entered base pairs preserving case"""
        return self.input_field.text() 