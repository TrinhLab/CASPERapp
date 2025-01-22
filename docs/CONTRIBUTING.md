## Project Overview

- MVC Pattern: The application follows the Model-View-Controller (MVC) design pattern, which separates the application logic (Model), user interface (View), and user interaction (Controller).
- PyQt6 Framework: The application uses PyQt6 for building the GUI.
- Logging: The application uses logging to track events and errors, which helps in debugging and monitoring the application's behavior.

## Codebase Overview

#### 1. **Controllers**
Controllers are responsible for handling user interactions and updating the views and models accordingly. They act as intermediaries between the user interface and the data logic.

- **CoTargetingController**: Manages co-targeting functionalities, handling user inputs and updating the view.
- **ExportSelectedgRNAsController**: Handles exporting selected gRNA sequences to a file.
- **FindTargetsController**: Manages the process of finding targets based on user input.
- **GenerateLibraryController**: Handles the generation of a library of targets.
- **HomeWindowController**: Manages the home tab and its interactions.
- **MainWindowController**: Controls the main application window, managing tabs and user interactions.
- **MultitargetingWindowController**: Handles multitargeting analysis functionalities.
- **NCBIWindowController**: Manages interactions with the NCBI database for downloading genomic data.
- **NewEndonucleaseController**: Handles the creation and management of new endonucleases.
- **NewGenomeWindowController**: Manages the process of adding new genomes to the application.
- **OffTargetController**: Handles off-target analysis functionalities.
- **PopulationAnalysisWindowController**: Manages population analysis features.
- **ScoringOptionsController**: Handles scoring options and configurations.
- **StartupWindowController**: Manages the startup window and initial configurations.

#### 2. **Models**
Models are responsible for handling data and business logic. They interact with databases, perform calculations, and manage the state of the application.

- **AnnotationParser**: Parses and manages genomic annotation data.
- **BaseModel**: Provides a base class for other models, handling common functionalities.
- **Other Models**: Each controller typically has a corresponding model that handles specific data operations.

#### 3. **Views**
Views are responsible for displaying the user interface. They define the layout and appearance of the application windows and dialogs.

- **Each controller has a corresponding view** that defines the UI components and layout for that part of the application.

#### 4. **Utilities**
Utility modules provide helper functions and classes that are used throughout the application.

- **ui.py**: Contains utility functions for UI operations, such as showing messages and errors.
- **web.py**: Provides functions for web-related operations, like opening URLs.

#### 5. **Main Application Entry Points**
- **main.py**: Contains the main function that initializes the application, checks dependencies, and starts the PyQt6 application loop.

## Advices for New Developers
   - Use try-catch blocks for file operations and data processing
   - Log errors appropriately using the logger
   -  Check the app.log file for error information.
   - Show error messages through the GUI

## Running CASPER
1) Clone the repository using the command: `git clone https://github.com/TrinhLab/CASPERapp`
2) Ensure that you have Conda installed on your computer.
3) Create a new Conda environment with Python 3.11: `conda create --name casper_env python=3.11`
4) Activate the Conda environment: `conda activate casper_env`
5) cd into the CASPERapp directory: `cd CASPERapp`
7) Run CASPER! `python3 main.py`