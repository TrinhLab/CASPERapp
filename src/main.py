import sys
import os
import platform
import subprocess
import importlib
import importlib.metadata
import concurrent.futures

def check_dependencies():
    """Check and install all required dependencies"""
    try:
        # Read requirements from requirements.txt
        requirements_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'requirements.txt')
        
        if not os.path.exists(requirements_path):
            return True  # Skip if requirements.txt doesn't exist
            
        with open(requirements_path, 'r') as f:
            requirements = [line.strip().split('==')[0] for line in f.readlines() 
                          if line.strip() and not line.startswith('#')]
        
        # Check which packages are missing
        missing = []
        for package in requirements:
            try:
                importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError:
                missing.append(package)
        
        if missing:
            print("\nThe following packages are required but not installed:")
            print(", ".join(missing))
            response = input("\nWould you like to install them now? (y/n): ").lower().strip()
            
            if response != 'y':
                print("\nWarning: The application requires these packages to run.")
                print("You can install them manually using:")
                print(f"pip install {' '.join(missing)}")
                return False
            
            print("\nInstalling required packages...")
            try:
                for package in missing:
                    print(f"\nInstalling {package}...")
                    subprocess.check_call([sys.executable, "-m", "pip", "install", package])
                print("\nAll packages installed successfully!")
                return True
            except subprocess.CalledProcessError as e:
                print(f"\nError: Failed to install packages.")
                print(f"Please install them manually using:")
                print(f"pip install {' '.join(missing)}")
                return False
                
        return True
        
    except Exception as e:
        print(f"\nError during dependency check: {str(e)}")
        return False

def main():
    # Check all dependencies first
    if not check_dependencies():
        sys.exit(1)
    
    # Now we can safely import PyQt6 and other modules
    from PyQt6.QtWidgets import QApplication, QMessageBox
    from PyQt6.QtCore import Qt, QCoreApplication
    
    # Enable Qt's built-in caching mechanisms BEFORE creating QApplication
    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    
    RESTART_CODE = 1000
    
    while True:
        # Create a single QApplication instance that will be used throughout
        app = QApplication(sys.argv)
        app.setOrganizationName("TrinhLab-UTK")
        app.setApplicationName("CASPER")

        # Enable high DPI scaling
        if hasattr(Qt.ApplicationAttribute, 'AA_UseHighDpiPixmaps'):
            app.setAttribute(Qt.ApplicationAttribute.AA_UseHighDpiPixmaps)
            app.setAttribute(Qt.ApplicationAttribute.AA_EnableHighDpiScaling, True)
        
        # Now we can safely import other modules
        from models.GlobalSettings import GlobalSettings
        from utils.ui import show_error
        
        # Preload modules in parallel
        preload_modules()
        
        # Get the application directory
        app_dir_path = get_app_directory()

        try:
            global_settings = GlobalSettings(app_dir_path)
            
            # Import here after preloading
            from controllers.MainWindowController import MainWindowController
            main_window_controller = MainWindowController(global_settings)
            global_settings.set_main_window(main_window_controller)
            main_window_controller.show()

            exit_code = app.exec()

            if exit_code != RESTART_CODE:
                sys.exit(exit_code)
                break

            main_window_controller = None
            global_settings = None
            app = None
        except Exception as e:
            show_error(global_settings, "An error occurred during application initialization", e)

def preload_modules():
    """Preload commonly used modules in parallel"""
    modules_to_load = [
        'PyQt6.QtWidgets',
        'PyQt6.QtCore',
        'PyQt6.QtGui',
        'controllers.MainWindowController',
        'views.MainWindowView',
        'models.MainWindowModel',
        'models.DatabaseManager',
        'models.ConfigManager'
    ]
    
    def import_module(module_name):
        try:
            importlib.import_module(module_name)
            return True, module_name
        except Exception as e:
            return False, f"Failed to load {module_name}: {str(e)}"
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(import_module, modules_to_load)

def get_app_directory():
    """Determine the application root directory based on whether we're frozen or not"""
    if hasattr(sys, 'frozen'):
        if platform.system() == 'Darwin':  # macOS
            bundle_dir = os.path.abspath(os.path.dirname(sys.executable))
            return os.path.join(os.path.dirname(os.path.dirname(bundle_dir)), 'Contents', 'Resources')
        else:
            return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if __name__ == '__main__':
    main()

