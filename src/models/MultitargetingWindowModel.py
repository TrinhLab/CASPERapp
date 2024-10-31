import os

class MultitargetingWindowModel:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        
        self.cspr_file = ""
        self.db_file = ""
        self.row_limit = 1000
        
        # Get organism and endo mappings from DatabaseManager
        self.organisms_to_files, self.organisms_to_endos = self.settings.db_manager.get_organisms_and_endos()

    def get_organisms(self):
        """Get list of available organisms"""
        return list(self.organisms_to_endos.keys())

    def get_endos_for_organism(self, organism):
        """Get available endonucleases for given organism"""
        return self.organisms_to_endos.get(organism, [])

    def set_files(self, organism, endo):
        """Set the CSPR and DB files for analysis"""
        if not organism or not endo:
            self.logger.error("Organism or endonuclease not specified")
            raise ValueError("Organism and endonuclease must be specified")

        self.cspr_file = self._get_cspr_file_path(organism, endo)
        self.db_file = self._get_db_file_path(organism, endo)
        
        if not self.cspr_file or not self.db_file:
            raise FileNotFoundError("Required files not found")

    def get_repeats_data(self):
        """Get repeats data for the seeds table"""
        if not self.db_file:
            raise ValueError("Database file not set. Please select an organism and endonuclease first.")
        return self.settings.db_manager.get_repeats_data(self.db_file, self.row_limit)

    def get_seed_data(self, seed):
        """Get detailed data for a specific seed"""
        return self.settings.db_manager.get_seed_data(self.db_file, seed)

    def get_chro_bar_data(self, seed):
        """Get chromosome distribution data for a seed"""
        return self.settings.db_manager.get_chro_bar_data(self.db_file, seed)

    def get_seeds_vs_repeats_data(self):
        """Get data for seeds vs repeats plot"""
        return self.settings.db_manager.get_seeds_vs_repeats_data(self.db_file)

    def get_repeats_vs_seeds_data(self):
        """Get data for repeats vs seeds plot"""
        return self.settings.db_manager.get_repeats_vs_seeds_data(self.db_file)

    def calculate_statistics(self):
        """Calculate global statistics"""
        return self.settings.db_manager.calculate_statistics(self.db_file)

    def get_sql_settings(self):
        """Get current SQL query settings"""
        return {
            'row_limit': self.row_limit
        }

    def update_sql_settings(self, settings):
        """Update SQL query settings"""
        if 'row_limit' in settings:
            self.row_limit = settings['row_limit']

    # Keep the file path methods as they are specific to this model
    def _get_cspr_file_path(self, organism, endo):
        """Get path to CSPR file for organism/endo combination"""
        try:
            if organism not in self.organisms_to_files or endo not in self.organisms_to_files[organism]:
                self.logger.error(f"No CSPR file mapping found for {organism} with {endo}")
                return None
                
            file_name = self.organisms_to_files[organism][endo][0]
            cspr_path = os.path.join(self.settings.get_db_path(), file_name)
            
            if not os.path.exists(cspr_path):
                self.logger.error(f"CSPR file not found: {cspr_path}")
                return None
                
            return cspr_path
        except Exception as e:
            self.logger.error(f"Error getting CSPR file path: {str(e)}")
            return None

    def _get_db_file_path(self, organism, endo):
        """Get path to DB file for organism/endo combination"""
        try:
            if organism not in self.organisms_to_files or endo not in self.organisms_to_files[organism]:
                self.logger.error(f"No DB file mapping found for {organism} with {endo}")
                return None
                
            file_name = self.organisms_to_files[organism][endo][1]
            db_path = os.path.join(self.settings.get_db_path(), file_name)
            
            if not os.path.exists(db_path):
                self.logger.error(f"Database file not found: {db_path}")
                return None
                
            return db_path
        except Exception as e:
            self.logger.error(f"Error getting database file path: {str(e)}")
            return None

    def get_kstats(self):
        """Get kstats from CSPR file"""
        try:
            if not self.cspr_file:
                self.logger.error("CSPR file not set")
                raise ValueError("CSPR file not set. Please select an organism and endonuclease first.")
                
            kstats = []
            with open(self.cspr_file, "r") as f:
                for line in f:
                    if "KARYSTATS" in line:
                        kstats = line.replace("KARYSTATS: ", "").strip().split(',')[:-1]
                        break
                        
            if not kstats:
                raise ValueError("No KARYSTATS found in CSPR file")
                
            return kstats
            
        except Exception as e:
            self.logger.error(f"Error getting kstats: {str(e)}")
            raise
