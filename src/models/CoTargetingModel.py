class CoTargetingModel:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        self.endo_data = {}
        self.get_endo_data()

    def get_endo_data(self):
        """Load endonuclease data from CASPERinfo file"""
        try:
            f = open(self.settings.get_casper_info_path())
            while True:
                line = f.readline()
                if line.startswith('ENDONUCLEASES'):
                    while True:
                        line = f.readline()
                        if line[0] == "-":
                            break
                        line_tokened = line.split(";")
                        endo = line_tokened[0]
                        self.endo_data[endo] = ([line_tokened[2], line_tokened[3], line_tokened[4]], line_tokened[5])
                    break
            f.close()
        except Exception as e:
            self.logger.error(f"Error loading endonuclease data: {str(e)}")
            raise

    def validate_endonucleases(self, endo_list):
        """Validate that selected endonucleases are compatible"""
        try:
            for endo1 in endo_list:
                for endo2 in endo_list:
                    if endo1 == endo2:
                        continue
                    # Check gRNA length compatibility
                    endo1_len = sum([int(x) for x in self.endo_data[endo1][0]])
                    endo2_len = sum([int(x) for x in self.endo_data[endo2][0]])
                    
                    # Check directionality compatibility
                    if (endo1_len != endo2_len or 
                        self.endo_data[endo1][1] != self.endo_data[endo2][1]):
                        return False
            return True
        except Exception as e:
            self.logger.error(f"Error validating endonucleases: {str(e)}")
            return False

    def format_endo_combination(self, endo_list):
        """Format selected endonucleases into combined string"""
        return "|".join(endo_list)
