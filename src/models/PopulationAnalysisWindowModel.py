import os
import sqlite3
import itertools
from utils.ui import show_error

class PopulationAnalysisWindowModel:
    def __init__(self, global_settings):
        self.settings = global_settings
        self.logger = self.settings.get_logger()
        self.app_dir = self.settings.get_app_dir_path()
        self.cspr_files = []
        self.db_files = []
        self.org_names = {}
        self.seeds = []
        self.counts = []
        self.index_to_cspr = {}
        self.index_to_db = {}

    def load_endonucleases(self):
        """Load endonucleases from GlobalSettings"""
        try:
            self.logger.info("Starting load_endonucleases()")
            
            # Get endonucleases from global settings
            endos = self.settings.get_endonucleases()
            self.logger.debug(f"Raw endonucleases from settings: {endos}")
            
            if not endos:
                self.logger.warning("No endonucleases returned from settings")
                return {}
            
            # Format the endonucleases for display
            formatted_endos = {}
            for endo, data in endos.items():
                self.logger.debug(f"Processing endo: {endo}, data: {data}")
                pam = data.get('pam', '').strip()
                # Remove any extra "PAM:" text that might be in the PAM string
                pam = pam.replace('PAM:', '').strip()
                # Create display name without duplicate "PAM:" text
                display_name = f"{endo}"
                
                formatted_endos[display_name] = (endo, pam, 
                                              data.get('default_five_length', ''),
                                              data.get('default_seed_length', ''),
                                              data.get('default_three_length', ''))
            
            self.logger.info(f"Successfully formatted {len(formatted_endos)} endonucleases")
            self.logger.debug(f"Formatted endonucleases: {formatted_endos}")
            return formatted_endos
            
        except Exception as e:
            self.logger.error(f"Error loading endonucleases: {str(e)}")
            self.logger.exception("Full traceback:")
            show_error(self.settings, "Error loading endonucleases", str(e))
            return {}

    def get_organism_files(self, endo_display_name):
        """Get organism files for selected endonuclease using DatabaseManager"""
        org_files = []
        try:
            # Extract just the endonuclease name from the display text (remove PAM)
            endo = endo_display_name.split(" - PAM:")[0].strip()
            self.logger.info(f"Getting organism files for endonuclease: {endo}")
            
            # Get organism mappings from database manager
            organisms_to_files, organisms_to_endos = self.settings.db_manager.get_organisms_and_endos()
            
            # Process each organism that has this endonuclease
            for organism, endos in organisms_to_endos.items():
                if endo in endos:
                    cspr_file = os.path.join(self.settings.CSPR_DB, organisms_to_files[organism][endo][0])
                    db_file = os.path.join(self.settings.CSPR_DB, organisms_to_files[organism][endo][1])
                    
                    if not os.path.exists(db_file):
                        self.logger.warning(f"Database file not found: {db_file}")
                        continue
                    
                    org_files.append((organism, cspr_file, db_file))
                    
                    # Store the mapping for later use
                    index = len(org_files) - 1
                    self.index_to_cspr[index] = cspr_file
                    self.index_to_db[index] = db_file
            
            self.logger.info(f"Found {len(org_files)} organism files")
        except Exception as e:
            self.logger.error(f"Error getting organism files: {str(e)}")
            show_error(self.settings, "Error getting organism files", str(e))
        return org_files

    def get_shared_seeds(self, db_files, limit=False):
        try:
            aliases = [f"main{i}" for i in range(1, len(db_files) + 1)]
            
            new_conn = sqlite3.connect(os.path.join(self.app_dir, "temp_join.db"))
            new_c = new_conn.cursor()
            new_c.execute("PRAGMA synchronous = OFF;")
            new_c.execute("PRAGMA journal_mode = OFF;")
            new_c.execute("PRAGMA locking_mode = EXCLUSIVE;")
            new_c.execute("DROP TABLE IF EXISTS repeats;")
            new_c.execute("VACUUM;")
            new_c.execute("DROP TABLE IF EXISTS join_results;")
            new_c.execute("CREATE table join_results (seed TEXT PRIMARY KEY);")

            for i, db_file in enumerate(db_files):
                new_c.execute(f"ATTACH DATABASE '{db_file}' AS {aliases[i]};")

            new_c.execute("BEGIN TRANSACTION;")

            sql_inner_join = "INSERT into main.join_results select main1.repeats.seed from main1.repeats "
            for i in range(len(aliases[:-1])):
                sql_inner_join += f"inner join {aliases[i + 1]}.repeats on {aliases[i]}.repeats.seed = {aliases[i + 1]}.repeats.seed "

            new_c.execute(sql_inner_join)

            if limit:
                shared_seeds = new_c.execute("select * from join_results limit 0,1000").fetchall()
            else:
                shared_seeds = new_c.execute("select count(*) from join_results").fetchall()

            new_c.execute("END TRANSACTION;")
            new_c.close()
            new_conn.close()

            return [seed[0] for seed in shared_seeds] if limit else shared_seeds[0][0]

        except Exception as e:
            show_error(self.global_settings, "Error in get_shared_seeds()", str(e))
            return [] if limit else 0

    def get_seed_data(self, seed, db_files):
        data = {'total_count': 0, 'org_count': 0, 'threes': [], 'fives': [], 'scores': [], 'pams': [], 'locs': []}
        try:
            for db_file in db_files:
                with sqlite3.connect(db_file) as conn:
                    c = conn.cursor()
                    result = c.execute("SELECT count, three, five, pam, score, location FROM repeats WHERE seed = ?", (seed,)).fetchone()
                    if result:
                        data['org_count'] += 1
                        data['total_count'] += int(result[0])
                        data['threes'].extend(result[1].split(","))
                        data['fives'].extend(result[2].split(","))
                        data['pams'].extend(result[3].split(","))
                        data['scores'].extend(result[4].split(","))
                        data['locs'].extend(result[5].split(","))
        except Exception as e:
            show_error(self.global_settings, f"Error getting data for seed {seed}", str(e))
        return data

    def get_heatmap_data(self, db_files):
        try:
            size = len(db_files)
            arr = [[0 for _ in range(size)] for _ in range(size)]

            for i, j in itertools.combinations(range(size), 2):
                shared_seeds = self.get_shared_seeds([db_files[i], db_files[j]])
                arr[i][j] = arr[j][i] = shared_seeds

            for i in range(size):
                with sqlite3.connect(db_files[i]) as conn:
                    c = conn.cursor()
                    arr[i][i] = c.execute("SELECT COUNT(*) FROM repeats").fetchone()[0]

            return arr
        except Exception as e:
            show_error(self.global_settings, "Error generating heatmap data", str(e))
            return []

    def get_seed_locations(self, seeds, db_files):
        locations = []
        try:
            for db_file in db_files:
                # Get organism name from CSPR file
                cspr_file = db_file.replace("_repeats.db", ".cspr")
                with open(cspr_file, 'r') as f:
                    organism_name = f.readline().split(":")[-1].strip()

                with sqlite3.connect(db_file) as conn:
                    c = conn.cursor()
                    for seed in seeds:
                        data = c.execute("SELECT chromosome, location, five, three FROM repeats WHERE seed = ?", (seed,)).fetchone()
                        if data:
                            chroms = data[0].split(',')
                            locs = data[1].split(',')
                            fives = data[2].split(',') if data[2] else []
                            threes = data[3].split(',') if data[3] else []
                            for i, chrom in enumerate(chroms):
                                sequence = (fives[i] if fives else '') + seed + (threes[i] if threes else '')
                                locations.append({
                                    'seed': seed,
                                    'sequence': sequence,
                                    'organism': organism_name, 
                                    'chromosome': chrom,
                                    'location': abs(int(locs[i]))
                                })
            self.logger.debug(f"Found {len(locations)} locations")
        except Exception as e:
            self.logger.error(f"Error getting seed locations: {str(e)}")
            show_error(self.global_settings, "Error getting seed locations", str(e))
        return locations

    def get_org_names(self):
        try:
            self.org_names = {}
            for i, cspr_file in self.index_to_cspr.items():
                with open(cspr_file, 'r') as f:
                    org_name = f.readline().split(":")[-1].strip()
                    self.org_names[i] = org_name
            self.logger.info(f"Loaded {len(self.org_names)} organism names")
        except Exception as e:
            self.logger.error(f"Error getting organism names: {str(e)}")
            show_error(self.settings, "Error getting organism names", str(e))

