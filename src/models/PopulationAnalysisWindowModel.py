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
            # Get endonucleases from global settings
            endos = self.settings.get_endonucleases()
            
            if not endos:
                self.logger.warning("No endonucleases returned from settings")
                return {}
            
            # Format the endonucleases for display
            formatted_endos = {}
            for endo, data in endos.items():
                pam = data.get('pam', '').strip()
                # Remove any extra "PAM:" text that might be in the PAM string
                pam = pam.replace('PAM:', '').strip()
                # Create display name without duplicate "PAM:" text
                display_name = f"{endo}"
                
                formatted_endos[display_name] = (endo, pam, 
                                              data.get('default_five_length', ''),
                                              data.get('default_seed_length', ''),
                                              data.get('default_three_length', ''))
            
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
            
            # Create a list to sort alphabetically
            sorted_organisms = []
            
            # Process each organism that has this endonuclease
            for organism, endos in organisms_to_endos.items():
                if endo in endos:
                    cspr_file = os.path.join(self.settings.CSPR_DB, organisms_to_files[organism][endo][0])
                    db_file = os.path.join(self.settings.CSPR_DB, organisms_to_files[organism][endo][1])
                    
                    if not os.path.exists(db_file):
                        self.logger.warning(f"Database file not found: {db_file}")
                        continue
                    
                    sorted_organisms.append((organism, cspr_file, db_file))
            
            # Sort organisms alphabetically by organism name
            sorted_organisms.sort(key=lambda x: x[0].lower())
            
            # Store the sorted results
            for index, (organism, cspr_file, db_file) in enumerate(sorted_organisms):
                org_files.append((organism, cspr_file, db_file))
                
                # Store the mapping for later use
                self.index_to_cspr[index] = cspr_file
                self.index_to_db[index] = db_file
            
            self.logger.info(f"Found {len(org_files)} organism files")
        except Exception as e:
            self.logger.error(f"Error getting organism files: {str(e)}")
            show_error(self.settings, "Error getting organism files", str(e))
        return org_files

    def get_shared_seeds(self, db_files, limit=False):
        """Get shared seeds between organisms"""
        try:
            self.logger.debug(f"Getting shared seeds for {len(db_files)} organisms")
            
            # Create temporary database for join operations
            temp_db_path = os.path.join(self.app_dir, "temp_join.db")
            new_conn = sqlite3.connect(temp_db_path)
            new_c = new_conn.cursor()
            
            # Set pragmas for better performance
            new_c.execute("PRAGMA synchronous = OFF")
            new_c.execute("PRAGMA journal_mode = OFF")
            new_c.execute("PRAGMA locking_mode = EXCLUSIVE")
            
            # Clean up any existing tables
            new_c.execute("DROP TABLE IF EXISTS repeats")
            new_c.execute("DROP TABLE IF EXISTS join_results")
            new_c.execute("VACUUM")
            
            # Create results table
            new_c.execute("CREATE TABLE join_results (seed TEXT)")
            
            try:
                # Attach all databases
                for i, db_file in enumerate(db_files, 1):
                    new_c.execute(f"ATTACH DATABASE '{db_file}' AS main{i}")
                
                # Start transaction
                new_c.execute("BEGIN TRANSACTION")
                
                # Build query to find seeds shared across all organisms
                base_sql = """INSERT into main.join_results 
                            select main1.repeats.seed from main1.repeats"""
                
                joins = []
                for i in range(2, len(db_files) + 1):
                    joins.append(
                        f"inner join main{i}.repeats on "
                        f"main{i-1}.repeats.seed = main{i}.repeats.seed"
                    )
                
                full_sql = base_sql + " " + " ".join(joins)
                self.logger.debug(f"Executing SQL for shared seeds: {full_sql}")
                new_c.execute(full_sql)
                
                # Get results based on limit parameter
                if limit:
                    shared_seeds = new_c.execute("SELECT DISTINCT seed FROM join_results LIMIT 1000").fetchall()
                    result = [seed[0] for seed in shared_seeds]
                else:
                    result = new_c.execute("SELECT COUNT(DISTINCT seed) FROM join_results").fetchone()[0]
                
                # Commit and cleanup
                new_c.execute("END TRANSACTION")
                new_c.close()
                new_conn.close()
                
                try:
                    os.remove(temp_db_path)
                except:
                    self.logger.warning("Could not remove temporary database file")
                
                self.logger.debug(f"Found {len(result) if limit else result} shared seeds")
                return result
                
            except Exception as e:
                new_c.execute("ROLLBACK")
                raise e
                
        except Exception as e:
            self.logger.error(f"Error getting shared seeds: {str(e)}")
            self.logger.exception("Full traceback:")
            show_error(self.settings, "Error getting shared seeds", str(e))
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
        """Get data for heatmap visualization"""
        try:
            size = len(db_files)
            arr = [[0 for _ in range(size)] for _ in range(size)]

            # Get shared seeds between pairs of organisms
            for i, j in itertools.combinations(range(size), 2):
                # Create temporary database for join operations
                temp_db_path = os.path.join(self.app_dir, "temp_join.db")
                new_conn = sqlite3.connect(temp_db_path)
                new_c = new_conn.cursor()
                
                try:
                    # Attach databases
                    new_c.execute(f"ATTACH DATABASE '{db_files[i]}' AS main1")
                    new_c.execute(f"ATTACH DATABASE '{db_files[j]}' AS main2")
                    
                    # Count shared seeds
                    sql = """SELECT COUNT(DISTINCT main1.repeats.seed) 
                            FROM main1.repeats 
                            INNER JOIN main2.repeats 
                            ON main1.repeats.seed = main2.repeats.seed"""
                            
                    shared_count = new_c.execute(sql).fetchone()[0]
                    arr[i][j] = arr[j][i] = shared_count
                    
                finally:
                    new_c.close()
                    new_conn.close()
                    try:
                        os.remove(temp_db_path)
                    except:
                        pass

            # Get individual organism seed counts
            for i in range(size):
                with sqlite3.connect(db_files[i]) as conn:
                    c = conn.cursor()
                    arr[i][i] = c.execute("SELECT COUNT(*) FROM repeats").fetchone()[0]

            return arr
            
        except Exception as e:
            self.logger.error(f"Error generating heatmap data: {str(e)}")
            self.logger.exception("Full traceback:")
            show_error(self.settings, "Error generating heatmap data", str(e))
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

