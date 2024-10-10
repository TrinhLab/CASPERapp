import os
import sqlite3
import itertools
from utils.ui import show_error

class PopulationAnalysisWindowModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.app_dir = global_settings.get_app_dir()
        self.cspr_files = []
        self.db_files = []
        self.org_names = {}
        self.seeds = []
        self.counts = []
        self.index_to_cspr = {}
        self.index_to_db = {}

    def load_endonucleases(self):
        endos = {}
        try:
            with open(self.global_settings.get_casper_info_path(), 'r') as f:
                for line in f:
                    if line.startswith('ENDONUCLEASES'):
                        for line in f:
                            if line.startswith('-'):
                                break
                            line_tokens = line.strip().split(';')
                            endo = line_tokens[0]
                            pam = line_tokens[1].split(',')[0] if ',' in line_tokens[1] else line_tokens[1]
                            default_five_length = line_tokens[2]
                            default_seed_length = line_tokens[3]
                            default_three_length = line_tokens[4]
                            endos[f"{endo} PAM: {pam}"] = (endo, pam, default_five_length, default_seed_length, default_three_length)
        except Exception as e:
            show_error(self.global_settings, "Error loading endonucleases", str(e))
        return endos

    def get_organism_files(self, endo):
        org_files = []
        try:
            for file in os.listdir(self.global_settings.CSPR_DB):
                if file.endswith('.cspr') and file[file.rfind('_') + 1:file.find('.cspr')] == endo:
                    cspr_file = os.path.join(self.global_settings.CSPR_DB, file)
                    db_file = cspr_file.replace(".cspr", "_repeats.db")
                    with open(cspr_file, 'r') as f:
                        org_name = f.readline().split(":")[-1].strip()
                    org_files.append((org_name, cspr_file, db_file))
        except Exception as e:
            show_error(self.global_settings, "Error getting organism files", str(e))
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
                                    'organism': os.path.basename(db_file).split('_')[0],
                                    'chromosome': chrom,
                                    'location': abs(int(locs[i]))
                                })
        except Exception as e:
            show_error(self.global_settings, "Error getting seed locations", str(e))
        return locations

    def get_org_names(self):
        # Implement this method to populate self.org_names
        pass

