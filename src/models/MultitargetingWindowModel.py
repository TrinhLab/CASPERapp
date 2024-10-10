import os
import sqlite3
from collections import Counter
import statistics
import models.GlobalSettings as GlobalSettings
from models.CSPRparser import CSPRparser

class MultitargetingWindowModel:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.cspr_file = ""
        self.db_file = ""
        self.organisms_to_files = {}
        self.organisms_to_endos = {}
        self.chromo_length = []
        self.max_repeats = 1
        self.average = 0
        self.median = 0
        self.mode = 0
        self.average_unique = 0
        self.average_rep = 0
        self.repeat_count = 0
        self.row_limit = 1000
        self.parser = CSPRparser("", self.global_settings.get_casper_info_path())
        self.load_organisms_and_endos()

    def load_organisms_and_endos(self):
        # This method should populate the organisms_to_endos dictionary
        # You'll need to implement the logic to load this data from your database or files
        # For example:
        # self.organisms_to_endos = {
        #     "E. coli": ["Cas9", "Cas12a"],
        #     "S. cerevisiae": ["Cas9"],
        #     # ... other organisms and their associated endonucleases
        # }
        pass

    def get_organisms(self):
        # Return the list of organisms
        return list(self.organisms_to_endos.keys())

    def get_endos_for_organism(self, organism):
        # Return the list of endonucleases for a given organism
        return self.organisms_to_endos.get(organism, [])

    def set_files(self, organism, endo):
        self.cspr_file = self.organisms_to_files[organism][endo][0]
        self.db_file = self.organisms_to_files[organism][endo][1]

    def get_kstats(self):
        kstats = []
        with open(self.cspr_file, "r") as f:
            for line in f:
                if "KARYSTATS" in line:
                    kstats = line.replace("KARYSTATS: ", "").strip().split(',')[:-1]
                    break
        return kstats

    def get_repeats_data(self):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        
        if self.row_limit == -1:
            sql_query = "SELECT * FROM repeats ORDER BY count DESC;"
        else:
            sql_query = f"SELECT * FROM repeats ORDER BY count DESC LIMIT 0, {self.row_limit};"
        
        data = c.execute(sql_query).fetchall()
        c.close()
        conn.close()
        return data

    def get_seed_data(self, seed):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT chromosome, location, pam, score, five, three FROM repeats WHERE seed = ?", (seed,)).fetchone()
        c.close()
        conn.close()
        return data

    def process_seed_data(self, seed_data, kstats):
        chromo, pos, pam, score, five, three = seed_data
        chromo = chromo.split(',')
        pos = pos.split(',')
        pam = pam.split(',')
        score = score.split(',')
        five = five.split(',')
        three = three.split(',')

        seed_data = {}
        event_data = {}
        for i in range(len(chromo)):
            curr_chromo = int(chromo[i])
            dir = "+" if int(pos[i]) >= 0 else "-"
            normalized_location = abs(float(pos[i]) / float(kstats[curr_chromo - 1]))
            if curr_chromo in seed_data:
                seed_data[curr_chromo].append(normalized_location)
                event_data[curr_chromo].append([normalized_location, pos[i], five[i] + seed + three[i], pam[i], score[i], dir])
            else:
                seed_data[curr_chromo] = [normalized_location]
                event_data[curr_chromo] = [[normalized_location, pos[i], five[i] + seed + three[i], pam[i], score[i], dir]]

        return seed_data, event_data

    def get_chro_bar_data(self, seed):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT chromosome FROM repeats WHERE seed = ?", (seed,)).fetchone()
        c.close()
        conn.close()
        
        data = [int(x) for x in data[0].split(',')]
        bar_data = Counter(data)
        return bar_data

    def get_seeds_vs_repeats_data(self):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("select count, COUNT(count) as cnt from repeats group by count order by cnt DESC;").fetchall()
        c.close()
        conn.close()
        return data

    def get_repeats_vs_seeds_data(self):
        conn = sqlite3.connect(self.db_file)
        c = conn.cursor()
        data = c.execute("SELECT count from repeats;").fetchall()
        c.close()
        conn.close()
        
        y1 = [row[0] for row in data]
        self.average = statistics.mean(y1)
        self.mode = statistics.mode(y1)
        self.median = statistics.median(y1)
        self.repeat_count = len(y1)
        
        return y1

    def set_row_limit(self, limit):
        self.row_limit = limit
