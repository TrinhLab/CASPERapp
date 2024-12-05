import os
import sqlite3
import statistics
from functools import lru_cache
import time

class MultitargetingWindowModel:
    def __init__(self, global_settings):
        start_time = time.time()
        self.settings = global_settings
        self.logger = global_settings.get_logger()
        
        self.cspr_file = ""
        self.db_file = ""
        self.row_limit = 1000
        
        # Get organism and endo mappings from DatabaseManager
        db_start = time.time()
        self.organisms_to_files, self.organisms_to_endos = self.settings.db_manager.get_organisms_and_endos()
        self.logger.debug(f"Getting DB mappings took: {time.time() - db_start:.2f} seconds")
        
        self._cache = {}
        self.logger.debug(f"Model initialization took: {time.time() - start_time:.2f} seconds")

    @lru_cache(maxsize=32)
    def get_organisms(self):
        """Get list of available organisms with caching"""
        return list(self.organisms_to_endos.keys())

    @lru_cache(maxsize=32)
    def get_endos_for_organism(self, organism):
        """Get available endonucleases for given organism with caching"""
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
        start_time = time.time()
        if not self.db_file:
            raise ValueError("Database file not set. Please select an organism and endonuclease first.")
        
        try:
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            
            query_start = time.time()
            # Use row limit in query
            if self.row_limit == -1:  # No limit
                query = "SELECT * FROM repeats ORDER BY count DESC;"
            else:
                query = f"SELECT * FROM repeats ORDER BY count DESC LIMIT 0, {self.row_limit};"
            
            results = []
            for repeat in c.execute(query):
                # Extract repeat info
                seed = repeat[0]
                chroms = repeat[1].split(",")
                locs = repeat[2].split(",")
                threes = repeat[3].split(",")
                fives = repeat[4].split(",")
                pams = repeat[5].split(",")
                scores = repeat[6].split(",")
                count = repeat[7]

                # Calculate average repeats per scaffold
                location_repeat_counts = {}
                for chrom in chroms:
                    location_repeat_counts[chrom] = location_repeat_counts.get(chrom, 0) + 1
                avg_per_scaffold = sum(location_repeat_counts.values()) / len(location_repeat_counts)
                avg_per_scaffold = float("%.2f" % avg_per_scaffold)

                # Find majority sequence
                majority_index = 0
                sequences = ""
                consensus_percent = 0
                
                if threes[0] == '':
                    majority = max(set(fives), key=fives.count)
                    majority_index = fives.index(majority)
                    sequences = fives[majority_index] + seed
                    consensus_percent = (fives.count(fives[majority_index]) / len(fives)) * 100
                elif fives[0] == '':
                    majority = max(set(threes), key=threes.count)
                    majority_index = threes.index(majority)
                    sequences = seed + threes[majority_index]
                    consensus_percent = (threes.count(threes[majority_index]) / len(threes)) * 100
                else:
                    # Both 3' and 5' present
                    combined_seqs = [f"{fives[i]}{threes[i]}" for i in range(len(threes))]
                    majority = max(set(combined_seqs), key=combined_seqs.count)
                    majority_index = combined_seqs.index(majority)
                    sequences = fives[majority_index] + seed + threes[majority_index]
                    consensus_percent = (combined_seqs.count(majority) / len(combined_seqs)) * 100

                # Determine strand
                strand = "+" if int(locs[majority_index]) >= 0 else "-"
                
                # Format consensus percent
                consensus_percent = float("%.1f" % consensus_percent)
                
                results.append((
                    seed,               # Seed sequence
                    count,             # Total repeats
                    avg_per_scaffold,  # Average repeats per scaffold
                    sequences,         # Consensus sequence
                    consensus_percent, # Consensus percentage
                    scores[majority_index],  # Score
                    pams[majority_index],    # PAM
                    strand             # Strand
                ))
            
            self.logger.debug(f"Query and processing took: {time.time() - query_start:.2f} seconds")
            
            c.close()
            conn.close()
            
            self.logger.debug(f"Total get_repeats_data took: {time.time() - start_time:.2f} seconds")
            return results
            
        except Exception as e:
            self.logger.error(f"Error getting repeats data: {str(e)}")
            raise

    def get_seed_data(self, seed):
        """Get detailed data for a specific seed"""
        return self.settings.db_manager.get_seed_data(self.db_file, seed)

    def get_chro_bar_data(self, seed):
        """Get chromosome distribution data for a seed"""
        return self.settings.db_manager.get_chro_bar_data(self.db_file, seed)

    def get_seeds_vs_repeats_data(self):
        """Get data for seeds vs repeats plot"""
        try:
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            
            # Query to get count of sequences for each repeat count, ordered by count DESC
            # This matches the original implementation
            query = """
                SELECT count, COUNT(count) as cnt 
                FROM repeats 
                GROUP BY count 
                ORDER BY cnt DESC;
            """
            
            x_vals = []  # Number of repeats
            y_vals = []  # Number of sequences
            
            for row in c.execute(query):
                x_vals.append(row[0])  # count
                y_vals.append(row[1])  # cnt
                
            # Sort x_vals after collecting all data, just like in original
            x_vals = sorted(x_vals)
                
            c.close()
            conn.close()
            
            return {
                'x_vals': x_vals,
                'y_vals': y_vals
            }
            
        except Exception as e:
            self.logger.error(f"Error getting seeds vs repeats data: {str(e)}")
            return None

    def get_repeats_vs_seeds_data(self):
        """Get data for repeats vs seeds plot"""
        return self.settings.db_manager.get_repeats_vs_seeds_data(self.db_file)

    def calculate_statistics(self):
        """Calculate statistics for the repeats data"""
        try:
            conn = sqlite3.connect(self.db_file)
            c = conn.cursor()
            
            # Get all repeat counts
            counts = []
            for obj in c.execute("SELECT count FROM repeats;"):
                counts.append(obj[0])
                
            if not counts:
                return None
                
            # Calculate statistics
            stats = {
                'average': statistics.mean(counts),
                'mode': statistics.mode(counts),
                'median': statistics.median(counts),
                'repeat_count': len(counts)
            }
            
            c.close()
            conn.close()
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Error calculating statistics: {str(e)}")
            raise

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

    def set_row_limit(self, limit):
        """Set the maximum number of rows to return"""
        try:
            self.row_limit = int(limit)
        except ValueError:
            self.logger.error(f"Invalid row limit value: {limit}")
            self.row_limit = 1000  # Reset to default if invalid

    def get_row_limit(self):
        """Get current row limit setting"""
        return self.row_limit

    def _clear_cache(self):
        """Clear the internal cache"""
        self._cache.clear()
        self.get_organisms.cache_clear()
        self.get_endos_for_organism.cache_clear()