from utils.sequence_utils import SeqTranslate
import logging
from multiprocessing import Pool, cpu_count
from functools import partial
import pickle
import os
import traceback

class CSPRparser:
    def __init__(self, inputFileName, casper_info_path):
        self.fileName = inputFileName
        self.seqTrans = SeqTranslate(casper_info_path)
        self.logger = logging.getLogger(__name__)
        self._cached_results = {}
        self.index_file = f"{inputFileName}.index"

    def _create_index(self):
        try:
            self.logger.debug("Creating CSPR index file...")
            
            # Initialize index structure
            index_data = {}
            
            with open(self.fileName, 'rb') as f:
                # Skip header lines
                for _ in range(3):
                    f.readline()
                
                current_chrom = None
                chrom_data = []
                
                # Process file line by line
                for line in f:
                    if line.startswith(b'>'):
                        # Save previous chromosome data if exists
                        if current_chrom and chrom_data:
                            index_data[current_chrom] = chrom_data
                            
                        # Start new chromosome
                        current_chrom = line.decode().split()[0][1:]  # Remove '>' and get chromosome id
                        chrom_data = []
                        continue
                    
                    if not line.strip():
                        continue
                        
                    try:
                        # Parse position and store line offset
                        first_comma = line.find(b',')
                        if first_comma != -1:
                            pos = int(line[:first_comma])
                            abs_pos = abs(pos)
                            chrom_data.append((abs_pos, line))
                    except ValueError:
                        continue
                
                # Save last chromosome data
                if current_chrom and chrom_data:
                    index_data[current_chrom] = chrom_data

            # Save index to file
            with open(self.index_file, 'wb') as f:
                pickle.dump(index_data, f)

            self._index = index_data
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating index: {str(e)}")
            return False

    def _load_index(self):
        try:
            if not os.path.exists(self.index_file):
                return False
                
            if os.path.getmtime(self.index_file) < os.path.getmtime(self.fileName):
                return False
                
            with open(self.index_file, 'rb') as f:
                self._index = pickle.load(f)
            return True
            
        except Exception as e:
            self.logger.error(f"Error loading index: {str(e)}")
            return False

    def read_targets_batch(self, chromosome, targets, endonuclease):
        try:


            print(f"Reading targets for chromosome: {chromosome}")
            # Load or create index
            if not hasattr(self, '_index'):
                if not self._load_index():
                    self._create_index()
            
            # Sort targets by start position
            sorted_targets = sorted(targets, key=lambda x: x['start'])
            min_start = sorted_targets[0]['start']
            max_end = max(t['end'] for t in sorted_targets)
            
            self.logger.debug(f"Processing targets from {min_start} to {max_end}")
            self.logger.debug(f"Looking for chromosome: {chromosome}")
            
            results = []
            lines_processed = 0
            lines_skipped = 0
            
            # Find chromosome by full ID
            found_chrom = None
            for chrom_id in self._index:
                # Decode bytes to string if necessary
                chrom_str = chrom_id.decode() if isinstance(chrom_id, bytes) else chrom_id
                
                # Match the full chromosome ID
                if chrom_str == chromosome:
                    found_chrom = chrom_id
                    self.logger.debug(f"Found matching chromosome: {chrom_str}")
                    break
                    
            if found_chrom:
                chrom_data = self._index[found_chrom]
                
                # Binary search for start position
                start_idx = 0
                end_idx = len(chrom_data)
                
                for target in sorted_targets:
                    target_start = target['start']
                    target_end = target['end']
                    feature_id = target.get('feature_id', '')
                    feature_name = target.get('feature_name', '')
                    
                    # Find relevant positions for this target
                    while start_idx < end_idx and chrom_data[start_idx][0] < target_start:
                        start_idx += 1
                    
                    current_idx = start_idx
                    while current_idx < end_idx and chrom_data[current_idx][0] < target_end:
                        try:
                            pos, line = chrom_data[current_idx]
                            parts = line.split(b',')
                            
                            if len(parts) >= 4:
                                pos = int(parts[0])
                                results.append({
                                    'feature_name': feature_name,
                                    'feature_id': feature_id,
                                    'chromosome': found_chrom,
                                    'position': abs(pos),
                                    'location': f"{abs(pos)}-{abs(pos) + 23}",
                                    'sequence': parts[1].decode(),
                                    'pam': parts[2].decode(),
                                    'strand': "-" if pos < 0 else "+",
                                    'score': float(parts[3]),
                                    'endonuclease': endonuclease
                                })
                                lines_processed += 1
                            
                        except (ValueError, IndexError) as e:
                            self.logger.error(f"Error processing line: {str(e)}")
                            lines_skipped += 1
                            
                        current_idx += 1
            else:
                self.logger.error(f"Chromosome {chromosome} not found in index")
                self.logger.debug(f"Available chromosomes: {list(self._index.keys())}")
            
            return results
            
        except Exception as e:
            self.logger.error(f"Error in read_targets_batch: {str(e)}")
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return []

    def parse_targets(self, file_path, region):
        """Parse targets with parallel processing and caching"""
        cache_key = f"{file_path}:{region}"
        if cache_key in self._cached_results:
            return self._cached_results[cache_key]
            
        # Split the region into chunks for parallel processing
        chunks = self._split_region(region)
        
        with Pool() as pool:
            results = pool.map(partial(self._parse_chunk, file_path), chunks)
            
        # Combine results
        combined_targets = []
        for chunk_result in results:
            combined_targets.extend(chunk_result)
            
        self._cached_results[cache_key] = combined_targets
        return combined_targets
        
    def _split_region(self, region):
        """Split a region into chunks for parallel processing"""
        start, end = region
        chunk_size = (end - start) // cpu_count()
        chunks = []
        for i in range(cpu_count()):
            chunk_start = start + (i * chunk_size)
            chunk_end = chunk_start + chunk_size if i < cpu_count()-1 else end
            chunks.append((chunk_start, chunk_end))
        return chunks

