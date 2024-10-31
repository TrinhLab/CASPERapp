from utils.sequence_utils import SeqTranslate
import logging
from multiprocessing import Pool, cpu_count
from functools import partial

class CSPRparser:
    def __init__(self, inputFileName, casper_info_path):
        self.fileName = inputFileName
        self.filename = inputFileName
        self.seqTrans = SeqTranslate(casper_info_path)
        self.logger = logging.getLogger(__name__)
        self._line_buffer = []  # Pre-allocate buffer for lines
        self._cached_results = {}

    def read_targets_batch(self, chromosome, targets, endonuclease):
        """Ultra-fast target reading using direct tuple creation"""
        try:
            # Pre-process targets into a sorted list of ranges for faster lookup
            target_ranges = []
            for t in targets:
                start = t['start']
                end = t['end']
                target_ranges.append((start, end, t['feature_name']))
            target_ranges.sort()  # Sort by start position
            
            # Pre-allocate results list
            results = []
            results_append = results.append
            
            # Read file in binary mode for speed
            with open(self.fileName, 'rb') as f:
                # Skip header
                for _ in range(3):
                    f.readline()
                
                # Find chromosome section
                header = False
                for line in f:
                    if b'>' in line and str(chromosome).encode() in line:
                        header = True
                        break
                
                # Read targets
                if header:
                    current_range_idx = 0
                    max_ranges = len(target_ranges)
                    
                    while current_range_idx < max_ranges:
                        line = f.readline()
                        if not line or line.startswith(b'>'):
                            break
                            
                        if not line.strip():
                            continue
                            
                        # Fast string splitting without decode
                        parts = line.strip().split(b',')
                        if not parts:
                            continue
                            
                        try:
                            pos = int(parts[0])
                            abs_pos = abs(pos)
                            
                            # Get current target range
                            start, end, feature_name = target_ranges[current_range_idx]
                            
                            # Skip if position is past current range
                            if abs_pos >= end:
                                current_range_idx += 1
                                continue
                                
                            # Check if position is in range
                            if start <= abs_pos < end:
                                sequence = parts[1].decode()
                                pam = sequence[-3:]
                                target_seq = sequence[:-3]
                                
                                results_append({
                                    'feature_name': feature_name,
                                    'chromosome': chromosome,
                                    'position': abs_pos,
                                    'location': f"{abs_pos}-{abs_pos + 23}",
                                    'sequence': target_seq,
                                    'pam': pam,
                                    'strand': "-" if pos < 0 else "+",
                                    'score': float(parts[3]) if len(parts) > 3 else 0.0,
                                    'endonuclease': endonuclease
                                })
                                
                        except (ValueError, IndexError):
                            continue
            
            return results
                
        except Exception as e:
            self.logger.error(f"Error in read_targets_batch: {str(e)}")
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
