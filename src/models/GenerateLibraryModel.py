from models.CSPRparser import CSPRparser
from models.HomeWindowModel import HomeWindowModel
import os
import re
import traceback

class GenerateLibraryModel(HomeWindowModel):
    def __init__(self, global_settings):
        super().__init__(global_settings)
        self.logger = global_settings.logger
        self.parser = None
        self.targets_data = {}
        self._deleted_targets = {}
        
    def initialize_parser(self, cspr_file):
        """Initialize CSPR parser"""
        self.parser = CSPRparser(cspr_file, self.global_settings.get_casper_info_path())
        
    def generate_library(self, selected_targets, settings):
        """Generate library with given settings"""
        try:
            self.logger.debug(f"Generating library with settings: {settings}")
            
            # Process targets based on settings
            processed_targets = self._process_targets(
                selected_targets,
                settings['min_score'],
                settings['five_prime_seq'],
                settings['target_range_start'],
                settings['target_range_end']
            )

            # Generate output for each target
            output_data = self._generate_output(
                processed_targets,
                settings['guides_per_gene'],
                settings['space_between_guides']
            )

            self.logger.debug(f"Output data: {output_data}")
            
            # Write output to file
            self._write_output(output_data, settings)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error generating library: {str(e)}")
            self.logger.error(traceback.format_exc())
            raise
            
    def _process_targets(self, targets, min_score, five_prime_seq, start_range, end_range):
        """Process and filter targets based on criteria"""
        processed = {}
        self._deleted_targets = {}  # Store deleted targets for modify parameters option
        
        self.logger.debug(f"Processing {len(targets)} targets with min_score={min_score}")
        
        for target in targets:
            # Use feature_name (gene name) as key instead of feature_id
            gene_name = target.get('feature_name', target.get('feature_id'))
            
            if gene_name not in processed:
                processed[gene_name] = []
                self._deleted_targets[gene_name] = []
                
            # Create a copy of the target data
            target_data = target.copy()
            
            # Check if target passes filters
            if self._passes_filters(target_data, min_score, five_prime_seq, start_range, end_range):
                processed[gene_name].append(target_data)
            else:
                self._deleted_targets[gene_name].append(target_data)
                
        # Sort targets for each gene
        for gene in processed:
            # First sort by score (ascending)
            processed[gene].sort(key=lambda x: float(x['score']))
            
            # Then sort by position (ascending)
            processed[gene].sort(key=lambda x: abs(int(x['position'])))
            
            # Reverse list if gene is on negative strand
            if processed[gene] and processed[gene][0].get('strand', '+') == '-':
                processed[gene].reverse()
                
        return processed
        
    def _passes_filters(self, target, min_score, five_prime_seq, start_range, end_range):
        """Check if target passes all filters"""
        try:
            # Score filter - convert score to float and compare
            target_score = float(target.get('score', 0))
            if target_score < min_score:
                self.logger.debug(f"Target failed score filter: {target_score} < {min_score}")
                return False
            
            # Poly-T filter
            if re.search("T{5,10}", target['sequence']):
                self.logger.debug(f"Target failed poly-T filter: {target['sequence']}")
                return False
            
            # 5' sequence filter
            if five_prime_seq and not target['sequence'].startswith(five_prime_seq.upper()):
                self.logger.debug(f"Target failed 5' sequence filter")
                return False
            
            # Range filter
            if start_range != 0 or end_range != 100:
                position_ratio = self._calculate_position_ratio(target)
                if not (start_range/100 <= position_ratio <= end_range/100):
                    self.logger.debug(f"Target failed range filter: {position_ratio}")
                    return False
                
            return True
            
        except Exception as e:
            self.logger.error(f"Error in _passes_filters: {str(e)}")
            self.logger.error(f"Target data: {target}")
            return False
        
    def _calculate_position_ratio(self, target):
        """Calculate relative position ratio in gene"""
        try:
            # Get gene length from feature info
            gene_length = abs(int(target.get('end', 0)) - int(target.get('start', 0)))
            if gene_length == 0:
                return 0
            
            # Calculate position relative to gene start
            target_pos = abs(int(target['position']))
            gene_start = int(target.get('start', 0))
            
            # If gene is on reverse strand, flip the position calculation
            if target.get('strand', '+') == '-':
                ratio = (target.get('end', 0) - target_pos) / gene_length
            else:
                ratio = (target_pos - gene_start) / gene_length
            
            self.logger.debug(f"Position ratio: {ratio} (pos={target_pos}, start={gene_start}, len={gene_length}, strand={target.get('strand', '+')}")
            return ratio
            
        except Exception as e:
            self.logger.error(f"Error calculating position ratio: {str(e)}")
            return 0
        
    def _generate_output(self, processed_targets, guides_per_gene, space_between):
        output = {}
        
        for gene_id, targets in processed_targets.items():
            output[gene_id] = []
            i = 0
            vec_index = 0
            prev_target = None
            
            while i < guides_per_gene:
                if len(targets) == 0 or vec_index >= len(targets):
                    break
                    
                current = targets[vec_index]
                
                # Check spacing from previous target
                if prev_target is None or abs(int(current['position']) - int(prev_target['position'])) >= space_between:
                    # If current target has better score than previous
                    if (prev_target and float(current['score']) > float(prev_target['score'])):
                        output[gene_id].pop()
                        output[gene_id].append(current)
                    else:
                        output[gene_id].append(current)
                    prev_target = current
                    i += 1
                
                vec_index += 1
                if vec_index >= len(targets):
                    break
                    
            # Add deleted targets if needed
            if len(output[gene_id]) < guides_per_gene:
                deleted_sorted = sorted(
                    self._deleted_targets.get(gene_id, []),
                    key=lambda x: (float(x['score']), abs(int(x['position'])))
                )
                
                for deleted_target in deleted_sorted:
                    if len(output[gene_id]) >= guides_per_gene:
                        break
                    deleted_target['modified'] = True
                    output[gene_id].append(deleted_target)
        
        return output
        
    def _write_output(self, output_data, settings):
        """Write library to output file"""
        output_file = settings['output_file']
        if not output_file.endswith('.csv'):
            output_file += '.csv'
            
        with open(output_file, 'w') as f:
            # Write header
            headers = ['Gene Name', 'Sequence', 'On-Target Score']
            if settings.get('find_off_targets'):
                headers.append('Off-Target Score')
            headers.extend(['Location', 'PAM', 'Strand'])
            f.write(','.join(headers) + '\n')
            
            # Write data
            for gene_id, targets in output_data.items():
                if not targets:  # Skip genes with no targets
                    continue
                    
                for i, target in enumerate(targets, 1):
                    # Use feature_name (gene name) instead of feature_id
                    gene_name = target.get('feature_name', gene_id)
                    
                    # Format gene name with index
                    tag_id = f"{gene_name}-{i}"
                    tag_id = tag_id.replace(',', '')  # Remove any commas
                    
                    # If target was modified by parameters, add asterisks
                    if settings.get('modify_params') and target.get('modified'):
                        tag_id = "**" + tag_id
                    
                    row = [
                        tag_id,
                        target['sequence'],
                        str(target['score'])
                    ]
                    
                    if settings.get('find_off_targets'):
                        row.append(str(target.get('off_target_score', '')))
                        
                    row.extend([
                        str(abs(int(target['position']))),
                        target['pam'],
                        target['strand'][0]  # Just take first character of strand
                    ])
                    
                    f.write(','.join(row) + '\n')
