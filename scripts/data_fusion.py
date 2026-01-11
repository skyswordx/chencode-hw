import os
import csv
import re
import glob
from collections import defaultdict

# Configuration
INPUT_DIR = r'd:\Musii-SnapShot\GithubRepo\chencode-hw\overnight_output'
OUTPUT_DIR = r'd:\Musii-SnapShot\GithubRepo\chencode-hw\output\fused_results'
FILE_PATTERN = r'^(.*)_\d{8}_\d{6}(?:_merged)?\.csv$'

# Integer limits for risk analysis
INT32_MAX = 2**31 - 1
UINT32_MAX = 2**32 - 1
INT64_MAX = 2**63 - 1
DOUBLE_PRECISION_LIMIT = 2**53  # Loss of precision for integers in double

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def parse_csv_file(filepath):
    """
    Reads a CSV file and returns a list of dictionaries (rows).
    Ensures numerical columns are converted to floats/ints.
    """
    data = []
    try:
        with open(filepath, 'r', newline='', encoding='utf-8-sig', errors='replace') as f:
            # Filter comments and empty lines
            filtered_lines = (line for line in f if line.strip() and not line.strip().startswith('#'))
            reader = csv.DictReader(filtered_lines)
            # Normalize headers (strip whitespace)
            reader.fieldnames = [name.strip() for name in reader.fieldnames] if reader.fieldnames else []
            
            for row in reader:
                # Clean keys and values
                clean_row = {}
                for k, v in row.items():
                    if k is None: continue
                    clean_key = k.strip()
                    try:
                        clean_row[clean_key] = float(v)
                    except ValueError:
                        clean_row[clean_key] = v
                data.append(clean_row)
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return data

def main():
    ensure_dir(OUTPUT_DIR)
    
    # 1. Group files
    all_files = glob.glob(os.path.join(INPUT_DIR, '*.csv'))
    groups = defaultdict(list)
    
    print(f"Scanning {INPUT_DIR}...")
    
    for filepath in all_files:
        filename = os.path.basename(filepath)
        
        # Skip 'parts' if we are looking for 'merged' files to combine
        if filename.startswith('part_'):
            continue
            
        match = re.match(FILE_PATTERN, filename)
        if match:
            group_key = match.group(1)
            groups[group_key].append(filepath)
            
    print(f"Found {len(groups)} distinct simulation groups.")
    
    # 2. Process each group
    for group_name, file_list in groups.items():
        print(f"\nProcessing Group: {group_name}")
        print(f"  -> Found {len(file_list)} files to merge.")
        
        # Aggregation Dictionary: stats[Eb_N0] = {Bit_Errors: 0, Total_Bits: 0, ...}
        aggregated_stats = defaultdict(lambda: {
            'Bit_Errors': 0.0, 
            'Total_Bits': 0.0, 
            'Frame_Errors': 0.0, 
            'Total_Frames': 0.0
        })
        
        # Track if we found valid data
        has_data = False
        
        for filepath in file_list:
            rows = parse_csv_file(filepath)
            for row in rows:
                if 'Eb_N0' not in row:
                    continue
                
                try:
                    snr = float(row['Eb_N0'])
                except (ValueError, TypeError):
                    # print(f"Skipping invalid SNR: {row['Eb_N0']} in {filepath}")
                    continue
                    
                stats = aggregated_stats[snr]
                
                # safely add values
                stats['Bit_Errors'] += row.get('Bit_Errors', 0)
                stats['Total_Bits'] += row.get('Total_Bits', 0)
                stats['Frame_Errors'] += row.get('Frame_Errors', 0)
                stats['Total_Frames'] += row.get('Total_Frames', 0)
                has_data = True

        if not has_data:
            print("  -> No valid data found in this group.")
            continue
            
        # 3. Recalculate Metrics and Check Risk
        final_rows = []
        max_bits = 0
        sorted_snrs = sorted(aggregated_stats.keys())
        
        for snr in sorted_snrs:
            stats = aggregated_stats[snr]
            
            bits = stats['Total_Bits']
            frames = stats['Total_Frames']
            bit_err = stats['Bit_Errors']
            frame_err = stats['Frame_Errors']
            
            if bits > max_bits:
                max_bits = bits
            
            # Recalculate Rates
            ber = bit_err / bits if bits > 0 else 0.0
            fer = frame_err / frames if frames > 0 else 0.0
            
            final_rows.append({
                'Eb_N0': snr,
                'Bit_Errors': int(bit_err),
                'Total_Bits': int(bits),
                'BER': ber,
                'Frame_Errors': int(frame_err),
                'Total_Frames': int(frames),
                'FER': fer
            })
            
        # 4. Overflow Risk Analysis
        risk_level = "None"
        risk_msg = ""
        
        if max_bits > DOUBLE_PRECISION_LIMIT:
            risk_level = "CRITICAL"
            risk_msg = f"Total bits ({max_bits:.2e}) exceeds double precision exactness (9e15). BER calculation might lose precision."
        elif max_bits > INT32_MAX:
             risk_level = "Notice"
             risk_msg = f"Total bits ({max_bits:.2e}) exceeds 32-bit integer limit (2e9)."
        
        print(f"  -> Merged {len(file_list)} files into {len(final_rows)} SNR points.")
        if risk_msg:
            print(f"  -> [Overflow Check] {risk_level}: {risk_msg}")
        else:
            print(f"  -> [Overflow Check] OK (Max bits: {max_bits:.2e} < 2^31)")

        # 5. Write Output
        output_filename = f"{group_name}_FUSED.csv"
        output_path = os.path.join(OUTPUT_DIR, output_filename)
        
        try:
            with open(output_path, 'w', newline='') as f:
                fieldnames = ['Eb_N0', 'Bit_Errors', 'Total_Bits', 'BER', 'Frame_Errors', 'Total_Frames', 'FER']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                
                # Write header with metadata
                f.write(f"# Fused Result for {group_name}\n")
                f.write(f"# Source Files: {len(file_list)}\n")
                f.write(f"# Max Total Bits: {max_bits}\n")
                writer.writeheader()
                writer.writerows(final_rows)
            print(f"  -> Saved to {output_path}")
        except Exception as e:
            print(f"  -> Error writing output: {e}")

if __name__ == "__main__":
    main()
