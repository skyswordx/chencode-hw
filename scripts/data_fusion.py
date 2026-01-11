import os
import csv
import re
import glob
from collections import defaultdict

# ============================================================================
# Configuration
# ============================================================================
# Input: Folder containing the simulation records to fuse
INPUT_DIR = r'd:\Musii-SnapShot\GithubRepo\chencode-hw\output\10_error_test'
# Output: Folder where fused results will be saved
OUTPUT_DIR = r'd:\Musii-SnapShot\GithubRepo\chencode-hw\output\10_error_test\fused_results'

# File pattern to match: captures the base name (e.g., "ber_Turbo_LogMAP_K1784_R1-3")
# Expects format: <base_name>_YYYYMMDD_HHMMSS[_merged].csv
FILE_PATTERN = r'^(.*)_\d{8}_\d{6}(?:_merged)?\.csv$'

# Integer limits for risk analysis
INT32_MAX = 2**31 - 1
UINT32_MAX = 2**32 - 1
INT64_MAX = 2**63 - 1
DOUBLE_PRECISION_LIMIT = 2**53  # Loss of precision for integers in double

# ============================================================================
# Utility Functions
# ============================================================================

def ensure_dir(path):
    """Create directory if it doesn't exist."""
    if not os.path.exists(path):
        os.makedirs(path)

def parse_csv_file(filepath):
    """
    Reads a CSV file and returns a list of dictionaries (rows).
    Ensures numerical columns are converted to floats/ints.
    Handles comment lines (starting with #) and empty lines.
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

def format_snr(snr):
    """Format SNR for consistent display (handles -0.0 -> 0.0)."""
    if snr == 0.0:
        return "0.0"
    return f"{snr:.1f}"

# ============================================================================
# Main Fusion Logic
# ============================================================================

def main():
    ensure_dir(OUTPUT_DIR)
    
    # 1. Group files by their base pattern (e.g., by K value and Rate)
    all_files = glob.glob(os.path.join(INPUT_DIR, '*.csv'))
    groups = defaultdict(list)
    
    print("=" * 70)
    print("Data Fusion Script for Different SNR Range Files")
    print("=" * 70)
    print(f"Input Directory:  {INPUT_DIR}")
    print(f"Output Directory: {OUTPUT_DIR}")
    print()
    
    for filepath in all_files:
        filename = os.path.basename(filepath)
        
        # Skip 'part_' prefix files (intermediate parallel simulation files)
        if filename.startswith('part_'):
            continue
        
        # Skip already fused files
        if '_FUSED' in filename:
            continue
            
        match = re.match(FILE_PATTERN, filename)
        if match:
            group_key = match.group(1)
            groups[group_key].append(filepath)
            
    print(f"Found {len(groups)} distinct simulation groups to fuse.\n")
    
    # 2. Process each group
    for group_name, file_list in sorted(groups.items()):
        print("-" * 70)
        print(f"Group: {group_name}")
        print(f"Files: {len(file_list)} to merge")
        
        # ====================================================================
        # KEY: Aggregation using SNR as key
        # This correctly handles files with different SNR ranges!
        # - File A has SNR [0.0, 0.1, 0.2]
        # - File B has SNR [0.1, 0.2, 0.3]
        # Result: SNR [0.0, 0.1, 0.2, 0.3] with correct cumulative stats
        # ====================================================================
        aggregated_stats = defaultdict(lambda: {
            'Bit_Errors': 0.0, 
            'Total_Bits': 0.0, 
            'Frame_Errors': 0.0, 
            'Total_Frames': 0.0,
            'Source_Files': 0  # Track how many files contributed to this SNR point
        })
        
        # Track info about the source files
        has_data = False
        all_snr_ranges = []
        
        for filepath in file_list:
            filename = os.path.basename(filepath)
            rows = parse_csv_file(filepath)
            
            if not rows:
                print(f"  ! Warning: No data in {filename}")
                continue
            
            # Track SNR range of this file
            file_snrs = [float(row['Eb_N0']) for row in rows if 'Eb_N0' in row]
            if file_snrs:
                snr_min, snr_max = min(file_snrs), max(file_snrs)
                all_snr_ranges.append((snr_min, snr_max, filename))
            
            for row in rows:
                if 'Eb_N0' not in row:
                    continue
                
                try:
                    snr = float(row['Eb_N0'])
                except (ValueError, TypeError):
                    continue
                    
                stats = aggregated_stats[snr]
                
                # Accumulate raw counts (not rates!)
                stats['Bit_Errors'] += row.get('Bit_Errors', 0)
                stats['Total_Bits'] += row.get('Total_Bits', 0)
                stats['Frame_Errors'] += row.get('Frame_Errors', 0)
                stats['Total_Frames'] += row.get('Total_Frames', 0)
                stats['Source_Files'] += 1
                has_data = True

        if not has_data:
            print("  -> No valid data found in this group. Skipping.")
            continue
        
        # 3. Show SNR range coverage info
        print("\n  SNR Coverage Analysis:")
        if all_snr_ranges:
            # Overall coverage
            overall_min = min(r[0] for r in all_snr_ranges)
            overall_max = max(r[0] for r in all_snr_ranges)
            print(f"    Overall SNR range: [{format_snr(overall_min)}, {format_snr(overall_max)}] dB")
            
            # Group files by their SNR range
            range_groups = defaultdict(list)
            for snr_min, snr_max, fname in all_snr_ranges:
                range_key = (snr_min, snr_max)
                range_groups[range_key].append(fname)
            
            print(f"    Different SNR ranges found: {len(range_groups)}")
            for (rmin, rmax), fnames in sorted(range_groups.items()):
                print(f"      [{format_snr(rmin)}, {format_snr(rmax)}]: {len(fnames)} files")
            
        # 4. Recalculate Metrics and Check Risk
        final_rows = []
        max_bits = 0
        sorted_snrs = sorted(aggregated_stats.keys())
        
        print(f"\n  Fused SNR points: {len(sorted_snrs)}")
        print(f"  SNR range: [{format_snr(min(sorted_snrs))}, {format_snr(max(sorted_snrs))}] dB")
        
        # Statistics table header
        print("\n  Per-SNR Statistics:")
        print("  " + "-" * 66)
        print(f"  {'SNR':>6} | {'Sources':>7} | {'Frames':>12} | {'BER':>12} | {'FER':>12}")
        print("  " + "-" * 66)
        
        for snr in sorted_snrs:
            stats = aggregated_stats[snr]
            
            bits = stats['Total_Bits']
            frames = stats['Total_Frames']
            bit_err = stats['Bit_Errors']
            frame_err = stats['Frame_Errors']
            source_files = stats['Source_Files']
            
            if bits > max_bits:
                max_bits = bits
            
            # Recalculate rates from fused counts
            ber = bit_err / bits if bits > 0 else 0.0
            fer = frame_err / frames if frames > 0 else 0.0
            
            # Print per-SNR summary
            ber_str = f"{ber:.4e}" if ber > 0 else "0"
            fer_str = f"{fer:.4e}" if fer > 0 else "0"
            print(f"  {format_snr(snr):>6} | {source_files:>7} | {int(frames):>12,} | {ber_str:>12} | {fer_str:>12}")
            
            final_rows.append({
                'Eb_N0': snr,
                'Bit_Errors': int(bit_err),
                'Total_Bits': int(bits),
                'BER': ber,
                'Frame_Errors': int(frame_err),
                'Total_Frames': int(frames),
                'FER': fer
            })
        
        print("  " + "-" * 66)
            
        # 5. Overflow Risk Analysis
        risk_level = "None"
        risk_msg = ""
        
        if max_bits > DOUBLE_PRECISION_LIMIT:
            risk_level = "CRITICAL"
            risk_msg = f"Total bits ({max_bits:.2e}) exceeds double precision exactness (9e15). BER calculation might lose precision."
        elif max_bits > INT32_MAX:
             risk_level = "Notice"
             risk_msg = f"Total bits ({max_bits:.2e}) exceeds 32-bit integer limit (2e9)."
        
        if risk_msg:
            print(f"\n  [Overflow Check] {risk_level}: {risk_msg}")
        else:
            print(f"\n  [Overflow Check] OK (Max bits: {max_bits:.2e} < 2^31)")

        # 6. Write Output
        output_filename = f"{group_name}_FUSED.csv"
        output_path = os.path.join(OUTPUT_DIR, output_filename)
        
        try:
            with open(output_path, 'w', newline='') as f:
                fieldnames = ['Eb_N0', 'Bit_Errors', 'Total_Bits', 'BER', 'Frame_Errors', 'Total_Frames', 'FER']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                
                # Write header with metadata
                f.write(f"# Fused Result for {group_name}\n")
                f.write(f"# Source Files: {len(file_list)}\n")
                f.write(f"# SNR Points: {len(final_rows)}\n")
                f.write(f"# SNR Range: [{format_snr(min(sorted_snrs))}, {format_snr(max(sorted_snrs))}] dB\n")
                f.write(f"# Max Total Bits: {max_bits}\n")
                writer.writeheader()
                writer.writerows(final_rows)
            print(f"\n  -> Output: {output_path}")
        except Exception as e:
            print(f"\n  -> Error writing output: {e}")

    print("\n" + "=" * 70)
    print("Data Fusion Complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
