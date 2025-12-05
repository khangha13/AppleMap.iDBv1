#!/usr/bin/env python3
"""
Beagle-exact VCF validator
Mimics Beagle's VcfRecGTParser.ninthTabPos() logic exactly
"""
import gzip
import sys

def validate_vcf_for_beagle(vcf_path):
    """
    Validates VCF exactly as Beagle does in VcfRecGTParser.ninthTabPos()
    
    Beagle's logic:
    1. Reads header to get sample count
    2. For EACH line (including header), counts tabs to find 9th tab position
    3. Fails if any line doesn't have at least 9 tabs
    """
    
    errors = []
    
    try:
        with gzip.open(vcf_path, 'rt', encoding='utf-8', errors='replace') as fh:
            header_line = None
            line_num = 0
            
            # Find #CHROM header
            for raw_line in fh:
                line_num += 1
                if raw_line.startswith('#CHROM'):
                    header_line = raw_line.rstrip('\r\n')
                    break
            
            if not header_line:
                print(f"[FATAL] No #CHROM header found", file=sys.stderr)
                return False
            
            # Count tabs in header using Beagle's method
            tab_count = 0
            ninth_tab_pos = -1
            for i, char in enumerate(header_line):
                if char == '\t':
                    tab_count += 1
                    if tab_count == 9:
                        ninth_tab_pos = i
                        break
            
            print(f"Header line {line_num}: found {tab_count} tabs before 9th position")
            if ninth_tab_pos == -1:
                errors.append(f"Header line {line_num}: ninthTabPos failed - only {tab_count} tabs found (need >=9)")
                print(f"[FAILED] Header doesn't have 9 tabs", file=sys.stderr)
                return False
            
            # Get expected sample count
            header_fields = header_line.split('\t')
            sample_count = len(header_fields) - 9
            expected_tabs = len(header_fields) - 1
            
            print(f"Header has {len(header_fields)} fields, {sample_count} samples")
            print(f"Each data line must have exactly {expected_tabs} tabs")
            
            # Check each data line
            data_lines_checked = 0
            for raw_line in fh:
                line_num += 1
                if raw_line.startswith('#'):
                    continue
                
                line = raw_line.rstrip('\r\n')
                
                # Count tabs exactly as Beagle does
                tab_count = 0
                ninth_tab_pos = -1
                for i, char in enumerate(line):
                    if char == '\t':
                        tab_count += 1
                        if tab_count == 9:
                            ninth_tab_pos = i
                            break
                
                # Check if we found 9th tab
                if ninth_tab_pos == -1:
                    errors.append(f"Data line {line_num}: ninthTabPos failed - only {tab_count} tabs found")
                    if data_lines_checked < 10:  # Only show first 10 errors
                        print(f"[ERROR] Line {line_num}: only {tab_count} tabs (need >=9)", file=sys.stderr)
                        # Show the line (truncated)
                        print(f"  Content: {line[:200]}", file=sys.stderr)
                
                # Also check total tab count
                total_tabs = line.count('\t')
                if total_tabs != expected_tabs:
                    errors.append(f"Data line {line_num}: has {total_tabs} tabs (expected {expected_tabs})")
                    if data_lines_checked < 10:
                        print(f"[ERROR] Line {line_num}: {total_tabs} tabs, expected {expected_tabs}", file=sys.stderr)
                
                data_lines_checked += 1
                if data_lines_checked >= 100:  # Check first 100 lines
                    print(f"Checked first {data_lines_checked} data lines...")
                    break
            
            print(f"Validated {data_lines_checked} data lines")
            
            if errors:
                print(f"\n[FAILED] Found {len(errors)} errors", file=sys.stderr)
                return False
            else:
                print("[PASSED] VCF should work with Beagle")
                return True
    
    except Exception as e:
        print(f"[FATAL] Error reading VCF: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: beagle_exact_validator.py <input.vcf.gz>", file=sys.stderr)
        sys.exit(2)
    
    vcf_path = sys.argv[1]
    success = validate_vcf_for_beagle(vcf_path)
    sys.exit(0 if success else 1)

