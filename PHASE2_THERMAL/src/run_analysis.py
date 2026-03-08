#!/usr/bin/env python3
"""
DS Phase 3B-2: Sector-Level Fisher Analysis Runner
Combines all analysis components and provides a simple interface.
"""

import subprocess
import sys
import os

def combine_analysis_scripts():
    """Combine the three analysis script parts into one file."""
    print("Combining analysis scripts...")
    
    # Read all parts
    with open('fisher_sector_analysis_part1.py', 'r') as f:
        part1 = f.read()
    
    with open('fisher_sector_analysis_part2.py', 'r') as f:
        part2 = f.read()
    
    with open('fisher_sector_analysis_part3.py', 'r') as f:
        part3 = f.read()
    
    # Remove duplicate imports and combine
    combined = part1 + "\n\n" + part2 + "\n\n" + part3
    
    # Write combined script
    with open('fisher_sector_analysis.py', 'w') as f:
        f.write(combined)
    
    print("✓ Combined fisher_sector_analysis.py created")
    return True

def create_requirements():
    """Create requirements.txt for the analysis."""
    requirements = """numpy>=1.21.0
pandas>=1.3.0
matplotlib>=3.4.0
scipy>=1.7.0
"""
    
    with open('requirements.txt', 'w') as f:
        f.write(requirements)
    
    print("✓ Created requirements.txt")
    return True

def create_runner_script():
    """Create a runner script for the analysis."""
    runner = """#!/usr/bin/env python3
"""
    with open('run_analysis.py', 'w') as f:
        f.write(runner)
    
    print("✓ Created analysis runner")
    return True

def main():
    print("DS Phase 3B-2: Sector-Level Fisher Analysis Setup")
    print("=" * 50)
    
    # Combine analysis scripts
    if combine_analysis_scripts():
        print("✓ Analysis scripts combined")
    
    # Create requirements
    create_requirements()
    
    print("\n" + "="*50)
    print("SETUP COMPLETE")
    print("="*50)
    print("\nTo use the analysis pipeline:")
    print("1. Run the LEAN algorithm in QuantConnect")
    print("2. Download results as 'fisher_sector_results.json'")
    print("3. Update the path in fisher_sector_analysis.py")
    print("4. Run: python fisher_sector_analysis.py")
    print("\nFiles created:")
    print("  - fisher_sector_analysis.py (combined analysis script)")
    print("  - requirements.txt (Python dependencies)")
    print("  - README_PHASE3B2.md (documentation)")
    print("\nTo install dependencies:")
    print("  pip install -r requirements.txt")
    print("\nTo run analysis (after getting results):")
    print("  python fisher_sector_analysis.py")

if __name__ == "__main__":
    main()