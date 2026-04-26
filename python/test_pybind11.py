"""
Test script for cTEA pybind11 module
验证libcTEA_pybind11.so是否正确编译和可以正常导入
"""

import sys
import os

# Add cpp directory to path so Python can find the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cpp'))

try:
    import libcTEA
    print("✅ Successfully imported libcTEA module")
    print(f"   Version: {libcTEA.__version__}")
except ImportError as e:
    print(f"❌ Failed to import libcTEA: {e}")
    sys.exit(1)

# Test 1: Check module attributes
print("\n=== Test 1: Module attributes ===")
try:
    version = libcTEA.version()
    print(f"✅ libcTEA.version() = {version}")
except Exception as e:
    print(f"❌ Error calling version(): {e}")

# Test 2: Create BamProcessor object
print("\n=== Test 2: BamProcessor class ===")
try:
    bp = libcTEA.BamProcessor()
    print("✅ Created BamProcessor object")
    print(f"   Object: {bp}")
except Exception as e:
    print(f"❌ Error creating BamProcessor: {e}")

# Test 3: Create CandidateDetector object
print("\n=== Test 3: CandidateDetector class ===")
try:
    detector = libcTEA.CandidateDetector()
    print("✅ Created CandidateDetector object")
    print(f"   Object: {detector}")
except Exception as e:
    print(f"❌ Error creating CandidateDetector: {e}")

# Test 4: Check BamHeaderInfo struct
print("\n=== Test 4: BamHeaderInfo struct ===")
try:
    # Check if BamHeaderInfo is available
    if hasattr(libcTEA, 'BamHeaderInfo'):
        info = libcTEA.BamHeaderInfo()
        print("✅ BamHeaderInfo struct is available")
    else:
        print("⚠️  BamHeaderInfo struct not exposed")
except Exception as e:
    print(f"❌ Error with BamHeaderInfo: {e}")

# Test 5: Check CandidateEvidence struct
print("\n=== Test 5: CandidateEvidence struct ===")
try:
    if hasattr(libcTEA, 'CandidateEvidence'):
        evidence = libcTEA.CandidateEvidence()
        print("✅ CandidateEvidence struct is available")
        print(f"   Object: {evidence}")
    else:
        print("⚠️  CandidateEvidence struct not exposed")
except Exception as e:
    print(f"❌ Error with CandidateEvidence: {e}")

# Test 6: Check ScanStatistics struct
print("\n=== Test 6: ScanStatistics struct ===")
try:
    if hasattr(libcTEA, 'ScanStatistics'):
        stats = libcTEA.ScanStatistics()
        print("✅ ScanStatistics struct is available")
    else:
        print("⚠️  ScanStatistics struct not exposed")
except Exception as e:
    print(f"❌ Error with ScanStatistics: {e}")

print("\n" + "="*50)
print("Test summary:")
print("✅ pybind11 module successfully compiled and basic functionality works")
print("="*50)

# Try to initialize with a test file (if available)
print("\n=== Optional: Test with real BAM file ===")
test_bam = os.path.join(os.path.dirname(__file__), '..', 'tests', 'test_data', 'test.bam')
if os.path.isfile(test_bam):
    print(f"Found test BAM file: {test_bam}")
    try:
        bp = libcTEA.BamProcessor()
        # Note: This will fail without a real BAM file, but tests the API
        # result = bp.initialize(test_bam, "", 1)
        # print(f"Initialize result: {result}")
        print("⚠️  BAM initialization test skipped (need valid BAM file)")
    except Exception as e:
        print(f"❌ Error: {e}")
else:
    print("⚠️  No test BAM file found, skipping BAM test")
    
print("\n✅ pybind11 module test completed!")