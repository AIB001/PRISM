#!/bin/bash

# PMF Remodeling Test Runner
# This script specifically tests the PMF remodeling process

set -e  # Exit on error

# Default parameters
INPUT_DIR="./gromacssim"
OUTPUT_DIR="./pmf_remodeling_test"
PULLING_DISTANCE="2.5"
VERBOSE=true

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Help function
show_help() {
    cat << EOF
PMF Remodeling Test Runner

This script tests the PMF system remodeling process specifically,
helping to diagnose and validate the remodeling workflow.

Usage: $0 [OPTIONS]

Options:
    -i, --input DIR         Input system directory (default: ./gromacssim)
    -o, --output DIR        Output directory (default: ./pmf_remodeling_test)
    -p, --pulling DIST      Pulling distance in nm (default: 2.5)
    -q, --quiet             Quiet mode (less verbose output)
    -h, --help              Show this help message

Examples:
    $0                                          # Use defaults
    $0 -i ./my_system -o ./remodel_test        # Custom directories
    $0 -p 3.0                                  # 3nm pulling distance
    $0 -i ./gromacssim -q                      # Quiet mode

Test Coverage:
    1. Input MD system validation
    2. Remodeling script execution
    3. Remodeled system structure validation
    4. PMF workflow compatibility check

Requirements:
    - Input directory must contain GROMACS MD system files
    - Python 3.6+ with PRISM PMF module installed
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--pulling)
            PULLING_DISTANCE="$2"
            shift 2
            ;;
        -q|--quiet)
            VERBOSE=false
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Main execution
main() {
    print_status "Starting PMF Remodeling Test"
    echo "======================================================"
    echo "Input directory:    $INPUT_DIR"
    echo "Output directory:   $OUTPUT_DIR"
    echo "Pulling distance:   $PULLING_DISTANCE nm"
    echo "Verbose mode:       $VERBOSE"
    echo "======================================================"

    # Validate input directory
    print_status "Validating input system..."

    if [[ ! -d "$INPUT_DIR" ]]; then
        print_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi

    # Check for MD files
    if ! find "$INPUT_DIR" -name "solv_ions.gro" | head -1 | read; then
        print_error "No solv_ions.gro found in input directory"
        exit 1
    fi

    if ! find "$INPUT_DIR" -name "topol.top" | head -1 | read; then
        print_error "No topol.top found in input directory"
        exit 1
    fi

    print_success "Input validation passed"

    # Prepare output directory
    print_status "Preparing output directory..."
    mkdir -p "$OUTPUT_DIR"
    print_success "Output directory ready: $OUTPUT_DIR"

    # Find Python script
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PYTHON_SCRIPT="$SCRIPT_DIR/prism/pmf/examples/test_pmf_remodeling.py"

    if [[ ! -f "$PYTHON_SCRIPT" ]]; then
        print_error "PMF remodeling test script not found: $PYTHON_SCRIPT"
        exit 1
    fi

    print_status "Using PMF remodeling test script"

    # Build Python command
    PYTHON_CMD="python $PYTHON_SCRIPT"
    PYTHON_CMD="$PYTHON_CMD --input '$INPUT_DIR'"
    PYTHON_CMD="$PYTHON_CMD --output '$OUTPUT_DIR'"
    PYTHON_CMD="$PYTHON_CMD --pulling-distance '$PULLING_DISTANCE'"

    if [[ "$VERBOSE" == "false" ]]; then
        PYTHON_CMD="$PYTHON_CMD --quiet"
    fi

    print_status "Executing PMF remodeling test..."
    echo "Command: $PYTHON_CMD"
    echo ""

    # Run the Python script
    START_TIME=$(date +%s)

    if eval "$PYTHON_CMD"; then
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))

        print_success "PMF remodeling test completed successfully!"
        print_success "Total runtime: ${DURATION} seconds"

        # Show results summary
        if [[ -f "$OUTPUT_DIR/pmf_remodeling_test_report.json" ]]; then
            print_status "Test report generated:"
            echo "  📄 $OUTPUT_DIR/pmf_remodeling_test_report.json"
        fi

        if [[ -d "$OUTPUT_DIR/GMX_PMF_SYSTEM" ]]; then
            print_status "Remodeled PMF system:"
            echo "  📁 $OUTPUT_DIR/GMX_PMF_SYSTEM/"
            ls -la "$OUTPUT_DIR/GMX_PMF_SYSTEM/" | sed 's/^/    /'
        fi

        echo ""
        echo "======================================================"
        echo -e "${GREEN}[SUCCESS] PMF REMODELING TEST COMPLETED${NC}"
        echo "======================================================"
        echo "Results are available in: $OUTPUT_DIR"
        echo ""
        echo "Next steps:"
        echo "1. Review the test report: $OUTPUT_DIR/pmf_remodeling_test_report.json"
        echo "2. Check the remodeled system: $OUTPUT_DIR/GMX_PMF_SYSTEM/"
        echo "3. If tests passed, the remodeling process is working correctly"
        echo "4. If tests failed, check the error messages for debugging"

    else
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))

        print_error "PMF remodeling test failed!"
        print_error "Runtime: ${DURATION} seconds"

        echo ""
        echo "======================================================"
        echo -e "${RED}[FAILED] PMF REMODELING TEST FAILED${NC}"
        echo "======================================================"
        echo "Check the output above for error details"
        echo ""
        echo "Debugging steps:"
        echo "1. Check input system structure and files"
        echo "2. Verify PRISM PMF installation"
        echo "3. Check Python environment and dependencies"
        echo "4. Review detailed error messages above"

        if [[ -f "$OUTPUT_DIR/pmf_remodeling_test_report.json" ]]; then
            echo "5. Review test report: $OUTPUT_DIR/pmf_remodeling_test_report.json"
        fi

        exit 1
    fi
}

# Run main function
main