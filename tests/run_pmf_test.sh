#!/bin/bash

# PMF End-to-End Test Runner
# This script runs complete PMF workflow testing on a server

set -e  # Exit on error

# Default parameters
INPUT_DIR="./gromacssim"
OUTPUT_DIR="./pmf_test_results"
CONFIG="fast"
PULLING_DISTANCE="2.5"
FORCE_REMODEL=false
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
PMF End-to-End Test Runner

Usage: $0 [OPTIONS]

Options:
    -i, --input DIR         Input system directory (default: ./gromacssim)
    -o, --output DIR        Output directory (default: ./pmf_test_results)
    -c, --config TYPE       Configuration type: fast|default|accurate (default: fast)
    -p, --pulling DIST      Pulling distance in nm (default: 2.5)
    -f, --force-remodel     Force remodeling even if PMF system exists
    -q, --quiet             Quiet mode (less verbose output)
    -h, --help              Show this help message

Examples:
    $0                                          # Use defaults
    $0 -i ./my_system -o ./my_results          # Custom directories
    $0 -c accurate -p 3.0                      # High accuracy, 3nm pulling
    $0 -i ./gromacssim -c fast -q              # Fast test, quiet mode

Requirements:
    - Input directory must contain GROMACS MD system files
    - Python 3.6+ with PRISM PMF module installed
    - GROMACS installation (for actual simulations)
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
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        -p|--pulling)
            PULLING_DISTANCE="$2"
            shift 2
            ;;
        -f|--force-remodel)
            FORCE_REMODEL=true
            shift
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
    print_status "Starting PMF End-to-End Test"
    echo "======================================================"
    echo "Input directory:    $INPUT_DIR"
    echo "Output directory:   $OUTPUT_DIR"
    echo "Configuration:      $CONFIG"
    echo "Pulling distance:   ${PULLING_DISTANCE} nm"
    echo "Verbose mode:       $VERBOSE"
    echo "======================================================"

    # Check if input directory exists
    if [[ ! -d "$INPUT_DIR" ]]; then
        print_error "Input directory not found: $INPUT_DIR"
        echo "Please provide a valid input directory with -i option"
        exit 1
    fi

    # Check for required files in input directory
    print_status "Validating input system..."

    REQUIRED_FILES=(
        "*/solv_ions.gro"
        "*/topol.top"
    )

    for pattern in "${REQUIRED_FILES[@]}"; do
        if ! ls $INPUT_DIR/$pattern 1> /dev/null 2>&1; then
            print_warning "Required file pattern not found: $pattern"
            print_warning "Make sure your input directory contains a complete GROMACS MD system"
        else
            print_success "Found files matching: $pattern"
        fi
    done

    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    print_success "Output directory ready: $OUTPUT_DIR"

    # Find Python script (try robust version first)
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    ROBUST_SCRIPT="$SCRIPT_DIR/prism/pmf/examples/robust_pmf_test.py"
    ORIGINAL_SCRIPT="$SCRIPT_DIR/prism/pmf/examples/end_to_end_pmf_test.py"

    if [[ -f "$ROBUST_SCRIPT" ]]; then
        PYTHON_SCRIPT="$ROBUST_SCRIPT"
        print_status "Using robust PMF test script"
    elif [[ -f "$ORIGINAL_SCRIPT" ]]; then
        PYTHON_SCRIPT="$ORIGINAL_SCRIPT"
        print_status "Using original PMF test script"
    else
        print_error "PMF test script not found in: $SCRIPT_DIR/prism/pmf/examples/"
        exit 1
    fi

    # Build Python command
    PYTHON_CMD="python $PYTHON_SCRIPT"
    PYTHON_CMD="$PYTHON_CMD --input '$INPUT_DIR'"
    PYTHON_CMD="$PYTHON_CMD --output '$OUTPUT_DIR'"
    PYTHON_CMD="$PYTHON_CMD --config '$CONFIG'"
    PYTHON_CMD="$PYTHON_CMD --pulling-distance '$PULLING_DISTANCE'"

    if [[ "$FORCE_REMODEL" == "true" ]]; then
        PYTHON_CMD="$PYTHON_CMD --force-remodel"
    fi

    if [[ "$VERBOSE" == "false" ]]; then
        PYTHON_CMD="$PYTHON_CMD --quiet"
    fi

    print_status "Executing PMF workflow test..."
    echo "Command: $PYTHON_CMD"
    echo ""

    # Run the Python script
    START_TIME=$(date +%s)

    if eval "$PYTHON_CMD"; then
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))

        print_success "PMF test completed successfully!"
        print_success "Total runtime: ${DURATION} seconds"

        # Show results summary
        if [[ -f "$OUTPUT_DIR/pmf_test_report.json" ]]; then
            print_status "Test report generated: $OUTPUT_DIR/pmf_test_report.json"
        fi

        if [[ -f "$OUTPUT_DIR/pmf_test.log" ]]; then
            print_status "Detailed log: $OUTPUT_DIR/pmf_test.log"
        fi

        echo ""
        echo "======================================================"
        print_success "PMF WORKFLOW TEST COMPLETED"
        echo "======================================================"
        echo "Results are available in: $OUTPUT_DIR"
        echo ""
        echo "Next steps:"
        echo "1. Review the test report: $OUTPUT_DIR/pmf_test_report.json"
        echo "2. Check the detailed log: $OUTPUT_DIR/pmf_test.log"
        echo "3. For actual calculations, run the generated scripts manually"
        echo ""

    else
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))

        print_error "PMF test failed after ${DURATION} seconds"

        if [[ -f "$OUTPUT_DIR/pmf_test.log" ]]; then
            print_status "Check error details in: $OUTPUT_DIR/pmf_test.log"
            echo ""
            echo "Last 10 lines of log:"
            tail -10 "$OUTPUT_DIR/pmf_test.log"
        fi

        exit 1
    fi
}

# Trap to handle interruptions
trap 'print_error "Test interrupted by user"; exit 130' INT TERM

# Run main function
main "$@"