#!/bin/bash
#------------------------------------------------------------------------------
# RGP-13 Real Fluid Models Installer for OpenFOAM-13
#
# Installs SRK EOS, Chung transport, Takahashi diffusion models
# into an existing OpenFOAM-13 source tree, patches Make/files,
# and compiles the modified libraries.
#
# Usage:
#   ./install.sh                          # auto-detect FOAM_SRC
#   ./install.sh /path/to/OpenFOAM-13     # specify OpenFOAM root
#   ./install.sh --dry-run                # show what would be done
#   ./install.sh --uninstall              # remove installed files
#------------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/src"

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

#------------------------------------------------------------------------------
# Parse arguments
#------------------------------------------------------------------------------
DRY_RUN=false
UNINSTALL=false
FOAM_ROOT=""

for arg in "$@"; do
    case "$arg" in
        --dry-run)   DRY_RUN=true ;;
        --uninstall) UNINSTALL=true ;;
        --help|-h)
            echo "Usage: $0 [--dry-run|--uninstall] [/path/to/OpenFOAM-13]"
            exit 0 ;;
        *)           FOAM_ROOT="$arg" ;;
    esac
done

#------------------------------------------------------------------------------
# Locate OpenFOAM-13
#------------------------------------------------------------------------------
if [ -z "$FOAM_ROOT" ]; then
    if [ -n "${FOAM_SRC:-}" ]; then
        FOAM_ROOT="$(dirname "$FOAM_SRC")"
        info "Using FOAM_SRC environment: $FOAM_ROOT"
    elif [ -d "$SCRIPT_DIR/../OpenFOAM-13" ]; then
        FOAM_ROOT="$SCRIPT_DIR/../OpenFOAM-13"
        info "Auto-detected: $FOAM_ROOT"
    else
        error "Cannot find OpenFOAM-13. Provide path or set FOAM_SRC."
        exit 1
    fi
fi

# Resolve to absolute path
FOAM_ROOT="$(cd "$FOAM_ROOT" && pwd)"
SPECIE="$FOAM_ROOT/src/thermophysicalModels/specie"
MCTHERMO="$FOAM_ROOT/src/thermophysicalModels/multicomponentThermo"

# Verify
if [ ! -f "$SPECIE/Make/files" ]; then
    error "Not a valid OpenFOAM-13 tree: $SPECIE/Make/files not found"
    exit 1
fi

info "OpenFOAM root: $FOAM_ROOT"

#------------------------------------------------------------------------------
# Target paths
#------------------------------------------------------------------------------
declare -A TARGETS=(
    ["rfSpecie"]="$SPECIE/rfSpecie"
    ["SRKGas"]="$SPECIE/equationOfState/SRKGas"
    ["chungTransport"]="$SPECIE/transport/chungTransport"
    ["SRKchungTakaMixture"]="$MCTHERMO/mixtures/SRKchungTakaMixture"
    ["include"]="$SPECIE/include"
)

# Files to copy per module
declare -A FILES=(
    ["rfSpecie"]="rfSpecie.H rfSpecieI.H rfSpecie.C"
    ["SRKGas"]="SRKGas.H SRKGasI.H SRKGas.C"
    ["chungTransport"]="chungTransport.H chungTransportI.H chungTransport.C"
    ["SRKchungTakaMixture"]="SRKchungTakaMixture.H SRKchungTakaMixture.C"
    ["include"]="forRealFluidGases.H"
)

# Lines to add to Make/files
SPECIE_MAKE_LINES=(
    "rfSpecie/rfSpecie.C"
    "equationOfState/SRKGas/SRKGas.C"
    "transport/chungTransport/chungTransport.C"
)

#------------------------------------------------------------------------------
# Uninstall
#------------------------------------------------------------------------------
if $UNINSTALL; then
    info "=== Uninstalling RGP-13 Real Fluid Models ==="

    for module in "${!TARGETS[@]}"; do
        target="${TARGETS[$module]}"
        if [ "$module" = "include" ]; then
            # Only remove the specific file, not the directory
            f="$target/forRealFluidGases.H"
            if [ -f "$f" ]; then
                info "Removing $f"
                $DRY_RUN || rm -f "$f"
            fi
        else
            if [ -d "$target" ]; then
                info "Removing directory $target"
                $DRY_RUN || rm -rf "$target"
            fi
        fi
    done

    # Restore Make/files from backup
    for bak in "$SPECIE/Make/files.bak.rgp13" "$MCTHERMO/Make/files.bak.rgp13"; do
        if [ -f "$bak" ]; then
            info "Restoring $(dirname "$bak")/files from backup"
            $DRY_RUN || mv "$bak" "${bak%.bak.rgp13}"
        fi
    done

    info "=== Uninstall complete ==="
    exit 0
fi

#------------------------------------------------------------------------------
# Install: Copy source files
#------------------------------------------------------------------------------
info "=== Installing RGP-13 Real Fluid Models ==="

for module in "${!TARGETS[@]}"; do
    target="${TARGETS[$module]}"
    src_sub="$SRC_DIR/$module"

    if [ ! -d "$src_sub" ] && [ "$module" != "include" ]; then
        error "Source directory not found: $src_sub"
        exit 1
    fi

    $DRY_RUN || mkdir -p "$target"

    for f in ${FILES[$module]}; do
        src_file="$src_sub/$f"
        dst_file="$target/$f"
        if [ ! -f "$src_file" ]; then
            error "Source file not found: $src_file"
            exit 1
        fi
        if $DRY_RUN; then
            info "[DRY-RUN] cp $src_file -> $dst_file"
        else
            cp "$src_file" "$dst_file"
            info "Copied: $module/$f"
        fi
    done
done

#------------------------------------------------------------------------------
# Patch Make/files: specie
#------------------------------------------------------------------------------
SPECIE_MAKEFILE="$SPECIE/Make/files"
MARKER="# >>> RGP-13 real fluid models >>>"
END_MARKER="# <<< RGP-13 real fluid models <<<"

patch_makefile() {
    local makefile="$1"
    shift
    local lines=("$@")
    local marker="$MARKER"
    local end_marker="$END_MARKER"

    # Already patched?
    if grep -q "$marker" "$makefile" 2>/dev/null; then
        warn "Already patched: $makefile (skipping)"
        return 0
    fi

    # Backup
    if [ ! -f "${makefile}.bak.rgp13" ]; then
        cp "$makefile" "${makefile}.bak.rgp13"
        info "Backup: ${makefile}.bak.rgp13"
    fi

    # Insert before the LIB = line
    local tmpfile
    tmpfile=$(mktemp)
    local inserted=false

    while IFS= read -r line; do
        if [[ "$line" == LIB* ]] && ! $inserted; then
            echo "" >> "$tmpfile"
            echo "$marker" >> "$tmpfile"
            for entry in "${lines[@]}"; do
                echo "$entry" >> "$tmpfile"
            done
            echo "$end_marker" >> "$tmpfile"
            echo "" >> "$tmpfile"
            inserted=true
        fi
        echo "$line" >> "$tmpfile"
    done < "$makefile"

    mv "$tmpfile" "$makefile"
    info "Patched: $makefile"
}

if $DRY_RUN; then
    info "[DRY-RUN] Would patch $SPECIE_MAKEFILE with:"
    for l in "${SPECIE_MAKE_LINES[@]}"; do echo "  $l"; done
else
    patch_makefile "$SPECIE_MAKEFILE" "${SPECIE_MAKE_LINES[@]}"
fi

# multicomponentThermo Make/files — currently the mixture is header-only template,
# so it gets compiled through the existing psiMulticomponentThermos.C instantiation.
# No Make/files change needed for multicomponentThermo unless explicit instantiation is added.

#------------------------------------------------------------------------------
# Compile
#------------------------------------------------------------------------------
if $DRY_RUN; then
    info "[DRY-RUN] Would compile:"
    info "  cd $SPECIE && wmake libso"
    info "  cd $MCTHERMO && wmake libso"
    info "=== Dry run complete ==="
    exit 0
fi

info ""
info "=== Compiling specie library ==="
cd "$SPECIE" && wmake libso 2>&1 | tail -20
SPECIE_RC=${PIPESTATUS[0]}

if [ $SPECIE_RC -ne 0 ]; then
    error "specie compilation failed (exit code $SPECIE_RC)"
    error "Run 'cd $SPECIE && wmake libso' manually to see full errors."
    exit 1
fi
info "specie library compiled successfully"

info ""
info "=== Compiling multicomponentThermo library ==="
cd "$MCTHERMO" && wmake libso 2>&1 | tail -20
MCTHERMO_RC=${PIPESTATUS[0]}

if [ $MCTHERMO_RC -ne 0 ]; then
    error "multicomponentThermo compilation failed (exit code $MCTHERMO_RC)"
    error "Run 'cd $MCTHERMO && wmake libso' manually to see full errors."
    exit 1
fi
info "multicomponentThermo library compiled successfully"

info ""
info "=== Installation complete ==="
info "To use in a case, set thermoType in constant/physicalProperties:"
info "  type            hePsiThermo;"
info "  mixture         SRKchungTakaMixture;"
info "  transport       chung;"
info "  thermo          janaf;"
info "  equationOfState SRKGas;"
info "  specie          rfSpecie;"
info "  energy          sensibleEnthalpy;"
