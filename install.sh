#!/usr/bin/env bash
# PanMiner Installation Script
# Usage:
#   bash install.sh          — build and install PanMiner
#   bash install.sh --dev     — build in debug mode
#   bash install.sh --uninstall — remove PanMiner

set -euo pipefail

BINARY_NAME="panminer"
INSTALL_DIR="${HOME}/.local/bin"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }

# ── Check dependencies ──────────────────────────────────────────────

check_rust() {
    if command -v rustc &>/dev/null && command -v cargo &>/dev/null; then
        local version=$(rustc --version | awk '{print $2}')
        info "Found Rust $version"
        return 0
    fi
    return 1
}

install_rust() {
    if check_rust; then return 0; fi

    warn "Rust not found."

    # Check if conda is available
    if command -v conda &>/dev/null; then
        info "Conda detected. Installing Rust via conda..."
        conda install -y -c conda-forge rust cargo
        return 0
    fi

    # Fall back to rustup
    info "Installing Rust via rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "${HOME}/.cargo/env"
    info "Rust installed via rustup."
}

# ── GPU Detection ────────────────────────────────────────────────────

check_gpu() {
    # Check if nvidia-smi is available and reports a GPU
    if command -v nvidia-smi &>/dev/null; then
        # Query GPU names - returns CSV with header "name" and GPU names below
        if nvidia-smi --query-gpu=name --format=csv 2>/dev/null | grep -q .; then
            # Check if output has more than just the header
            local gpu_count=$(nvidia-smi --query-gpu=name --format=csv 2>/dev/null | tail -n +2 | grep -c .)
            if [[ $gpu_count -gt 0 ]]; then
                info "NVIDIA GPU detected"
                return 0
            fi
        fi
    fi
    return 1
}

check_mmseqs2() {
    if command -v mmseqs &>/dev/null; then
        local version=$(mmseqs version 2>/dev/null | head -1 || echo "unknown")
        info "Found MMseqs2 ($version)"
        return 0
    fi
    return 1
}

install_mmseqs2() {
    if check_mmseqs2; then return 0; fi

    # Check if we should skip MMseqs2 or GPU detection
    if [[ "${PANMINER_NO_MMSEQS2:-}" == "1" ]]; then
        info "MMseqs2 skipped (PANMINER_NO_MMSEQS2=1). Built-in CPU clustering will be used."
        return 0
    fi

    # Check for GPU and decide if we should prompt for MMseqs2
    local gpu_detected=false
    if [[ "${PANMINER_NO_GPU:-}" != "1" ]] && check_gpu 2>/dev/null; then
        gpu_detected=true
    fi

    # If GPU detected, offer GPU version; otherwise offer CPU version or skip
    if [[ "$gpu_detected" == "true" ]]; then
        info "NVIDIA GPU detected. GPU-accelerated clustering will be available."
        if command -v conda &>/dev/null; then
            read -p "Install MMseqs2 with GPU support via conda? [y/N] " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                conda install -y -c bioconda mmseqs2
                info "MMseqs2 with GPU support installed."
                return 0
            fi
        fi
    else
        if command -v conda &>/dev/null; then
            read -p "Install MMseqs2 (CPU version) via conda? [y/N] " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                conda install -y -c bioconda mmseqs2
                info "MMseqs2 installed."
                return 0
            fi
        fi
    fi

    info "MMseqs2 not installed. Built-in CPU clustering will be used."
}

# ── Build ───────────────────────────────────────────────────────────

build_release() {
    info "Building PanMiner (release mode)..."
    cargo build --release
    info "Build complete: target/release/${BINARY_NAME}"
}

build_dev() {
    info "Building PanMiner (debug mode)..."
    cargo build
    info "Build complete: target/debug/${BINARY_NAME}"
}

# ── Install ──────────────────────────────────────────────────────────

install_binary() {
    local src="target/release/${BINARY_NAME}"
    [[ "${1:-}" == "--dev" ]] && src="target/debug/${BINARY_NAME}"

    if [[ ! -f "$src" ]]; then
        error "Binary not found at $src. Run build first."
    fi

    mkdir -p "$INSTALL_DIR"

    # Check if already in PATH
    if echo ":$PATH:" | grep -q ":$INSTALL_DIR:"; then
        : # already in PATH
    else
        warn "$INSTALL_DIR is not in your PATH."
        info "Add it with:  echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> ~/.bashrc && source ~/.bashrc"
    fi

    cp "$src" "${INSTALL_DIR}/${BINARY_NAME}"
    chmod +x "${INSTALL_DIR}/${BINARY_NAME}"
    info "Installed ${BINARY_NAME} to ${INSTALL_DIR}/"

    # Verify
    if command -v "$BINARY_NAME" &>/dev/null; then
        info "Verification: $(${BINARY_NAME} --version 2>/dev/null || echo 'installed successfully')"
    else
        info "Restart your shell or run: source ~/.bashrc"
    fi
}

# ── Uninstall ────────────────────────────────────────────────────────

uninstall() {
    if [[ -f "${INSTALL_DIR}/${BINARY_NAME}" ]]; then
        rm "${INSTALL_DIR}/${BINARY_NAME}"
        info "Uninstalled ${BINARY_NAME} from ${INSTALL_DIR}/"
    else
        warn "${BINARY_NAME} not found in ${INSTALL_DIR}/"
    fi
}

# ── Main ─────────────────────────────────────────────────────────────

main() {
    # Parse arguments for flags first
    local no_gpu=false
    local no_mmseqs2=false
    local mode="install"

    for arg in "$@"; do
        case "$arg" in
            --no-gpu)
                no_gpu=true
                ;;
            --no-mmseqs2)
                no_mmseqs2=true
                ;;
            --dev)
                mode="dev"
                ;;
            --uninstall)
                mode="uninstall"
                ;;
            --help|-h)
                mode="help"
                ;;
        esac
    done

    # Set environment variables for subfunctions based on flags
    if [[ "$no_gpu" == "true" ]]; then
        export PANMINER_NO_GPU=1
    fi
    if [[ "$no_mmseqs2" == "true" ]]; then
        export PANMINER_NO_MMSEQS2=1
    fi

    case "$mode" in
        --dev)
            install_rust
            install_mmseqs2
            build_dev
            install_binary --dev
            ;;
        --uninstall)
            uninstall
            ;;
        --help|-h)
            echo "Usage: bash install.sh [OPTION]"
            echo ""
            echo "Options:"
            echo "  (no option)    Build release and install"
            echo "  --dev           Build debug and install"
            echo "  --uninstall     Remove PanMiner"
            echo "  --no-gpu        Skip GPU detection and MMseqs2 GPU installation"
            echo "  --no-mmseqs2    Skip MMseqs2 installation entirely"
            echo "  --help          Show this help"
            ;;
        install|"")
            install_rust
            install_mmseqs2
            build_release
            install_binary
            ;;
        *)
            error "Unknown option: $mode. Use --help for usage."
            ;;
    esac
}

main "$@"