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

    warn "MMseqs2 not found (optional — built-in CPU clustering will be used)."

    if command -v conda &>/dev/null; then
        read -p "Install MMseqs2 via conda? [y/N] " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            conda install -y -c bioconda mmseqs2
            info "MMseqs2 installed."
        fi
    else
        info "To install MMseqs2 later:"
        info "  conda install -c bioconda mmseqs2"
        info "  or visit: https://github.com/soedinglab/MMseqs2"
    fi
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
    local mode="${1:-install}"

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