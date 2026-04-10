---
name: development_approach
description: User's approach to PanMiner development
type: user
---

User wants to build PanMiner (Panaroo replacement in Rust) to production-ready status with:
1. Systematic analysis of current state
2. Brainstorming for design decisions
3. Code implementation with cargo test/cargo check iterations
4. Web search for best tools/algorithms when needed
5. Memory maintenance of progress (done, in-progress, todo)

Key constraints:
- Use modern alternatives where Panaroo uses Python/CD-HIT/NetworkX
- GPU acceleration via MMseqs2 (already implemented)
- Rust implementation with memory-mapped I/O, DashMap, Rayon
- Maintain Panaroo/Roary output compatibility for drop-in replacement
