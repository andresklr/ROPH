# ROPH: A Robust, Optimized, and Parallelized Harris Detector

**ROPH** (Robust, Optimized, and Parallelized Harris) is a high-performance feature detection framework designed for multi-core architectures. It introduces a configurable FAST-based pruning mechanism that restores the Harris detector's ability to detect both corners and edges while achieving significant speedups over traditional implementations.

## 🚀 Key Features

* **FAST-Based Pruning:** A flexible first-stage filter that significantly reduces the computational burden of the Harris detector.
* **High Performance:** Achieves up to **9.24x speedup** over baseline implementations through efficient OpenMP parallelization and cache-friendly data structures.
* **Robustness:** Maintained accuracy across diverse real-world scenes, including urban landscapes, industrial environments, and traffic conditions.
* **Configurable Detection:** Unlike other optimized variants, ROPH allows for the precise detection of corners, edges, or both, depending on the application needs.

## 📄 Publication

This work was accepted for presentation and publication at the:
**16th Workshop on Applications for Multi-Core Architectures (WAMCA 2025)**
*Held in conjunction with the 37th IEEE/SBC International Symposium on Computer Architecture and High Performance Computing (SBAC-PAD 2025).*

## 👥 Authors & Credits

This research was developed by a collaborative team from the following institutions:

| Author | Institution | Email |
| :--- | :--- | :--- |
| **Andres Giraldo Morales** | Rio de Janeiro State University (UERJ) | andres1@outlook.com |
| **Cristiana Bentes** | Rio de Janeiro State University (UERJ) | cristiana@ime.uerj.br |
| **Maria Clicia Castro** | Rio de Janeiro State University (UERJ) | clicia@ime.uerj.br |
| **Claude Tadonki** | Mines Paris - PSL | claude.tadonki@minesparis.psl.eu |
| **Gilson Costa** | Rio de Janeiro State University (UERJ) | gilson@ime.uerj.br |

### Affiliated Institutions
* **UERJ** - Universidade do Estado do Rio de Janeiro, Brazil.
* **Mines Paris - PSL** - Center for Mathematical Morphology (CMM), France.

## 🛠️ Requirements & Installation

* C++ Compiler (GCC 9+ or MSVC)
* OpenMP 4.5+
* OpenCV (for image I/O)

```bash
git clone [https://github.com/yourusername/ROPH.git](https://github.com/yourusername/ROPH.git)
cd ROPH
# Add your specific build commands here (e.g., make, cmake, etc.)