# a2_MSD Pipeline App - Optimized Version

A comprehensive **molecular dynamics analysis pipeline** for generating **Mean‑Square‑Displacement (MSD)** and **non‑Gaussian parameter** analyses from NAMD DCD trajectories—now with **3-10x performance improvements**, parallel processing, and enhanced reliability.

**Available Calculations:**
- **α₂(t) and MSD**: Standard non-Gaussian parameter α₂(t) = 3⟨Δr⁴⟩/(5⟨Δr²⟩²) - 1 with mean square displacement
- **α_xz(t)**: Directional correlation parameter α_xz(t) = ⟨Δx²·Δz²⟩/(⟨Δx²⟩·⟨Δz²⟩) - 1 for anisotropic analysis

## 🚀 **New in Optimized Version v2.0**

- **🖥️ Windows Support**: Now available for Linux, macOS (Intel & ARM), and Windows platforms
- **3-10x Performance Boost**: Highly optimized pipeline functions with parallel processing
- **Scrollable GUI**: Responsive interface that works on any screen size
- **Parallel Processing**: Multi-core support for all pipeline steps
- **Enhanced Reliability**: Comprehensive error handling and data validation
- **Memory Optimization**: Efficient processing for large trajectories
- **Performance Benchmarking**: Built-in testing and comparison tools

---

## Available Formats

| Artifact                                                | Description                                                       | Platform                           |
| ------------------------------------------------------- | ----------------------------------------------------------------- | ---------------------------------- |
| `a2_MSD_pipeline.py` **(GUI + Compiled Modules)**      | **Main GUI with optimized compiled modules (.so files)**         | **Linux • macOS • Windows**       |
| `alpha2_MSD_pip_Linux`                                  | Stand‑alone Linux executable                                     | Linux                              |
| `alpha2_MSD_pip_mac_x86_64`                            | Stand‑alone macOS executable (Intel)                             | macOS Intel                        |
| `alpha2_MSD_pip_mac_arm64`                             | Stand‑alone macOS executable (Apple Silicon)                     | macOS ARM                          |
| `alpha2_MSD_pip_Windows.exe`                           | Stand‑alone Windows executable                                   | Windows                            |

---

## Quick Start

### **Option 1: Source Version**

```bash
# 1 – Clone the repository
$ git clone https://github.com/AdamPirnia/a2_MSD_pipeline.git
$ cd a2_MSD_pipeline

# 2 – Install dependencies (if needed)
$ sudo apt update && sudo apt install python3-tk -y  # Ubuntu/Debian
# or
$ brew install python-tk  # macOS with Homebrew

# 3 – Run the optimized GUI
$ python3 a2_MSD_pipeline.py
```

### **Option 2: Standalone Executable**

**🖱️ Double-Click to Run (All Platforms):**
Simply double-click the appropriate executable file for your operating system - no additional setup required!

**Command Line Alternative:**
```bash
# Linux
$ chmod +x alpha2_MSD_pip_Linux
$ ./alpha2_MSD_pip_Linux

# macOS 
$ xattr -dr com.apple.quarantine alpha2_MSD_pip_mac_arm64  # ARM CPUs
$ chmod +x alpha2_MSD_pip_mac_arm64 
$ ./alpha2_MSD_pip_mac_arm64 

# or for Intel CPUs
$ xattr -dr com.apple.quarantine alpha2_MSD_pip_mac_x86_64
$ chmod +x alpha2_MSD_pip_mac_x86_64
$ ./alpha2_MSD_pip_mac_x86_64

# Windows
# Simply double-click alpha2_MSD_pip_Windows.exe
# or run from command prompt:
> alpha2_MSD_pip_Windows.exe
```

---

## 🎯 **Enhanced GUI & Workflow**

![An example](example.png)

### **Basic Setup**
1. **Common Parameters** – Base directory, number of DCDs, **max workers** (auto-detected)
2. **Calculation Type Selection** – Choose between α₂(t)/MSD or α_xz(t) analysis
3. **Step-by-Step Optimization** – Configure parallel processing for each pipeline step
4. **SLURM Configuration** – Partition, wall‑time, CPUs (auto-populated), email
5. **Generate & Benchmark** – Create optimized scripts and performance tests

### **Advanced Features**
- **Scrollable interface**: Works perfectly on laptops and small screens
- **Parallel Processing**: Enable/disable for each step individually with auto-detected CPU cores
- **Memory Management**: Chunked processing and memory mapping for large datasets
- **Data Validation**: Quality checks and comprehensive error reporting
- **Progress Monitoring**: Detailed execution summaries with timing information

---

## 🚀 **Running the Generated Pipeline**

After generating your pipeline files through the GUI, you have two execution options:

### **Option 1: Local Execution**
Run the pipeline directly on your current machine:

```bash
# Navigate to the generated folder
$ cd your_output_folder_name

# Execute the pipeline
$ python your_main_script.py
```

### **Option 2: Supercomputer Execution**
Copy the entire generated folder to your target supercomputer and submit to SLURM:

```bash
# 1. Copy the entire folder to your supercomputer
$ scp -r your_output_folder_name username@supercomputer.edu:~/

# 2. Log into the supercomputer
$ ssh username@supercomputer.edu

# 3. Navigate to the copied folder
$ cd your_output_folder_name

# 4. Submit the SLURM job
$ sbatch your_submit_script.sh
```

**🎯 Key Points:**
- **Self-contained**: The generated folder includes everything needed (`main_functions`, scripts, etc.)
- **Portable**: Simply copy the entire folder - no additional setup required
- **Ready to run**: Both local Python execution and SLURM submission work immediately

---

## What the Optimized Pipeline Does

| Step | Task                                            | **New Optimizations**                                    | Output             |
| ---- | ----------------------------------------------- | -------------------------------------------------------- | ------------------ |
|  1   | Extract raw coordinates (`coordinates_extract`) | **Parallel VMD execution, timeout protection** | User‑chosen OUTdir |
|  2   | Unwrap PBC (`unwrap_coords`)                    | **Chunked processing, parallel files, auto-scaling**     | User‑chosen OUTdir |
|  3   | Center‑of‑Mass calc (`COM_calc`)                | **Vectorized NumPy operations, memory mapping, 10x faster** | User‑chosen OUTdir |
|  4   | **Statistical Analysis** (`alpha2_MSD` or `alpha_xz`) | **Dual calculation modes, numerical stability, enhanced validation** | User‑chosen OUTdir |

### **Step 4 Calculation Options:**
- **α₂(t) and MSD**: Computes standard non-Gaussian parameter and mean square displacement
- **α_xz(t)**: Computes directional correlation parameter for anisotropic diffusion analysis

### **Performance Improvements**
- **coordinates_extract**: Parallel processing with timeout protection, 3-5x faster
- **unwrap_coords**: Memory-efficient chunked processing, 3-5x faster
- **COM_calc**: Highly optimized vectorized operations, ~10x performance boost
- **alpha2_MSD**: Enhanced numerical stability for reliable calculations
- **alpha_xz**: New directional correlation analysis with optimized framework

---

## 📊 **Performance Benchmarking**

The optimized pipeline includes built-in benchmarking tools:

```bash
# Generate benchmark script from GUI
1. Open a2_MSD_pipeline.py
2. Fill in base directory
3. Click "Benchmark Performance"
4. Run the generated benchmark script

# Manual benchmarking
$ python3 performance_benchmark.py
```

**Expected Performance Gains:**
- Small systems (< 1000 particles): **3-5x faster**
- Medium systems (1000-10000 particles): **5-8x faster**  
- Large systems (> 10000 particles): **8-10x faster**

---



## Requirements

### **GUI + Compiled Modules Version**
- **Python 3.6+** with standard library
- **tkinter**: GUI framework (`python3-tk` package on Linux, included with Python on Windows/macOS)
- **NumPy**: For optimized calculations (auto-installed with most Python distributions)
- **VMD**: Required for coordinate extraction (Step 1)
- **Compiled modules**: Pre-compiled .so files included for optimal performance

### **Executable Version**
- No additional requirements—everything is bundled
- Available for **Linux**, **macOS (Intel & ARM)**, and **Windows**

---

## 🛠 **Advanced Configuration**

<dl>
  <dt><strong>Parallel Workers</strong></dt>
  <dd>Auto-detected based on CPU cores. Adjust in GUI for optimal performance on your system.</dd>

  <dt><strong>Chunk Processing</strong></dt>
  <dd>
    Automatically manages memory for large trajectories:<br>
    &emsp;• <code>Chunk Size: auto</code> ← automatically optimizes chunk size<br>
    &emsp;• <code>Chunk Size: 1000</code> ← process 1000 frames at a time
  </dd>

  <dt><strong>Memory Mapping</strong></dt>
  <dd>
    For very large coordinate files:<br>
    &emsp;• <code>Use Memory Map: ✓</code> ← reduces RAM usage for massive datasets
  </dd>

  <dt><strong>VMD Path</strong></dt>
  <dd>
    Specify VMD executable location:<br>
    &emsp;• <code>VMD path: /usr/local/bin/vmd</code><br>
    &emsp;• Use "Browse..." button to select automatically
  </dd>

  <dt><strong>Data Validation</strong></dt>
  <dd>
    Quality checks and error detection:<br>
    &emsp;• <code>Validate Data: ✓</code> ← enables comprehensive data quality reporting
  </dd>
</dl>

---



## 📈 **Generated Script Features**

The optimized GUI generates production-ready scripts with:

- **Comprehensive Error Handling**: Graceful failure recovery and detailed error messages
- **Progress Reporting**: Real-time status updates and timing information
- **Result Validation**: Automatic verification of output quality
- **Performance Metrics**: Execution time tracking for each pipeline step
- **SLURM Integration**: Optimized resource allocation and job management

---

## FAQ

1. **Q: Which version should I use?**  
   **A:** Use the **optimized source version** (`a2_MSD_pipeline.py`) for best performance and latest features.

2. **Q: How much faster is the optimized version?**  
   **A:** Typically **3-10x faster** depending on system size, with additional memory efficiency improvements.

3. **Q: Does it work on small screens?**  
   **A:** Yes! The new scrollable GUI works perfectly on laptops and small displays.

4. **Q: Can I benchmark the performance?**  
   **A:** Absolutely! Use the "Benchmark Performance" button to generate test scripts and measure improvements.

5. **Q: Are the results scientifically equivalent?**  
   **A:** Yes, and **more accurate** due to enhanced numerical stability and optimized calculations.

6. **Q: Can I still use my old configuration?**  
   **A:** Yes, the optimized version is backward compatible with existing setups.

7. **Q: How do I run the generated pipeline files?**  
   **A:** You can either run `python your_main_script.py` locally, or copy the entire folder to a supercomputer and submit with `sbatch your_submit_script.sh`. The generated folder is completely self-contained.

---

## 🏆 **Optimization Summary**

| Component | Original | Optimized | Improvement |
|-----------|----------|-----------|-------------|
| **coordinates_extract** | Serial VMD processing | Parallel VMD execution | **3-5x + Enhanced Reliability** |
| **unwrap_coords** | Memory intensive | Chunked processing | **3-5x + Memory Efficient** |
| **COM_calc** | Loop-based calculations | Vectorized NumPy operations | **~10x faster** |
| **alpha2_MSD** | Unstable numerics | Enhanced stability | **Reliable + Faster** |
| **alpha_xz** | New functionality | Optimized directional analysis | **New + Efficient** |
| **Overall Pipeline** | Sequential processing | Parallel + dual modes | **3-10x end-to-end** |

---

## 🔗 **Access to Source Code**

This repository contains compiled modules (.so files) for optimal performance. If you require access to the complete source code for research, development, or educational purposes, please contact the author directly:

**📧 Contact:** Please email the repository owner to request access to the private source repository. Include a brief description of your intended use case and institutional affiliation (if applicable).

Access to the source code may be granted for:
- Academic research and publications
- Educational use in courses and workshops  
- Collaborative development projects
- Custom modifications and extensions

---

© 2025 Adam Pirnia — All rights reserved.

**Optimized Version Contributors**: Enhanced performance, parallel processing, and reliability improvements.

