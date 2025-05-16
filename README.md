# IST-BubbleSort-Parallel

**A Parallel Algorithm for Constructing Multiple Independent Spanning Trees in Bubble-Sort Networks**

## Overview

This project implements a novel **parallel algorithm** for constructing **n-1 independent spanning trees (ISTs)** in bubble-sort networks. It solves an open problem from previous research by Kao et al. by enabling each vertex to determine its parent in **constant time** without using recursion.

Our implementation uses a **hybrid parallel approach**:
- **MPI** for distributed computing across multiple nodes
- **OpenMP** for thread-level parallelism within each node

---

## Background

**Bubble-sort networks** are a class of interconnection networks represented as graphs where:
- **Vertices** are permutations of numbers (e.g., `1234`, `4231`)
- **Edges** connect permutations that differ by one adjacent swap

**Independent spanning trees** have significant applications in:
- Fault-tolerant network communication
- Secure message distribution
- Enhancing reliability of interconnection networks

---

## Algorithm

The algorithm allows each vertex to independently compute its parent in each spanning tree using a set of rules based on:
- The vertex's permutation representation
- The last digits of the permutation
- The tree index being constructed

### Key Features
- **Non-recursive approach**
- **Constant-time** parent determination per vertex
- **Optimal total time complexity**: `O(n · n!)`
- **Tree height**: at most `D(Bn) + n - 1`

---

## Implementation

Implemented in **C++** using:
- **MPI** (e.g., MPICH or OpenMPI) for distributed memory parallelism
- **OpenMP** for shared memory parallelism

### Key Components
- **Permutation Class**: Represents vertices in the bubble-sort network and supports permutation operations
- **Parent Finding Algorithm**: Implements Algorithm 1 from the reference paper
- **Hybrid Parallelism**:
  - MPI distributes permutations across processes
  - OpenMP computes parents in parallel within each process
- **Result Collection**: Master process gathers and formats all outputs

---

## Prerequisites

- C++ compiler with OpenMP support
- MPI library (MPICH or OpenMPI)
- CMake (for build system)

---

## Building

```bash
mkdir build
cd build
cmake ..
make
```

---

## Running

```bash
mpirun -np <num_processes> ./ist_bubble_sort -n <dimension>
```

- `<num_processes>`: Number of MPI processes
- `<dimension>`: Dimension of the bubble-sort network (n)

---

## Performance

- Efficient scaling with constant-time parent determination per vertex
- Distributed workload for large-scale computation
- OpenMP-based parallel processing within each MPI process

**Note:** For large networks (n ≥ 7), distributed computing is crucial due to factorial growth (`n!`) in vertex count.

---

## Output Files

For `n < 10`, the program generates the following files:
- `bubble_sort_ist_table_n<dimension>.txt`: Table of parents for each vertex in each tree
- `bubble_sort_ist_vertices_n<dimension>.txt`: Mapping of vertices to MPI processes
- `bubble_sort_ist_log_n<dimension>.txt`: Performance statistics

---

## Implementation Details

### Parallel Workflow

#### Permutation Distribution
- Each MPI process receives roughly `n!/p` permutations
- Ensures balanced workload distribution

#### OpenMP Parallelization
- Parent computations run in parallel within MPI processes
- Uses thread-local storage to reduce contention

#### Result Collection
- Results sent from all processes to the master
- Master combines and formats outputs

### Key Optimizations
- Thread-local intermediate results (no synchronization)
- Dynamic scheduling for thread load balancing
- Efficient, non-recursive permutation processing
- Minimized MPI communication overhead

---

## Contributors

- **Sana Mir**
- **Ayesha Kiani**
- **Abdulrehman**

---

## Citation

This work is based on the paper:

> Kao, Shih-Shun, et al.  
> "A parallel algorithm for constructing multiple independent spanning trees in bubble-sort networks."  
> *Journal of Parallel and Distributed Computing*, vol. 117, 2023, pp. 139–145.

```bibtex
@article{kao2023parallel,
  title={A parallel algorithm for constructing multiple independent spanning trees in bubble-sort networks},
  author={Kao, Shih-Shun and Klasing, Ralf and Hung, Ling-Ju and Lee, Chia-Wei and Hsieh, Sun-Yuan},
  journal={Journal of Parallel and Distributed Computing},
  volume={117},
  pages={139--145},
  year={2023}
}
```

---

## License

This project is licensed under the [MIT License](./LICENSE).
