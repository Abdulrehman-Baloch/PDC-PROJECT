#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <map>
#include <set>
#include <mpi.h>
#include <omp.h>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <iostream>

// Function to format and print the table
void formatAndPrintTable(const std::vector<struct ParentEdge>& edges, int n, int numTrees, std::ofstream& outFile) {
    outFile << "Parent Table for Bubble Sort Network ISTs (n=" << n << "):" << std::endl;
    outFile << "==============================================" << std::endl;
    outFile << std::setw(15) << "Vertex" << " | " << std::setw(5) << "Tree" << " | " 
            << std::setw(10) << "Last Symbol" << " | " << std::setw(10) << "Rule" << " | " 
            << std::setw(15) << "Parent" << std::endl;
    outFile << "----------------------------------------------" << std::endl;
    
    // Group by vertex for better readability
    std::map<std::string, std::vector<struct ParentEdge>> vertexGroups;
    for (const auto& edge : edges) {
        vertexGroups[edge.vertex].push_back(edge);
    }
    
    for (const auto& group : vertexGroups) {
        for (const auto& edge : group.second) {
            outFile << std::setw(15) << edge.vertex << " | " << std::setw(5) << edge.treeId << " | " 
                    << std::setw(10) << edge.lastSymbol << " | " << std::setw(10) << edge.rule << " | " 
                    << std::setw(15) << edge.parent << std::endl;
        }
        outFile << "----------------------------------------------" << std::endl;
    }
}

// Permutation class to represent vertices in bubble-sort network
class Permutation {
private:
    std::vector<int> data;
    
public:
    Permutation() {}
    
    Permutation(int n) {
        data.resize(n);
        for (int i = 0; i < n; i++) {
            data[i] = i + 1;
        }
    }
    
    Permutation(const std::vector<int>& values) : data(values) {}
    
    Permutation(const std::string& str) {
        for (char c : str) {
            data.push_back(c - '0');
        }
    }
    
    int at(int i) const {
        if (i < 0 || i >= static_cast<int>(data.size())) {
            return -1;
        }
        return data[i];
    }
    
    int size() const {
        return static_cast<int>(data.size());
    }
    
    int findPosition(int x) const {
        for (int i = 0; i < static_cast<int>(data.size()); i++) {
            if (data[i] == x) return i;
        }
        return -1;
    }
    
    int findRightmostWrongPosition() const {
        for (int i = size() - 1; i >= 0; i--) {
            if (data[i] != i + 1) return i;
        }
        return -1;
    }
    
    Permutation swap(int i) const {
        if (i < 0 || i >= static_cast<int>(data.size()) - 1) {
            return *this;
        }
        
        Permutation result = *this;
        std::swap(result.data[i], result.data[i+1]);
        return result;
    }
    
    Permutation swapSymbol(int x) const {
        int pos = findPosition(x);
        if (pos == -1 || pos >= size() - 1) {
            return *this;
        }
        return swap(pos);
    }
    
    bool isIdentity() const {
        for (int i = 0; i < static_cast<int>(data.size()); i++) {
            if (data[i] != i + 1) return false;
        }
        return true;
    }
    
    std::string toString() const {
        std::stringstream ss;
        for (int val : data) {
            ss << val;
        }
        return ss.str();
    }
    
    bool operator==(const Permutation& other) const {
        return data == other.data;
    }
    
    bool operator<(const Permutation& other) const {
        return data < other.data;
    }
    
    const std::vector<int>& getData() const {
        return data;
    }
};

long long factorial(int n) {
    if (n <= 1) return 1;
    long long result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

// Implementation of Algorithm 1 from the paper
Permutation findParent(const Permutation& v, int t, int n) {
    // Handle edge cases
    if (v.isIdentity() || t < 1 || t >= n) {
        return v;
    }
    
    // Implementation of Algorithm 1
    if (v.at(n-1) == n) {
        // Case A: vn = n
        if (t != n-1) {
            // Case A.1: t != n-1
            // Function FindPosition(v)
            if (t == 2 && v.swapSymbol(t).isIdentity()) {
                // Rule (1.1)
                return v.swapSymbol(t-1);
            } else if (v.at(n-2) == t || v.at(n-2) == n-1) {
                // Rule (1.2)
                int j = v.findRightmostWrongPosition();
                return v.swap(j);
            } else {
                // Rule (1.3)
                return v.swapSymbol(t);
            }
        } else {
            // Case A.2: t = n-1
            // Rule (2)
            return v.swapSymbol(v.at(n-2));
        }
    } else if (v.at(n-1) == n-1 && v.at(n-2) == n && !v.swapSymbol(n).isIdentity()) {
        // Case B.2: vn = n-1, vn-1 = n, Swap(v,n) != 1n
        if (t == 1) {
            // Rule (3)
            return v.swapSymbol(n);
        } else {
            // Rule (4)
            return v.swapSymbol(t-1);
        }
    } else if (v.at(n-1) == t) {
        // Case B.1.1 or C.1: vn = t
        // Rule (5)
        return v.swapSymbol(n);
    } else {
        // Case B.1.2 or C.2: vn != t
        // Rule (6)
        return v.swapSymbol(t);
    }
}

struct ParentEdge {
    std::string vertex;
    int treeId;
    int lastSymbol;
    std::string rule;
    std::string parent;
};

std::string toString(const ParentEdge& edge) {
    std::stringstream ss;
    ss << edge.vertex << "," << edge.treeId << "," << edge.lastSymbol << "," 
       << edge.rule << "," << edge.parent;
    return ss.str();
}

std::string determineRule(const Permutation& v, int t, int n) {
    int vn = v.at(n-1);
    
    if (vn == n) {
        if (t != n - 1) {
            if (t == 2 && v.swapSymbol(t).isIdentity()) {
                return "(1.1)";
            } else if (v.at(n-2) == t || v.at(n-2) == n-1) {
                return "(1.2)";
            } else {
                return "(1.3)";
            }
        } else {
            return "(2)";
        }
    } else if (vn == n-1 && v.at(n-2) == n && !v.swapSymbol(n).isIdentity()) {
        if (t == 1) {
            return "(3)";
        } else {
            return "(4)";
        }
    } else if (vn == t) {
        return "(5)";
    } else {
        return "(6)";
    }
    return "(unknown)";
}

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n = 4;
    bool storeFiles = true;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-n" && i+1 < argc) {
            n = std::stoi(argv[i+1]);
            i++;
            if (n == 10 || n == 11 || n == 12) {
                storeFiles = false;
            }
        }
    }
    
    // Set up OpenMP threads based on node capabilities
    int max_threads = omp_get_max_threads();
    
    // Use all available cores on each machine for optimal performance
    // You mentioned 1 machine with 4 cores and 2 with 8 cores
    omp_set_num_threads(max_threads);
    
    if (rank == 0) {
        std::cout << "Running with " << size << " MPI processes" << std::endl;
        std::cout << "Each process using up to " << max_threads << " OpenMP threads" << std::endl;
        std::cout << "Constructing ISTs for bubble-sort network of dimension " << n << std::endl;
    }
    
    long long totalPerms = factorial(n);
    long long permsPerProcess = totalPerms / size;
    long long startIdx = rank * permsPerProcess;
    long long endIdx = (rank == size - 1) ? totalPerms : startIdx + permsPerProcess;
    
    std::vector<int> perm(n);
    for (int i = 0; i < n; i++) {
        perm[i] = i + 1;
    }
    
    // Skip to starting permutation for this process
    long long currentIdx = 0;
    while (currentIdx < startIdx && std::next_permutation(perm.begin(), perm.end())) {
        currentIdx++;
    }
    
    // Local storage for all permutations this process will handle
    std::vector<std::vector<int>> localPerms;
    
    // First collect all permutations this process will handle
    while (currentIdx < endIdx && std::next_permutation(perm.begin(), perm.end())) {
        Permutation vertex(perm);
        if (!vertex.isIdentity()) {
            localPerms.push_back(perm);
        }
        currentIdx++;
    }
    
    // Thread-local storage preparations
    int numTrees = n - 1;
    struct VertexProcessing {
        std::string vertex;
        int rank;
    };
    
    std::vector<ParentEdge> parentEdges;
    std::vector<VertexProcessing> myVertexProcessings;
    
    // Reserve space to avoid reallocation during parallel execution
    int estimated_size = localPerms.size() * numTrees;
    parentEdges.reserve(estimated_size);
    myVertexProcessings.reserve(localPerms.size());
    
    double constructionStartTime = MPI_Wtime();
    double parentComputationTime = 0.0;
    
    // Process permutations in parallel using OpenMP
    #pragma omp parallel
    {
        // Thread-local vectors to avoid contention
        std::vector<ParentEdge> thread_edges;
        std::vector<VertexProcessing> thread_vertices;
        
        // Pre-allocate to avoid reallocations
        thread_edges.reserve(estimated_size / omp_get_num_threads());
        thread_vertices.reserve(localPerms.size() / omp_get_num_threads());
        
        // Parallel loop over all permutations assigned to this MPI process
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < localPerms.size(); i++) {
            Permutation vertex(localPerms[i]);
            std::string vertexStr = vertex.toString();
            
            VertexProcessing proc;
            proc.vertex = vertexStr;
            proc.rank = rank;
            thread_vertices.push_back(proc);
            
            // For each tree, compute the parent
            for (int t = 1; t <= numTrees; t++) {
                double parentStart = omp_get_wtime();
                Permutation parent = findParent(vertex, t, n);
                double elapsed = omp_get_wtime() - parentStart;
                
                #pragma omp atomic
                parentComputationTime += elapsed;
                
                std::string parentStr = parent.toString();
                std::string rule = determineRule(vertex, t, n);
                
                ParentEdge edge;
                edge.vertex = vertexStr;
                edge.treeId = t;
                edge.lastSymbol = vertex.at(n-1);
                edge.rule = rule;
                edge.parent = parentStr;
                
                thread_edges.push_back(edge);
            }
        }
        
        // Merge thread-local results into shared vectors
        #pragma omp critical
        {
            parentEdges.insert(parentEdges.end(), thread_edges.begin(), thread_edges.end());
            myVertexProcessings.insert(myVertexProcessings.end(), thread_vertices.begin(), thread_vertices.end());
        }
    }
    
    double constructionTime = MPI_Wtime() - constructionStartTime;
    
    if (rank == 0) {
        std::vector<ParentEdge> allEdges = parentEdges;
        std::vector<VertexProcessing> allVertexProcessings = myVertexProcessings;
        
        // Collect results from all other MPI processes
        for (int src = 1; src < size; src++) {
            int edgeCount;
            MPI_Recv(&edgeCount, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i < edgeCount; i++) {
                char buffer[1024];
                MPI_Recv(buffer, sizeof(buffer), MPI_CHAR, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::string edgeStr(buffer);
                std::stringstream ss(edgeStr);
                std::string token;
                
                ParentEdge edge;
                std::getline(ss, token, ','); edge.vertex = token;
                std::getline(ss, token, ','); edge.treeId = std::stoi(token);
                std::getline(ss, token, ','); edge.lastSymbol = std::stoi(token);
                std::getline(ss, token, ','); edge.rule = token;
                std::getline(ss, token, ','); edge.parent = token;
                
                allEdges.push_back(edge);
            }
            
            int vertexCount;
            MPI_Recv(&vertexCount, 1, MPI_INT, src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i < vertexCount; i++) {
                char buffer[1024];
                MPI_Recv(buffer, sizeof(buffer), MPI_CHAR, src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                VertexProcessing vp;
                vp.vertex = std::string(buffer);
                vp.rank = src;
                allVertexProcessings.push_back(vp);
            }
        }
        
        if (storeFiles) {
            std::ofstream tableFile("bubble_sort_ist_table_n" + std::to_string(n) + ".txt");
            std::ofstream vertexFile("bubble_sort_ist_vertices_n" + std::to_string(n) + ".txt");
            std::ofstream logFile("bubble_sort_ist_log_n" + std::to_string(n) + ".txt");
            
            if (tableFile.is_open() && vertexFile.is_open() && logFile.is_open()) {
                // Write performance data to log file
                logFile << "Performance Statistics:" << std::endl;
                logFile << "======================" << std::endl;
                logFile << "Network dimension (n): " << n << std::endl;
                logFile << "Total vertices: " << factorial(n) << std::endl;
                logFile << "Number of ISTs constructed: " << numTrees << std::endl;
                logFile << "MPI processes: " << size << std::endl;
                logFile << "OpenMP threads per process: " << omp_get_max_threads() << std::endl;
                logFile << "Total construction time: " << constructionTime << " seconds" << std::endl;
                logFile << "Parent computation time: " << parentComputationTime << " seconds" << std::endl;
                
                // Write vertex-to-process mapping to vertex file
                vertexFile << "Vertex-to-Process Mapping:" << std::endl;
                vertexFile << "============================" << std::endl;
                vertexFile << std::setw(15) << "Vertex" << " | " << std::setw(10) << "Process" << std::endl;
                vertexFile << "----------------------------" << std::endl;
                
                std::sort(allVertexProcessings.begin(), allVertexProcessings.end(), 
                    [](const VertexProcessing& a, const VertexProcessing& b) {
                        return a.vertex < b.vertex;
                    });
                    
                for (const auto& vp : allVertexProcessings) {
                    vertexFile << std::setw(15) << vp.vertex << " | " << std::setw(10) << vp.rank << std::endl;
                }
                vertexFile << "============================" << std::endl;
                
                // Format and print the table to table file
                formatAndPrintTable(allEdges, n, numTrees, tableFile);
                
                tableFile.close();
                vertexFile.close();
                logFile.close();
            }
        }
        
        // Print summary statistics
        std::cout << "Total time to construct " << numTrees << " ISTs for n = " << n << ": " << constructionTime << " seconds" << std::endl;
        std::cout << "Average time per vertex-tree pair: " << (parentComputationTime / (factorial(n) * numTrees)) << " seconds" << std::endl;
    } else {
        int edgeCount = parentEdges.size();
        MPI_Send(&edgeCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
        for (const auto& edge : parentEdges) {
            std::string edgeStr = toString(edge);
            MPI_Send(edgeStr.c_str(), edgeStr.length()+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
        
        int vertexCount = myVertexProcessings.size();
        MPI_Send(&vertexCount, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        
        for (const auto& vp : myVertexProcessings) {
            MPI_Send(vp.vertex.c_str(), vp.vertex.length()+1, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
        }
    }
    
    MPI_Finalize();
    return 0;
}