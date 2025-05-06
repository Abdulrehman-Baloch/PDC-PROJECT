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

Permutation findParent(const Permutation& v, int t, int n) {
    // Handle edge cases
    if (v.isIdentity() || t < 1 || t >= n) {
        return v;
    }
    
    // Implementation of Algorithm 1 with corrected rule (1.2)
    if (v.at(n-1) == n) {
        // Case A: vn = n
        if (t != n-1) {
            // Case A.1: t != n-1
            if (t == 2 && v.swapSymbol(t).isIdentity()) {
                // Case A.1.2: t = 2 and Swap(v,t) = 1n
                return v.swapSymbol(t-1);
            } else if (v.at(n-2) == t || v.at(n-2) == n-1) {
                // Case A.1.1.1: vn-1 ∈ {t, n-1}
                if (v.at(n-2) == t) {
                    // Focus on moving symbol t toward its correct position (t-1)
                    int pos = v.findPosition(t);
                    if (pos == t-1) {
                        // t is in its correct position; find the leftmost inversion
                        for (int i = 0; i < n-1; i++) {
                            if (v.at(i) > v.at(i+1)) {
                                return v.swap(i);
                            }
                        }
                        return v; // Should not reach here
                    } else if (pos > t-1) {
                        // t is to the right of its correct position; swap left
                        return v.swap(pos-1);
                    } else {
                        // t is to the left of its correct position; swap right
                        return v.swap(pos);
                    }
                } else {
                    // v[n-2] == n-1; find the leftmost inversion to progress toward identity
                    for (int i = 0; i < n-1; i++) {
                        if (v.at(i) > v.at(i+1)) {
                            return v.swap(i);
                        }
                    }
                    return v; // Should not reach here
                }
            } else if (v.at(n-1) == t || v.at(n-1) == n-1) {
                // Case A.1.1.2: vn-1 ∈ {t, n-1} (additional condition)
                int j = v.findRightmostWrongPosition();
                return v.swapSymbol(v.at(j));
            } else {
                // Case A.1.1.3: vn-1 ∉ {t, n-1}
                return v.swapSymbol(t);
            }
        } else {
            // Case A.2: t = n-1
            return v.swapSymbol(v.at(n-2));
        }
    } else if (v.at(n-1) == n-1 && v.at(n-2) == n && !v.swapSymbol(n).isIdentity()) {
        // Case B.2: vn = n-1, vn-1 = n, Swap(v,n) != 1n
        if (t == 1) {
            // Case B.2.2: t = 1
            return v.swapSymbol(n);
        } else {
            // Case B.2.1: t != 1
            return v.swapSymbol(t-1);
        }
    } else if (v.at(n-1) == t) {
        // Case B.1.1 or C.1: vn = t
        return v.swapSymbol(n);
    } else {
        // Case B.1.2 or C.2: vn != t
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

void formatAndPrintTable(const std::vector<ParentEdge>& allEdges, int n, int numTrees, std::ofstream& tableFile) {
    long long expectedPerms = factorial(n);
    long long expectedEdges = (expectedPerms - 1) * numTrees;
    
    std::set<std::string> uniqueVertices;
    for (const auto& edge : allEdges) {
        uniqueVertices.insert(edge.vertex);
    }

    std::stringstream header;
    header << "\nTable 1: The parent of every vertex v ∈ V(B" << n << ") in T" << n 
           << "t for t ∈ {1, 2, ..., " << numTrees << "}" << std::endl;
    header << "Total vertices (excluding root): " << uniqueVertices.size() << " (expected: " 
           << expectedPerms - 1 << ")" << std::endl;
    header << "Total edges: " << allEdges.size() << " (expected: " << expectedEdges << ")" << std::endl;

    if (uniqueVertices.size() != expectedPerms - 1 || allEdges.size() != expectedEdges) {
        header << "WARNING: Missing vertices or edges in the table!" << std::endl;
    }
    
    tableFile << header.str();
    
    std::vector<ParentEdge> sortedEdges = allEdges;
    std::sort(sortedEdges.begin(), sortedEdges.end(), 
        [](const ParentEdge& a, const ParentEdge& b) {
            if (a.vertex != b.vertex) return a.vertex < b.vertex;
            return a.treeId < b.treeId;
        });
    
    std::string tableDivider = "------------------------------------------------------------------";
    
    tableFile << tableDivider << std::endl;
    tableFile << std::setw(15) << "v" << " | " 
              << std::setw(3) << "t" << " | " 
              << std::setw(3) << "vn" << " | " 
              << std::setw(6) << "Rule" << " | " 
              << std::setw(15) << "p" << " || " 
              << std::setw(15) << "v" << " | " 
              << std::setw(3) << "t" << " | " 
              << std::setw(3) << "vn" << " | " 
              << std::setw(6) << "Rule" << " | " 
              << std::setw(15) << "p" << std::endl;
    tableFile << tableDivider << std::endl;
    
    size_t numRows = (sortedEdges.size() + 1) / 2;
    
    for (size_t i = 0; i < numRows; i++) {
        std::stringstream row;
        
        const auto& left = sortedEdges[i];
        row << std::setw(15) << left.vertex << " | " 
            << std::setw(3) << left.treeId << " | " 
            << std::setw(3) << left.lastSymbol << " | " 
            << std::setw(6) << left.rule << " | " 
            << std::setw(15) << left.parent << " || ";
        
        size_t rightIdx = i + numRows;
        if (rightIdx < sortedEdges.size()) {
            const auto& right = sortedEdges[rightIdx];
            row << std::setw(15) << right.vertex << " | " 
                << std::setw(3) << right.treeId << " | " 
                << std::setw(3) << right.lastSymbol << " | " 
                << std::setw(6) << right.rule << " | " 
                << std::setw(15) << right.parent;
        } else {
            row << std::setw(15) << "-" << " | " 
                << std::setw(3) << "-" << " | " 
                << std::setw(3) << "-" << " | " 
                << std::setw(6) << "-" << " | " 
                << std::setw(15) << "-";
        }
        
        tableFile << row.str() << std::endl;
    }
    
    tableFile << tableDivider << std::endl;
}

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n = 4;
    std::string tableFilename = "bubble_sort_ist_table.txt";
    std::string vertexFilename = "bubble_sort_ist_vertices.txt";
    std::string logFilename = "bubble_sort_ist_log.txt";
    
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-n" && i+1 < argc) {
            n = std::stoi(argv[i+1]);
            i++;
        } else if (std::string(argv[i]) == "-t" && i+1 < argc) {
            tableFilename = argv[i+1];
            i++;
        } else if (std::string(argv[i]) == "-v" && i+1 < argc) {
            vertexFilename = argv[i+1];
            i++;
        } else if (std::string(argv[i]) == "-l" && i+1 < argc) {
            logFilename = argv[i+1];
            i++;
        }
    }
    
    int max_threads = omp_get_max_threads();
    int optimal_threads = std::min(max_threads, 2);
    omp_set_num_threads(optimal_threads);
    
    long long numPermutations = factorial(n);
    double estimatedMemoryMB = numPermutations * n * sizeof(int) / (1024.0 * 1024.0) / size;
    
    std::ofstream tableFile, vertexFile, logFile;
    if (rank == 0) {
        tableFile.open(tableFilename);
        if (!tableFile.is_open()) {
            std::cerr << "ERROR: Could not open table file: " << tableFilename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        vertexFile.open(vertexFilename);
        if (!vertexFile.is_open()) {
            std::cerr << "ERROR: Could not open vertex file: " << vertexFilename << std::endl;
            tableFile.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        logFile.open(logFilename);
        if (!logFile.is_open()) {
            std::cerr << "ERROR: Could not open log file: " << logFilename << std::endl;
            tableFile.close();
            vertexFile.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        std::cout << "Parallel IST Construction in Bubble-Sort Networks" << std::endl;
        std::cout << "Dimension: " << n << std::endl;
        std::cout << "MPI processes: " << size << std::endl;
        std::cout << "OpenMP threads per process: " << optimal_threads << std::endl;
        std::cout << "Expected permutations: " << numPermutations << std::endl;
        std::cout << "Estimated memory usage per process: " << estimatedMemoryMB << " MB" << std::endl;
        std::cout << "Table output file: " << tableFilename << std::endl;
        std::cout << "Vertex mapping file: " << vertexFilename << std::endl;
        std::cout << "Log file: " << logFilename << std::endl;
        
        if (n > 8 || estimatedMemoryMB * size > 1000) {
            std::cout << "WARNING: Large dimension selected." << std::endl;
        }
    }
    
    double startTime = MPI_Wtime();
    
    long long totalPerms = factorial(n);
    long long permsPerProcess = totalPerms / size;
    long long startIdx = rank * permsPerProcess;
    long long endIdx = (rank == size - 1) ? totalPerms : startIdx + permsPerProcess;
    
    std::vector<int> perm(n);
    for (int i = 0; i < n; i++) {
        perm[i] = i + 1;
    }
    
    long long currentIdx = 0;
    std::vector<Permutation> myVertices;
    struct VertexProcessing {
        std::string vertex;
        int rank;
    };
    std::vector<VertexProcessing> myVertexProcessings;
    std::vector<ParentEdge> parentEdges;
    int numTrees = n - 1;
    
    double constructionStartTime = MPI_Wtime();
    double parentComputationTime = 0.0;
    
    if (rank == 0) {
        std::string constructMsg = "Constructing " + std::to_string(numTrees) + " independent spanning trees...";
        std::cout << constructMsg << std::endl;
        if (logFile.is_open()) logFile << constructMsg << std::endl;
    }
    
    do {
        if (currentIdx >= startIdx && currentIdx < endIdx) {
            Permutation vertex(perm);
            if (!vertex.isIdentity()) {
                myVertices.push_back(vertex);
                std::string vertexStr = vertex.toString();
                
                VertexProcessing proc;
                proc.vertex = vertexStr;
                proc.rank = rank;
                myVertexProcessings.push_back(proc);
                
                for (int t = 1; t <= numTrees; t++) {
                    double parentStart = MPI_Wtime();
                    Permutation parent = findParent(vertex, t, n);
                    std::string rule = determineRule(vertex, t, n);
                    parentComputationTime += MPI_Wtime() - parentStart;
                    
                    std::string parentStr = parent.toString();
                    
                    ParentEdge edge;
                    edge.vertex = vertexStr;
                    edge.treeId = t;
                    edge.lastSymbol = vertex.at(n-1);
                    edge.rule = rule;
                    edge.parent = parentStr;
                    
                    parentEdges.push_back(edge);
                }
            }
        }
        currentIdx++;
    } while (currentIdx < endIdx && std::next_permutation(perm.begin(), perm.end()));
    
    double constructionTime = MPI_Wtime() - constructionStartTime;
    
    std::stringstream rankStats;
    rankStats << "Process " << rank << " computed parents for " << myVertices.size() 
              << " vertices in " << parentComputationTime << " seconds (total processing: " 
              << constructionTime << " seconds)";
    std::cout << rankStats.str() << std::endl;
    
    if (rank != 0) {
        std::string statsStr = rankStats.str();
        MPI_Send(statsStr.c_str(), statsStr.length()+1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
        // Send parentComputationTime to rank 0
        MPI_Send(&parentComputationTime, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    } else {
        if (logFile.is_open()) logFile << rankStats.str() << std::endl;
        
        for (int src = 1; src < size; src++) {
            char buffer[1024];
            MPI_Recv(buffer, sizeof(buffer), MPI_CHAR, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (logFile.is_open()) logFile << buffer << std::endl;
        }
    }
    
    if (rank == 0) {
        std::vector<ParentEdge> allEdges = parentEdges;
        std::vector<VertexProcessing> allVertexProcessings = myVertexProcessings;
        double totalParentComputationTime = parentComputationTime;
        
        // Receive parentComputationTime from other processes
        for (int src = 1; src < size; src++) {
            double receivedTime;
            MPI_Recv(&receivedTime, 1, MPI_DOUBLE, src, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totalParentComputationTime += receivedTime;
        }
        
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
        
        formatAndPrintTable(allEdges, n, numTrees, tableFile);
        
        // Modified output to show total parent computation time
        std::cout << "\nPerformance Summary:" << std::endl;
        std::cout << "Total parent computation time: " << totalParentComputationTime << " seconds" << std::endl;
        std::cout << "Table results saved to: " << tableFilename << std::endl;
        std::cout << "Vertex mapping saved to: " << vertexFilename << std::endl;
        std::cout << "Log saved to: " << logFilename << std::endl;
        
        tableFile.close();
        vertexFile.close();
        logFile.close();
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
        
        // Send parentComputationTime to rank 0
        MPI_Send(&parentComputationTime, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}