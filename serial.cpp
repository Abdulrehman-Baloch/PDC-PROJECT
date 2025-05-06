#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <map>
#include <set>
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
    // Parse command-line argument for n
    int n = 0;
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <n>" << std::endl;
        std::cerr << "where <n> is the dimension of the bubble-sort network (n >= 2)" << std::endl;
        return 1;
    }

    try {
        n = std::stoi(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Invalid value for n. It must be an integer." << std::endl;
        return 1;
    }

    if (n < 2) {
        std::cerr << "ERROR: n must be at least 2 to construct at least one tree." << std::endl;
        return 1;
    }

    // Generate filenames with n as suffix
    std::string tableFilename = "bubble_sort_ist_table_n" + std::to_string(n) + ".txt";
    std::string vertexFilename = "bubble_sort_ist_vertices_n" + std::to_string(n) + ".txt";
    std::string logFilename = "bubble_sort_ist_log_n" + std::to_string(n) + ".txt";
    
    long long numPermutations = factorial(n);
    double estimatedMemoryMB = numPermutations * n * sizeof(int) / (1024.0 * 1024.0);
    
    std::ofstream tableFile, vertexFile, logFile;
    
    tableFile.open(tableFilename);
    if (!tableFile.is_open()) {
        std::cerr << "ERROR: Could not open table file: " << tableFilename << std::endl;
        return 1;
    }
    
    vertexFile.open(vertexFilename);
    if (!vertexFile.is_open()) {
        std::cerr << "ERROR: Could not open vertex file: " << vertexFilename << std::endl;
        tableFile.close();
        return 1;
    }
    
    logFile.open(logFilename);
    if (!logFile.is_open()) {
        std::cerr << "ERROR: Could not open log file: " << logFilename << std::endl;
        tableFile.close();
        vertexFile.close();
        return 1;
    }
    
    std::cout << "Serial IST Construction in Bubble-Sort Networks" << std::endl;
    std::cout << "Dimension: " << n << std::endl;
    std::cout << "Expected permutations: " << numPermutations << std::endl;
    std::cout << "Estimated memory usage: " << estimatedMemoryMB << " MB" << std::endl;
    std::cout << "Table output file: " << tableFilename << std::endl;
    std::cout << "Vertex mapping file: " << vertexFilename << std::endl;
    std::cout << "Log file: " << logFilename << std::endl;
    
    if (n > 8 || estimatedMemoryMB > 1000) {
        std::cout << "WARNING: Large dimension selected. This may consume significant time and memory." << std::endl;
    }
    
    double startTime = static_cast<double>(clock()) / CLOCKS_PER_SEC;
    
    std::vector<int> perm(n);
    for (int i = 0; i < n; i++) {
        perm[i] = i + 1;
    }
    
    std::vector<Permutation> myVertices;
    struct VertexProcessing {
        std::string vertex;
    };
    std::vector<VertexProcessing> myVertexProcessings;
    std::vector<ParentEdge> parentEdges;
    int numTrees = n - 1;
    
    double constructionStartTime = static_cast<double>(clock()) / CLOCKS_PER_SEC;
    double parentComputationTime = 0.0;
    
    std::string constructMsg = "Constructing " + std::to_string(numTrees) + " independent spanning trees...";
    std::cout << constructMsg << std::endl;
    logFile << constructMsg << std::endl;
    
    // Generate all permutations and process them serially
    do {
        Permutation vertex(perm);
        if (!vertex.isIdentity()) {
            myVertices.push_back(vertex);
            std::string vertexStr = vertex.toString();
            
            VertexProcessing proc;
            proc.vertex = vertexStr;
            myVertexProcessings.push_back(proc);
            
            for (int t = 1; t <= numTrees; t++) {
                double parentStart = static_cast<double>(clock()) / CLOCKS_PER_SEC;
                Permutation parent = findParent(vertex, t, n);
                std::string rule = determineRule(vertex, t, n);
                parentComputationTime += (static_cast<double>(clock()) / CLOCKS_PER_SEC) - parentStart;
                
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
    } while (std::next_permutation(perm.begin(), perm.end()));
    
    double constructionTime = (static_cast<double>(clock()) / CLOCKS_PER_SEC) - constructionStartTime;
    
    std::stringstream stats;
    stats << "Computed parents for " << myVertices.size() 
          << " vertices in " << parentComputationTime << " seconds (total processing: " 
          << constructionTime << " seconds)";
    std::cout << stats.str() << std::endl;
    logFile << stats.str() << std::endl;
    
    // Write vertex mapping (no process info in serial version)
    vertexFile << "Vertex List:" << std::endl;
    vertexFile << "============================" << std::endl;
    vertexFile << std::setw(15) << "Vertex" << std::endl;
    vertexFile << "----------------------------" << std::endl;
    
    std::sort(myVertexProcessings.begin(), myVertexProcessings.end(), 
        [](const VertexProcessing& a, const VertexProcessing& b) {
            return a.vertex < b.vertex;
        });
        
    for (const auto& vp : myVertexProcessings) {
        vertexFile << std::setw(15) << vp.vertex << std::endl;
    }
    vertexFile << "============================" << std::endl;
    
    formatAndPrintTable(parentEdges, n, numTrees, tableFile);
    
    double totalTime = (static_cast<double>(clock()) / CLOCKS_PER_SEC) - constructionStartTime;
    
    std::cout << "\nPerformance Summary:" << std::endl;
    std::cout << "Total execution time: " << totalTime << " seconds" << std::endl;
    std::cout << "Table results saved to: " << tableFilename << std::endl;
    std::cout << "Vertex mapping saved to: " << vertexFilename << std::endl;
    std::cout << "Log saved to: " << logFilename << std::endl;
    
    tableFile.close();
    vertexFile.close();
    logFile.close();
    
    std::cout << "Completed dimension n = " << n << std::endl;
    
    return 0;
}