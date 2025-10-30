import csv
import math
import argparse
from typing import List, Dict

class Instance:
    """
    Przechowuje dane instancji, w tym listę węzłów i macierz odległości.
    """
    def __init__(self, nodes: List[Dict]):
        self.nodes = nodes
        self.n = len(nodes)
        self.dist_mat = self._compute_distance_matrix()

    def _compute_distance_matrix(self) -> List[List[int]]:
        """
        Oblicza macierz odległości euklidesowych, zaokrąglonych do najbliższej
        liczby całkowitej, zgodnie z definicją problemu.
        """
        matrix = [[0] * self.n for _ in range(self.n)]
        for i in range(self.n):
            for j in range(i + 1, self.n):
                node_a = self.nodes[i]
                node_b = self.nodes[j]
                
                dx = node_a['x'] - node_b['x']
                dy = node_a['y'] - node_b['y']
                
                # Oblicz odległość euklidesową i zaokrąglij matematycznie
                dist = math.sqrt(dx*dx + dy*dy)
                rounded_dist = int(round(dist))
                
                matrix[i][j] = rounded_dist
                matrix[j][i] = rounded_dist
        return matrix

def read_csv_instance(filepath: str) -> Instance:
    """
    Wczytuje instancję z pliku CSV (delimiter ';').
    Format: x;y;cost
    """
    nodes = []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=';')
            for i, row in enumerate(reader):
                if len(row) != 3:
                    print(f"Warning: Skipping malformed row {i+1}: {row}")
                    continue
                nodes.append({
                    'id': i,
                    'x': int(row[0]),
                    'y': int(row[1]),
                    'cost': int(row[2])
                })
        if not nodes:
            raise ValueError("No valid nodes loaded from file.")
        return Instance(nodes)
    except FileNotFoundError:
        print(f"Error: File not found at '{filepath}'")
        exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        exit(1)

def evaluate_solution(instance: Instance, path_indices: List[int]) -> int:
    """
    Oblicza całkowitą wartość rozwiązania (koszt ścieżki + koszt węzłów).
    Zakłada, że `path_indices` to cykl (ostatni łączy się z pierwszym).
    """
    total_cost = 0
    total_distance = 0
    
    # 1. Oblicz całkowity koszt wybranych węzłów
    for node_idx in path_indices:
        if node_idx >= instance.n:
            print(f"Error: Node index {node_idx} is out of bounds (max is {instance.n - 1}).")
            exit(1)
        total_cost += instance.nodes[node_idx]['cost']
        
    # 2. Oblicz całkowitą długość ścieżki (cyklu)
    num_nodes_in_path = len(path_indices)
    if num_nodes_in_path < 2:
        return total_cost # Tylko koszt, jeśli ścieżka jest za krótka

    for i in range(num_nodes_in_path):
        node_a_idx = path_indices[i]
        
        # Weź następny węzeł, zawijając do pierwszego na końcu
        node_b_idx = path_indices[(i + 1) % num_nodes_in_path]
        
        total_distance += instance.dist_mat[node_a_idx][node_b_idx]
        
    # 3. Zwróć sumę
    return total_distance + total_cost

def main():
    """
    Główna funkcja do parsowania argumentów i uruchamiania ewaluacji.
    """
    parser = argparse.ArgumentParser(
        description="Calculate the total value of a solution for the TSP-cost problem.",
        epilog="Example: python check_solution.py TSPA.csv '1,5,3,2,8'"
    )
    
    parser.add_argument(
        "csv_file", 
        type=str,
        help="Path to the instance CSV file (e.g., TSPA.csv)"
    )
    parser.add_argument(
        "solution_indices", 
        type=str,
        help="A comma-separated string of node indices in the cycle (e.g., '1,5,3,2,8')"
    )
    
    args = parser.parse_args()
    
    # Przetwórz listę indeksów
    try:
        path_indices = [int(idx.strip()) for idx in args.solution_indices.split(', ')]
        if not path_indices:
            raise ValueError("Solution list cannot be empty.")
    except ValueError as e:
        print(f"Error: Invalid solution format. Must be comma-separated numbers. Details: {e}")
        exit(1)

    # Wczytaj instancję
    print(f"Loading instance from: {args.csv_file}")
    instance = read_csv_instance(args.csv_file)
    print(f"Loaded {instance.n} nodes.")
    print(f"Verifying solution with {len(path_indices)} nodes: {path_indices}")
    
    # Oblicz wartość
    total_value = evaluate_solution(instance, path_indices)
    
    print("\n--- CALCULATION COMPLETE ---")
    print(f"Total calculated value: {total_value}")

if __name__ == "__main__":
    main()