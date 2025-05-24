#pragma once

#include <vector>
#include <ostream>
#include <iomanip>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include "solver.hpp" // For KokkosVector

// Псевдонимы для типов Kokkos
using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::HostSpace;
using KokkosVector = Kokkos::View<double*, memory_space>;
using KokkosCrsMatrix = KokkosSparse::CrsMatrix<double, int, execution_space, void, int>;

class GridSystem
{
private:
    int n, m;                     // Размерность сетки
    double a, b, c, d;            // Границы области
    double x_step, y_step;        // Шаги сетки
    double A, x_k, y_k;           // Коэффициенты разностной схемы
    
    // Разреженная матрица системы
    KokkosCrsMatrix matrix;
    
    // Вектор правой части
    KokkosVector rhs;
    
    // Векторы координат узлов
    std::vector<double> node_x_coords;
    std::vector<double> node_y_coords;
    
    // Вспомогательные структуры для сборки разреженной матрицы
    std::vector<int> row_map;           // Смещения по строкам (CSR формат)
    std::vector<int> entries;           // Индексы столбцов ненулевых элементов
    std::vector<double> values;         // Значения ненулевых элементов
    std::vector<double> rhs_values;     // Значения правой части
    
    // Вспомогательные функции
    double function(double x, double y);
    double solution(double x, double y);
    bool is_left_boundary(int x, int y) const;
    bool is_right_boundary(int x, int y) const;
    bool is_top_boundary(int x, int y) const;
    bool is_bottom_boundary(int x, int y) const;
    bool is_boundary(int x, int y) const;
    double calculate_value(int x, int y, double x_k, double y_k);
    double calculate_x(int x) const;
    double calculate_y(int y) const;
    
    // Новая реализация сборки матрицы для разреженного формата
    void initiate_matrix();
    int calculate_position_in_template(int x, int y) const;
    int calculate_position_in_upper_area(int x, int y) const;
    int calculate_position_in_bottom_edge(int x, int y) const;
    
    // Методы для добавления элементов в разреженную матрицу
    void add_matrix_entry(int row, int col, double value);
    void finalize_matrix();

public:
    // Структура для хранения координат точки
    struct NodeCoordinates {
        double x;
        double y;
    };
    
    // Конструктор и деструктор
    GridSystem(int m, int n, double a, double b, double c, double d);
    ~GridSystem();
    
    // Получить матрицу и вектор правой части
    const KokkosCrsMatrix& get_matrix() const { return matrix; }
    const KokkosVector& get_rhs() const { return rhs; }
    KokkosVector get_true_solution_vector(); // New method declaration
    
    // Методы для доступа к векторам координат
    const std::vector<double>& get_x_coords() const { return node_x_coords; }
    const std::vector<double>& get_y_coords() const { return node_y_coords; }
    
    // Метод для получения координат точки по индексу в массиве решения
    NodeCoordinates get_node_coordinates(int solution_index) const;
    
    // Вывод информации
    friend std::ostream &operator<<(std::ostream &os, const GridSystem &grid);
};