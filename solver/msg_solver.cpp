#include "msg_solver.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <chrono>

// Вспомогательная функция для вычисления скалярного произведения
double MSGSolver::dot(const KokkosVector& v1, const KokkosVector& v2) const {
    double result = 0.0;
    
    auto v1_host = Kokkos::create_mirror_view(v1);
    auto v2_host = Kokkos::create_mirror_view(v2);
    
    Kokkos::deep_copy(v1_host, v1);
    Kokkos::deep_copy(v2_host, v2);
    
    for (int i = 0; i < v1.extent(0); ++i) {
        result += v1_host(i) * v2_host(i);
    }
    
    return result;
}

// Вспомогательная функция для умножения матрицы на вектор A*v
KokkosVector MSGSolver::multiply(const KokkosCrsMatrix& A, const KokkosVector& v) const {
    int n = v.extent(0);
    KokkosVector result("A*v", n);
    
    KokkosSparse::spmv("N", 1.0, A, v, 0.0, result);
    
    return result;
}

// Вычисление нормы вектора (Евклидова норма)
double MSGSolver::norm(const KokkosVector& v) const {
    return std::sqrt(dot(v, v));
}

// Вычисление максимальной нормы вектора (максимальный модуль элемента)
double MSGSolver::max_norm(const KokkosVector& v) const {
    double max_val = 0.0;
    Kokkos::parallel_reduce(v.extent(0), KOKKOS_LAMBDA(const int i, double& lmax) {
        double val = Kokkos::abs(v(i));
        if (val > lmax) {
            lmax = val;
        }
    }, Kokkos::Max<double>(max_val));
    return max_val;
}

// Implementation of the new utility function for max norm calculation
double calculate_max_norm(const KokkosVector& v) {
    double max_val = 0.0;
    Kokkos::parallel_reduce(v.extent(0), KOKKOS_LAMBDA(const int i, double& lmax) {
        double val = Kokkos::abs(v(i));
        if (val > lmax) {
            lmax = val;
        }
    }, Kokkos::Max<double>(max_val));
    return max_val;
}

// Функция для проверки условий завершения итераций
bool MSGSolver::checkTerminationConditions(double precision_max_norm, double r_max_norm,
                                          double error_max_norm, int iterationsDone,
                                          StopCriterion& reason) {
    // Проверка прерывания пользователем (имеет наивысший приоритет)
    if (stop_requested) {
        reason = StopCriterion::INTERRUPTED;
        return true;
    }
    
    // Проверка критериев остановки в порядке приоритета
    
    // 1. Проверка по разности последовательных приближений
    if (use_precision && !std::isnan(precision_max_norm) && precision_max_norm < eps_precision) {
        reason = StopCriterion::PRECISION;
        return true;
    }
    
    // 2. Проверка по норме невязки
    if (use_residual && !std::isnan(r_max_norm) && r_max_norm < eps_residual) {
        reason = StopCriterion::RESIDUAL;
        return true;
    }
    
    // 3. Проверка по сравнению с истинным решением
    if (use_exact_error && !std::isnan(error_max_norm) && error_max_norm < eps_exact_error) {
        reason = StopCriterion::EXACT_ERROR;
        return true;
    }
    
    // 4. Проверка по максимальному числу итераций
    if (use_max_iterations && iterationsDone >= maxIterations) {
        reason = StopCriterion::ITERATIONS;
        return true;
    }
    
    // Ни один из критериев не выполнен - продолжаем итерации
    return false;
}

// Метод для решения СЛАУ методом простых итераций
KokkosVector MSGSolver::solve() {
    std::cout << "Using tau: " << this->tau << std::endl;
    // Reset stop flag at the beginning of the solve process
    converged = false;
    stop_requested = false;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    // Размерность задачи
    int n = b.extent(0);
    
    // Инициализируем вектор решения нулями
    KokkosVector x("x", n);
    Kokkos::deep_copy(x, 0.0);
    
    KokkosVector x_prev("x_prev", n);
    KokkosVector r_current("r_current", n); // Текущая невязка r_k = b - A*x_k
    KokkosVector Ax_current("Ax_current", n); // Произведение A*x_k

    // Определяем, нужно ли вычислять различные нормы
    bool need_residual_norm_check = use_residual;
    bool need_precision_norm_check = use_precision;
    bool has_true_solution = use_exact_error && exact_solution.extent(0) > 0;
    bool need_error_norm_check = has_true_solution;

    // Нормы для критериев остановки (max-norms)
    double r_max_norm = std::numeric_limits<double>::quiet_NaN();
    double initial_r_max_norm = std::numeric_limits<double>::quiet_NaN(); 
    double initial_r_L2_norm = 0.0; 

    // Значения для callback
    double cb_precision_max_norm = std::numeric_limits<double>::quiet_NaN();
    double cb_r_max_norm = std::numeric_limits<double>::quiet_NaN();
    double cb_error_max_norm = std::numeric_limits<double>::quiet_NaN();
    
    // Вычисляем начальную невязку r_0 = b - A*x_0. Так как x_0 = 0, r_0 = b.
    // Для общности, если x_0 не ноль:
    // Ax_current = multiply(a, x);
    // Kokkos::parallel_for(Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
    //    r_current(i) = b(i) - Ax_current(i);
    // });
    Kokkos::deep_copy(r_current, b); // Если x_0 = 0

    if (need_residual_norm_check) {
        r_max_norm = max_norm(r_current);
        initial_r_max_norm = r_max_norm;
        cb_r_max_norm = r_max_norm;
    }
    initial_r_L2_norm = norm(r_current); 
    
    int iterationsDone = 0;
    stop_reason = StopCriterion::ITERATIONS; // По умолчанию
    double precision_max_norm = std::numeric_limits<double>::quiet_NaN();
    double error_max_norm = std::numeric_limits<double>::quiet_NaN();
    
    if (need_error_norm_check) {
        KokkosVector error_vec("error_vec", n);
        Kokkos::parallel_for(Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
            error_vec(i) = x(i) - exact_solution(i);
        });
        error_max_norm = max_norm(error_vec);
        cb_error_max_norm = error_max_norm;
    }
    
    if (iteration_callback) {
        iteration_callback(0, cb_precision_max_norm, cb_r_max_norm, cb_error_max_norm);
    }
    
    if (checkTerminationConditions(precision_max_norm, r_max_norm, error_max_norm, iterationsDone, stop_reason)) {
        converged = (stop_reason != StopCriterion::INTERRUPTED && stop_reason != StopCriterion::ITERATIONS && stop_reason != StopCriterion::NUMERICAL_ERROR);
    } else {
        // Основной цикл итераций МПИ
        while (true) {
            Kokkos::deep_copy(x_prev, x); // x_prev = x_k
            
            // r_current это r_k = b - A*x_k (вычислено на предыдущем шаге или инициализации)
            // Обновление решения: x_{k+1} = x_k + tau * r_k
            Kokkos::parallel_for(Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
                x(i) = x_prev(i) - tau * r_current(i);
            });
            
            iterationsDone++;
            
            // Вычисляем новую невязку r_{k+1} = b - A*x_{k+1}
            Ax_current = multiply(a, x);
            Kokkos::parallel_for(Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
               r_current(i) = b(i) - Ax_current(i);
            });

            // Сбрасываем значения callback перед пересчетом
            cb_precision_max_norm = std::numeric_limits<double>::quiet_NaN();
            cb_r_max_norm = std::numeric_limits<double>::quiet_NaN();
            cb_error_max_norm = std::numeric_limits<double>::quiet_NaN();
            
            if (need_residual_norm_check) {
                r_max_norm = max_norm(r_current);
                cb_r_max_norm = r_max_norm;
            }
            
            if (need_precision_norm_check) {
                KokkosVector diff_x("diff_x", n);
                Kokkos::parallel_for(Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
                    diff_x(i) = x(i) - x_prev(i); // x_{k+1} - x_k
                });
                precision_max_norm = max_norm(diff_x);
                cb_precision_max_norm = precision_max_norm;
            }
            
            if (need_error_norm_check) {
                KokkosVector error_vec("error_vec_loop", n); 
                Kokkos::parallel_for(Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int i) {
                    error_vec(i) = x(i) - exact_solution(i);
                });
                error_max_norm = max_norm(error_vec);
                cb_error_max_norm = error_max_norm;
            }
            
            if (checkTerminationConditions(precision_max_norm, r_max_norm, error_max_norm, iterationsDone, stop_reason)) {
                converged = (stop_reason != StopCriterion::ITERATIONS && stop_reason != StopCriterion::INTERRUPTED && stop_reason != StopCriterion::NUMERICAL_ERROR);
                break;
            }
            
            // Вывод информации в консоль и вызов колбэка
            if (iteration_callback || (iterationsDone % 100 == 0 || iterationsDone == 1)) {
                 if (iterationsDone % 100 == 0 || iterationsDone == 1) { // Console output condition
                    std::cout << "Итерация: " << iterationsDone << "\n";
                    if (need_precision_norm_check && !std::isnan(precision_max_norm)) {
                        std::cout << "Точность ||x(n)-x(n-1)||: max-норма = " << std::scientific << precision_max_norm << "\n";
                    }
                    if (need_residual_norm_check && !std::isnan(r_max_norm)) {
                        std::cout << "Невязка ||Ax-b||: max-норма = " << std::scientific << r_max_norm << "\n";
                    }
                    if (need_error_norm_check && !std::isnan(error_max_norm)) {
                        std::cout << "Ошибка ||u-x||: max-норма = " << std::scientific << error_max_norm << "\n";
                    } else if (need_error_norm_check) {
                        std::cout << "Ошибка ||u-x||: не вычисляется (нет точного решения или критерий отключен)\n";
                    }
                    std::cout << "\n";
                }
                if (iteration_callback) { 
                    iteration_callback(iterationsDone, cb_precision_max_norm, cb_r_max_norm, cb_error_max_norm);
                }
            }
        }
    }
    
    iterations = iterationsDone;
    final_residual_norm = need_residual_norm_check ? r_max_norm : std::numeric_limits<double>::quiet_NaN();
    final_precision = need_precision_norm_check ? precision_max_norm : std::numeric_limits<double>::quiet_NaN();
    final_error_norm = need_error_norm_check ? error_max_norm : std::numeric_limits<double>::quiet_NaN();
    
    if (iteration_callback) {
        iteration_callback(iterationsDone, 
                           final_precision, 
                           final_residual_norm, 
                           final_error_norm);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    
    std::cout << "Метод простых итераций (МПИ)\n" // Изменено название метода
              << "Итераций: " << iterationsDone << "\n"
              << "Время: " << duration << " мс\n";
    
    if (need_residual_norm_check) { 
        std::cout << "Начальная невязка (L2): " << initial_r_L2_norm << "\n"; 
        if (!std::isnan(final_residual_norm)) {
             std::cout << "Конечная невязка (max-norm): " << final_residual_norm << "\n";
        } else {
             std::cout << "Конечная невязка (max-norm): не вычислялась\n";
        }
    }
    
    std::cout << "Сходимость: " << (converged ? "Да" : "Нет") << "\n"
              << "Причина остановки: " << getStopReasonText() << std::endl;
    
    if (need_error_norm_check && !std::isnan(final_error_norm)) { 
        std::cout << "Ошибка сравнения с точным решением (max-norm): " << final_error_norm << std::endl;
    }
    
    return x;
}

// Генерация отчета о решении
std::string MSGSolver::generateReport(int n, int m, double a, double b, double c, double d) const {
    std::stringstream ss;
    
    ss << "ОТЧЕТ О РЕШЕНИИ ЗАДАЧИ ДИРИХЛЕ\n";
    ss << "===========================\n\n";
    
    // Параметры задачи
    ss << "ПАРАМЕТРЫ ЗАДАЧИ:\n";
    ss << "----------------\n";
    ss << "Размер сетки: " << n << "x" << m << " внутренних узлов\n";
    ss << "Область: [" << a << ", " << b << "] x [" << c << ", " << d << "]\n";
    ss << "Шаг по x: " << (b - a) / (n + 1) << "\n";
    ss << "Шаг по y: " << (d - c) / (m + 1) << "\n";
    ss << "Общее количество неизвестных: " << n * m << "\n\n";
    
    // Информация о методе решения
    ss << "МЕТОД РЕШЕНИЯ:\n";
    ss << "-------------\n";
    ss << "Название метода: " << name << "\n";
    
    // Выводим информацию только о включенных критериях остановки
    ss << "Критерии остановки:\n";
    
    if (use_max_iterations) {
        ss << "  - Максимальное число итераций: " << maxIterations << "\n";
    }
    
    if (use_precision) {
        ss << "  - Точность ||xn-x(n-1)||: " << eps_precision << "\n";
    }
    
    if (use_residual) {
        ss << "  - Норма невязки ||Ax-b||: " << eps_residual << "\n";
    }
    
    if (use_exact_error) {
        ss << "  - Норма ошибки ||u-x||: " << eps_exact_error << "\n";
    }
    
    ss << "\n";
    
    // Результаты решения
    ss << "РЕЗУЛЬТАТЫ РЕШЕНИЯ:\n";
    ss << "-----------------\n";
    ss << "Выполнено итераций: " << iterations << "\n";
    ss << "Сходимость: " << (converged ? "Да" : "Нет") << "\n";
    ss << "Причина остановки: " << getStopReasonText() << "\n";
    ss << "Достигнутые величины:\n";
    
    // Выводим результаты только для включенных критериев
    if (use_precision) {
        ss << "  - Точность ||xn-x(n-1)||: " << std::scientific << final_precision << "\n";
    }
    
    if (use_residual) {
        ss << "  - Норма невязки ||Ax-b||: " << std::scientific << final_residual_norm << "\n";
    }
    
    if (use_exact_error) {
        ss << "  - Норма ошибки ||u-x||: " << std::scientific << final_error_norm << "\n";
    }
    
    ss << "\n";
    
    // Примечания
    ss << "ПРИМЕЧАНИЯ:\n";
    ss << "----------\n";
    ss << "- Все нормы вычислены как maximum-norm (максимальный модуль элемента)\n";
    
    if (use_exact_error) {
        ss << "- Для сравнения с истинным решением используется аналитическое решение задачи\n";
    }
    
    return ss.str();
}
