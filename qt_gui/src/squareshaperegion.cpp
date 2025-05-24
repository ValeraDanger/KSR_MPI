#include "squareshaperegion.h"
#include <QDebug>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QValue3DAxis>
#include <algorithm>
#include <set>
#include <cmath>
#include <limits>

// Constructor
SquareShapeRegion::SquareShapeRegion(Q3DSurface* graph3D)
    : ShapeRegion(graph3D) {
}

// Destructor
SquareShapeRegion::~SquareShapeRegion() {
    clearAllSurfaces();
}

bool SquareShapeRegion::createSurfaces(
    const std::vector<double>& numericalSolution,
    const std::vector<double>& trueSolution,
    const std::vector<double>& errorValues,
    const std::vector<double>& xCoords,
    const std::vector<double>& yCoords,
    double domainXMin, double domainXMax,
    double domainYMin, double domainYMax,
    int decimationFactor,
    int connectorRows,
    const std::vector<double>& initialApproximation)
{
    // Очистка предыдущих поверхностей
    clearAllSurfaces();

    // Проверка входных данных
    if (numericalSolution.empty() || xCoords.empty() || yCoords.empty() ||
        numericalSolution.size() != xCoords.size() ||
        numericalSolution.size() != yCoords.size()) {
        qDebug() << "SquareShapeRegion::createSurfaces: Incorrect data provided";
        return false;
    }

    // Сохраняем границы домена
    m_currentDomainXMin = domainXMin;
    m_currentDomainXMax = domainXMax;
    m_currentDomainYMin = domainYMin;
    m_currentDomainYMax = domainYMax;

    // Устанавливаем диапазоны осей для графика
    m_graph3D->axisX()->setRange(domainXMin, domainXMax);
    m_graph3D->axisZ()->setRange(domainYMin, domainYMax);

    // Сохраняем данные для последующего использования
    m_lastNumericalSolutionData = numericalSolution;
    
    // Создаем вектор нулевых значений для начального приближения или используем предоставленное
    if (initialApproximation.empty()) {
        m_lastInitialApproximationData = std::vector<double>(numericalSolution.size(), 0.0);
    } else {
        m_lastInitialApproximationData = initialApproximation;
    }

    // Обновляем диапазоны осей
    updateAxesRanges(numericalSolution);

    // Создание поверхности численного решения
    m_solutionSeries = createSquareSurface(
        numericalSolution,
        xCoords,
        yCoords,
        QColor(0, 0, 255, 255), // Синий цвет
        "Численное решение",
        decimationFactor
    );
    
    // Создание поверхности начального приближения
    m_initialApproximationSeries = createSquareSurface(
        m_lastInitialApproximationData,
        xCoords,
        yCoords,
        QColor(150, 150, 150, 200), // Серый цвет с прозрачностью
        "Начальное приближение",
        decimationFactor
    );
    
    // Делаем серию начального приближения изначально невидимой
    if (m_initialApproximationSeries) {
        m_initialApproximationSeries->setVisible(false);
    }

    // Создание поверхности точного решения (если данные предоставлены)
    if (!trueSolution.empty() && trueSolution.size() == numericalSolution.size()) {
        m_lastTrueSolutionData = trueSolution;
        m_trueSolutionSeries = createSquareSurface(
            trueSolution,
            xCoords,
            yCoords,
            QColor(0, 255, 0, 200), // Зеленый цвет с прозрачностью
            "Точное решение",
            decimationFactor
        );
        
        // Делаем серию точного решения изначально невидимой
        if (m_trueSolutionSeries) {
            m_trueSolutionSeries->setVisible(false);
        }
    }

    // Создание поверхности ошибки (если данные предоставлены)
    if (!errorValues.empty() && errorValues.size() == numericalSolution.size()) {
        m_lastErrorData = errorValues;
        m_errorSeries = createSquareSurface(
            errorValues,
            xCoords,
            yCoords,
            QColor(255, 0, 0, 200), // Красный цвет с прозрачностью
            "Ошибка",
            decimationFactor
        );
        
        // Делаем серию ошибки изначально невидимой
        if (m_errorSeries) {
            m_errorSeries->setVisible(false);
        }
    }

    return true;
}

void SquareShapeRegion::setNumericalSolutionVisible(bool visible) {
    if (m_solutionSeries) {
        m_solutionSeries->setVisible(visible);
        updateDynamicAxesRanges();
    }
}

void SquareShapeRegion::setTrueSolutionVisible(bool visible) {
    if (m_trueSolutionSeries) {
        m_trueSolutionSeries->setVisible(visible);
        updateDynamicAxesRanges();
    }
}

void SquareShapeRegion::setErrorSurfaceVisible(bool visible) {
    if (m_errorSeries) {
        m_errorSeries->setVisible(visible);
        updateDynamicAxesRanges();
    }
}

// Методы для работы с решением на более мелкой сетке

void SquareShapeRegion::setRefinedGridSolutionVisible(bool visible) {
    if (m_refinedGridSolutionSeries) {
        m_refinedGridSolutionSeries->setVisible(visible);
        updateDynamicAxesRanges();
    }
}

void SquareShapeRegion::setSolutionRefinedDiffVisible(bool visible) {
    if (m_solutionRefinedDiffSeries) {
        m_solutionRefinedDiffSeries->setVisible(visible);
        updateDynamicAxesRanges();
    }
}

void SquareShapeRegion::setInitialApproximationVisible(bool visible) {
    if (m_initialApproximationSeries) {
        m_initialApproximationSeries->setVisible(visible);
        updateDynamicAxesRanges();
    }
}

bool SquareShapeRegion::createSurfacesWithRefinedGrid(
    const std::vector<double>& numericalSolution,
    const std::vector<double>& refinedGridSolution,
    const std::vector<double>& solutionRefinedDiff,
    const std::vector<double>& xCoords,
    const std::vector<double>& yCoords,
    const std::vector<double>& refinedGridXCoords,
    const std::vector<double>& refinedGridYCoords,
    double domainXMin, double domainXMax,
    double domainYMin, double domainYMax,
    int decimationFactor)
{
    // Очистка предыдущих поверхностей
    clearAllSurfaces();

    // Проверка входных данных
    if (numericalSolution.empty() || xCoords.empty() || yCoords.empty() ||
        numericalSolution.size() != xCoords.size() || numericalSolution.size() != yCoords.size()) {
        qDebug() << "SquareShapeRegion::createSurfacesWithRefinedGrid: Incorrect data provided";
        return false;
    }

    // Сохраняем границы домена
    m_currentDomainXMin = domainXMin;
    m_currentDomainXMax = domainXMax;
    m_currentDomainYMin = domainYMin;
    m_currentDomainYMax = domainYMax;

    // Устанавливаем диапазоны осей для графика
    m_graph3D->axisX()->setRange(domainXMin, domainXMax);
    m_graph3D->axisZ()->setRange(domainYMin, domainYMax);

    // Сохраняем данные для последующего использования
    m_lastNumericalSolutionData = numericalSolution;
    
    // Создаем вектор нулевых значений для начального приближения
    m_lastInitialApproximationData = std::vector<double>(numericalSolution.size(), 0.0);

    // Обновляем диапазоны осей
    updateAxesRanges(numericalSolution);

    // Создание поверхности численного решения
    m_solutionSeries = createSquareSurface(
        numericalSolution,
        xCoords,
        yCoords,
        QColor(0, 0, 255, 255), // Синий цвет
        "Численное решение v(N)",
        decimationFactor
    );

    // Создание поверхности начального приближения
    m_initialApproximationSeries = createSquareSurface(
        m_lastInitialApproximationData,
        xCoords,
        yCoords,
        QColor(150, 150, 150, 200), // Серый цвет с прозрачностью
        "Начальное приближение",
        decimationFactor
    );
    
    // Делаем серию начального приближения изначально невидимой
    if (m_initialApproximationSeries) {
        m_initialApproximationSeries->setVisible(false);
    }

    // Создание поверхности решения на более мелкой сетке (если данные предоставлены)
    if (!refinedGridSolution.empty() && !refinedGridXCoords.empty() && !refinedGridYCoords.empty()) {
        m_lastRefinedGridSolutionData = refinedGridSolution;
        m_refinedGridSolutionSeries = createSquareSurface(
            refinedGridSolution,
            refinedGridXCoords,
            refinedGridYCoords,
            QColor(0, 255, 0, 200), // Зеленый цвет с прозрачностью
            "Решение на мелкой сетке v2(2N)",
            decimationFactor
        );
        
        // Делаем серию решения на мелкой сетке изначально невидимой
        if (m_refinedGridSolutionSeries) {
            m_refinedGridSolutionSeries->setVisible(false);
        }
    }

    // Создание поверхности разницы между решениями (если данные предоставлены)
    if (!solutionRefinedDiff.empty() && solutionRefinedDiff.size() == numericalSolution.size()) {
        m_lastSolutionRefinedDiffData = solutionRefinedDiff;
        m_solutionRefinedDiffSeries = createSquareSurface(
            solutionRefinedDiff,
            xCoords,
            yCoords,
            QColor(255, 0, 0, 200), // Красный цвет с прозрачностью
            "Разница между решениями",
            decimationFactor
        );
        
        // Делаем серию разницы изначально невидимой
        if (m_solutionRefinedDiffSeries) {
            m_solutionRefinedDiffSeries->setVisible(false);
        }
    }

    return true;
}

void SquareShapeRegion::clearAllSurfaces() {
    // Проверяем, что у нас есть объект графика
    if (!m_graph3D) {
        return;
    }

    // Удаляем поверхности из графика
    QList<QSurface3DSeries*> series = m_graph3D->seriesList();
    for (QSurface3DSeries* s : series) {
        m_graph3D->removeSeries(s);
    }

    // Сбрасываем указатели на серии
    m_solutionSeries = nullptr;
    m_trueSolutionSeries = nullptr;
    m_errorSeries = nullptr;
    m_refinedGridSolutionSeries = nullptr;
    m_solutionRefinedDiffSeries = nullptr;
    m_initialApproximationSeries = nullptr;

    // Очищаем сохраненные данные
    m_lastNumericalSolutionData.clear();
    m_lastTrueSolutionData.clear();
    m_lastErrorData.clear();
    m_lastRefinedGridSolutionData.clear();
    m_lastSolutionRefinedDiffData.clear();
    m_lastInitialApproximationData.clear();
}

void SquareShapeRegion::updateAxesRanges(const std::vector<double>& values) {
    if (values.empty()) {
        return;
    }

    // Находим минимальное и максимальное значения
    auto [minIt, maxIt] = std::minmax_element(values.begin(), values.end());
    m_currentValueMin = *minIt;
    m_currentValueMax = *maxIt;

    // Добавляем небольшой отступ для лучшей визуализации
    double padding = (m_currentValueMax - m_currentValueMin) * 0.05;
    
    // Если min и max совпадают, добавляем стандартный отступ
    if (std::abs(m_currentValueMax - m_currentValueMin) < 1e-10) {
        padding = 0.5;
    }

    // Устанавливаем диапазон оси Y (значение)
    if (m_graph3D && m_graph3D->axisY()) {
        m_graph3D->axisY()->setRange(m_currentValueMin - padding, m_currentValueMax + padding);
    }
}

void SquareShapeRegion::updateDynamicAxesRanges() {
    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();
    bool dataAvailable = false;

    // Проверяем видимые серии и обновляем мин/макс значения
    if (m_solutionSeries && m_solutionSeries->isVisible() && !m_lastNumericalSolutionData.empty()) {
        auto [minIt, maxIt] = std::minmax_element(m_lastNumericalSolutionData.begin(), m_lastNumericalSolutionData.end());
        minVal = std::min(minVal, *minIt);
        maxVal = std::max(maxVal, *maxIt);
        dataAvailable = true;
    }

    if (m_trueSolutionSeries && m_trueSolutionSeries->isVisible() && !m_lastTrueSolutionData.empty()) {
        auto [minIt, maxIt] = std::minmax_element(m_lastTrueSolutionData.begin(), m_lastTrueSolutionData.end());
        minVal = std::min(minVal, *minIt);
        maxVal = std::max(maxVal, *maxIt);
        dataAvailable = true;
    }

    if (m_errorSeries && m_errorSeries->isVisible() && !m_lastErrorData.empty()) {
        auto [minIt, maxIt] = std::minmax_element(m_lastErrorData.begin(), m_lastErrorData.end());
        minVal = std::min(minVal, *minIt);
        maxVal = std::max(maxVal, *maxIt);
        dataAvailable = true;
    }

    if (m_refinedGridSolutionSeries && m_refinedGridSolutionSeries->isVisible() && !m_lastRefinedGridSolutionData.empty()) {
        auto [minIt, maxIt] = std::minmax_element(m_lastRefinedGridSolutionData.begin(), m_lastRefinedGridSolutionData.end());
        minVal = std::min(minVal, *minIt);
        maxVal = std::max(maxVal, *maxIt);
        dataAvailable = true;
    }

    if (m_solutionRefinedDiffSeries && m_solutionRefinedDiffSeries->isVisible() && !m_lastSolutionRefinedDiffData.empty()) {
        auto [minIt, maxIt] = std::minmax_element(m_lastSolutionRefinedDiffData.begin(), m_lastSolutionRefinedDiffData.end());
        minVal = std::min(minVal, *minIt);
        maxVal = std::max(maxVal, *maxIt);
        dataAvailable = true;
    }

    // Добавляем проверку для начального приближения
    if (m_initialApproximationSeries && m_initialApproximationSeries->isVisible() && !m_lastInitialApproximationData.empty()) {
        auto [minIt, maxIt] = std::minmax_element(m_lastInitialApproximationData.begin(), m_lastInitialApproximationData.end());
        minVal = std::min(minVal, *minIt);
        maxVal = std::max(maxVal, *maxIt);
        dataAvailable = true;
    }

    // Если есть видимые данные, обновляем диапазон оси Y
    if (dataAvailable && m_graph3D && m_graph3D->axisY()) {
        // Добавляем небольшой отступ для лучшей визуализации
        double padding = (maxVal - minVal) * 0.05;
        
        // Если min и max совпадают, добавляем стандартный отступ
        if (std::abs(maxVal - minVal) < 1e-10) {
            padding = 0.5;
        }

        m_graph3D->axisY()->setRange(minVal - padding, maxVal + padding);
    }
}

QSurface3DSeries* SquareShapeRegion::createSquareSurface(
    const std::vector<double>& data,
    const std::vector<double>& xCoords,
    const std::vector<double>& yCoords,
    const QColor& color,
    const QString& name,
    int decimationFactor)
{
    // Проверяем входные данные
    if (data.empty() || xCoords.empty() || yCoords.empty() ||
        data.size() != xCoords.size() || data.size() != yCoords.size()) {
        qDebug() << "SquareShapeRegion::createSquareSurface: Incorrect data provided";
        return nullptr;
    }

    // Создаем множества для уникальных координат x и y
    std::set<double> uniqueX(xCoords.begin(), xCoords.end());
    std::set<double> uniqueY(yCoords.begin(), yCoords.end());

    // Преобразуем в вектора для индексирования
    std::vector<double> uniqueXVec(uniqueX.begin(), uniqueX.end());
    std::vector<double> uniqueYVec(uniqueY.begin(), uniqueY.end());

    // Определяем необходимый шаг для прореживания сетки до максимум 100x100 узлов
    const int MAX_GRID_SIZE = 100;
    int stepX = 1;
    int stepY = 1;

    // Если размер сетки превышает 100x100, вычисляем шаг прореживания
    if (uniqueX.size() > MAX_GRID_SIZE || uniqueY.size() > MAX_GRID_SIZE) {
        stepX = std::ceil(static_cast<double>(uniqueX.size()) / MAX_GRID_SIZE);
        stepY = std::ceil(static_cast<double>(uniqueY.size()) / MAX_GRID_SIZE);
        
        // Учитываем пользовательский коэффициент прореживания, если он больше полученного
        if (decimationFactor > stepX) stepX = decimationFactor;
        if (decimationFactor > stepY) stepY = decimationFactor;
        
        qDebug() << "Grid size reduced: " << uniqueX.size() << "x" << uniqueY.size() 
                 << " -> ~" << (uniqueX.size() / stepX) << "x" << (uniqueY.size() / stepY)
                 << " (steps: " << stepX << "x" << stepY << ")";
    } else if (decimationFactor > 1) {
        // Если сетка меньше порога, но задан пользовательский decimationFactor
        stepX = decimationFactor;
        stepY = decimationFactor;
    }

    // Создаем массив для данных поверхности
    QSurfaceDataArray* dataArray = new QSurfaceDataArray();
    
    // Индексы для первого и последнего X, Y (гарантируем их включение)
    const size_t firstXIndex = 0;
    const size_t firstYIndex = 0;
    const size_t lastXIndex = uniqueXVec.size() - 1;
    const size_t lastYIndex = uniqueYVec.size() - 1;
    
    // Отслеживаем индексы для Y
    std::vector<size_t> yIndices;
    
    // Добавляем первый индекс
    yIndices.push_back(firstYIndex);
    
    // Добавляем промежуточные индексы с шагом stepY
    for (size_t yIndex = firstYIndex + stepY; yIndex < lastYIndex; yIndex += stepY) {
        yIndices.push_back(yIndex);
    }
    
    // Добавляем последний индекс, если он еще не добавлен
    if (yIndices.empty() || yIndices.back() != lastYIndex) {
        yIndices.push_back(lastYIndex);
    }

    // Проходим по выбранным координатам Y
    for (size_t yIndex : yIndices) {
        double yValue = uniqueYVec[yIndex];
        QSurfaceDataRow* dataRow = new QSurfaceDataRow();
        
        // Отслеживаем индексы для X
        std::vector<size_t> xIndices;
        
        // Добавляем первый индекс
        xIndices.push_back(firstXIndex);
        
        // Добавляем промежуточные индексы с шагом stepX
        for (size_t xIndex = firstXIndex + stepX; xIndex < lastXIndex; xIndex += stepX) {
            xIndices.push_back(xIndex);
        }
        
        // Добавляем последний индекс, если он еще не добавлен
        if (xIndices.empty() || xIndices.back() != lastXIndex) {
            xIndices.push_back(lastXIndex);
        }
        
        // Проходим по выбранным координатам X
        for (size_t xIndex : xIndices) {
            double xValue = uniqueXVec[xIndex];
            
            // Находим значение функции для данной точки (x, y)
            double functionValue = 0.0;
            
            // Находим индекс в оригинальных данных
            for (size_t i = 0; i < xCoords.size(); ++i) {
                if (std::abs(xCoords[i] - xValue) < 1e-9 && 
                    std::abs(yCoords[i] - yValue) < 1e-9) {
                    functionValue = data[i];
                    break;
                }
            }
            
            // Добавляем точку в ряд данных
            *dataRow << QVector3D(xValue, functionValue, yValue);
        }
        
        // Добавляем ряд в массив данных
        *dataArray << dataRow;
    }
    
    // Создаем прокси и серию
    QSurfaceDataProxy* proxy = new QSurfaceDataProxy();
    proxy->resetArray(dataArray);
    
    QSurface3DSeries* series = new QSurface3DSeries(proxy);
    series->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
    series->setFlatShadingEnabled(true);
    series->setBaseColor(color);
    series->setName(name);
    
    // Добавляем серию в график
    if (m_graph3D) {
        m_graph3D->addSeries(series);
    }
    
    return series;
}