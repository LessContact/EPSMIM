#ifndef GNUPLOT_REALTIME_H
#define GNUPLOT_REALTIME_H

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
#include <stdexcept>

class GnuplotRealTime {
private:
    FILE* pipe;
    const size_t nx;
    const size_t ny;
    bool is_active;
    std::string terminal_type;

public:
    GnuplotRealTime(size_t nx_, size_t ny_, const std::string& term="wxt");
    ~GnuplotRealTime();

    void setupPlot();
    void updatePlot(const std::vector<std::vector<double>>& data, const std::string& filename, bool save_png = false);

private:
    void sendCommand(const std::string& cmd);
};

#endif // GNUPLOT_REALTIME_H