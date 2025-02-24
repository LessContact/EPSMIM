#include "GnuplotRealTime.h"

GnuplotRealTime::GnuplotRealTime(size_t nx_, size_t ny_, const std::string& term)
    : nx(nx_), ny(ny_), is_active(true), terminal_type(term) {
    pipe = popen("gnuplot", "w");
    if (!pipe) {
        throw std::runtime_error("Could not open pipe to gnuplot");
    }
    // Initialize gnuplot settings
    setupPlot();
}

GnuplotRealTime::~GnuplotRealTime() {
    is_active = false;
    if (pipe) {
        pclose(pipe);
    }
}

void GnuplotRealTime::setupPlot() {
    if (terminal_type == "png") {
        sendCommand("set terminal png size 700,600");
    } else {
        sendCommand("set terminal " + terminal_type + " size 700,600");
    }

    // Fix the ranges explicitly
    sendCommand("unset autoscale"); // Disable autoscaling
    sendCommand("set xrange [0:" + std::to_string(nx) + "]");
    sendCommand("set yrange [0:" + std::to_string(ny) + "]");

    sendCommand("set zrange[-4e-5:4e-5]");
    sendCommand("set cbrange[-4e-5:4e-5]");

    sendCommand("set palette gray");

    // Add these settings to reduce visual noise during updates
    sendCommand("unset grid");
    sendCommand("unset mouse");
    sendCommand("set border 15"); // This ensures consistent border display
}

void GnuplotRealTime::updatePlot(const std::vector<std::vector<double>>& data, const std::string& filename, bool save_png) {
    if (data.size() != ny || data[0].size() != nx) {
        throw std::runtime_error("Data size doesn't match dimensions");
    }

    // Write binary data to temporary file
    std::string tempfile = "temp.bin";
    FILE* fp = fopen(tempfile.c_str(), "wb+");
    if (!fp) {
        throw std::runtime_error("Could not open temporary file");
    }
    for (const auto& row : data) {
        size_t written = fwrite(row.data(), sizeof(double), row.size(), fp);
        if (written != row.size()) {
            fclose(fp);
            throw std::runtime_error("Failed to write the expected amount of data to temporary file");
        }
    }
    fclose(fp);

    if (save_png) {
        sendCommand("set output '" + filename + ".png'");
    } else {
        sendCommand("unset output");
    }

    sendCommand("set title '" + filename + "'");

    // Use 'every' to reduce update frequency if needed
    std::string plot_cmd = "plot '" + tempfile + "' binary array=(" +
                          std::to_string(ny) + "," + std::to_string(nx) +
                          ") format='%double' with image notitle";
    sendCommand(plot_cmd);

    if (!save_png) {
        // Add a small delay to reduce flicker
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        sendCommand("refresh");
    }
}


void GnuplotRealTime::sendCommand(const std::string& cmd) {
    fprintf(pipe, "%s\n", cmd.c_str());
    fflush(pipe);
}