#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class Particle {
public:
    double x, height;
    int type;

    Particle(double x, double height, int type) : x(x), height(height), type(type) {}
};

void readParticlesFromProf(const std::string& filename, std::vector<Particle>& particles) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "ファイルを開けませんでした: " << filename << std::endl;
        return;
    }

    std::string line;
    double time, posX, posY, height, velX, velY, volume, initialHeight, initialDensity2D;
    int particleCount, id, type;

    std::getline(ifs, line); // 先頭行は時間
    std::getline(ifs, line); // 次の行は粒子数
    while (std::getline(ifs, line)) {
        std::stringstream ss(line);
        ss >> id >> type >> posX >> posY >> height >> velX >> velY >> volume >> initialHeight >> initialDensity2D;
        if (type == 1) {
            particles.emplace_back(posX, height, type);
        }
    }
}

int determineInterval(double x) {
    if (x <= 0.025) {
        return 0;
    } else {
        return static_cast<int>((x - 0.025) / 0.05) + 1;
    }
}

void writeAverageHeightToCSV(const std::string& filename, const std::vector<Particle>& particles) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "CSVファイルを開けませんでした: " << filename << std::endl;
        return;
    }

    std::map<int, std::pair<double, int>> heightSums;

    for (const auto& particle : particles) {
        int interval = determineInterval(particle.x);
        heightSums[interval].first += particle.height;
        heightSums[interval].second++;
    }

    ofs << "Interval, Average Height\n";
    for (const auto& pair : heightSums) {
        double averageHeight = pair.second.first / pair.second.second;
        ofs << pair.first << ", " << averageHeight << std::endl;
    }
}

int main() {
    std::vector<Particle> particles;

    readParticlesFromProf("../result/prof/output_0050.prof", particles);
    writeAverageHeightToCSV("average_particle_heights.csv", particles);

    return 0;
}