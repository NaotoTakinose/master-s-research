#pragma once

#include <fstream>
#include <iostream>

#include "common.hpp"
#include "output.hpp"
#include "particle.hpp"

class System {
 private:
 public:
  // 入力変数
  double finishTime, outputPeriod;
  bool   debugModeFlag;

  double  time, initialTime, courant;
  int     fileNum, timestep;
  clock_t simStartTime, timestepStartTime, timestepFinishTime;
  FILE*   logFile;

  System(double finishTime, double outputPeriod, bool debugModeFlag);

  void startSimulation();
  void endSimulation();

  void dataOutput(const std::vector<Particle>& particles);
  void logOutput(double dt);

  void calSecondMinuteHour(int& second, int& minute, int& hour);
};

System::System(double finishTime, double outputPeriod, bool debugModeFlag) {
  this->finishTime    = finishTime;
  this->outputPeriod  = outputPeriod;
  this->debugModeFlag = debugModeFlag;
  this->time          = 0.0;
  this->fileNum       = 0;
  this->timestep      = 0;
  char filename[256];
  sprintf(filename, "../result/result.log");
  logFile = fopen(filename, "w");
  if (!logFile) {
    std::cerr << "cannnot write result/result.log" << std::endl;
    exit(-1);
  }
}

void System::startSimulation() {
  std::cout << std::endl
            << "*** START SIMULATION ***" << std::endl;
  simStartTime = std::clock();

  return;
}

void System::endSimulation() {
  clock_t simFinishTime = std::clock();

  int second, minute, hour;
  second = (double)(simFinishTime - simStartTime) / CLOCKS_PER_SEC;
  calSecondMinuteHour(second, minute, hour);
  printf("\nTotal Simulation Time = %dh %02dm %02ds\n", hour, minute, second);

  fclose(logFile);

  std::cout << std::endl
            << "*** END SIMULATION ***" << std::endl;
}

void System::dataOutput(const std::vector<Particle>& particles) {
  if (debugModeFlag || (time - initialTime) >= outputPeriod * double(fileNum)) {
    std::stringstream ss;
    ss << "../result/prof/output_";
    ss << std::setfill('0') << std::setw(4) << fileNum << ".prof";
    writeData::writeProf(ss, time, particles);

    ss.str("");
    ss << "../result/vtu/output_";
    ss << std::setfill('0') << std::setw(4) << fileNum << ".vtu";
    writeData::writeVtu(ss, time, particles);

    fileNum++;
  }
}

void System::logOutput(double dt) {
  int second, minute, hour;

  char elapsed[256];
  second = (double)(timestepFinishTime - simStartTime) / CLOCKS_PER_SEC;
  calSecondMinuteHour(second, minute, hour);
  sprintf(elapsed, "elapsed=%dh %02dm %02ds", hour, minute, second);

  double ave = ((double)(timestepFinishTime - simStartTime) / CLOCKS_PER_SEC) / timestep;
  if (timestep == 0) ave = 0.0;

  char remain[256];
  second = ((finishTime - time) / time) * ave * timestep;
  calSecondMinuteHour(second, minute, hour);
  if (timestep == 0)
    sprintf(remain, "remain=-h --m --s");
  else
    sprintf(remain, "remain=%dh %02dm %02ds", hour, minute, second);

  double last = (double)(timestepFinishTime - timestepStartTime) / CLOCKS_PER_SEC;

  // ターミナル出力
  printf("%d: dt=%.2es   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs   last=%.3lfs   out=%dfiles   Courant=%.2lf\n",
         timestep, dt, time, finishTime, elapsed, remain, ave, last, fileNum, courant);

  // ログ出力
  fprintf(logFile, "%d: dt=%.2es   t=%.3lfs   fin=%.1lfs   %s   %s   ave=%.3lfs   last=%.3lfs   out=%dfiles   Courant=%.2lf\n",
          timestep, dt, time, finishTime, elapsed, remain, ave, last, fileNum, courant);

  // エラー出力
  fprintf(stderr, "%4d: t=%.3lfs\n", timestep, time);
}

void System::calSecondMinuteHour(int& second, int& minute, int& hour) {
  hour = second / 3600;
  second %= 3600;
  minute = second / 60;
  second %= 60;
}