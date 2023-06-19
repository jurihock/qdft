#include <qdft/qdft.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <iostream>
#include <vector>

int main()
{
  const auto samplerate = 44100;
  const auto bandwidth = std::make_pair(50.0, samplerate / 2.0);
  const auto resolution = 24;
  const auto latency = 0;
  const auto window = std::make_pair(+0.5, -0.5);

  const auto ta0 = std::chrono::high_resolution_clock::now();
  QDFT qdft(samplerate, bandwidth, resolution, latency, window);
  const auto tb0 = std::chrono::high_resolution_clock::now();
  const auto e0 = std::chrono::duration_cast<std::chrono::microseconds>(tb0 - ta0).count();

  std::cout << "PREP\t" << "CPP\t" << e0 << " us" << std::endl;

  const auto n = 1 * samplerate;
  const auto m = qdft.size();

  std::vector<float> x(n);
  std::vector<std::complex<double>> y(n * m);

  const auto runs = 10;

  for (auto run = 1; run <= runs; ++run)
  {
    std::cout << "RUN\t" << run << "/" << runs << std::endl;

    std::fill(x.begin(), x.end(), 0);
    std::fill(y.begin(), y.end(), 0);

    const auto ta1 = std::chrono::high_resolution_clock::now();
    qdft.qdft(n, x.data(), y.data());
    const auto tb1 = std::chrono::high_resolution_clock::now();
    const auto e1 = std::chrono::duration_cast<std::chrono::microseconds>(tb1 - ta1).count();

    const auto ta2 = std::chrono::high_resolution_clock::now();
    qdft.iqdft(n, y.data(), x.data());
    const auto tb2 = std::chrono::high_resolution_clock::now();
    const auto e2 = std::chrono::duration_cast<std::chrono::microseconds>(tb2 - ta2).count();

    std::cout << "\tQDFT\t" << e1 << " us" << std::endl;
    std::cout << "\tIQDFT\t" << e2 << " us" << std::endl;
  }

  return 0;
}
