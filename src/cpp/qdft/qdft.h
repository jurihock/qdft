/**
 * Copyright (c) 2023 Juergen Hock
 *
 * SPDX-License-Identifier: MIT
 *
 * Constant Q Sliding DFT implementation according to [1] and [2].
 *
 * [1] Russell Bradford and Richard Dobson and John ffitch
 *     Sliding with a Constant Q
 *     International Conference on Digital Audio Effects (2008)
 *     https://www.dafx.de/paper-archive/2008/papers/dafx08_63.pdf
 *
 * [2] Benjamin Blankertz
 *     The Constant Q Transform
 *     https://doc.ml.tu-berlin.de/bbci/material/publications/Bla_constQ.pdf
 *
 * Source: https://github.com/jurihock/qdft
 **/

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <map>
#include <optional>
#include <utility>
#include <vector>

template <typename T = float, typename F = double>
class QDFT
{

public:

  QDFT(const double samplerate,
       const std::pair<double, double> bandwidth,
       const double resolution = 24,
       const double latency = 0,
       const std::optional<std::pair<double, double>> window = std::make_pair(+0.5,-0.5))
  {
    const F pi = F(2) * std::acos(F(-1));

    config.samplerate = samplerate;
    config.bandwidth = bandwidth;
    config.resolution = resolution;
    config.latency = latency;
    config.window = window;
    config.quality = std::pow(std::pow(2.0, 1.0 / resolution) - 1.0, -1.0);
    config.size = std::ceil(resolution * std::log2(bandwidth.second / bandwidth.first));

    config.frequencies.resize(config.size);

    for (size_t i = 0; i < config.size; ++i)
    {
      const double frequency = config.bandwidth.first * std::pow(2.0, i / config.resolution);

      config.frequencies[i] = frequency;
    }

    for (const auto k : { -1, 0, +1 })
    {
      auto& kernel = kernels[k];

      kernel.data.resize(config.size + 1);

      for (size_t i = 0; i < config.size; ++i)
      {
        const double period = std::ceil(config.quality * config.samplerate / config.frequencies[i]);

        kernel.data[i].period = static_cast<size_t>(period);

        const double offset = (kernel.data.front().period - period)
          * std::clamp(config.latency * 0.5 + 0.5, 0.0, 1.0);

        kernel.data[i].offset = static_cast<size_t>(offset);

        const F weight = F(1) / period;

        kernel.data[i].weight = weight;

        const std::complex<F> twiddle = std::polar(F(1), +pi * (config.quality + k) / period);

        kernel.data[i].twiddle = twiddle;
      }

      kernel.data.back().twiddle = std::polar(F(1), -pi * (config.quality + k));

      kernel.buffer.input.resize(kernel.data.front().period + 1);
      kernel.buffer.output.resize(config.size);
    }
  }

  size_t size() const
  {
    return config.size;
  }

  double quality() const
  {
    return config.quality;
  }

  double latency() const
  {
    return config.latency;
  }

  const std::vector<double>& frequencies() const
  {
    return config.frequencies;
  }

  void qdft(const T sample, std::complex<F>* const dft)
  {
    const auto core = [&](const int k)
    {
      const auto& data = kernels[k].data;

      auto& input = kernels[k].buffer.input;
      auto& output = kernels[k].buffer.output;

      std::rotate(input.begin(), input.begin() + 1, input.end());

      input.back() = sample;

      const std::complex<F> fiddle = data.back().twiddle;

      for (size_t i = 0; i < config.size; ++i)
      {
        const size_t period = data[i].period;
        const size_t offset = data[i].offset;

        const F weight = data[i].weight;
        const std::complex<F> twiddle = data[i].twiddle;

        const F left = input[offset + period];
        const F right = input[offset];

        const std::complex<F> delta = fiddle * left - right;

        output[i] = twiddle * (output[i] + delta * weight);
      }
    };

    if (config.window)
    {
      const auto k = { -1, 0, +1 };

      std::for_each(k.begin(), k.end(), core);

      const auto& left = kernels[-1].buffer.output;
      const auto& middle = kernels[0].buffer.output;
      const auto& right = kernels[+1].buffer.output;

      const F a = config.window.value().first;
      const F b = config.window.value().second / 2;

      for (size_t i = 0; i < config.size; ++i)
      {
        dft[i] = a * middle[i] + b * (left[i] + right[i]);
      }
    }
    else
    {
      core(0);

      const auto& middle = kernels[0].buffer.output;

      for (size_t i = 0; i < config.size; ++i)
      {
        dft[i] = middle[i];
      }
    }
  }

  void qdft(const size_t nsamples, const T* samples, std::complex<F>* const dfts)
  {
    for (size_t i = 0; i < nsamples; ++i)
    {
      qdft(samples[i], &dfts[i * config.size]);
    }
  }

  T iqdft(const std::complex<F>* dft)
  {
    const auto& data = kernels[0].data;

    F sample = F(0);

    for (size_t i = 0; i < config.size; ++i)
    {
      sample += (dft[i] * data[i].twiddle).real();
    }

    return static_cast<T>(sample);
  }

  void iqdft(const size_t nsamples, const std::complex<F>* dfts, T* const samples)
  {
    for (size_t i = 0; i < nsamples; ++i)
    {
      samples[i] = iqdft(&dfts[i * config.size]);
    }
  }

private:

  struct qdft_config_t
  {
    double samplerate;
    std::pair<double, double> bandwidth;
    double resolution;
    double latency;
    double quality;
    size_t size;
    std::vector<double> frequencies;
    std::optional<std::pair<double, double>> window;
  };

  struct qdft_data_t
  {
    size_t period;
    size_t offset;
    F weight;
    std::complex<F> twiddle;
  };

  struct qdft_buffer_t
  {
    std::vector<T> input;
    std::vector<std::complex<F>> output;
  };

  struct qdft_kernel_t
  {
    std::vector<qdft_data_t> data;
    qdft_buffer_t buffer;
  };

  qdft_config_t config;
  std::map<int, qdft_kernel_t> kernels;

};
