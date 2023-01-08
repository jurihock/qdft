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
    const F pi = std::acos(F(-1));

    config.samplerate = samplerate;
    config.bandwidth = bandwidth;
    config.resolution = resolution;
    config.latency = latency;
    config.window = window;
    config.quality = std::pow(std::pow(2.0, 1.0 / resolution) - 1.0, -1.0);
    config.size = std::ceil(resolution * std::log2(bandwidth.second / bandwidth.first));

    config.frequencies.resize(config.size);
    common.periods.resize(config.size);
    common.offsets.resize(config.size);
    common.weights.resize(config.size);

    for (size_t i = 0; i < config.size; ++i)
    {
      const double frequency = config.bandwidth.first * std::pow(2.0, i / config.resolution);

      config.frequencies[i] = frequency;

      const double period = std::ceil(config.quality * config.samplerate / frequency);

      common.periods[i] = static_cast<size_t>(period);

      const double offset = std::ceil((common.periods.front() - period)
        * std::clamp(config.latency * 0.5 + 0.5, 0.0, 1.0));

      common.offsets[i] = static_cast<size_t>(offset);

      const F weight = F(1) / period;

      common.weights[i] = weight;
    }

    common.inputs.resize(common.periods.front() + 1);

    for (const auto k : { -1, 0, +1 })
    {
      auto& kernel = kernels[k];

      kernel.outputs.resize(config.size);
      kernel.twiddles.resize(config.size);

      for (size_t i = 0; i < config.size; ++i)
      {
        const std::complex<F> twiddle = std::polar(F(1), F(+2) * pi * (config.quality + k) / common.periods[i]);

        kernel.twiddles[i] = twiddle;
      }

      kernel.fiddle = std::polar(F(1), F(-2) * pi * (config.quality + k));
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
    std::rotate(common.inputs.begin(), common.inputs.begin() + 1, common.inputs.end());

    common.inputs.back() = sample;

    const auto core = [&](const int k)
    {
      const auto& kernel = kernels[k];

      const auto& inputs = common.inputs;
      auto& outputs = kernels[k].outputs;

      for (size_t i = 0; i < config.size; ++i)
      {
        const size_t period = common.periods[i];
        const size_t offset = common.offsets[i];
        const F weight = common.weights[i];
        const std::complex<F> twiddle = kernel.twiddles[i];
        const std::complex<F> fiddle = kernel.fiddle;

        const std::complex<F> delta = (fiddle * inputs[offset + period] - inputs[offset]) * weight;

        outputs[i] = twiddle * (outputs[i] + delta);
      }
    };

    if (config.window)
    {
      const auto k = { -1, 0, +1 };

      std::for_each(k.begin(), k.end(), core);

      const auto& left = kernels[-1].outputs;
      const auto& middle = kernels[0].outputs;
      const auto& right = kernels[+1].outputs;

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

      const auto& middle = kernels[0].outputs;

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
    const auto& twiddles = kernels[0].twiddles;

    F sample = F(0);

    for (size_t i = 0; i < config.size; ++i)
    {
      sample += (dft[i] * twiddles[i]).real();
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

  qdft_config_t config;

  struct qdft_common_t
  {
    std::vector<size_t> periods;
    std::vector<size_t> offsets;
    std::vector<F> weights;
    std::vector<T> inputs;
  };

  qdft_common_t common;

  struct qdft_kernel_t
  {
    std::complex<F> fiddle;
    std::vector<std::complex<F>> twiddles;
    std::vector<std::complex<F>> outputs;
  };

  std::map<int, qdft_kernel_t> kernels;

};
