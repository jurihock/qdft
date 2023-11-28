/**
 * Copyright (c) 2023 Juergen Hock
 *
 * SPDX-License-Identifier: MIT
 *
 * Constant-Q Sliding DFT implementation according to [1], [2], and [3].
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
 * [3] Christian Schörkhuber and Anssi Klapuri and Nicki Holighaus and Monika Dörfler
 *     A Matlab Toolbox for Efficient Perfect Reconstruction
 *     Time-Frequency Transforms with Log-Frequency Resolution
 *     http://www.aes.org/e-lib/browse.cfm?elib=17112
 *
 * Source: https://github.com/jurihock/qdft
 **/

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <deque>
#include <map>
#include <optional>
#include <utility>
#include <vector>

namespace qdft
{
  template <typename T = float, typename F = double>
  class QDFT
  {

  public:

    QDFT(const double samplerate,
         const std::pair<double, double> bandwidth,
         const double resolution = 24,
         const double quality = 0,
         const double latency = 0,
         const std::optional<std::pair<double, double>> window = std::make_pair(+0.5,-0.5))
    {
      const F pi = std::acos(F(-1));

      config.samplerate = samplerate;
      config.bandwidth = bandwidth;
      config.resolution = resolution;
      config.quality = quality;
      config.latency = latency;
      config.window = window;
      config.size = static_cast<size_t>(std::ceil(resolution * std::log2(bandwidth.second / bandwidth.first)));

      data.frequencies.resize(config.size);
      data.qualities.resize(config.size);
      data.latencies.resize(config.size);
      data.periods.resize(config.size);
      data.offsets.resize(config.size);
      data.weights.resize(config.size);

      const double alpha = std::pow(2.0, 1.0 / config.resolution) - 1.0;
      const double beta = (config.quality < 0) ? (alpha * 24.7 / 0.108) : config.quality;

      for (size_t i = 0; i < config.size; ++i)
      {
        const double frequency = config.bandwidth.first * std::pow(2.0, i / config.resolution);

        data.frequencies[i] = frequency;

        const double quality = frequency / (alpha * frequency + beta);

        data.qualities[i] = quality;

        const double period = std::ceil(quality * config.samplerate / frequency);

        data.periods[i] = static_cast<size_t>(period);

        const double offset = std::ceil((data.periods.front() - period)
          * std::clamp(config.latency * 0.5 + 0.5, 0.0, 1.0));

        data.offsets[i] = static_cast<size_t>(offset);

        const double latency = (data.periods.front() - offset) / config.samplerate;

        data.latencies[i] = latency;

        const F weight = F(1) / period;

        data.weights[i] = weight;
      }

      data.fiddles.resize(config.size * 3);
      data.twiddles.resize(config.size * 3);

      for (const int k : { -1, 0, +1 })
      {
        for (size_t i = 0, j = 1; i < config.size; ++i, j+=3)
        {
          const std::complex<F> fiddle = std::polar(F(1), F(-2) * pi * (data.qualities[i] + k));

          data.fiddles[j + k] = fiddle;

          const std::complex<F> twiddle = std::polar(F(1), F(+2) * pi * (data.qualities[i] + k) / data.periods[i]);

          data.twiddles[j + k] = twiddle;
        }
      }

      data.inputs.resize(data.periods.front() + 1);
      data.outputs.resize(config.size * 3);

      #ifndef NDEBUG
      for (size_t i = 0, j = 1; i < config.size; ++i, j+=3)
      {
        assert((data.offsets[i]) < data.inputs.size());
        assert((data.offsets[i] + data.periods[i]) < data.inputs.size());

        assert((-1 + 1) < data.fiddles.size());
        assert(( 0 + 1) < data.fiddles.size());
        assert((+1 + 1) < data.fiddles.size());

        assert((-1 + j) < data.twiddles.size());
        assert(( 0 + j) < data.twiddles.size());
        assert((+1 + j) < data.twiddles.size());

        assert((-1 + j) < data.outputs.size());
        assert(( 0 + j) < data.outputs.size());
        assert((+1 + j) < data.outputs.size());
      }
      #endif
    }

    size_t size() const
    {
      return config.size;
    }

    double samplerate() const
    {
      return config.samplerate;
    }

    const std::pair<double, double>& bandwidth() const
    {
      return config.bandwidth;
    }

    double resolution() const
    {
      return config.resolution;
    }

    double quality() const
    {
      return config.quality;
    }

    double latency() const
    {
      return config.latency;
    }

    const std::optional<std::pair<double, double>>& window() const
    {
      return config.window;
    }

    const std::vector<double>& frequencies() const
    {
      return data.frequencies;
    }

    const std::vector<double>& qualities() const
    {
      return data.qualities;
    }

    const std::vector<double>& latencies() const
    {
      return data.latencies;
    }

    void qdft(const T sample, std::complex<F>* const dft)
    {
      std::deque<T>& inputs = data.inputs;

      inputs.pop_front();
      inputs.push_back(sample);

      if (config.window)
      {
        const std::pair<double, double> w = config.window.value();
        const F a = w.first;
        const F b = w.second / 2;

        for (size_t i = 0, j = 1; i < config.size; ++i, j+=3)
        {
          const size_t period = data.periods[i];
          const size_t offset = data.offsets[i];
          const F weight = data.weights[i];

          const std::complex<F>* fiddles = data.fiddles.data() + j;
          const std::complex<F>* twiddles = data.twiddles.data() + j;

          const F left = inputs[offset + period];
          const F right = inputs[offset];

          const std::complex<F> delta1 = (fiddles[-1] * left - right) * weight;
          const std::complex<F> delta2 = (fiddles[ 0] * left - right) * weight;
          const std::complex<F> delta3 = (fiddles[+1] * left - right) * weight;

          std::complex<F>* const outputs = data.outputs.data() + j;

          outputs[-1] = twiddles[-1] * (outputs[-1] + delta1);
          outputs[ 0] = twiddles[ 0] * (outputs[ 0] + delta2);
          outputs[+1] = twiddles[+1] * (outputs[+1] + delta3);

          dft[i] = outputs[0] * a + (outputs[-1] + outputs[+1]) * b;
        }
      }
      else
      {
        for (size_t i = 0, j = 1; i < config.size; ++i, j+=3)
        {
          const size_t period = data.periods[i];
          const size_t offset = data.offsets[i];
          const F weight = data.weights[i];

          const std::complex<F>& fiddle = data.fiddles[j];
          const std::complex<F>& twiddle = data.twiddles[j];

          const F left = inputs[offset + period];
          const F right = inputs[offset];

          const std::complex<F> delta = (fiddle * left - right) * weight;

          std::complex<F>& output = data.outputs[j];

          output = twiddle * (output + delta);

          dft[i] = output;
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
      F sample = F(0);

      for (size_t i = 0, j = 1; i < config.size; ++i, j+=3)
      {
        const std::complex<F>& twiddle = data.twiddles[j];

        sample += (dft[i] * twiddle).real();
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
      double quality;
      double latency;
      std::optional<std::pair<double, double>> window;
      size_t size;
    };

    qdft_config_t config;

    struct qdft_data_t
    {
      std::vector<double> frequencies;
      std::vector<double> qualities;
      std::vector<double> latencies;
      std::vector<size_t> periods;
      std::vector<size_t> offsets;
      std::vector<F> weights;

      std::vector<std::complex<F>> fiddles;
      std::vector<std::complex<F>> twiddles;

      std::deque<T> inputs;
      std::vector<std::complex<F>> outputs;
    };

    qdft_data_t data;

  };
}
