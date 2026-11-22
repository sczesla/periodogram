import 'dart:math';

import 'frequency_grid.dart';
import 'numerics.dart';
import 'sine_fit.dart';

/// Normalizations of the periodogram power; see [Periodogram.normalized].
///
/// All are conversions of the power `p` of Zechmeister & Kürster (2009),
/// the default, which is the relative reduction of chi² by the sine:
/// `p = (chi²_0 - chi²(f)) / chi²_0`, where `chi²_0` is the chi² of the
/// constant model (the weighted mean).
enum Normalization {
  /// `p` itself, between 0 and 1 (default of [Periodogram.power]).
  zk,

  /// `(N - 1) / 2 * p` (Horne & Baliunas 1986).
  horneBaliunas,

  /// `(N - 3) / 2 * p / (1 - p_max)` (Cumming et al. 1999), with the
  /// highest power `p_max` of the periodogram.
  cumming,

  /// chi² of the sine fit, `chi²_0 * (1 - p)`. Without errors, unit
  /// variance is assumed (the sum of the squared residuals).
  chiSquared,

  /// Weighted rms of the residuals of the sine fit, `sqrt(YY * (1 - p))`,
  /// with the weighted variance `YY` of the data.
  wrms,

  /// Gaussian log-likelihood of the sine fit,
  /// `-chi²/2 - sum_i ln(2 pi err_i²) / 2`. Without errors, unit variance
  /// is assumed.
  lnL,

  /// Log-likelihood improvement over the constant model, `chi²_0 * p / 2`.
  /// For correct errors, this is the normalization of Scargle (1982),
  /// generalized to weights.
  dlnL,
}

/// Generalized Lomb-Scargle periodogram (GLS) of the data ([t], [y]) with
/// optional errors [err], as described by Zechmeister & Kürster (2009,
/// A&A 496, 577) and earlier by Ferraz-Mello (1981, AJ 86, 619).
///
/// At each frequency f, the sine `a cos(2 pi f t) + b sin(2 pi f t) + c`
/// is fitted by weighted least squares (weights `1 / err²`, all equal
/// without [err]). The power is the relative improvement of chi² over the
/// constant model, between 0 (no improvement) and 1 (perfect fit); see
/// [Periodogram.power].
///
/// [fitMean] fits the offset `c` together with the sine (the "floating
/// mean"). If false, the weighted mean is subtracted beforehand and only
/// `a` and `b` are fitted; without errors, this is the classical
/// Lomb-Scargle periodogram, see [lombScargle].
///
/// The frequencies (in units of 1 / unit of [t]) are either given
/// explicitly by [frequencies], or form an evenly spaced grid from [fMin]
/// to [fMax] (exclusive) with step `1 / (T * oversampling)`, where `T` is
/// the time span of the data. By default, the grid starts at one step and
/// ends at [nyquistFactor] times the average Nyquist frequency
/// `N / (2 T)`. For uneven sampling, periods shorter than the average
/// Nyquist period can be detectable; increase [nyquistFactor] to search
/// them.
///
/// Throws an [ArgumentError] if there are fewer than 4 points, the lengths
/// differ, values are not finite, errors are not positive, all times are
/// equal, frequencies are not positive, or [frequencies] is combined with
/// [fMin] or [fMax].
///
/// The computation takes `O(N * number of frequencies)` time. For large
/// data sets in a Flutter app, run it in a background isolate, e.g. with
/// `Isolate.run(() => gls(t, y))`.
Periodogram gls(
  List<double> t,
  List<double> y, {
  List<double>? err,
  bool fitMean = true,
  List<double>? frequencies,
  double? fMin,
  double? fMax,
  double oversampling = 10,
  double nyquistFactor = 1,
}) {
  _checkData(t, y, err);
  final n = t.length;
  final tMin = t.reduce(min);
  final tBase = t.reduce(max) - tMin;
  if (tBase == 0) {
    throw ArgumentError("All times are equal; the time span must be > 0");
  }
  final grid = frequencyGrid(
    n: n,
    timeSpan: tBase,
    frequencies: frequencies,
    fMin: fMin,
    fMax: fMax,
    oversampling: oversampling,
    nyquistFactor: nyquistFactor,
  );
  final freq = grid.frequencies;

  // Normalized weights (sum 1) and times relative to the first one
  final w = [
    for (int i = 0; i < n; i++) err == null ? 1.0 : 1 / (err[i] * err[i]),
  ];
  final wSum = w.fold(0.0, (s, v) => s + v);
  for (int i = 0; i < n; i++) {
    w[i] /= wSum;
  }
  final th = [for (final ti in t) ti - tMin];

  // Weighted mean (Eq. 7) and weighted variance about it (Eq. 10)
  double mean = 0;
  for (int i = 0; i < n; i++) {
    mean += w[i] * y[i];
  }
  final wy = [for (int i = 0; i < n; i++) w[i] * (y[i] - mean)];
  double yy = 0;
  for (int i = 0; i < n; i++) {
    yy += wy[i] * (y[i] - mean);
  }
  if (yy == 0) {
    throw ArgumentError("All y values are equal; nothing to search");
  }

  final nf = freq.length;
  final power = List<double>.filled(nf, 0);
  final a = List<double>.filled(nf, 0);
  final b = List<double>.filled(nf, 0);
  final c = List<double>.filled(nf, 0);
  final s = List<double>.filled(nf, 0);

  for (int k = 0; k < nf; k++) {
    final omega = 2 * pi * freq[k];
    double sumC = 0, sumS = 0, yc = 0, ys = 0, cc = 0, cs = 0;
    for (int i = 0; i < n; i++) {
      final x = omega * th[i];
      final cosx = cos(x), sinx = sin(x);
      sumC += w[i] * cosx; // Eq. 8
      sumS += w[i] * sinx; // Eq. 9
      yc += wy[i] * cosx; // Eq. 11
      ys += wy[i] * sinx; // Eq. 12
      cc += w[i] * cosx * cosx; // Eq. 13
      cs += w[i] * cosx * sinx; // Eq. 15
    }
    double ss = 1 - cc;
    if (fitMean) {
      cc -= sumC * sumC; // Eq. 13
      ss -= sumS * sumS; // Eq. 14
      cs -= sumC * sumS; // Eq. 15
    }
    final r = _solve(yc, ys, cc, ss, cs, yy);
    power[k] = r.p;
    a[k] = r.a;
    b[k] = r.b;
    c[k] = sumC;
    s[k] = sumS;
  }

  double sumLnVar = 0;
  for (int i = 0; i < n; i++) {
    sumLnVar += log(2 * pi * (err == null ? 1.0 : err[i] * err[i]));
  }

  return Periodogram._(
    frequencies: List.unmodifiable(freq),
    power: List.unmodifiable(power),
    fitMean: fitMean,
    n: n,
    weightedMean: mean,
    weightedVariance: yy,
    weightSum: wSum,
    sumLn2piVariance: sumLnVar,
    independentFrequencies: grid.independentFrequencies,
    timeOrigin: tMin,
    cosCoefficients: a,
    sinCoefficients: b,
    meanCos: c,
    meanSin: s,
  );
}

/// Classical Lomb-Scargle periodogram (Lomb 1976, Scargle 1982) of the
/// data ([t], [y]): the mean is subtracted, and a sine without offset is
/// fitted with equal weights at each frequency.
///
/// Same as `gls(t, y, fitMean: false)`; see [gls] for the frequency grid
/// and the exceptions. The power is in the normalization of Zechmeister &
/// Kürster (2009); the classical normalization by the data variance is
/// `N / 2` times it (`(N - 1) / 2` times for the sample variance, see
/// [Normalization.horneBaliunas]).
Periodogram lombScargle(
  List<double> t,
  List<double> y, {
  List<double>? frequencies,
  double? fMin,
  double? fMax,
  double oversampling = 10,
  double nyquistFactor = 1,
}) =>
    gls(
      t,
      y,
      fitMean: false,
      frequencies: frequencies,
      fMin: fMin,
      fMax: fMax,
      oversampling: oversampling,
      nyquistFactor: nyquistFactor,
    );

void _checkData(List<double> t, List<double> y, List<double>? err) {
  if (t.length != y.length || (err != null && err.length != t.length)) {
    throw ArgumentError(
      "t, y, and err must have the same length "
      "(got ${t.length}, ${y.length}, ${err?.length})",
    );
  }
  if (t.length < 4) {
    throw ArgumentError("At least 4 data points are needed (got ${t.length})");
  }
  for (int i = 0; i < t.length; i++) {
    if (!t[i].isFinite || !y[i].isFinite) {
      throw ArgumentError(
        "t and y must be finite (index $i: t = ${t[i]}, "
        "y = ${y[i]})",
      );
    }
    if (err != null && !(err[i] > 0 && err[i].isFinite)) {
      throw ArgumentError(
        "err must be finite and > 0 (index $i: err = ${err[i]})",
      );
    }
  }
}

/// Least-squares solution for the sine coefficients and the power (Eqs. 5
/// and 6). If cosine and sine are (nearly) collinear at the data points,
/// e.g. for even sampling at multiples of the Nyquist frequency, the fit
/// uses their common direction only.
({double a, double b, double p}) _solve(
  double yc,
  double ys,
  double cc,
  double ss,
  double cs,
  double yy,
) {
  final d = cc * ss - cs * cs; // Eq. 6
  final trace = cc + ss;
  if (trace <= 1e-15) return (a: 0, b: 0, p: 0);
  if (d > 1e-12 * trace * trace) {
    return (
      a: (yc * ss - ys * cs) / d,
      b: (ys * cc - yc * cs) / d,
      p: (ss * yc * yc + cc * ys * ys - 2 * cs * yc * ys) / (yy * d), // Eq. 5
    );
  }
  // Rank one: [[cc, cs], [cs, ss]] = trace * v v^T
  var vx = cc >= ss ? cc : cs, vy = cc >= ss ? cs : ss;
  final norm = sqrt(vx * vx + vy * vy);
  vx /= norm;
  vy /= norm;
  final proj = vx * yc + vy * ys;
  final alpha = proj / trace;
  return (a: alpha * vx, b: alpha * vy, p: proj * proj / (trace * yy));
}

/// Result of [gls] or [lombScargle]: the power at each frequency, the
/// best-fit sine, and the false-alarm probability.
class Periodogram {
  Periodogram._({
    required this.frequencies,
    required this.power,
    required this.fitMean,
    required this.n,
    required this.weightedMean,
    required double weightedVariance,
    required double weightSum,
    required double sumLn2piVariance,
    required this.independentFrequencies,
    required double timeOrigin,
    required List<double> cosCoefficients,
    required List<double> sinCoefficients,
    required List<double> meanCos,
    required List<double> meanSin,
  })  : _yy = weightedVariance,
        _wSum = weightSum,
        _sumLnVar = sumLn2piVariance,
        _tMin = timeOrigin,
        _a = cosCoefficients,
        _b = sinCoefficients,
        _c = meanCos,
        _s = meanSin {
    int k = 0;
    for (int i = 1; i < power.length; i++) {
      if (power[i] > power[k]) k = i;
    }
    bestIndex = k;
  }

  /// The frequencies, in units of 1 / unit of the times.
  final List<double> frequencies;

  /// Power at each frequency in the normalization of Zechmeister & Kürster
  /// (2009): the relative reduction of chi² by the sine fit, between 0 and
  /// 1. Other normalizations: [normalized].
  final List<double> power;

  /// Whether the offset was fitted together with the sine ([gls]) or the
  /// mean was subtracted beforehand ([lombScargle]).
  final bool fitMean;

  /// Number of data points.
  final int n;

  /// Weighted mean of the data (all weights equal without errors).
  final double weightedMean;

  /// Ad hoc estimate of the number of independent frequencies in the
  /// searched band, `(f_max - f_min) * T` with the time span `T`
  /// (Zechmeister & Kürster 2009, Eq. 24), but at least 1. It enters the
  /// false-alarm probability [fap].
  final double independentFrequencies;

  /// Index of the highest power in [frequencies] and [power] (the first
  /// one if several are equal).
  late final int bestIndex;

  final double _yy, _wSum, _sumLnVar, _tMin;
  final List<double> _a, _b, _c, _s;

  /// Weighted rms of the data about [weightedMean].
  double get weightedRms => sqrt(_yy);

  /// The sine fit at the highest power, see [sineFitAt].
  SineFit get best => sineFitAt(bestIndex);

  /// The sine fitted at `frequencies[index]`.
  ///
  /// The errors are the approximations of Zechmeister & Kürster (2009),
  /// which assume well sampled data. The frequency error follows from the
  /// curvature of the power at [index] and is only defined at a local
  /// maximum that is not at the edge of the grid; otherwise it is NaN.
  /// Throws a [RangeError] if [index] is out of range.
  SineFit sineFitAt(int index) {
    RangeError.checkValidIndex(index, frequencies, "index");
    final f = frequencies[index], p = power[index];
    final a = _a[index], b = _b[index];
    final amp = sqrt(a * a + b * b);
    final phase = atan2(a, b) / (2 * pi);
    final offset =
        fitMean ? weightedMean - a * _c[index] - b * _s[index] : weightedMean;
    final rms = sqrt(_yy * (1 - p));
    final ePhase = sqrt(2 / n) * rms / amp / (2 * pi);

    double eF = double.nan;
    if (index > 0 && index < frequencies.length - 1) {
      // Parabola y = k x² through the neighbors, origin at the peak
      double sxy = 0, sxx = 0;
      for (final j in [index - 1, index + 1]) {
        final dx = frequencies[j] - f;
        final x = dx * dx, yv = power[j] - p;
        sxy += x * yv;
        sxx += x * x;
      }
      final curvature = sxy / sxx;
      if (curvature < 0) eF = sqrt(-2 / n / curvature * (1 - p));
    }

    return SineFit(
      frequency: f,
      amplitude: amp,
      phase: phase,
      t0: _tMin - phase / f,
      offset: offset,
      power: p,
      residualRms: rms,
      amplitudeError: sqrt(2 / n) * rms,
      phaseError: ePhase,
      t0Error: ePhase / f,
      offsetError: sqrt(1 / n) * rms,
      frequencyError: eF,
    );
  }

  /// [power] in the normalization [norm].
  List<double> normalized(Normalization norm) {
    final chi0 = _yy * _wSum;
    final pMax = power[bestIndex];
    return List.unmodifiable([
      for (final p in power)
        switch (norm) {
          Normalization.zk => p,
          Normalization.horneBaliunas => (n - 1) / 2 * p,
          Normalization.cumming => (n - 3) / 2 * p / (1 - pMax),
          Normalization.chiSquared => chi0 * (1 - p),
          Normalization.wrms => sqrt(_yy * (1 - p)),
          Normalization.lnL => -0.5 * chi0 * (1 - p) - 0.5 * _sumLnVar,
          Normalization.dlnL => 0.5 * chi0 * p,
        },
    ]);
  }

  /// Probability that pure Gaussian noise gives a power of at least [p]
  /// ([power] normalization) at one given frequency:
  /// `(1 - p)^((N - 3) / 2)`.
  double probability(double p) => exp((n - 3) / 2 * log1p(-p));

  /// False-alarm probability of the power [p] ([power] normalization): the
  /// probability that pure Gaussian noise gives a power of at least [p]
  /// anywhere in the searched band, `1 - (1 - probability(p))^M` with
  /// M = [independentFrequencies].
  ///
  /// Small values indicate a significant signal. Computed without the
  /// linear approximation for small values used in some implementations.
  double fap(double p) =>
      -expm1(independentFrequencies * log1p(-probability(p)));

  /// False-alarm probability of the highest peak, `fap(best.power)`.
  double get bestFap => fap(power[bestIndex]);

  /// The power ([power] normalization) at which the false-alarm
  /// probability is [fap]; powers above it are more significant. The
  /// inverse of [fap]. Throws an [ArgumentError] unless `0 < fap < 1`.
  double powerLevel(double fap) {
    if (!(fap > 0 && fap < 1)) {
      throw ArgumentError.value(fap, "fap", "must be > 0 and < 1");
    }
    final prob = -expm1(log1p(-fap) / independentFrequencies);
    return -expm1(2 / (n - 3) * log(prob));
  }
}
