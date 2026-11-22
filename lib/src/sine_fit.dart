import 'dart:math';

/// A sine fitted at one frequency of a periodogram:
///
///     y(t) = amplitude * sin(2 pi frequency (t - t0)) + offset
///
/// Obtained from `Periodogram.best` or `Periodogram.sineFitAt`. The errors
/// are approximations that assume well sampled data (Zechmeister & Kürster
/// 2009).
class SineFit {
  /// Creates a sine fit; usually obtained from a periodogram instead.
  const SineFit({
    required this.frequency,
    required this.amplitude,
    required this.phase,
    required this.t0,
    required this.offset,
    required this.power,
    required this.residualRms,
    required this.amplitudeError,
    required this.phaseError,
    required this.t0Error,
    required this.offsetError,
    required this.frequencyError,
  });

  /// Frequency, in units of 1 / unit of the times.
  final double frequency;

  /// Semi-amplitude (>= 0).
  final double amplitude;

  /// Phase in cycles (between -0.5 and 0.5) at the first time of the data.
  final double phase;

  /// A time of zero phase (ascending zero crossing of the sine).
  final double t0;

  /// Offset: fitted (floating mean) or the weighted mean of the data.
  final double offset;

  /// Power at [frequency] (normalization of Zechmeister & Kürster 2009).
  final double power;

  /// Weighted rms of the residuals.
  final double residualRms;

  /// Approximate error of [amplitude].
  final double amplitudeError;

  /// Approximate error of [phase], in cycles.
  final double phaseError;

  /// Approximate error of [t0].
  final double t0Error;

  /// Approximate error of [offset].
  final double offsetError;

  /// Approximate error of [frequency] from the curvature of the power
  /// peak; NaN if not at a local maximum inside the frequency grid.
  final double frequencyError;

  /// Period, `1 / frequency`.
  double get period => 1 / frequency;

  /// Approximate error of [period]; NaN if [frequencyError] is NaN.
  double get periodError => frequencyError / (frequency * frequency);

  /// The sine at the times [t].
  List<double> evaluate(List<double> t) => [
        for (final ti in t)
          amplitude * sin(2 * pi * frequency * (ti - t0)) + offset,
      ];

  @override
  String toString() => "SineFit(frequency: $frequency, period: $period, "
      "amplitude: $amplitude, phase: $phase, t0: $t0, offset: $offset, "
      "power: $power)";
}
