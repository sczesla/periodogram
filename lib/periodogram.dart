/// Period search in unevenly sampled data: the generalized Lomb-Scargle
/// periodogram (`gls()`) with errors and floating mean, the classical
/// Lomb-Scargle periodogram (`lombScargle()`), false-alarm probabilities,
/// the best-fit sine (`SineFit`), and the spectral window
/// (`spectralWindow()`).
library;

export 'src/gls.dart' show gls, lombScargle, Periodogram, Normalization;
export 'src/sine_fit.dart';
export 'src/spectral_window.dart';
