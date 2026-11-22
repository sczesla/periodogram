import 'package:test/test.dart';

import '../example/periodogram_example.dart' as example;

void main() {
  test('example runs', () {
    expect(() => example.main(), prints(contains("Second period: 2.93 days")));
  });

  test('simulated data: 40 points in nights outside the gap', () {
    final d = example.simulatedData();
    expect(d.t, hasLength(40));
    expect(d.t.every((t) => t < 5 || t >= 7), isTrue);
    expect(d.err.every((e) => e >= 0.1 && e <= 0.5), isTrue);
  });
}
