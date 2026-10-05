// Reference checks for calc.js. Run with: node test/reference.js
// Each case reproduces a textbook value or a hazard distance published by a manufacturer.
const C = require('../calc.js');
let fail = 0;
function check(name, got, want, tol) {
  const ok = Math.abs(got - want) <= tol * want;
  if (!ok) fail++;
  console.log((ok ? 'ok   ' : 'FAIL ') + name + ': ' + got.toFixed(2) + ' (expected ' + want + ' ± ' + tol * 100 + ' %)');
}
const K = C.DEF_GAUSS;

// Laser: 3 × 1 W, 5 mm, 1 mrad (1/e²), Gaussian, 0.25 s. Closed form sqrt(4P/(pi·MPE) - a²)/phi with d63 = 547.8 m.
const laser = { lines: [{ nm: 638, p: 1 }, { nm: 520, p: 1 }, { nm: 450, p: 1 }], a: 0.005, phi: 0.001, kA: K.e2, kPhi: K.e2,
  fA: 1, fPhi: 1, dpt: 0, ed: '2014', profile: 'gauss', scan: null };
check('NOHD, Gaussian 3 W 1 mrad', C.laserNOHD(0.25, laser), 547.8, 0.002);
check('MPE at 0.25 s (W/m²)', C.hThermal(0.25, '2014') / 0.25, 25.46, 0.002);

// Black-body spectral constant (ICNIRP 2013 R-weighting, 380-1400 nm): 187.7 lm/W at 6500 K.
check('Black body K_R at 6500 K', C.blackBody(6500).KR13, 187.7, 0.005);

// Moving lights against published hazard distances (0.25 s, RG2/RG3 boundary).
function light(lux, d, beamDeg, lensMm, engine, cct, std) {
  const I = lux * d * d, e = C.engine(engine, cct);
  return C.lightHazardDistance(0.25, Object.assign({ I, Ifull: I, dsrc: lensMm / 1000, xo: 0, dlens: lensMm / 1000,
    cap: Math.tan(beamDeg * Math.PI / 360) * 1.2, peak: 1.5, std, pwm: 0, duty: 1 }, e));
}
check('Ayrton Cobra, IEC 62471-5 (45 m)', light(1544000, 10, 0.6, 170, 'laserph', 6500, '2013'), 45, 0.06);
check('Claypaky Skylos, IEC 62471-5 (65 m)', light(28332, 100, 0.6, 300, 'laserph', 10000, '2013'), 65, 0.06);
check('Ayrton Kyalami, IEC 62471:2006 (76 m, US)', light(412000, 10, 1, 126, 'laserph', 9000, '2006'), 76, 0.06);
check('Martin MAC Viper XIP, EN 62471 (4.5 m)', light(99678, 5, 5.4, 150, 'led', 5800, '2006'), 4.5, 0.06);

process.exit(fail ? 1 : 0);
