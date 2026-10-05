/* MPE maths. Pure functions, SI units throughout (m, rad, W, s) unless a name says otherwise.
 *
 * Lasers: IEC 60825-1:2014 (= ICNIRP 2013, Health Phys. 105(3):271-295) with ISH1:2017, or the 2007 edition
 *         (= EU Directive 2006/25/EC Annex II). 400-700 nm, point source (C6 = 1), 7 mm limiting aperture.
 * Lights: ICNIRP 2013 (Health Phys. 105(1):74-96) as used by IEC 62471-5:2015, or IEC 62471:2006 / EN 62471 /
 *         EU 2006/25/EC Annex I (ICNIRP 1997 weighting). Retinal thermal and blue light.
 */
(function (root) {
  'use strict';

  const AP = 0.007, R_AP = AP / 2;       // limiting aperture (dark-adapted pupil) and its radius
  const A_AP = Math.PI * AP * AP / 4;
  const EDITIONS = {
    '2014': { ti: 5e-6, hShort: 2e-3 },    // Ti and the short-pulse MPE (J/m²) below it
    '2007': { ti: 1.8e-5, hShort: 5e-3 },
  };
  const T2 = 10;                         // s, point source: pulses for C5 are counted within T2

  // ---------- generic solvers ----------
  // Largest x in [lo, hi] where the (non-increasing) f still exceeds 1; log bisection.
  function solveFalling(f, lo, hi) {
    if (f(lo) <= 1) return 0;
    if (f(hi) > 1) return Infinity;
    for (let i = 0; i < 70; i++) {
      const mid = Math.sqrt(lo * hi);
      if (f(mid) > 1) lo = mid; else hi = mid;
    }
    return hi;
  }
  // Largest x in [lo, hi] where the (non-decreasing) f is still <= 1.
  function solveRising(f, lo, hi) {
    if (f(lo) > 1) return 0;
    if (f(hi) <= 1) return Infinity;
    for (let i = 0; i < 70; i++) {
      const mid = Math.sqrt(lo * hi);
      if (f(mid) > 1) hi = mid; else lo = mid;
    }
    return lo;
  }
  // Farthest x in [lo, hi] where f exceeds 1, for an f that may jump: walk a log grid in from hi, then bisect.
  function lastAbove(f, lo, hi, perDecade, extra) {
    if (f(hi) > 1) return Infinity;
    const n = Math.ceil(Math.log10(hi / lo) * (perDecade || 40)), xs = [];
    for (let i = 0; i < n; i++) xs.push(lo * Math.pow(hi / lo, i / n));
    (extra || []).forEach(x => { if (x > lo && x < hi) xs.push(x * (1 - 1e-9), x * (1 + 1e-9)); });
    xs.sort((a, b) => a - b);
    let prev = hi;
    for (let i = xs.length - 1; i >= 0; i--) {
      const x = xs[i];
      if (f(x) > 1) {
        let a = x, b = prev;
        for (let k = 0; k < 60; k++) { const m = Math.sqrt(a * b); if (f(m) > 1) a = m; else b = m; }
        return b;
      }
      prev = x;
    }
    return 0;
  }
  // Largest x such that f stays <= 1 on all of [lo, x]: walk a log grid out from lo (plus known breakpoints), then bisect.
  function firstAbove(f, lo, hi, breaks, perDecade) {
    const n = Math.ceil(Math.log10(hi / lo) * (perDecade || 60));
    const xs = [];
    for (let i = 0; i <= n; i++) xs.push(lo * Math.pow(hi / lo, i / n));
    (breaks || []).forEach(b => { if (b > lo && b < hi) xs.push(b * (1 - 1e-9), b); });
    xs.sort((a, b) => a - b);
    if (f(xs[0]) > 1) return 0;
    for (let i = 1; i < xs.length; i++) if (f(xs[i]) > 1) {
      let a = xs[i - 1], b = xs[i];
      for (let k = 0; k < 60; k++) { const m = Math.sqrt(a * b); if (f(m) > 1) b = m; else a = m; }
      return a;
    }
    return Infinity;
  }

  // ---------- laser limits (IEC 60825-1, 400-700 nm, point source) ----------
  // Thermal MPE as radiant exposure, J/m².
  function hThermal(t, ed) {
    const e = EDITIONS[ed] || EDITIONS['2014'];
    if (t < e.ti) return e.hShort;
    if (t < 10) return 18 * Math.pow(t, 0.75);
    return 10 * t;
  }
  function c3(nm) { return nm <= 450 ? 1 : Math.pow(10, 0.02 * (nm - 450)); }   // C_B
  // Photochemical MPE as radiant exposure, J/m² (only 400-600 nm and t >= 10 s).
  function hPhoto(t, nm) {
    if (t < 10 || nm > 600) return Infinity;
    return t <= 100 ? 100 * c3(nm) : c3(nm) * t;
  }
  // Exposure / MPE for irradiances `items` [{e (W/m², 7 mm average), nm}] lasting t seconds. Thermal adds over
  // all colours; photochemical adds each colour against its own limit; the lower limit counts.
  function mpeRatio(t, items, ed) {
    let eTot = 0, photo = 0;
    for (const it of items) {
      eTot += it.e;
      const hp = hPhoto(t, it.nm);
      if (hp !== Infinity) photo += it.e / (hp / t);
    }
    const thermal = eTot / (hThermal(t, ed) / t);
    return { hr: Math.max(thermal, photo), thermal, photo, e: eTot };
  }

  // ---------- laser beam geometry ----------
  // b = { a, phi, kA, kPhi, fA, fPhi, dpt, gauss }. kA/kPhi convert the data-sheet definitions to d63 (1/e, 63 %
  // power) for a Gaussian: 1/e² -> 1/√2, FWHM -> 1/√ln2. fA/fPhi convert them to the full width of a flat-top
  // beam: 1/e² and FWHM -> 1, 1/e -> 1/√(1-1/e). dpt = power of a diverging lens at the aperture (waist there).
  const DEF_GAUSS = { e2: Math.SQRT1_2, e1: 1, fwhm: 1 / Math.sqrt(Math.LN2) };
  const DEF_FLAT = { e2: 1, e1: 1 / Math.sqrt(1 - Math.exp(-1)), fwhm: 1 };
  // d63 diameter at range r (second-moment growth with the lens: near field a(1 + rD), far field phi·r).
  function d63(r, b) {
    return Math.hypot(b.kA * b.a * (1 + r * (b.dpt || 0)), b.kPhi * b.phi * r);
  }
  // Area of the overlap of two discs with radii R1, R2 whose centres are s apart.
  function overlap(s, R1, R2) {
    if (s >= R1 + R2) return 0;
    if (s <= Math.abs(R1 - R2)) return Math.PI * Math.min(R1, R2) ** 2;
    const c1 = Math.min(1, Math.max(-1, (s * s + R1 * R1 - R2 * R2) / (2 * s * R1)));
    const c2 = Math.min(1, Math.max(-1, (s * s + R2 * R2 - R1 * R1) / (2 * s * R2)));
    return R1 * R1 * Math.acos(c1) + R2 * R2 * Math.acos(c2) -
      0.5 * Math.sqrt(Math.max(0, (-s + R1 + R2) * (s + R1 - R2) * (s - R1 + R2) * (s + R1 + R2)));
  }
  // Fraction of a flat-top beam (an even disc of diameter dNear spread by an even angular disc that has grown
  // to dFar) passing the centred 7 mm aperture. For any even profile shape, the equal-area round discs are the
  // worst case (Riesz rearrangement), so this bounds square, rectangular and elliptical flat-top beams too.
  function flatFrac(dNear, dFar) {
    const R1 = Math.max(dNear, 1e-9) / 2, R2 = Math.max(dFar, 1e-9) / 2, top = Math.min(R_AP, R1 + R2);
    const n = 128, ds = top / n;
    let f = 0;
    for (let i = 0; i < n; i++) { const s = (i + 0.5) * ds; f += overlap(s, R1, R2) * 2 * Math.PI * s * ds; }
    return Math.min(1, f / (Math.PI * R1 * R1 * Math.PI * R2 * R2));
  }
  // Fraction of the beam's power through the 7 mm aperture, by profile (b.profile):
  //  'gauss' single-mode Gaussian; 'flat' round flat-top, the worst case of any flat profile;
  //  'multimode' (default) a typical multimode diode beam, flat along the slow axis and near-Gaussian along the
  //  fast axis: halfway (geometric mean) between Gaussian and flat-top, which matches the 1.1-1.4× over Gaussian
  //  such beams give in a 7 mm aperture.
  function apertureFraction(r, b) {
    const g = -Math.expm1(-((AP / d63(r, b)) ** 2));
    if (b.profile === 'gauss') return g;
    const f = Math.max(g, flatFrac(b.fA * b.a * (1 + r * (b.dpt || 0)), b.fPhi * b.phi * r));
    return b.profile === 'flat' ? f : Math.sqrt(g * f);
  }
  // Geometry for one colour: its own diameter/divergence if given, else the common one.
  function lineGeom(L, ln) {
    if (ln.a === undefined && ln.phi === undefined) return L;
    return Object.assign({}, L, { a: ln.a !== undefined ? ln.a : L.a, phi: ln.phi !== undefined ? ln.phi : L.phi });
  }
  // Per colour: irradiance averaged over the 7 mm aperture (W/m²) and d63 at range r.
  function lineData(r, L) {
    let shared = null;
    return L.lines.map(ln => {
      const g = lineGeom(L, ln);
      if (g === L) { if (!shared) shared = { f: apertureFraction(r, L), d: d63(r, L) }; return { nm: ln.nm, e: ln.p * shared.f / A_AP, d: shared.d }; }
      return { nm: ln.nm, e: ln.p * apertureFraction(r, g) / A_AP, d: d63(r, g) };
    });
  }
  // Exposure / MPE for a beam resting on the eye for t seconds.
  function staticRatio(r, t, L) { return mpeRatio(t, lineData(r, L), L.ed); }

  // ---------- one pass of a Gaussian beam across the aperture ----------
  // A Gaussian beam (d63) crossing the 7 mm aperture through its centre at speed v delivers the energy
  // e·LEFF/v (e = peak 7 mm-averaged irradiance) in a pulse whose FWHM is FWHM/v. Both lengths come from the
  // exact aperture integral, tabulated once against u = d63 / 7 mm (IEC: energy over the FWHM duration).
  function erf(x) {                       // Abramowitz & Stegun 7.1.26, |error| < 1.5e-7
    const s = x < 0 ? -1 : 1, z = Math.abs(x), t = 1 / (1 + 0.3275911 * z);
    return s * (1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z));
  }
  const PASS = (function () {
    const N = 64, us = [], le = [], fw = [];
    // power through the aperture with the beam centre at offset x (up to a constant), y = R sin(phi)
    const pAt = (x, w) => {
      let s = 0;
      for (let i = 0; i < N; i++) {
        const ph = -Math.PI / 2 + (i + 0.5) * Math.PI / N, c = R_AP * Math.cos(ph), y = R_AP * Math.sin(ph);
        s += c * Math.exp(-(y * y) / (w * w)) * (erf((c - x) / w) + erf((c + x) / w));
      }
      return s;
    };
    const total = w => {                  // integral of pAt over all x, using ∫(erf+erf)dx = 4c
      let s = 0;
      for (let i = 0; i < N; i++) {
        const ph = -Math.PI / 2 + (i + 0.5) * Math.PI / N, c = R_AP * Math.cos(ph), y = R_AP * Math.sin(ph);
        s += c * Math.exp(-(y * y) / (w * w)) * 4 * c;
      }
      return s;
    };
    for (let i = 0; i <= 60; i++) {
      const u = Math.pow(10, -1.5 + 3 * i / 60), w = u * AP / 2;   // w = d63 / 2 = 1/e radius
      const p0 = pAt(0, w);
      let a = 0, b = R_AP + 6 * w;
      for (let k = 0; k < 40; k++) { const m = (a + b) / 2; if (pAt(m, w) > p0 / 2) a = m; else b = m; }
      us.push(Math.log(u)); le.push(Math.log(total(w) / p0)); fw.push(Math.log(a + b));
    }
    return { us, le, fw };
  })();
  function passLengths(d) {               // { leff, fwhm } in metres for a Gaussian beam of d63 = d
    const u = d / AP, lu = Math.log(u), n = PASS.us.length;
    if (lu <= PASS.us[0]) return { leff: AP, fwhm: AP };
    if (lu >= PASS.us[n - 1]) return { leff: Math.sqrt(Math.PI) * d / 2, fwhm: Math.sqrt(Math.LN2) * d };
    const i = Math.min(n - 2, Math.floor((lu - PASS.us[0]) / (PASS.us[1] - PASS.us[0])));
    const f = (lu - PASS.us[i]) / (PASS.us[i + 1] - PASS.us[i]);
    return { leff: Math.exp(PASS.le[i] + f * (PASS.le[i + 1] - PASS.le[i])), fwhm: Math.exp(PASS.fw[i] + f * (PASS.fw[i + 1] - PASS.fw[i])) };
  }

  // ---------- scanned beams (repetitive pulses) ----------
  // C5 for N pulses (point source). 2014: pulses longer than Ti give 1 (alpha <= 5 mrad); shorter ones give
  // 5·N^-0.25 (at least 0.4) when the time base exceeds 0.25 s and N > 600. 2007 / EU: N^-0.25.
  function c5(n, tp, t, ed) {
    if (tp > 0.25) return 1;
    if (ed === '2007') return Math.min(1, Math.pow(n, -0.25));
    if (tp > EDITIONS['2014'].ti) return 1;
    if (t <= 0.25 || n <= 600) return 1;
    return Math.min(1, Math.max(0.4, 5 * Math.pow(n, -0.25)));
  }
  // Ratio for one spot on the pattern: `items` per colour [{e, nm, tp (FWHM, s), teff (s at peak irradiance)}]
  // hit `rate` times per second. Rule 1: each pulse (with C5); rule 2: the average over t (covers every group).
  function spotRatio(t, items, rate, ed, nRate) {
    const nr = nRate || rate, nC5 = Math.max(1, ed === '2007' ? nr * t : nr * Math.min(t, T2));
    const evalWith = its => {
      let single = 0;
      for (const it of its) single += it.e * it.teff / (hThermal(it.tp, ed) * c5(nC5, it.tp, t, ed));
      const avg = mpeRatio(t, its.map(it => ({ e: it.e * Math.min(1, it.teff * rate + (it.extraDuty || 0)), nm: it.nm })), ed).hr;
      return { hr: Math.max(single, avg), single, avg };
    };
    let res = evalWith(items);
    // IEC 60825-1:2014 ISH1 6 a): a pulse shorter than Ti may be assessed as a Ti-long pulse of the same peak power.
    const ti = EDITIONS['2014'].ti;
    if (ed !== '2007' && items.some(it => it.tp <= ti)) {
      const alt = evalWith(items.map(it => it.tp <= ti ? Object.assign({}, it, { tp: ti * 1.0000001, teff: Math.max(it.teff, ti) }) : it));
      if (alt.hr < res.hr) res = Object.assign(alt, { ish1: true });
    }
    return res;
  }

  /* L.scan, one of:
   *  { kind: 'pattern', pattern: 'line2' | 'line1' | 'circle' | 'dots', angle (rad, pattern size), pps (points/s),
   *    points (per frame), anchors (points held at each line end), k (repeats drawn over the same spot per frame),
   *    g (beams following each other along the path: grating, splitter), dotPoints }
   *  { kind: 'measured', tp (s, FWHM through 7 mm at the audience), rate (pulses/s at the eye) }
   * Returns the spots to check: [{ name, rate, items }]. */
  function scanSpots(r, L, data) {
    const s = L.scan;
    if (s.kind === 'measured') {
      return [{ name: 'measured pulse', rate: s.rate, items: data.map(d => ({ e: d.e, nm: d.nm, tp: s.tp, teff: 1.064 * s.tp })) }];
    }
    const f = s.pps / s.points, k = s.k || 1, g = s.g || 1, td = (s.anchors || 0) / s.pps;
    if (s.pattern === 'dots') {
      const tDot = (s.dotPoints || 1) / s.pps;
      return [{ name: 'beam dot', rate: k * f, items: data.map(d => ({ e: d.e, nm: d.nm, tp: tDot, teff: tDot })) }];
    }
    const omega = { line2: 2, line1: 1, circle: Math.PI }[s.pattern] * s.angle * k * f;   // rad/s along the path
    const v = omega * r;
    const crossRate = (s.pattern === 'line2' ? 2 : 1) * k * f;
    const cross = data.map(d => { const p = passLengths(d.d); return { e: d.e, nm: d.nm, tp: p.fwhm / v, teff: p.leff / v }; });
    const spots = [{ name: 'line crossing', rate: crossRate * g, items: cross }];
    if (s.pattern !== 'circle') {
      // Line ends: the beam stops and turns. It rests for the anchor points plus about two points of galvo lag, and
      // the passes in and out merge into one pulse.
      const fc = k * f;
      spots.push({
        name: 'line end', rate: fc, nRate: fc + (g - 1) * crossRate,   // other beams cross the end spot too
        items: data.map(d => {
          const w = Math.max(d.d, AP), tc = w / v;
          const tE = td + 2 / s.pps + 2 * tc;
          return { e: d.e, nm: d.nm, tp: tE, teff: tE, extraDuty: (g - 1) * crossRate * cross[data.indexOf(d)].teff };
        }),
      });
    }
    return spots;
  }

  /* L = { lines: [{nm, p, a?, phi?}], a, phi, kA, kPhi, fA, fPhi, gauss, dpt, ed, scan }
   * Exposure / MPE at range r for exposure time t. With scanning, every spot of the pattern is checked and the
   * worst one counts, never more than a beam resting on the eye. */
  function laserRatio(r, t, L) {
    const data = lineData(r, L);
    const st = mpeRatio(t, data, L.ed);
    const out = { hr: st.hr, e: st.e, static: st.hr, limit: st.photo > st.thermal ? 'photochemical' : 'thermal', spots: [] };
    if (!L.scan) return out;
    let worst = null;
    for (const sp of scanSpots(r, L, data)) {
      const q = spotRatio(t, sp.items, sp.rate, L.ed, sp.nRate);
      const rec = Object.assign({ name: sp.name, rate: sp.rate, tp: Math.max(...sp.items.map(i => i.tp)), n: sp.rate * t,
        c5: c5(Math.max(1, L.ed === '2007' ? (sp.nRate || sp.rate) * t : (sp.nRate || sp.rate) * Math.min(t, T2)), Math.min(...sp.items.map(i => i.tp)), t, L.ed) }, q);
      if (q.ish1) { rec.c5 = 1; rec.tp = Math.max(rec.tp, EDITIONS['2014'].ti); }
      out.spots.push(rec);
      if (!worst || rec.hr > worst.hr) worst = rec;
    }
    out.worst = worst;
    if (worst.hr < st.hr) { out.hr = worst.hr; out.limit = worst.name + (worst.single >= worst.avg ? ', single pulse' : ', average'); }
    else out.limit = 'as a resting beam';
    return out;
  }

  function laserNOHD(t, L) {
    const f = r => laserRatio(r, t, L).hr;
    return L.scan ? lastAbove(f, 0.01, 1e6, 40) : solveFalling(f, 0.01, 1e6);
  }
  function laserMaxTime(r, L) {
    return solveRising(t => staticRatio(r, t, L).hr, 1e-9, 3e4);
  }
  // Power (dioptres, as a positive number) of a diverging lens at the aperture that meets the MPE at r.
  function laserNeededLens(r, t, L) {
    const f = dpt => laserRatio(r, t, Object.assign({}, L, { dpt })).hr;
    if (f(0) <= 1) return 0;
    return L.scan ? lastAbove(f, 1e-4, 1e4, 40) : solveFalling(f, 1e-4, 1e4);
  }

  // ---------- laser dazzle (Williamson & McLin: maximum dazzle exposure, MDE) ----------
  // Simplified MDE limits at 555 nm in W/m² (published in µW/cm²: 1 µW/cm² = 0.01 W/m²), by the width
  // of the field of view that is blocked and the ambient light: night 0.1, dusk 10, day 1000 cd/m².
  const MDE = [
    { deg: 2, name: 'Very low', night: 0.001e-2, dusk: 0.6e-2, day: 40e-2 },
    { deg: 10, name: 'Low', night: 0.04e-2, dusk: 30e-2, day: 2000e-2 },
    { deg: 20, name: 'Medium', night: 0.16e-2, dusk: 120e-2, day: 8000e-2 },
    { deg: 40, name: 'High', night: 0.6e-2, dusk: 450e-2, day: 30000e-2 },
  ];
  // Eye sensitivity the dazzle framework specifies: CIE 2008 2° luminous efficiency (CIE 170-2:2015,
  // CVRL linCIE2008v2e_5), 5 nm steps from 400 to 700 nm. It sees blue about 1.7 times brighter than the
  // 1924 curve that lux is defined with (V(450) 0.0647 against 0.038).
  const V8_TAB = [.002452, .004972, .00908, .01429, .02027, .02612, .03319, .04158, .05034, .05743, .06472, .07238, .08515, .106, .1299, .1535,
    .1788, .2065, .2379, .2851, .3484, .4278, .5205, .6206, .7181, .7946, .8576, .9071, .9545, .9814, .989, .9995,
    .9968, .9903, .9733, .9425, .8964, .8587, .8116, .7545, .6919, .627, .5584, .4896, .423, .3609, .2981, .2417,
    .1943, .1547, .1193, .0898, .06671, .049, .0356, .02554, .01808, .01262, .008661, .006028, .004196];
  function vLambda2008(nm) {
    const x = Math.min(Math.max(nm, 400), 700), i = Math.min(V8_TAB.length - 2, Math.floor((x - 400) / 5)), f = (x - 400) / 5 - i;
    return V8_TAB[i] + (V8_TAB[i + 1] - V8_TAB[i]) * f;
  }
  // Irradiance of a resting beam weighted by how bright each colour looks, W/m² (555 nm equivalent).
  function visualIrradiance(r, L) {
    return lineData(r, L).reduce((sum, d) => sum + d.e * vLambda2008(d.nm), 0);
  }
  // Nominal ocular dazzle distance: beyond it the beam no longer blocks that much of the view.
  function laserNODD(L, mde) { return solveFalling(r => visualIrradiance(r, L) / mde, 0.01, 1e7); }

  // ---------- spectral weightings ----------
  // CIE 1924 photopic V(lambda), 5 nm steps from 380 to 780 nm.
  const V_TAB = [3.9e-5, 6.4e-5, 1.2e-4, 2.17e-4, 3.96e-4, 6.4e-4, 1.21e-3, 2.18e-3, 4e-3, 7.3e-3, .0116, .01684, .023, .0298, .038, .048,
    .06, .0739, .09098, .1126, .13902, .1693, .20802, .2586, .323, .4073, .503, .6082, .71, .7932, .862, .91485, .954, .9803, .99495, 1,
    .995, .9786, .952, .9154, .87, .8163, .757, .6949, .631, .5668, .503, .4412, .381, .321, .265, .217, .175, .1382, .107, .0816, .061,
    .04458, .032, .0232, .017, .01192, .00821, .005723, .004102, .002929, .002091, .001484, .001047, .00074, .00052, .0003611, .0002492,
    .0001719, .00012, 8.48e-5, 6e-5, 4.24e-5, 3e-5, 2.12e-5, 1.499e-5];
  function vLambda(nm) {
    if (nm < 380 || nm > 780) return 0;
    const i = Math.min(V_TAB.length - 2, Math.floor((nm - 380) / 5)), f = (nm - 380) / 5 - i;
    return V_TAB[i] + (V_TAB[i + 1] - V_TAB[i]) * f;
  }
  // Blue-light hazard weighting B(lambda), ICNIRP 2013 Table 2 (0 above 700 nm).
  const B_TAB = [[380, .01], [385, .0125], [390, .025], [395, .05], [400, .1], [405, .2], [410, .4], [415, .8],
    [420, .9], [425, .95], [430, .98], [435, 1], [440, 1], [445, .97], [450, .94], [455, .9], [460, .8],
    [465, .7], [470, .62], [475, .55], [480, .45], [485, .4], [490, .22], [495, .16], [500, .1]];
  function bLambda(nm) {
    if (nm < 380) return 0.01;
    if (nm > 700) return 0;
    if (nm >= 600) return 0.001;
    if (nm >= 500) return Math.pow(10, (450 - nm) / 50);
    const i = Math.min(B_TAB.length - 2, Math.floor((nm - 380) / 5));
    const [x0, y0] = B_TAB[i], [x1, y1] = B_TAB[i + 1];
    return y0 + (y1 - y0) * (nm - x0) / (x1 - x0);
  }
  // Retinal thermal weighting R(lambda): ICNIRP 2013 Table 2 (380-1400 nm).
  function r13(nm) {
    if (nm < 380 || nm > 1400) return 0;
    if (nm < 435) return bLambda(nm);
    if (nm <= 700) return 1;
    if (nm <= 1050) return Math.pow(10, (700 - nm) / 500);
    if (nm <= 1150) return 0.2;
    if (nm <= 1200) return 0.2 * Math.pow(10, 0.02 * (1150 - nm));
    return 0.02;
  }
  // The older weighting of ICNIRP 1997 / IEC 62471:2006 / EU 2006/25/EC Table 1.3: 10 × B in the blue.
  function r97(nm) {
    if (nm < 380 || nm > 1400) return 0;
    if (nm < 500) return 10 * bLambda(nm);
    if (nm <= 700) return 1;
    if (nm <= 1050) return Math.pow(10, (700 - nm) / 500);
    return 0.2;
  }

  // ---------- light engines ----------
  // Spectral constants per lumen: KR13 / KR97 = lumens per R-weighted watt (ICNIRP 2013 / 1997 weighting),
  // KB = blue-light-weighted watts per lumen, V8 = CIE 2008 / CIE 1924 brightness ratio (for dazzle).
  function blackBody(T) {                 // full black body 380-1400 nm, including its infrared
    let sv = 0, sb = 0, s13 = 0, s97 = 0, s8 = 0;
    for (let nm = 380; nm <= 1400; nm += 2) {
      const l = nm * 1e-9, p = 1 / (Math.pow(l, 5) * (Math.exp(0.014387769 / (l * T)) - 1));
      sv += p * vLambda(nm); sb += p * bLambda(nm); s13 += p * r13(nm); s97 += p * r97(nm);
      if (nm >= 400 && nm <= 700) s8 += p * vLambda2008(nm);
    }
    return { KR13: 683 * sv / s13, KR97: 683 * sv / s97, KB: sb / (683 * sv), V8: s8 / sv };
  }
  function interp(tab, x) {               // linear, extended with the end slopes
    let i = 0;
    while (i < tab.length - 2 && x > tab[i + 1][0]) i++;
    const [x0, y0] = tab[i], [x1, y1] = tab[i + 1];
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }
  // Typical values from CIE 015:2018 LED-B1..B5, 13 measured LED / discharge moving heads (DGUV BoSS
  // spectroradiometer tables) and blue-laser + garnet phosphor models; see the notes on the page.
  const LED_KR13 = [[2700, 314], [6500, 314], [8000, 300], [9000, 292], [10000, 285], [12000, 273]];
  const LED_KR97 = [[2700, 190], [4100, 132], [5100, 111], [6500, 93], [9000, 76], [10000, 72], [12000, 64]];
  const ENGINES = {
    led: { name: 'White LED (blue-pumped phosphor)' },
    ledv: { name: 'White LED, violet-pumped, high-CRI or unknown' },
    laserph: { name: 'Laser-phosphor' },
    discharge: { name: 'Discharge lamp in a moving head (cold-light reflector)' },
    hmi: { name: 'Open-face HMI / MSR, no IR filter' },
    tungsten: { name: 'Tungsten / halogen' },
    xenon: { name: 'Xenon short arc, no IR filter (estimate)' },
    rgb: { name: 'RGB / RGBW LED or RGB laser (colour mix)' },
  };
  function engine(type, cct) {
    const T = Math.min(Math.max(cct || 6500, 1800), 20000), bb = blackBody(T);
    const kr13 = Math.max(250, Math.min(314, interp(LED_KR13, T))), kr97 = Math.max(55, Math.min(190, interp(LED_KR97, T)));
    switch (type) {
      case 'led': return { KR13: kr13, KR97: kr97, KB: bb.KB * 0.9, V8: bb.V8 };
      case 'ledv': return { KR13: 245, KR97: kr97, KB: bb.KB * 1.45, V8: bb.V8 };
      case 'laserph': return { KR13: Math.min(300, kr13), KR97: kr97, KB: bb.KB * 0.9, V8: bb.V8 };
      case 'discharge': return { KR13: 270, KR97: 78, KB: bb.KB, V8: bb.V8 };
      case 'hmi': return { KR13: Math.min(200, bb.KR13), KR97: Math.min(70, bb.KR97), KB: bb.KB * 1.4, V8: bb.V8 };
      case 'tungsten': return bb;
      case 'xenon': return { KR13: 130, KR97: Math.min(60, bb.KR97), KB: bb.KB * 1.4, V8: bb.V8 };   // unverified placeholder
      case 'rgb': return { KR13: 140, KR97: 24, KB: 3.7e-3, V8: 1.12 };   // measured Claypaky Xtylos at full RGB, blue-heavy bound
      default: return { KR13: kr13, KR97: kr97, KB: bb.KB, V8: bb.V8 };
    }
  }

  // ---------- broadband lights: geometry ----------
  function intensityFromLux(lux, dist) { return lux * dist * dist; }
  // Peak intensity at zoom angle z from known points [[angle, cd], ...] (the narrowest zoom first), fitted to the
  // published intensity-vs-zoom data of 36 fixtures (Ayrton, Martin, JB-Lighting, Robe, Claypaky, Elation, Chauvet).
  // With only the narrowest point: flat to 1.2× its angle, then angle^-1.7 (median error 8 %, on the high side).
  // Between two known points the curve bows above a straight log-log line; the midpoint between that line and
  // the bounding curve (flat to 1.2×, then angle^-1.45 from below, angle^-2.2 from above) hits the published
  // middle points at a median of 1.01.
  const ZOOM = { c: 1.2, n1: 1.7, n: 1.45, m: 2.2 };
  function zoomIntensity(points, z) {
    const pts = points.slice().sort((a, b) => a[0] - b[0]);
    const lo = pts.filter(p => p[0] <= z).pop() || pts[0], hi = pts.find(p => p[0] >= z);
    if (!hi || hi === lo || hi[0] === lo[0]) return lo[1] * Math.min(1, Math.pow(ZOOM.c * lo[0] / z, ZOOM.n1));
    const line = lo[1] * Math.pow(lo[0] / z, Math.log(lo[1] / hi[1]) / Math.log(hi[0] / lo[0]));
    const bound = Math.min(lo[1] * Math.min(1, Math.pow(ZOOM.c * lo[0] / z, ZOOM.n)), hi[1] * Math.pow(hi[0] / z, ZOOM.m));
    return Math.sqrt(line * Math.max(line, bound));
  }
  // Luminous flux spread evenly over a cone with the given full angle (rad).
  function intensityFromLumens(lm, angle) { return lm / (2 * Math.PI * (1 - Math.cos(angle / 2))); }
  // Distance from the lens at which a converging beam (negative beam angle) is narrowest.
  // dlens = front lens diameter, w = beam diameter at that crossover, angle = full beam angle (rad).
  function crossover(dlens, w, angle) {
    return Math.sqrt(Math.max(dlens * dlens - w * w, 0)) / (2 * Math.tan(Math.abs(angle) / 2));
  }
  // Tangent of the half angle the apparent source fills at range r. Normal beam: the lit front lens,
  // but each point of the lens only radiates into the beam, so close to the fixture the eye sees just
  // a patch of the lens and the angle cannot exceed the beam's spread (M.cap = tan of half of it).
  // Converging beam (M.xo > 0): the crossover, a real image of diameter M.dsrc at M.xo in front of the
  // lens, seen through the lens before the crossover and directly after it; never larger than the lens.
  function halfTan(r, M) {
    if (!M.xo) return Math.min(M.dsrc / (2 * r), M.cap || Infinity);
    return Math.min(M.dlens / (2 * r), M.dsrc / (2 * Math.abs(r - M.xo)));
  }

  /* M = { I (cd, peak, after dimming), Ifull (cd, at full output), dsrc, xo, dlens, cap,
   *       KR13, KR97, KB, V8, peak (hot-spot factor of the lit lens, 1 = even), std ('2013' | '2006'),
   *       pwm (Hz, PWM dimming; 0 = current dimming), duty (0-1) }
   * Returns exposure / limit for the retinal thermal and the blue-light hazard. */
  function alphaMax(t) { return t < 625e-6 ? 0.005 : t < 0.25 ? 0.2 * Math.sqrt(t) : 0.1; }
  function lightRatio(r, t, M) {
    const area = Math.PI * M.dsrc * M.dsrc / 4, lv = M.I / area;   // luminance, cd/m²
    const h = halfTan(r, M), alpha = 2 * Math.atan(h), p = M.peak || 1;
    // luminance averaged over the acceptance angle g; a lens with hot spots is brighter than its mean by up to p,
    // but no average can exceed all the light concentrated inside g
    const avgOver = (l, g) => l * Math.min(p, (alpha / g) ** 2);
    let thermal, blue, pulse = 0;
    if (M.std === '2006') {
      // IEC 62471:2006 / EU 2006/25/EC: L_R <= 5e4 / (C t^0.25) up to 10 s, 2.8e4 / C after; C = alpha clamped to
      // 1.7-100 mrad; field of view 1.7 mrad up to 0.25 s, 11·√(t/10) mrad to 10 s, then 11 mrad.
      const g = t <= 0.25 ? 0.0017 : t <= 10 ? 0.011 * Math.sqrt(t / 10) : 0.011;
      const Cs = Math.min(Math.max(alpha, g, 0.0017), 0.1);
      thermal = avgOver(lv, g) / M.KR97 / (t <= 10 ? 5e4 / (Cs * Math.pow(t, 0.25)) : 2.8e4 / Cs);
      if (alpha >= 0.011) blue = avgOver(lv, 0.011) * M.KB / (t <= 1e4 ? 1e6 / t : 100);
      else {                                  // small source: irradiance, 100 J/m² (normal viewing, at most 100 s)
        const eB = lv * M.KB * Math.PI * alpha * alpha / 4;
        blue = eB * Math.min(t, 100) / 100;
      }
    } else {
      // ICNIRP 2013 / IEC 62471-5: 2e4 / (alpha t^0.25) up to 0.25 s, 2.8e4 / alpha after; continuous light is
      // averaged over 11 mrad and alpha is not taken below that (IEC 62471-5, 62471-7 ISH1).
      const aEff = Math.min(Math.max(alpha, 0.011), alphaMax(t));
      thermal = avgOver(lv, 0.011) / M.KR13 / (t >= 0.25 ? 2.8e4 / aEff : 2e4 / (aEff * Math.pow(t, 0.25)));
      const gph = t <= 100 ? 0.011 : t <= 1e4 ? 0.0011 * Math.sqrt(t) : 0.11;
      blue = avgOver(lv, gph) * M.KB / (t <= 1e4 ? 1e6 / t : 100);
    }
    // PWM dimming: every pulse runs at full peak radiance, so each pulse must meet the single-pulse limit. For
    // extended sources (alpha > 5 mrad) at more than 5 Hz, ICNIRP 2013 lowers that limit by n^-0.25 (at least
    // 0.4, or 0.2 above alpha_max, none from 100 mrad).
    // PWM dimming: every pulse runs at full peak radiance, so each pulse must meet the single-pulse limit.
    // Pulses of 0.25 s or longer are continuous light at full brightness. Above 5 Hz and for alpha > 5 mrad,
    // ICNIRP 2013 also lowers the single-pulse limit by n^-0.25 (at least 0.4, or 0.2 above alpha_max, none
    // from 100 mrad); at 5 Hz or less there is no such reduction.
    if (M.pwm > 0 && M.duty < 1) {
      const tp = M.duty / M.pwm, lFull = M.Ifull / area;
      if (tp >= 0.25) pulse = lightRatio(r, Math.min(tp, t), Object.assign({}, M, { I: M.Ifull, duty: 1, pwm: 0 })).hr;
      else {
        const n = M.pwm * Math.min(Math.max(t, 0.25), 100), aMaxP = alphaMax(tp);
        const k = M.pwm <= 5 || alpha <= 0.005 || alpha >= 0.1 ? 1 : Math.max(Math.pow(n, -0.25), alpha <= aMaxP ? 0.4 : 0.2);
        if (M.std === '2006') {
          const Cs = Math.min(Math.max(alpha, 0.0017), 0.1);
          pulse = avgOver(lFull, 0.0017) / M.KR97 / (5e4 / (Cs * Math.pow(tp, 0.25)) * k);
        } else {
          const aEff = Math.min(Math.max(alpha, 0.005), aMaxP);
          pulse = avgOver(lFull, 0.005) / M.KR13 / (2e4 / (aEff * Math.pow(tp, 0.25)) * k);
        }
      }
    }
    const ev = Math.PI * lv * h * h / (1 + h * h);      // on-axis illuminance of an evenly lit disc
    const hr = Math.max(thermal, blue, pulse);
    return { hr, thermal, blue, pulse, alpha, ev, lv,
      limit: hr === pulse ? 'PWM pulse' : thermal >= blue ? 'retinal thermal' : 'blue light' };
  }
  // Range at which the apparent source looks largest, so exposure and illuminance peak: at the lens for a
  // normal beam (search from 1 cm), and for a crossover where the lens edge and the crossover look the
  // same size. Every exposure here grows with that size, so it falls on both sides of this point.
  function peakRange(M) { return M.xo ? M.xo / (1 + M.dsrc / M.dlens) : 0.01; }
  // Farthest range up to `far` at which f (exposure / limit) is above 1.
  function farthestAbove(f, M, far) {
    const lo = peakRange(M);
    if (f(far) > 1) return Infinity;
    if (f(lo) <= 1) return 0;
    return solveFalling(f, lo, far);
  }
  // Ranges where the apparent source crosses an angle at which a limit changes (alpha_max, 100, 11, 5, 1.7 mrad).
  function alphaBreaks(M) {
    const a = [0.1, 0.011, 0.005, 0.0017];
    if (M.pwm > 0 && M.duty < 1) a.push(alphaMax(M.duty / M.pwm));
    const xs = a.map(x => M.dsrc / (2 * Math.tan(x / 2)));
    if (M.xo) a.forEach(x => xs.push(M.xo + M.dsrc / (2 * Math.tan(x / 2)), M.dlens / (2 * Math.tan(x / 2))));
    xs.push(peakRange(M));
    return xs;
  }
  // The PWM pulse rule steps with the source size, so walk a grid in from far away instead of bisecting.
  function lightHazardDistance(t, M) { return lastAbove(r => lightRatio(r, t, M).hr, 0.01, 1e6, 100, alphaBreaks(M)); }
  // Dazzle distance for a light. Lux is the light weighted by the 1924 eye-sensitivity curve (683 lm per W
  // at 555 nm); the dazzle framework uses the 2008 curve, so lux / 683 is scaled by V8 (about 1.05 for white).
  function lightNODD(M, mde) { return farthestAbove(r => lightRatio(r, 1, M).ev / 683 * (M.V8 || 1) / mde, M, 1e7); }
  // Longest stare at r that stays within the limit (the older standard is not monotonic in t, so search outward).
  function lightMaxTime(r, M) { return firstAbove(t => lightRatio(r, t, M).hr, 0.01, 3e4, [0.25, 10, 100, 1e4]); }
  // Largest duty (dimmer fraction of full output) that keeps the exposure at r within the limit for time t.
  function lightMaxDuty(r, t, M) {
    const f = d => lightRatio(r, t, Object.assign({}, M, { I: M.Ifull * d, duty: d })).hr;
    if (f(1) <= 1) return 1;
    return solveRising(f, 1e-9, 1);
  }

  const api = {
    ZOOM, zoomIntensity, AP, A_AP, EDITIONS, DEF_GAUSS, DEF_FLAT, hThermal, hPhoto, c3, c5, d63, flatFrac, apertureFraction, lineData,
    staticRatio, passLengths, scanSpots, laserRatio, laserNOHD, laserMaxTime, laserNeededLens, MDE, vLambda2008,
    visualIrradiance, laserNODD, vLambda, bLambda, r13, r97, blackBody, ENGINES, engine, intensityFromLux,
    intensityFromLumens, crossover, lightRatio, peakRange, lightHazardDistance, lightNODD, lightMaxTime, lightMaxDuty,
    solveFalling, solveRising, lastAbove, firstAbove,
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.MPE = api;
})(typeof self !== 'undefined' ? self : this);
