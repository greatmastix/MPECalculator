# MPE Calculator

**▶ Open the calculator: [greatmastix.github.io/MPECalculator](https://greatmastix.github.io/MPECalculator/)**

Eye-safety calculator for **show lasers**, **high-power moving lights** and **projectors**. It works out the maximum permissible exposure, hazard distances and the reduction needed at the audience, from data-sheet or measured values. It runs entirely in the browser: no server, no cookies, no tracking. The whole calculation lives in the page address, so a link reproduces it exactly.

## What it calculates

**Show lasers** (IEC 60825-1:2007 / EU Directive 2006/25/EC, or IEC 60825-1:2014 / ICNIRP 2013; 400–700 nm)
- Exposure at the audience as a multiple of the MPE, averaged over a 7 mm pupil
- Nominal ocular hazard distance (NOHD), static and scanned
- Power reduction (% of full, mW, optical density) or diverging lens (dioptres) needed at a given distance
- Scanned patterns described the way Pangolin BEYOND builds them: scan rate (kpps), points per frame, pattern size, anchor points, repeats and grating beams. Every spot of the pattern (line crossings, line ends, beam dots) is checked against the single-pulse and average limits. Measured pulse widths can be entered instead.
- What a stopped beam delivers within the scan-fail reaction time
- Per-colour beams, diameter and divergence in 1/e², 1/e or FWHM, and a beam-profile choice (multimode diode, flat-top or Gaussian)

**Moving lights** (IEC 62471:2006 / EN 62471 / EU 2006/25/EC, or IEC 62471-5:2015 / ICNIRP 2013)
- Risk group at the audience for a 0.25 s look (RG2 or lower is acceptable, RG3 is not) and the RG3 hazard distance
- Dimming needed, including PWM dimming, where each pulse runs at full brightness
- Zoom slider: peak intensity across the zoom from the photometric data
- Retinal thermal and blue-light hazard with spectral constants per light engine (white LED, laser-phosphor, discharge, HMI, tungsten, xenon, RGB mixes), negative beam angles (crossover beams), and how long someone can stare into the fixture

**Projectors** (same standards as moving lights)
- Hazard distance and risk group from ANSI lumens, throw ratio (zoom slider over the lens range) and native chip aspect ratio
- Panasonic PT-RQ / PT-RZ presets with their lens lists; the colour temperature setting is taken into account

**All**
- Dazzle distances (NODD) after Williamson & McLin (2018)
- Presets with data-sheet values for common lasers, fixtures and projectors. Values a preset provides are greyed out; **Unlock** lets you replace them with measured values
- A compact text of the results for a risk assessment

## Using it

Open `index.html` in a browser, or serve the folder with any static web server. There is no build step and nothing to install. To share a calculation, use **Share link** or copy the address.

## How results relate to reality

The limits follow the named standards. Where a data sheet does not pin something down, the calculator uses typical real-world values from published measurements: beam profile of multimode diodes, galvo dwell at line ends, light-engine spectra, hot spots on the front lens and intensity fall-off with zoom. The *How this is calculated* section on the page lists every model and its sources.

With these typical values, the calculator reproduces the hazard distances that manufacturers publish within a few percent (Ayrton Cobra and Kyalami, Claypaky Skylos, Martin MAC Viper XIP; see `test/reference.js`).

Results are only as good as the inputs. Data-sheet power and divergence are often module ratings or averages. Measured values of the unit in use (power through a 7 mm aperture, beam size, pulse width, lux at a distance) give the most reliable answer. Audience scanning additionally needs a working scan-fail safeguard and whatever your local rules require.

## Tests

```bash
node test/reference.js
```

The script checks textbook laser values and the published hazard distances of four moving lights.

## Files

| File | What |
|---|---|
| `index.html` | the page: inputs, results, chart, tables, notes |
| `calc.js` | all maths as pure functions (also loads in Node) |
| `presets.js` | data-sheet values of the presets, with source links |
| `test/reference.js` | reference checks |

## License

MIT, see [LICENSE](LICENSE). Provided without warranty; using the results for a safety assessment is your responsibility.
