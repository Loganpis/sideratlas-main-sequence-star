# Validation and remaining uncertainty

The numerical reports accompany the final pipeline. A successful numerical check does not make the assembled sphere a globally coupled MHD solution.

| Check | Result | Meaning |
|---|---:|---|
| Native GOL states retained | 71, 970–1630 s | Actual irregular timestamps bracket the film |
| C IV emissivity volume integral after stored averaging | maximum relative error 3.4 × 10⁻⁹ | Spatial averaging conserves calculated emission |
| UV depth samples, 48 vs 96 across eleven views/times | maximum 2.282% projected L1 | Oblique-view quadrature checked; original identical-geometry report retained with hash |
| UV emissivity averaging, 4×4 vs 2×2 | 1.90% raw L1, 0.526% after the declared Gaussian | Retained averaging is adequate at presentation resolution in the tested state |
| Photospheric sampling, 1 vs 4 regular rays per output pixel | 1.533% post-imaging RMS relative to disk mean | Native 1080 delivery sampling has been checked |
| Golden grade | maximum absolute linear-Y error 1.8 × 10⁻¹⁵ | Deliberate chromaticity change preserves pre-tone-map luminance |
| Cached source/radiation reproduction | completed successfully | Retrieved/cached inputs and calculated UV arrays match the manifest |
| Native UV alternate-output holdout | median 116% L1 | **Subcadence emissivity motion is not independently validated** |
| Four UV projected holdouts | 21.9–39.0% raw L1; 9.8–19.4% with the declared Gaussian | Fine UV changes cannot be described as precision continuous fluid trajectories |

The temporal holdouts intentionally omit real recorded outputs. They demonstrate the importance of retaining every selected state; they do not provide a proven bound on the remaining ten-second interpolation errors. The final animation linearly interpolates calculated emissivity and radiation, not a newly integrated velocity field.

The snapshot-325 source check finds real cool upflows above 2 Mm at roughly 30–45 km/s in the selected diagnostic class. The bright C IV phase is generally around 0.07–0.20 MK, with typical electron densities near 10⁹–10¹⁰ cm⁻³. Its velocity distribution includes downward and slowly moving material. The diagnostic finds **no cells meeting its additional >50 km/s, T<0.2 MK, density>10⁻¹³ kg m⁻³ cut at that state**. This is not a survey of the complete simulation, and no larger jet is fabricated to compensate. A bright UV strand must not automatically be called an ejected mass parcel.

Inherited optical tests include 0.035% RMS change from 17 to 33 transfer wavelengths, 1.58% RMS vertical formal-integration refinement, 2.28% RMS at an independent oblique angle, and 1.14% RMS for a 40-second optical midpoint interpolation. MURaM photospheric temporal holdouts average 1.29% RMS of the fixed quiet reference. These concern optical radiation; they must not be substituted for the larger UV errors.

The source's isolated ancillary Bifrost TEFF header anomaly is retained and flagged in the inherited source report. It never sets the emission temperature or exposure; local thermodynamic arrays supply the radiation. No native MHD mesh convergence or global energy transport test is claimed from these archive subsets.

The final export report verifies 960 frames, exactly 32 seconds, 1080 × 1920, H.264 High, Rec.709 limited YUV420p, complete decoding, nonidentical adjacent frames and the normalized PSF checks. Its recorded SHA256 identifies the delivered movie. Encoding follows one high-depth-to-H.264 conversion; source physics and display grading remain separate in the reproducible pipeline.

## Four-orbit revision

The 960 frame centers sample a continuous unwrapped azimuth at exactly 45°/s, or 1.5° per frame. The boundaries at 8, 16 and 24 seconds are ordinary consecutive frame pairs, not cuts or animation resets. Each eight-second interval contains exactly 240 frames. The analytic start/end camera transformation agrees to 1.5 × 10⁻¹⁵; actual quiet-surface and active-map samples return identically when physical time is held fixed. The source clocks and independent stellar rotation are not frozen in the film.

The new fixed chromatic LUT preserves pre-tone CIE Y to 1.8 × 10⁻¹⁵ in the tested full-resolution frame. Its amber lanes and warm-white highlights follow existing calculated intensity. One-versus-four-ray photospheric sampling changes the post-imaging result by 1.533% RMS relative to disk mean. This is presentation convergence, not source MHD convergence.

`validation/reused-inputs.json` records SHA256 hashes of the existing processed radiation inputs used here. `validation/scientific-checks.json` records the configuration-matched inherited UV quadrature and fresh camera/anchoring/grade tests. The final encoded-camera report checks the 240-frame orbit boundaries and four nonidentical evolving rotations. None of these checks removes the disclosed 2.5D geometry or native UV cadence limitation.

The original report behind the reused whole-orbit UV quadrature is retained as `validation/inherited-UV/whole-orbit-projection-original.json`; its SHA256 matches the fresh report's provenance record. Only its physical UV projection results are reused. Its unrelated historical palette test is not the golden film's grading or imaging model. A clean reproduction of the new validator recomputes the projection test directly.
