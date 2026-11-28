# Fitted Enceladus plume source distributions (end-to-end nightly tests)

Input of the AMPS Table tests `Enceladus--E3-INMS`, `Enceladus--UVIS-2010`
and their long-tier variants (`srcEnceladus/test/end_to_end/README.md` in
AMPS). These are byte-identical copies; the AMPS test decks read the `.in`
files directly with `source_distribution_file`.

| File | Use |
|---|---|
| `tenishev2014_hybrid_8__independent_each_observation__E3__H2O__absolute.in` | E3 INMS test: H2O source rates fitted to the E3 INMS profile alone |
| `tenishev2014_hybrid_8__independent_each_observation__uvis_2010_e10_solar__H2O__absolute.in` | UVIS 2010 test: H2O source rates fitted to the 2010 E10 solar occultation alone |
| `*.manifest.json` | the exporter's manifest of each `.in` (per-source rates, fit provenance) |
| `tenishev2014_hybrid_8__resolved_geometry.json`, `tenishev2014_hybrid_8__stripe_geometry_authority.json` | source geometry named (by SHA-256) in the `.in` headers |
| `SHA256SUMS` | checksums of every file here (`sha256sum -c SHA256SUMS`) |

## Provenance

Source: `/data/vtenishe/cassini_validation_data/studies/fitting_compain_111126/h2o/tenishev2014_hybrid_8/amps_sources/`
(standalone Stage-0 fitting campaign of 2026-11-11), files
`independent_each_observation__{E3,uvis_2010_e10_solar}__H2O__absolute.{in,manifest.json}`,
`resolved_geometry.json` and `stripe_geometry_authority.json`, copied
unchanged on 2026-11-28 and renamed only by the `tenishev2014_hybrid_8__`
prefix. The original absolute paths in the `.in` header comments are
provenance; AMPS does not read them.

Source basis `tenishev2014_hybrid_8`: 8 point jets and 7 line stripes;
`mode=absolute` (physical rates in molecules/s, no `total_source_rate`
key allowed), `source_law=tenishev2010_point`.

| Fit | Total jet+stripe H2O rate | Fit statistic (campaign `fitting_study.md`) |
|---|---|---|
| E3 INMS, independent | 3.142482e28 s^-1 | reduced chi-square 6.440 (local density) |
| UVIS 2010 E10 solar, independent | 4.576082e27 s^-1 | shape reduced chi-square 3.876 (optical depth) |

## Limitations (from the fitting campaign)

- The fits use a steady, collisionless, gravity-free model; AMPS adds
  gravity and parent-loss photolysis, so its profiles differ from the fit.
- E3 is a digitized, sticking-corrected INMS profile, not a raw-count
  inversion.
- The UVIS fit constrains optical-depth shape; its absolute scale comes from
  one published peak column (9.0e19 +/- 2.3e19 m^-2), not an absolute
  calibration of every sample.
- The scaled design matrices are ill-conditioned; individual source rates are
  not uniquely determined.

Replacing a file changes the test's regression reference; regenerate the
references (AMPS `srcEnceladus/test/end_to_end/README.md`) in the same change.
