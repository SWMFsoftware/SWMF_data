# Cassini observations for the Enceladus end-to-end nightly tests

Input of the AMPS Table tests `Enceladus--E3-INMS`, `Enceladus--UVIS-2010`
and their long-tier variants (`srcEnceladus/test/end_to_end/README.md` in
AMPS). The test decks set `observation_data_root` to this directory and
compare the AMPS model with these products after every sampling window.

| File | Content | Used by |
|---|---|---|
| `e3_h2o_sticking_corrected.csv` (+ `.metadata.json`) | INMS E3 H2O number density along the body-fixed (IAU_ENCELADUS) spacecraft track: `time_ca_s`, `x_m`, `y_m`, `z_m`, `density_m3`, `sigma_m3` | E3 INMS test |
| `uvis_2010_e10_solar_optical_depth.csv` (+ `.metadata.json`) | UVIS 2010 E10 solar occultation: tangent positions, time from minimum tangent altitude, optical-depth proxy and its uncertainty | UVIS 2010 test |
| `uvis_2010_e10_solar_observer_range.csv` (+ `.metadata.json`) | row-aligned ray table (signed observer coordinate, body-fixed ray direction) for `ray_geometry=finite_observer_to_sun` | UVIS 2010 test |
| `SHA256SUMS` | checksums of every file here (`sha256sum -c SHA256SUMS`) | |

## Provenance

Copied unchanged on 2026-11-28 from the processed Cassini validation release
`/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/`. The
`.metadata.json` file next to each CSV records its source, units, frame and
processing. Do not edit these files: a new reduction belongs in the
validation release first and is then copied here.

## Limitations

- E3 is a digitized, sticking-corrected INMS profile, not an independent
  raw-count inversion.
- The UVIS optical depth is a shape proxy. The comparison fits one
  column-to-optical-depth coefficient and one baseline; the absolute scale is
  tested only against the published peak column (9.0e19 +/- 2.3e19 m^-2, set
  in the test deck).
- In the 2010 ray table, 308 of 1,133 rows have a nonnegative observer
  coordinate without the full Cassini position; their rays start at the
  tangent point (`srcEnceladus/RUNTIME_INPUT.md`, UVIS contract). The UVIS
  comparison is therefore diagnostic, not a publication-qualified column.
