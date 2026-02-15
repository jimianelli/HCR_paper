# HCR_paper

Concise applied‑fisheries manuscript built from the HCR_eval simulation, centered on a **risk‑calibrated HCR** under steepness uncertainty.

## Main file
- `HCR_paper.qmd` — manuscript source (HTML/PDF via Quarto)

## Supporting analysis
- `R/pollock_simulation.R` — model and helper functions
- `ref.bib` — bibliography

## Render
```sh
quarto render HCR_paper.qmd
```

## Notes
- The paper uses a single, risk‑calibrated Ftarget applied across steepness scenarios.
- Output is configured for a concise applied fisheries journal format.
