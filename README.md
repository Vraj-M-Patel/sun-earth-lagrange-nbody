# Sun–Earth Lagrange Points (N-body)

Computational physics demo: a massless satellite is placed at the Sun–Earth **L2** and **L4** points and integrated with a 3-body gravity solver. Trajectories are shown in the **rotating frame** that turns with Earth.

## Physics takeaway

| Case | Expected behavior |
|------|-------------------|
| L2, no kick | Unstable saddle — probe eventually leaves |
| L4, no kick | Stable for Earth–Sun mass ratio — stays put |
| L2 + small \(v_x\) kick | Leaves sooner |
| L4 + same kick | Short tadpole arc near L4 |

Units: AU, years, solar masses, with \(G = 4\pi^2\) so a 1 AU circular orbit has a period of ~1 year.

## Run

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter notebook lagrange_points.ipynb
```

Or from a Python shell / script using `nbodygravity.NBodyGravity` the same way as in the notebook.

## Repo contents

- `nbodygravity.py` — N-body integrator (leapfrog / RK4)
- `lagrange_points.ipynb` — L2 / L4 setup, integration, rotating-frame plots
- `requirements.txt` — Python dependencies

## Notes

This is the circular restricted three-body problem (Sun + Earth + test mass). Real spacecraft near L2 (e.g. JWST) need station-keeping; Trojan asteroids occupy the stable triangular points of a much heavier planet.
