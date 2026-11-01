import numpy as np


class NBodyGravity:
    """General N-body gravity solver in 3D with force softening."""

    def __init__(self, masses, g_const=4 * np.pi**2, eps=1e-6):
        self.m = np.array(masses, dtype=float)
        self.N = len(self.m)
        self.G = float(g_const) # g_const is in units where G=4*pi^2 for solar system simulations (At 1 AU, 1 yr, 1 Msun).
        self.eps = float(eps) # Softening length to avoid singularities in force calculations.

    def accelerations(self, r):
        """Compute accelerations for all bodies, shape (N, 3)."""
        r = np.array(r, dtype=float)
        a = np.zeros_like(r)

        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    continue
                rij = r[j] - r[i]
                dist2 = np.dot(rij, rij) + self.eps**2
                a[i] += self.G * self.m[j] * rij / (dist2 ** 1.5)
        return a

    def total_energy(self, r, v):
        """Total energy = kinetic + potential."""
        kinetic = 0.5 * np.sum(self.m[:, None] * v * v)

        potential = 0.0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                rij = r[j] - r[i]
                dist = np.sqrt(np.dot(rij, rij) + self.eps**2)
                potential += -self.G * self.m[i] * self.m[j] / dist

        return kinetic + potential

    def total_angular_momentum(self, r, v):
        """Total angular momentum vector."""
        total_L = np.zeros(3) #Lx, Ly, Lz
        for i in range(self.N):
            total_L += self.m[i] * np.cross(r[i], v[i])
        return total_L

    def f(self, r, v):
        """Return dr/dt and dv/dt for given positions and velocities."""
        a = self.accelerations(r)
        return v, a

    def step_leapfrog(self, r, v, h):
        """Kick-Drift-Kick leapfrog step."""
        a0 = self.f(r, v)[1]
        v_half = v + 0.5 * h * a0
        r_new = r + h * v_half
        a1 = self.f(r_new, v_half)[1]
        v_new = v_half + 0.5 * h * a1
        return r_new, v_new

    def step_rk4(self, r, v, h):
        """Classic 4th-order Runge-Kutta step."""
        k1_r, k1_v = self.f(r, v)
        k2_r, k2_v = self.f(r + 0.5 * h * k1_r, v + 0.5 * h * k1_v)
        k3_r, k3_v = self.f(r + 0.5 * h * k2_r, v + 0.5 * h * k2_v)
        k4_r, k4_v = self.f(r + h * k3_r, v + h * k3_v)

        r_new = r + (h / 6.0) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)
        v_new = v + (h / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)
        return r_new, v_new

    def integrate(self, r0, v0, times, method="leapfrog"):
        """Integrate for an arbitrary increasing time array."""
        t = np.array(times, dtype=float)
        r = np.array(r0, dtype=float)
        v = np.array(v0, dtype=float)

        r_hist = [r.copy()]
        v_hist = [v.copy()]

        for k in range(1, len(t)):
            h = t[k] - t[k - 1]
            if method == "leapfrog":
                r, v = self.step_leapfrog(r, v, h)
            elif method == "rk4":
                r, v = self.step_rk4(r, v, h)

            r_hist.append(r.copy())
            v_hist.append(v.copy())

        return t, np.array(r_hist), np.array(v_hist)
