"""Smoke tests that pin crnpy functions to values and identities taken from the source papers.

These tests were written during the September 2026 audit of crnpy.py. Each one pins a function to a
value, identity or range taken from the manuscript that the function implements, so that a future
change which departs from the published equations is caught.

Sources
-------
Schrön et al. 2017, HESS 21, 5009-5030, Appendix A and Table A1 (weighting functions, D86 in cm,
    Fp with fixed 1013 mbar).
Schrön et al. 2018, WRR 54, 6441-6459, Eq. 3-6 and Table 1 (road correction,
    F2' = p2 - p3*theta_road - (p4 + theta_road)/(p5 + theta(N)), N_corr = N / C_road).
McJannet & Desilets 2023, WRR 59, e2022WR033889, Eq. 1, 8 and 10 and Table B1 (incoming flux
    correction N_corr = N * [tau*I/I_ref + 1 - tau]^-1, location factor, neutron monitor Rc values).
Smart & Shea 2008, Proc. 30th ICRC (Merida) 1, 733-736: world grid of vertical cutoff rigidities,
    epoch 1995.0, every 5 degrees in latitude and 15 degrees in longitude.
Baatz et al. 2015, WRR 51, 2030-2046: 0.9 % count reduction per kg m-2 dry biomass, about
    0.5 % per kg m-2 biomass water equivalent.
Jakobi et al. 2020, Front. Water 2:10, Eq. 7: sigma_N = s * sqrt(N).
Köhli et al. 2015, WRR 51, 5772-5790, Eq. 8: D86 ranges 15-83 cm.
Franz et al. 2012, WRR 48, W08515: z* = 5.8 / (rho_bd * tau + theta + 0.0829) in cm.
"""

import numpy as np
import pandas as pd
import pytest

import crnpy


# ---------------------------------------------------------------------------
# 1. Cutoff rigidity: Smart & Shea (2008) epoch 1995.0 world grid, checked against neutron monitors
# ---------------------------------------------------------------------------

def test_cutoff_rigidity_longitude_wraparound():
    # A longitude and the same longitude expressed on the 0-360 grid must give identical results.
    lat, lon = 39.68, -75.75
    assert crnpy.cutoff_rigidity(lat, lon) == pytest.approx(crnpy.cutoff_rigidity(lat, lon + 360), abs=1e-6)


# Neutron monitors with latitude, longitude and vertical cutoff rigidity (GV) from Table B1 of
# McJannet & Desilets (2023), Water Resources Research 59, e2022WR033889.
NEUTRON_MONITORS = [
    ("ALERT", 82.50, -62.33, 0), ("BARENT", 78.06, 14.22, 0), ("RESOLU", 74.69, -94.91, 0),
    ("VOSTOK", -78.47, 106.87, 0), ("TERA", -66.65, 140.00, 0.01), ("WILKES", -66.27, 110.53, 0.01),
    ("MRNY", -66.55, 93.02, 0.03), ("HEISS", 80.62, 58.05, 0.1), ("NEU3", -70.63, 8.26, 0.1),
    ("SOPO", -90.00, 0.00, 0.1), ("CHURCH", 58.75, -94.08, 0.21), ("MAWSON", -67.60, 62.88, 0.22),
    ("FSMT", 60.02, -111.93, 0.3), ("INVK", 68.36, -133.72, 0.3), ("JBGO", -74.60, 164.20, 0.3),
    ("MCMU", -77.90, 166.60, 0.3), ("NAIN", 56.55, -61.68, 0.3), ("PWNK", 54.98, -85.44, 0.3),
    ("THUL", 76.50, -68.70, 0.3), ("TIXIE", 71.36, 128.54, 0.48), ("KIRUNA", 67.83, 20.43, 0.54),
    ("GOOSE", 53.27, -60.40, 0.64), ("APTY", 67.57, 33.40, 0.65), ("SNAE", -71.67, -2.85, 0.73),
    ("OULU", 65.05, 25.47, 0.81), ("CALGAR", 51.08, -114.13, 1.08), ("SULPHU", 51.20, -115.60, 1.09),
    ("KERG", -49.35, 70.25, 1.14), ("SULIGY", 51.20, -115.60, 1.14), ("OTTAWA", 45.44, -75.68, 1.22),
    ("UPPSAL", 59.86, 17.62, 1.43), ("DURHAM", 43.10, -70.83, 1.58), ("YAKUTS", 62.01, 129.43, 1.65),
    ("CHICAG", 41.83, -87.67, 1.72), ("KINGST", -42.99, 147.29, 1.81), ("VICTOR", 48.42, -123.32, 1.86),
    ("HOBART", -42.88, 147.33, 1.88), ("SWARTH", 39.90, -75.35, 1.92), ("MGDN", 60.04, 151.05, 2.1),
    ("LEEDS", 53.80, -1.55, 2.2), ("KIEL", 54.34, 10.12, 2.36), ("KIEL2", 54.34, 10.12, 2.36),
    ("NEWK", 39.68, -75.75, 2.4), ("MOSC", 55.47, 37.32, 2.43), ("MCRL", 55.47, 37.32, 2.46),
    ("LONDON", 51.53, -0.10, 2.73), ("UTRECH", 52.10, 5.12, 2.76), ("NOVOSI", 54.48, 83.00, 2.87),
    ("DENVER", 39.67, -104.97, 2.91), ("HERSTM", 50.87, 0.33, 2.92), ("CLMX", 39.37, -106.18, 3),
    ("LARC", -62.20, -58.96, 3), ("LINDAU", 51.65, 10.13, 3), ("LINIGY", 51.65, 10.13, 3),
    ("HALLE", 51.48, 11.97, 3.07), ("DRBS", 50.10, 4.60, 3.18), ("KIEV", 50.43, 30.18, 3.57),
    ("IRK2", 52.37, 100.55, 3.64), ("IRKT", 52.47, 104.03, 3.64), ("LMKS", 49.20, 20.22, 3.84),
    ("MUNCHE", 48.20, 11.60, 4.14), ("ZUGSPI", 47.42, 10.98, 4.24), ("DALLAS", 32.98, -96.73, 4.35),
    ("HAFELE", 47.31, 11.38, 4.38), ("JUNG", 46.55, 7.98, 4.49), ("JUNG1", 46.55, 7.98, 4.49),
    ("BERKEL", 37.87, -122.27, 4.54), ("HRMS", -34.43, 19.23, 4.58), ("BURE", 44.63, 5.91, 5),
    ("USHUAI", -54.80, -68.32, 5.68), ("BKSN", 43.28, 42.69, 5.7), ("ROME", 41.00, 12.52, 6.32),
    ("AATB", 43.13, 76.55, 6.61), ("TBILIS", 41.43, 44.48, 6.73), ("CALM", 40.56, -3.16, 6.95),
    ("PTFM", -26.68, 27.09, 6.98), ("NANM", 40.37, 44.25, 7.1), ("BRISBA", -27.43, 153.08, 7.21),
    ("TASHKE", 41.20, 69.37, 7.5), ("MXCO", 19.33, -99.18, 8.28), ("ATHN", 37.97, 23.78, 8.53),
    ("TSMB", -19.20, 17.58, 9.15), ("BEIJIN", 39.08, 116.26, 10), ("MORIOK", 39.70, 141.13, 10.16),
    ("FUKUSH", 37.75, 140.48, 10.46), ("BUENOS", -34.60, -58.48, 10.63), ("ESOI", 33.30, 35.80, 10.75),
    ("SEOUL", 37.53, 126.93, 10.79), ("SNTIAG", -33.48, -70.71, 11), ("DJON", 36.24, 127.22, 11.2),
    ("CORDOB", -31.42, -64.19, 11.45), ("TOKYO", 35.75, 139.72, 11.63), ("HALESM", 20.72, -156.27, 12.91),
    ("HALIGY", 20.72, -156.25, 12.91), ("HUAN", -12.03, -75.33, 12.92), ("CHACAL", -16.32, -68.15, 13.1),
    ("MAKAPU", 21.30, -157.65, 13.23), ("KULA", 20.73, -156.33, 13.3), ("TIBET", 30.11, 90.53, 13.71),
    ("DARWIN", -12.43, 130.87, 14.09), ("AHMD", 23.01, 72.61, 15.94), ("PSNM", 18.59, 98.49, 16.8),
]


def test_cutoff_rigidity_against_neutron_monitor_table():
    errors = np.array([float(crnpy.cutoff_rigidity(lat, lon)) - rc for _, lat, lon, rc in NEUTRON_MONITORS])
    western = np.array([e for e, (_, _, lon, _) in zip(errors, NEUTRON_MONITORS) if lon < 0])
    assert np.mean(np.abs(errors)) < 0.25    # bilinear check on the grid gave 0.18 GV over 102 monitors
    assert np.mean(np.abs(western)) < 0.40   # 0.29 GV for the 38 monitors in the western hemisphere
    assert np.max(np.abs(errors)) < 1.8      # largest deviations are Buenos Aires and Cordoba, about 1.5 GV


def test_cutoff_rigidity_reference_stations():
    # NMDB values for well-known monitors in both hemispheres and on both sides of Greenwich.
    for lat, lon, rc in [(39.68, -75.75, 2.40), (39.37, -106.18, 3.00), (46.55, 7.98, 4.49), (37.97, 23.78, 8.53),
                         (-34.43, 19.23, 4.58), (18.59, 98.49, 16.80), (20.72, -156.25, 12.91), (-49.35, 70.25, 1.14)]:
        assert abs(float(crnpy.cutoff_rigidity(lat, lon)) - rc) < 0.5


# ---------------------------------------------------------------------------
# 2. Incoming flux: McJannet & Desilets branch must reduce to I/Iref when tau = 1
# ---------------------------------------------------------------------------

def _incoming_series():
    return np.array([1150.0, 1200.0, 1250.0]), 1200.0


def test_incoming_flux_hawdon_reduces_to_ratio_when_rc_equal():
    I, Iref = _incoming_series()
    fi = crnpy.correction_incoming_flux(I, Iref, Rc_method="Hawdonetal2014", Rc_site=4.5, Rc_ref=4.5)
    np.testing.assert_allclose(fi, I / Iref)


def test_incoming_flux_mcjannet_reduces_to_ratio_when_site_equals_reference():
    I, Iref = _incoming_series()
    # Identical site and reference give tau = 1, so the factor must equal the plain ratio I/Iref.
    fi = crnpy.correction_incoming_flux(I, Iref, Rc_method="McJannetandDesilets2023",
                                        Rc_site=4.5, Rc_ref=4.5, site_atmdepth=665.18, ref_atmdepth=665.18)
    np.testing.assert_allclose(fi, I / Iref)


def test_incoming_flux_mcjannet_matches_eq10():
    # McJannet & Desilets 2023: Eq. 1 N_corr = N * F with Eq. 10 F = [tau I/Iref + 1 - tau]^-1.
    # With the crnpy convention N_corr = N / fi, fi must equal tau (I/Iref - 1) + 1.
    I, Iref = _incoming_series()
    tau = crnpy.location_factor(963.47, 10.08, 665.18, 4.50)
    expected = tau * (I / Iref - 1) + 1
    fi = crnpy.correction_incoming_flux(I, Iref, Rc_method="McJannetandDesilets2023",
                                        Rc_site=10.08, Rc_ref=4.50, site_atmdepth=963.47, ref_atmdepth=665.18)
    np.testing.assert_allclose(fi, expected)
    # Direction check: a stronger incoming flux (I > Iref) must lower the corrected counts.
    assert fi[-1] > 1 and fi[0] < 1


def test_location_factor_against_mcjannet_table_b1():
    # McJannet & Desilets 2023, Table B1 lists for each monitor its atmospheric depth x, cutoff rigidity Rc
    # and K = 1/tau(x, Rc) from Eq. 8 (epsilon = 1). crnpy.location_factor(site, ref) returns
    # tau(site)/tau(ref), which must therefore equal K_ref / K_site.
    table_b1 = {  # code: (x g/cm2, Rc GV, K)
        "CLMX": (680.5, 3.00, 1.07),
        "JUNG": (665.2, 4.49, 1.36),
        "HALIGY": (714.8, 12.91, 3.53),
        "PSNM": (758.0, 16.80, 4.67),
    }
    for site, ref in [("JUNG", "CLMX"), ("HALIGY", "JUNG"), ("PSNM", "JUNG"), ("PSNM", "CLMX")]:
        xs, rs, ks = table_b1[site]
        xr, rr, kr = table_b1[ref]
        assert crnpy.location_factor(xs, rs, xr, rr) == pytest.approx(kr / ks, rel=0.02)


# ---------------------------------------------------------------------------
# 3. Schrön 2017 weighting: horizontal weights must depend on air humidity
# ---------------------------------------------------------------------------

def _two_profile_setup():
    # Two profiles at 5 m and 150 m, three depths each, different soil moisture.
    profiles = np.array([1, 1, 1, 2, 2, 2])
    distances = np.array([5.0, 5.0, 5.0, 150.0, 150.0, 150.0])
    depth = np.array([0.05, 0.15, 0.25, 0.05, 0.15, 0.25])
    theta = np.array([0.10, 0.12, 0.14, 0.30, 0.32, 0.34])
    p = np.full(6, 1013.25)
    return profiles, distances, depth, theta, p


def test_nrad_weight_schron_depends_on_air_humidity():
    profiles, distances, depth, theta, p = _two_profile_setup()
    dry, _ = crnpy.nrad_weight(np.full(6, 2.0), theta, distances, depth, profiles=profiles, rhob=1.4,
                               p=p, Hveg=0)
    humid, _ = crnpy.nrad_weight(np.full(6, 25.0), theta, distances, depth, profiles=profiles, rhob=1.4,
                                 p=p, Hveg=0)
    assert dry != pytest.approx(humid, abs=1e-6)


def test_schron_2017_weighting_functions_match_excel_supplement():
    # Cached values from the authors' CRNS-weighting.xlsx (Schrön et al. 2017 supplement) with
    # inputs air humidity x = 10 g/m3, soil moisture y = 0.02 m3/m3, p = 1013.25 hPa, Hveg = 0.
    # The workbook evaluates W at r_xl * Fp (cell $M$2). crnpy evaluates W at the scaled distance
    # r / Fp, so a sample placed at r_xl * Fp**2 reproduces the workbook value exactly.
    # Only r <= 50 m rows are used: the workbook's WrB uses b25 = 33474 instead of the 3347.4
    # given in Table A1 of the paper (a typo in the supplement).
    x, y, p = 10.0, 0.02, 1013.25
    Fp = 0.4922 / (0.86 - np.exp(-p / 1013.25))  # = 1.000161426 (cell M2)
    r_xl = np.array([0.5, 2.0, 3.1])                                   # rows 14, 29, 41
    W_xl = np.array([332414.75704399636, 154190.34096952606, 82522.228017326124])
    # One sample per profile at the surface: Wd = 1, so theta_P = theta and, with uniform theta,
    # the field average is fixed at y and the returned horizontal weights are W_r(r*, x, y).
    profiles = np.array([0, 1, 2])
    theta_avg, w = crnpy.nrad_weight(x, np.full(3, y), r_xl * Fp ** 2, np.zeros(3),
                                     profiles=profiles, rhob=1.0, p=p, Hveg=0)
    theta_P, r_stars, Wrs = w
    assert theta_avg == pytest.approx(y)
    np.testing.assert_allclose(r_stars, r_xl * Fp, rtol=1e-12)
    np.testing.assert_allclose(Wrs, W_xl, rtol=1e-9)


# ---------------------------------------------------------------------------
# 4. Biomass correction default parameter (Baatz et al. 2015)
# ---------------------------------------------------------------------------

def test_correction_bwe_with_baatz_coefficient():
    # Baatz et al. 2015: about 0.5 % reduction per kg m-2 BWE, i.e. r2/N0 ~ 0.005.
    r2_N0 = 6.4 / 1210
    assert abs(crnpy.correction_bwe(420, 51.5, r2_N0) - 574) < 5


def test_correction_bwe_default_reproduces_baatz_table3():
    # Baatz et al. 2015, Table 3, site 1: BWE = 51.5 kg m-2, N_epih = 420 cph, N_epihv = 574 cph,
    # using the paper's r2/N0 = 6.4/1210 (Sect. 3.3) as the library default.
    assert abs(crnpy.correction_bwe(420, 51.5) - 574) < 5
    corrected = crnpy.correction_bwe(1000.0, 25.0)
    assert np.isfinite(corrected) and 1000.0 < corrected < 1300.0


# ---------------------------------------------------------------------------
# 5. Road correction (Schrön et al. 2018, Eq. 4-6 and Table 1)
# ---------------------------------------------------------------------------

def _c_road_schron2018(theta_N, theta_road, w, r):
    # Published Eq. 4-6 with Table 1 (F1, F2' and F3 rows) of Schrön et al. 2018.
    p0, p1 = 0.42, 0.50
    p2, p3, p4, p5 = 1.06, 4.00, 0.16, 0.39
    p6, p7, p8, p9 = 0.94, 1.10, 2.70, 0.01
    F1 = p0 * (1 - np.exp(-p1 * w))
    F2 = p2 - p3 * theta_road - (p4 + theta_road) / (p5 + theta_N)
    F3 = p6 * np.exp(-p7 * w ** (-p8) * r ** 4) + (1 - p6) * np.exp(-p9 * r)
    return 1 + F1 * F2 * F3


def test_correction_road_matches_schron2018():
    N, theta_N, theta_road, w, r = 1000.0, 0.30, 0.12, 3.0, 0.0
    expected = N / _c_road_schron2018(theta_N, theta_road, w, r)  # about 946 counts
    assert crnpy.correction_road(N, theta_N, w, road_distance=r, theta_road=theta_road) == pytest.approx(expected, rel=1e-6)


def test_correction_road_reduces_counts_for_dry_road():
    corrected = crnpy.correction_road(1000.0, 0.30, 3.0, road_distance=0.0, theta_road=0.05)
    assert corrected < 1000.0
    # Paper requirement 4: no correction for zero road width, and the effect vanishes far from the road.
    assert crnpy.correction_road(1000.0, 0.30, 0.0, road_distance=0.0) == pytest.approx(1000.0)
    assert crnpy.correction_road(1000.0, 0.30, 3.0, road_distance=200.0) == pytest.approx(1000.0, rel=1e-2)


# ---------------------------------------------------------------------------
# 6. Sensing depth units and pressure scaling
# ---------------------------------------------------------------------------

def test_sensing_depth_franz_2012_is_in_centimetres():
    # Franz et al. 2012: z* = 5.8 / (rho_bd * tau + theta + 0.0829), result in cm.
    vwc, bd, wlat = 0.20, 1.4, 0.03
    expected_cm = 5.8 / (bd * wlat + vwc + 0.0829)  # ~17.9 cm
    z = crnpy.sensing_depth(vwc, 1013.25, 1013.25, bd, wlat, method="Franz_2012")
    assert z == pytest.approx(expected_cm, rel=1e-9)
    assert 10 < z < 90  # cm, not m


def test_sensing_depth_schron_2017_is_in_centimetres():
    # Köhli et al. 2015 / Schrön et al. 2017: D86 between about 15 cm (wet) and 83 cm (dry).
    dry = crnpy.sensing_depth(0.02, 1013.25, 1013.25, 1.4, 0.0, dist=[0.0], method="Schron_2017")[0]
    wet = crnpy.sensing_depth(0.45, 1013.25, 1013.25, 1.4, 0.0, dist=[0.0], method="Schron_2017")[0]
    assert 40 < dry < 90
    assert 10 < wet < 25


def test_sensing_depth_schron_2017_with_lattice_water():
    # Schrön 2017 Eq. 2: theta is the volumetric water equivalent theta_sm + theta_lw, so lattice water
    # (g/g) enters as bulk_density * Wlat. Hand-computed D86 directly below the sensor for
    # vwc = 0.20, bd = 1.4, Wlat = 0.03 (theta = 0.242) at standard pressure: 20.16 cm.
    got = crnpy.sensing_depth(0.20, 1013.25, 1013.25, 1.4, 0.03, dist=[0.0], method="Schron_2017")[0]
    assert got == pytest.approx(20.16, abs=0.05)


def test_sensing_depth_methods_return_similar_values():
    # Both methods describe the depth containing 86 % of the detected neutrons, so they should agree
    # to within about 25 % and both must become shallower as the soil gets wetter.
    bd, wlat = 1.4, 0.03
    previous = (np.inf, np.inf)
    for vwc in [0.05, 0.10, 0.20, 0.30]:
        z_schron = crnpy.sensing_depth(vwc, 1013.25, 1013.25, bd, wlat, dist=[0.0], method="Schron_2017")[0]
        z_franz = crnpy.sensing_depth(vwc, 1013.25, 1013.25, bd, wlat, method="Franz_2012")
        assert 0.8 < z_schron / z_franz < 1.25
        assert z_schron < previous[0] and z_franz < previous[1]
        previous = (z_schron, z_franz)


# ---------------------------------------------------------------------------
# 7. Counting uncertainty (Jakobi et al. 2020, Eq. 7)
# ---------------------------------------------------------------------------

def test_uncertainty_counts_std_follows_jakobi_eq7():
    N, fp, fw, fi = 1000.0, 1.1, 0.95, 1.05
    s = fw / (fp * fi)
    assert crnpy.uncertainty_counts(N, "std", fp=fp, fw=fw, fi=fi) == pytest.approx(s * np.sqrt(N))


def test_uncertainty_counts_cv_is_independent_of_corrections():
    N = 1000.0
    cv_plain = crnpy.uncertainty_counts(N, "cv")
    cv_corr = crnpy.uncertainty_counts(N, "cv", fp=1.1, fw=0.95, fi=1.05)
    assert cv_plain == pytest.approx(1 / np.sqrt(N))
    assert cv_corr == pytest.approx(cv_plain)


def test_uncertainty_vwc_reduces_to_first_order_propagation():
    # Jakobi et al. 2020 Eq. 10 is a 3rd-order Taylor expansion; for large N it must approach the
    # first-order propagation |d theta/dN| * sigma_N * rho_bd (Eq. 8 and 11). Hand-computed: 0.008703 m3/m3.
    N, N0, bd = 10000.0, 15000.0, 1.4
    a0, a1 = 0.0808, 0.372
    sigma_N = np.sqrt(N)
    first_order = a0 * N0 / (N - a1 * N0) ** 2 * sigma_N * bd
    got = crnpy.uncertainty_vwc(N, N0, bd)
    assert got == pytest.approx(0.008703, abs=1e-5)
    assert got == pytest.approx(first_order, rel=5e-3)


# ---------------------------------------------------------------------------
# 8. Atmospheric corrections and conversions
# ---------------------------------------------------------------------------

def test_atmospheric_depth_against_mcjannet_table_b1():
    # McJannet & Desilets 2023, Table B1: elevation (m), latitude and atmospheric depth x (g/cm2).
    for elev, lat, x in [(57, 82.50, 1023.7), (3570, 46.55, 665.2), (3400, 39.37, 680.5),
                         (2565, 18.59, 758.0), (4300, 30.11, 606.4), (33, -49.35, 1028.8)]:
        assert crnpy.atmospheric_depth(elev, lat) == pytest.approx(x, abs=0.1)


def test_abs_humidity_textbook_values():
    # Saturation vapour pressure from Campbell & Norman (1998) Eq. 3.8; absolute humidity at 20 C and 100 % RH
    # is 17.3 g/m3 and at 30 C and 50 % RH is 15.2 g/m3 (standard psychrometric tables).
    assert crnpy.abs_humidity(100, 20) == pytest.approx(17.27, abs=0.05)
    assert crnpy.abs_humidity(50, 30) == pytest.approx(15.16, abs=0.05)


def test_correction_pressure_zreda_eq5():
    # Zreda et al. 2012 Eq. 5: fp = exp((P0 - P)/L) and counts are divided by fp, so a pressure above the
    # reference (more air mass, fewer neutrons) must give fp < 1 and raise the corrected counts.
    # Hawdon et al. 2014 use beta = 0.0075 mb-1, i.e. L = 1/beta = 133.3 mb.
    assert crnpy.correction_pressure(1000.0, 1000.0, 133.3) == pytest.approx(1.0)
    fp = crnpy.correction_pressure(1010.0, 1000.0, 133.3)
    assert fp == pytest.approx(0.9277, abs=1e-3)
    assert 1000.0 / fp > 1000.0


def test_correction_humidity_rosolem_2013():
    # Rosolem et al. 2013 (Andreasen et al. 2017 Eq. 3): fw = 1 + 0.0054 (A - Aref), counts multiplied by fw.
    assert crnpy.correction_humidity(10.0, 0.0) == pytest.approx(1.054)
    assert crnpy.correction_humidity(5.0, 5.0) == pytest.approx(1.0)


def test_counts_to_vwc_desilets_and_hawdon():
    # Desilets et al. 2010: gravimetric water content at N/N0 = 0.6 with a0, a1, a2 defaults is 0.2394 g/g.
    assert crnpy.counts_to_vwc(600.0, 1000.0, 0.0, 0.0, 1.0) == pytest.approx(0.2394, abs=1e-4)
    # Hawdon et al. 2014 Eq. 7: theta_v = (theta_g - w_lat - w_SOM) * rho_bd
    theta_g = crnpy.counts_to_vwc(600.0, 1000.0, 0.0, 0.0, 1.0)
    assert crnpy.counts_to_vwc(600.0, 1000.0, 0.03, 0.01, 1.4) == pytest.approx((theta_g - 0.03 - 0.01) * 1.4)


# ---------------------------------------------------------------------------
# 9. Exponential filter (Albergel et al. 2008, Eqs. 4 and 6)
# ---------------------------------------------------------------------------

def test_exp_filter_recursion_matches_albergel():
    # Hand-computed with K1 = 1, K_n = K_{n-1} / (K_{n-1} + exp(-dt/T)), SWI_n = SWI_{n-1} + K_n (ms_n - SWI_{n-1}).
    sm = np.array([0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3])
    expected = [0.1, 0.1, 0.1, 0.191011, 0.237730, 0.263515, 0.278318, 0.287008]
    np.testing.assert_allclose(crnpy.exp_filter(sm, T=2), expected, atol=1e-5)


def test_exp_filter_tracks_input_for_very_short_time_scale():
    sm = np.array([0.1, 0.2, 0.3, 0.25])
    np.testing.assert_allclose(crnpy.exp_filter(sm, T=1e-6), sm, atol=1e-9)


def test_exp_filter_continues_across_missing_values():
    # Albergel et al. 2008: the gain uses the time elapsed since the last available observation, so a gap must not
    # turn the rest of the series into NaN. Hand-computed with dt = 2 across the gap.
    result = crnpy.exp_filter(np.array([0.1, np.nan, 0.3, 0.3]), T=2)
    assert np.isnan(result[1])
    np.testing.assert_allclose(result[[0, 2, 3]], [0.1, 0.246212, 0.275610], atol=1e-5)


# ---------------------------------------------------------------------------
# 10. Data handling functions: shapes, NaN handling and documented options
# ---------------------------------------------------------------------------

def test_remove_incomplete_intervals_keeps_first_row_by_default():
    ts = pd.date_range("2021-01-01 00:00", periods=4, freq="h")
    df = pd.DataFrame({"timestamp": ts, "counts": [1, 2, 3, 4]})
    assert len(crnpy.remove_incomplete_intervals(df.copy(), "timestamp", 3600)) == 4
    assert len(crnpy.remove_incomplete_intervals(df.copy(), "timestamp", 3600, remove_first=True)) == 3
    # A row following a 2-hour gap is an incomplete interval and is removed
    df_gap = pd.DataFrame({"timestamp": ts[[0, 1, 3]], "counts": [1, 2, 4]})
    assert len(crnpy.remove_incomplete_intervals(df_gap, "timestamp", 3600)) == 2


def test_fill_missing_timestamps_adds_missing_rows():
    ts = pd.to_datetime(["2021-01-01 00:00", "2021-01-01 03:00"])
    df = pd.DataFrame({"timestamp": ts, "counts": [1.0, 4.0]})
    filled = crnpy.fill_missing_timestamps(df, timestamp_col="timestamp", freq="h")
    assert len(filled) == 4
    assert filled["counts"].isna().sum() == 2


def test_total_raw_counts_fills_missing_detector_with_mean_of_others():
    counts = pd.DataFrame({"N1": [100.0, 110.0], "N2": [np.nan, 90.0]})
    np.testing.assert_allclose(crnpy.total_raw_counts(counts).values, [200.0, 200.0])
    # A single detector is passed through unchanged, zeros become NaN
    single = pd.DataFrame({"N1": [100.0, 0.0]})
    result = crnpy.total_raw_counts(single)
    assert result.iloc[0] == 100.0 and np.isnan(result.iloc[1])


def test_is_outlier_methods():
    x = pd.Series([100.0, 102.0, 98.0, 101.0, 99.0, 100.0, 300.0, 101.0, 99.0, 5.0, 100.0])
    # scaled MAD is two-sided, as in MATLAB isoutlier: both the high and the low spike are flagged
    flagged = crnpy.is_outlier(x, method="scaled_mad")
    assert flagged.dtype == bool and flagged[6] and flagged[9] and not flagged[0]
    # the range check is always applied in addition to the selected method
    flagged = crnpy.is_outlier(x, method="iqr", min_val=50, max_val=200)
    assert flagged[6] and flagged[9]
    assert crnpy.is_outlier(x, method="range", min_val=50, max_val=200).sum() == 2
    with pytest.raises(ValueError):
        crnpy.is_outlier(x, method="range")


def test_smooth_1d_dataframe_and_series():
    df = pd.DataFrame({"a": np.arange(20.0), "b": np.arange(20.0) ** 2})
    original = df.copy()
    smoothed = crnpy.smooth_1d(df, window=5, order=2, method="savitzky_golay")
    assert smoothed.shape == df.shape
    pd.testing.assert_frame_equal(df, original)  # input must not be modified
    series = crnpy.smooth_1d(df["a"], window=5, method="moving_median")
    assert len(series) == 20


def test_spatial_average_and_idw():
    x = np.array([0.0, 10.0, 20.0, 1000.0])
    y = np.zeros(4)
    z = np.array([1.0, np.nan, 3.0, 100.0])
    smoothed = crnpy.spatial_average(x, y, z, buffer=50, min_neighbours=1, method="mean")
    assert len(smoothed) == 4 and np.isnan(smoothed[1])
    assert smoothed[0] == pytest.approx(2.0)      # mean of the available neighbours within the buffer
    assert smoothed[3] == pytest.approx(100.0)    # isolated point keeps its value
    # IDW returns the observed value at a coincident point and the distance-weighted mean elsewhere
    xo, yo, zo = np.array([0.0, 10.0]), np.array([0.0, 0.0]), np.array([1.0, 3.0])
    pred = crnpy.idw(xo, yo, zo, np.array([0.0, 5.0]), np.array([0.0, 0.0]), neighborhood=100, p=1)
    np.testing.assert_allclose(pred, [1.0, 2.0])


def test_interpolate_incoming_flux_keeps_length_and_nans():
    nmdb_ts = pd.date_range("2021-01-01 00:00", periods=4, freq="h")
    nmdb_counts = np.array([100.0, np.nan, 102.0, 103.0])
    crnp_ts = pd.date_range("2021-01-01 00:10", periods=4, freq="h")
    flux = crnpy.interpolate_incoming_flux(nmdb_ts, nmdb_counts, crnp_ts)
    assert len(flux) == 4 and np.isnan(flux[1]) and flux[0] == 100.0


def test_find_neutron_monitor_offline_listing():
    # Without dates no network access is needed; the ten closest stations by cutoff rigidity are returned
    result = crnpy.find_neutron_monitor(2.40)
    assert len(result) == 10
    assert "NEWK" in result["STID"].values


# NMDB draw_graph.php intermittently returns the HTML page without the ASCII data block; a minimal
# good body (parsed to two rows) and a bad body without the RCORR_E marker.
_NMDB_GOOD = ("<html>preamble RCORR_E\n2020-01-01 00:00:00;100.0\n2020-01-01 01:00:00;101.0\n"
              "\n</code></pre><br>Total time")
_NMDB_BAD = "<html>the following query returned 0 rows, no data block here</html>"


class _FakeResponse:
    def __init__(self, text):
        self.content = text.encode("utf-8")


def test_get_incoming_neutron_flux_retries_then_succeeds(monkeypatch):
    calls = {"n": 0}

    def fake_get(url):
        calls["n"] += 1
        return _FakeResponse(_NMDB_BAD if calls["n"] == 1 else _NMDB_GOOD)

    monkeypatch.setattr(crnpy.crnpy.requests, "get", fake_get)
    monkeypatch.setattr(crnpy.crnpy.time, "sleep", lambda s: None)
    df = crnpy.get_incoming_neutron_flux(pd.to_datetime("2020-01-01 00:00"),
                                         pd.to_datetime("2020-01-01 02:00"), station="NEWK")
    assert calls["n"] == 2          # first response was retried
    assert df is not None and len(df) == 2


def test_get_incoming_neutron_flux_gives_up_after_retries(monkeypatch):
    calls = {"n": 0}

    def fake_get(url):
        calls["n"] += 1
        return _FakeResponse(_NMDB_BAD)

    monkeypatch.setattr(crnpy.crnpy.requests, "get", fake_get)
    monkeypatch.setattr(crnpy.crnpy.time, "sleep", lambda s: None)
    df = crnpy.get_incoming_neutron_flux(pd.to_datetime("2020-01-01 00:00"),
                                         pd.to_datetime("2020-01-01 02:00"), station="NEWK")
    assert calls["n"] == 3          # at most three attempts in total
    assert df is None
