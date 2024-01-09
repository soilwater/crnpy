import crnpy


def test_correction_bwe():
    r2_N0 = 6.4/1210
    N = 420
    Nv = 574 # Eg. 1 from Table 3. in Baatz et al. 2015
    BWE = 51.5
    Nv_ = crnpy.correction_bwe(N, BWE, r2_N0)
    assert abs(Nv - Nv_) < 5 # 5 counts of difference is acceptable

def test_correction_humidity():
    # From Rosolem et al. 2013 Fig 2.
    fp = crnpy.correction_humidity(20, 0)
    assert round(fp,1) == 1.1

def test_correction_pressure():