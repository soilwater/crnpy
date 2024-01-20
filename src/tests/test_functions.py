import crnpy
import numpy as np
import pandas as pd

def test_correction_bwe():
    r2_N0 = 6.4/1210
    N = 420
    Nv = 574 # Eg. 1 from Table 3. in Baatz et al. 2015
    BWE = 51.5
    Nv_ = crnpy.correction_bwe(N, BWE, r2_N0)
    assert abs(Nv - Nv_) < 5 # 5 counts of difference is acceptable

def test_correction_humidity():
    # From Rosolem et al. 2013 Fig 2.
    fw = crnpy.correction_humidity(20, 0)
    assert round(fw,1) == 1.1

def test_correction_inncoming_neutrons():
    incomming_neutrons = crnpy.get_incoming_neutron_flux(pd.to_datetime('2013-01-01 12:00:00'), pd.to_datetime('2013-01-05 12:00:00'), station='IRKT')
    fi = crnpy.correction_incoming_flux(incomming_neutrons['counts'], 210)
    print(incomming_neutrons['counts'])
    # check that fi is between 0.9 and 1.1
    assert (fi > 0.9).all() and (fi < 1.1).all()
    #check that highest count has the highest correction factor
    min = np.argmin(incomming_neutrons['counts'])
    assert fi[min] == np.min(fi)
    #check that lowest count has the lowest correction factor
    max = np.argmax(incomming_neutrons['counts'])
    assert fi[max] == np.max(fi)

def test_correction_pressure():
    fp = crnpy.correction_pressure(1013, 1013, 130)
    assert round(fp,1) == 1.0

def test_count_to_vwc():
    # From Patrigani et al. 2021
    N0 = 3767
    N = 2580
    vwc = crnpy.counts_to_vwc(N, N0, 0.033, 0, 1.33)
    assert round(vwc,2) == 0.15

def test_weighting():
    # create vector from 0 to 350 m with .1 spacing
    x = np.arange(0, 350, 0.1)
    y = np.repeat(0, len(x))
    h = np.repeat(10, len(x))


    p = np.repeat(1013.25, len(x))
    Hveg = np.repeat(0, len(x))
    bd = 1.1
    SM = 0.02
    sm = np.repeat(SM, len(x))

    # calculate weighting function using Schrön et al. 2017 method
    weights = crnpy.nrad_weight(h, sm, x, y, bd, method="Schron_2017", p=p, Hveg=Hveg)
    # check that the sum of the weights is 1
    Max = 0.01147
    Avg = 0.00029
    Min = 0.00000
    assert round(np.sum(weights),5) == 1
    # Compare to the values in supplementary material of Schrön et al. 2017
    # check that the max weight is 0.01
    assert round(np.max(weights),3) == round(Max,3)
    # check that the avg weight is 0.00029
    assert round(np.mean(weights),5) == round(Avg,5)
    # check that the min weight is 0.00000
    assert round(np.min(weights),5) == round(Min,5)

def test_correction_inncoming_neutrons_RcMethods():
    # Test Hawdon 2014 Method using Baldry location from the original paper
    RcJUNG = 4.50 #4.49
    Rc_local = crnpy.cutoff_rigidity(-32.87, 148.54)
    reference_counts = 159
    # GEt JUNG for entire 2012
    reference_data = crnpy.get_incoming_neutron_flux(pd.to_datetime('2012-01-01 12:00:00'), pd.to_datetime('2013-01-01 12:00:00'), station='JUNG')
    reference_neutron_flux = reference_data['counts']
    assert np.abs(Rc_local - 4.7) < 0.5 # check that the calculated local cutoff rigidity is close to the measured value reported in the paper
    fi = crnpy.correction_incoming_flux(reference_neutron_flux, reference_counts, Rc_method='Hawdonetal2014', Rc_site=Rc_local, Rc_ref=RcJUNG)
    paper_fi_range = [0.87, 1.04]
    assert round(np.min(fi),2) == round(paper_fi_range[0],2)
    assert round(np.max(fi),2) == round(paper_fi_range[1],2)

    # Test McJannet and Desilets 2023 Method using data from the supplementary material
    CRNSelevation = 604 #m
    CRNSlatitude = -22.282 #deg
    CRNSlongitude = 133.251 #deg
    CutoffRigidity = 10.08 #GV
    # Simulate a 1:1 relationship between incoming neutron flux and counts to get Tau
    # This is not a realistic scenario but it is useful for testing
    reference_counts = np.array([2400])
    reference_neutron_flux = np.array([1200])
    site_atmdepth = crnpy.atmospheric_depth(CRNSelevation, CRNSlatitude)
    assert round(site_atmdepth,2) == 963.47
    JUNG_atmdepth = 665.18 # From supplementary material of McJannet and Desilets 2023
    fi = crnpy.correction_incoming_flux(reference_neutron_flux, reference_counts, Rc_method='McJannetandDesilets2023', Rc_site=CutoffRigidity, Rc_ref=RcJUNG, site_atmdepth=site_atmdepth, ref_atmdepth=JUNG_atmdepth)
    f0 = 0.5
    #if the formula for fi is fi = 1 / (tau * f0 + 1 - tau) get tau
    tau = (1/fi[0] - 1) / (f0 - 1)
    Tau_paper = 0.4753 # From supplementary material of McJannet and Desilets 2023
    assert round(tau,2) == round(Tau_paper,2)






