import numpy as np
'''
Dictionaries with sample selections
'''
def sample_cuts(tdata,cut):
   #several selections are presented below
    if cut == 'OFFICIAL_RED':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.6) & (tdata['DESDM_ZP'] <= 1.0)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_060_065':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.6) & (tdata['DESDM_ZP'] <= 0.65)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_065_070':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.65) & (tdata['DESDM_ZP'] <= 0.70)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_070_075':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.70) & (tdata['DESDM_ZP'] <= 0.75)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_075_080':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.75) & (tdata['DESDM_ZP'] <= 0.80)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_080_085':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.80) & (tdata['DESDM_ZP'] <= 0.85)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_085_090':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.85) & (tdata['DESDM_ZP'] <= 0.90)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_090_095':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.90) & (tdata['DESDM_ZP'] <= 0.95)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    elif cut == 'OFFICIAL_RED_095_100':
        mask_mag = (tdata['MAG_AUTO_I'] > 17.5) & (tdata['MAG_AUTO_I'] < 22)
        mask_gr_col = (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] > -1) & (tdata['MAG_AUTO_G'] - tdata['MAG_AUTO_R'] < 3)
        mask_ri_col = (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] > -1) & (tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I'] < 2.5)
        mask_iz_col = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] > -1) & (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z'] < 2)
        mask_red = (tdata['MAG_AUTO_I'] - tdata['MAG_AUTO_Z']) + 0.6*(tdata['MAG_AUTO_R'] - tdata['MAG_AUTO_I']) > 0.9
        mask_sg = tdata['SPREAD_MODEL_I'] + (5./3.)*(tdata['SPREADERR_MODEL_I']) > 0.007
        mask_zph = (tdata['DESDM_ZP'] > 0.95) & (tdata['DESDM_ZP'] <= 1.0)
        mask_qlt = (tdata['FLAGS_GOLD'] == 0) & (tdata['FLAGS_BADREGION'] == 0)
        mask = mask_mag & mask_gr_col & mask_ri_col & mask_iz_col & mask_red & mask_sg & mask_zph & mask_qlt
    else:
        print 'SAMPLE NOT FOUND, using default cut'
        mask = tdata['MAG_AUTO_I'] < 1000 ## any dummy cut not actually removing anything
    return mask
