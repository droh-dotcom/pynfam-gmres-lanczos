#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc

pynfam_inputs = {
 'directories': {
     'outputs' : 'TEST_GMRES',
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    :  (-2, (-0.2,0,0.2)),
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0,
     #'retry_nonconv': True,
     },

 'fam_mode': {
     'fam_contour': 'CONSTr',
     'beta_type'  : '-',
     'fam_ops'    : 'GAMOWTELLER'
     }
}

override_settings = {
 'hfb' : {'proton_number'    : 3, 
          'neutron_number'   : 3, 
          'number_of_shells' : 2,
          #'number_gauss'     : 20,
          #'number_laguerre'  : 20,
          #'number_legendre'  : 40,
          'functional'       : 'SKM*',
          'pairing_cutoff': 60.0,
          'pairing_feature': 0.5,
          #"vpair_n": 0.0,
     	#"vpair_p": 0.0,
     	"number_iterations": 600,
     	"accuracy": 1e-08,
          #'force_parity'     : False
          },

 'ctr' : {
	    'energy_min': 0.0,
	    'energy_max': 10.0, 
	    'half_width': 0.5, 
	    'de_hw_ratio': 1,
	    #'nr_points' : 1 
	   },

 'fam' : {
 	 'use_fam_storage': 0,
 	 'interaction_name' : 'SKM*',
 	 'force_j2_terms': False,
 	 'convergence_epsilon': 1e-07,
 	 'require_self_consistency': True,
 	 'require_gauge_invariance': True,
 	 #'vpair_t0': None,
 	 #'vpair_t1': None,
 	 #'override_cs0': None,
 	 #'override_csr': 0.0,
 	 #'override_cds': 0.0,
 	 'override_ct' : None,
 	 #'override_cgs' : 0.0,
 	 #'override_cf': 0.0,
 	 'max_iter': 200,
 	 'broyden_history_size': 200,
    },

 'psi' : {
 	#"GA": -1.1639079994315826,
 }

}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)

