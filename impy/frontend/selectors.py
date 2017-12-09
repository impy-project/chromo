'''
Created on 09.05.2014
This file contains a collection of strings for the histogram
classes. They can be individually defined in the plot definitions, too.
The prefered way is to recycle and generalize them here, because
the code still changes too rapidly.

@author: afedynitch
'''

#===============================================================================
# "Filler" are definitions of the observables which are filled in the histograms
#===============================================================================

# Fill histogram with feynman-x distribution
fillxf = 'lambda event, sel: event.xf[sel]'

# Fill histogram with Elab/Elabmax
fillxl = 'lambda event, sel: event.xl[sel]'

# Fill histogram with pt**2 distribution
fillpt2 = 'lambda event, sel: event.pt2[sel]'

# Fill histogram with pt**2 distribution
fillpt = 'lambda event, sel: event.pt[sel]'

# Fill histogram with rapidity (y) distribution
filly = 'lambda event, sel: event.y[sel]'

# Fill histogram with energy distribution
fillen = 'lambda event, sel: event.en[sel]'

# Fill histogram with rapidity (y) distribution
filleta = 'lambda event, sel: event.eta[sel]'

# Fill multiplicity
fillmult = 'lambda event, sel: np.count_nonzero(sel)'

#fillrapgap = ('lambda event, sel: np.max(np.sort(event.eta[sel])[1:] - '
#              + 'np.sort(event.eta[sel])[:-1]) if np.count_nonzero(sel) > 1 else -1.')

fillrapgap_ATLAS = (
    'lambda event, sel: max((4.9 - np.max(event.eta[sel]), ' +
    '4.9 + np.min(event.eta[sel]))) if np.count_nonzero(sel) > 1 else -1.')

#===============================================================================
# Weighted
#===============================================================================
# Fill histogram with invariant cross section as a function of xf
fillfxf = 'lambda event, sel: (event.xf[sel], event.fw[sel])'

# Fill histogram with invariant cross section as a function of xf
fille_flow = 'lambda event, sel: (event.eta[sel], event.en[sel])'

# Fill average pt as function of multiplicity
fillavpt_mult = 'lambda event, sel: (np.count_nonzero(sel), ' + \
                'np.mean(event.pt[sel]))'

# Fill average pt as function of multiplicity
fillmpi_mult = 'lambda event, sel: (np.count_nonzero(sel), ' + \
                'event.mpi)'

# Fill average pt as function of multiplicity
fillavp_eta = 'lambda event, sel: (event.eta[sel], ' + \
                                  'event.p_tot[sel])'

#===============================================================================
# Selectors are conditions about how to select the particles in the event
#===============================================================================

sel_pdgid = lambda pdgid: 'lambda event: event.p_ids == {0}'.format(pdgid)

sel_etarange = lambda eta: 'lambda event: np.where(np.abs(event.eta) < {0})'.format(eta)

sel_yrange = lambda y: 'lambda event: np.abs(event.y) < {0}'.format(y)

sel_all = 'lambda event: np.where(event.p_ids != 0)'

sel_charged = 'lambda event: event.charge != 0'

sel_charged_plus = 'lambda event: event.charge > 0'

sel_charged_minus = 'lambda event: event.charge < 0'

sel_pdgid_list = lambda pdgid_list: 'lambda event: np.logical_or.reduce( ' + \
                                    'np.vstack([(event.p_ids == pdgid) ' + \
                                    'for pdgid in {0}]), 0)'.format(pdgid_list)

sel_pdgid_list_yrange = lambda pdgid_list, yrange: ('lambda event: ' +
                        '(np.abs(event.y) < {1}) & (np.logical_or.reduce( ' +
                        'np.vstack([(event.p_ids == pdgid) ' +
                        'for pdgid in {0}]), 0))').format(pdgid_list, yrange)

sel_pdgid_list_ylims = lambda pdgid_list, ylims: ('lambda event: ' +
                        '(event.y > {1}) & (event.y < {2}) & ' +
                        '(np.logical_or.reduce( ' +
                        'np.vstack([(event.p_ids == pdgid) ' +
                        'for pdgid in {0}]), 0))').format(pdgid_list, *ylims)

sel_pdgid_ptrange = lambda pdgid, ptmin, ptmax: \
                    ('lambda event: (event.p_ids == {0}) & ' +
                     '((event.pt > {1}) & (event.pt < {2}))').format(pdgid,
                                                                     ptmin,
                                                                     ptmax)

sel_pdgid_yrange = lambda pdgid, yrange: \
                   ('lambda event: (event.p_ids == {0}) & ' +
                    '(np.abs(event.y) < {1})').format(pdgid, yrange)

sel_etalims_minpt_ptot = lambda etamin, etamax: \
                   ('lambda event: (event.eta >= {0}) & ' +
                    '(event.eta <= {1}) & (event.pt > 0.2) & ' +
                    '(event.p_tot > 2.0)').format(etamin, etamax)

sel_ptlims_ptot_etasel = lambda ptmin, ptmax: \
                   ('lambda event: (event.eta > {0}) & ' +
                    '(event.eta < {1}) & (event.pt > {2}) & ' +
                    '(event.pt < {3}) & (event.p_tot > 2.0)').format(2.0, 4.8,
                                                                     ptmin, ptmax)

sel_minpt = lambda minpt: 'lambda event: event.pt > {0}'.format(minpt)

sel_minpt_ptot = lambda minpt, ptot: ('lambda event: ((event.pt > {0}) & ' + \
                                        '(np.abs(event.p_tot) > {1}))').format(minpt, ptot)

sel_minch = lambda minch: ('lambda arr: arr if np.count_nonzero(arr) > {0} else ' + 
                            'np.zeros_like(arr,dtype="bool"'.format(minch))

sel_minpt_etarange = lambda minpt, eta: ('lambda event: (event.pt > {0}) & ' +
                                         '(np.abs(event.eta) < {1})').format(minpt, eta)
