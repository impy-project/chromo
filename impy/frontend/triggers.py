'''
Created on 20.04.2012

@author: afedynitch
'''
import numpy as np

#===============================================================================
# Trigger (base class)
#===============================================================================
class Trigger():
    def __init__(self, *args):
        self.histogram_list = {}

    def trigger_event(self, event):
        raise Exception('Warning(): Base class trig_eventtype() called. This \
                         method has to be overridden by the user.')
        
    def __ne__(self, other_trigger):
        if other_trigger.__class__.__name__ != self.__class__.__name__:
            print "Found unequal triggers."
            return False
    
    def __eq__(self, other_trigger):
        return not self.__ne__(other_trigger)
    
    def associate_hist(self, histObject):
        if histObject.unique_name not in self.histogram_list.keys():
            self.histogram_list[histObject.unique_name] = histObject
        else:
            raise Exception(("Error: The name of the histogram '{0}' " + 
                             "is already in the list.").format(histObject.unique_name))
            
    def reset(self):
        raise Exception("Warning: Base class reset() called. This method has to be overridden by the user.")
    
    def result(self):
        raise Exception("Warning: Base class result() called. This method has to be overridden by the user.")


#===============================================================================
# SPS Experiment Triggers
#===============================================================================          
class UA5_NSD(Trigger):
    def trigger_event(self, event):
        pos_counts = np.count_nonzero((event.eta >= 2.0) & \
                                     (event.eta <= 5.6))
        neg_counts = np.count_nonzero((event.eta <= -2.0) & \
                                     (event.eta >= -5.6))
        return pos_counts > 0 and neg_counts > 0   

class UA1_MinBias(Trigger):
    def trigger_event(self, event):
        pos_counts = np.count_nonzero((event.eta > 1.5) & \
                                     (event.eta < 5.5))
        neg_counts = np.count_nonzero((event.eta < -1.5) & \
                                     (event.eta > -5.5))
        return pos_counts > 0 and neg_counts > 0   

class UA5_NSD_53GeV_Ppbar(Trigger):
    def trigger_event(self, event):
        pos_counts = np.count_nonzero((event.eta >= 2.0) & \
                                     (event.eta <= 5.6))
        neg_counts = np.count_nonzero((event.eta <= -2.0) & \
                                     (event.eta >= -5.6))
        return pos_counts > 1 and neg_counts > 1 
        
#===============================================================================
# #CDF Triggers
#===============================================================================
class CDF_SD(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((np.abs(event.eta) > 3.2) & \
                                (np.abs(event.eta) < 5.9)) > 0

class CDF_NSD(Trigger):
    def trigger_event(self, event):
        return (np.count_nonzero((event.eta > 3.2) & (event.eta < 5.9)) > 0) & \
               (np.count_nonzero((event.eta < -3.2) & (event.eta > -5.9)) > 0) 



#===============================================================================
# LHCb
#===============================================================================
class LHCb_VELO(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((event.eta > 2.) & \
                                (event.eta < 5.)) > 0


class LHCb_MB_2014(Trigger):
    def trigger_event(self, event):
        return np.any((event.eta > 2.) 
                      & (event.eta < 4.8) 
                      & (event.p_tot > 2.0) 
                      & (event.pt > 0.2))
                                
class LHCb_flow_MB(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((event.eta > 1.9) & 
                                (event.eta < 4.9)) > 0

class LHCb_flow_HARD(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((event.eta > 1.9) & 
                                (event.eta < 4.9) & 
                                (event.pt > 3.0)) > 0
class LHCb_flow_DIFF(Trigger):
    def trigger_event(self, event):
        side1 = np.any((event.eta > 1.9) & 
                                (event.eta < 4.9))
        side2 = not np.any((event.eta > -3.5) & 
                                (event.eta < -1.5))
        return side1 and side2
    
class LHCb_flow_ND(Trigger):
    def trigger_event(self, event):
        side1 = np.count_nonzero((event.eta > 1.9) & 
                                (event.eta < 4.9)) > 0
        side2 = np.count_nonzero((event.eta > -3.5) & 
                                (event.eta < -1.5)) > 0
        return side1 and side2                                  
#===============================================================================
# CMS Triggers
#===============================================================================
class CMS_BEAM(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((np.abs(event.eta) > 3.23) & \
                                (np.abs(event.eta) < 4.65)) > 0
                                
class CMS_eflow_MB(Trigger):
    def trigger_event(self, event):
        side1 = np.any((event.eta > 3.23) & (event.eta < 4.65))
        side2 = np.any((event.eta > -4.65) & (event.eta < -3.23))
        return side1 and side2

class CMS_HFCAL(Trigger):
    def trigger_event(self, event):
        return np.any(event.en[(event.eta < -2.9) & (event.eta > -5.2)] > 3.0) & \
                                np.any(event.en[(event.eta > 2.9) & (event.eta < 5.2)] > 3.0)

class CMS_NSD(Trigger):
    def __init__(self):
        Trigger.__init__(self)
        self.beam_trigger = CMS_BEAM()
        self.hfcal = CMS_HFCAL()
        
    def trigger_event(self, event):
        return self.beam_trigger.trigger_event(event) and \
               self.hfcal.trigger_event(event)

class CMS_DS(Trigger):
    def trigger_event(self, event):
        pos_counts = np.count_nonzero((event.eta > 3) & \
                                     (event.eta < 5) & \
                                     (event.en > 3))
        neg_counts = np.count_nonzero((event.eta < -3) & \
                                     (event.eta > -5) & \
                                     (event.en > 3))
        return pos_counts > 1 and neg_counts > 1 
    
    
#===============================================================================
# ATLAS Triggers
#===============================================================================
class ATLAS_MBTS(Trigger):
    def trigger_event(self, event):
        try:
            return np.any((np.abs(event.eta) > 2.09) &
                      (np.abs(event.eta) < 3.84) &
                      event.charge != 0)
        except:
            return np.any((np.abs(event.eta) > 2.09) &
                      (np.abs(event.eta) < 3.84))
                                
                                
class ATLAS_MB_2011(Trigger):
    def __init__(self, nch_min, pt_min):
        self.nch_min = nch_min
        self.pt_min = pt_min
        Trigger.__init__(self)
        
    def trigger_event(self, event):
        return np.count_nonzero((np.abs(event.eta) < 2.5) & \
                                (event.pt > self.pt_min)) >= self.nch_min
#===============================================================================
# ALICE triggers
#===============================================================================
class ALICE_INEL(Trigger):
    def __init__(self):
        Trigger.__init__(self)
        self.spd = ALICE_SPD()
        self.vzero = ALICE_VZERO_OR()
        
    def trigger_event(self, event):
        return self.spd.trigger_event(event) or self.vzero.trigger_event(event)

class ALICE_NSD(Trigger):
    def __init__(self):
        Trigger.__init__(self)
        self.vzero = ALICE_VZERO_AND()
        
    def trigger_event(self, event):
        return self.vzero.trigger_event(event)


class ALICE_VZERO_OR(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((event.eta > -3.23) & \
                                (event.eta < -1.7)) > 0 or \
               np.count_nonzero((event.eta > 2.8) & \
                                (event.eta < 5.1)) > 0

class ALICE_VZERO_AND(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((event.eta > -3.23) & \
                                (event.eta < -1.7)) > 0 and \
               np.count_nonzero((event.eta > 2.8) & \
                                (event.eta < 5.1)) > 0


class ALICE_SPD(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero(np.abs(event.eta) < 2.0) > 0

class ALICE_INEL_GT_0(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero(abs(event.eta) < 0.5) > 0

class ALICE_INEL_GT_0_7TeV(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero(abs(event.eta) < 1.0) > 0

#===============================================================================
# Totem trigger
#===============================================================================
class Totem_Inclusive(Trigger):
    def trigger_event(self, event):
        return np.count_nonzero((np.abs(event.eta) > 5.3) & \
                                (np.abs(event.eta) < 6.5)) > 0

class Totem_NSD(Trigger):
    def trigger_event(self, event):
        return (np.count_nonzero((event.eta > -6.5) & 
                                 (event.eta < -5.3)) > 0) and \
               (np.count_nonzero((event.eta > 5.3) & 
                                 (event.eta < 6.5)) > 0)
class Totem_SD(Trigger):
    def trigger_event(self, event):
        return (np.count_nonzero((event.eta > -6.5) & 
                                 (event.eta < -5.3)) > 0) != \
               (np.count_nonzero((event.eta > 5.3) & 
                                 (event.eta < 6.5)) > 0)

class Totem_DISP(Trigger):
    def trigger_event(self, event):
        return (np.any((event.eta > -7) & (event.eta < -6)) or
               np.any((event.eta > 3.7) & (event.eta < 4.8)))

#===============================================================================
# Misc
#===============================================================================
class MT_trigger(Trigger):
    def trigger_event(self, event):
        return True
