from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import plot

jet_systs = {"central":"",
             "jes_up":"jes_up","jes_down":"jes_down",
        "jer_up":"jer_up","jer_down":"jer_down"}

muon_systs = {"central":"",
                "mutrig_up":"trig_up","mutrig_down":"trig_down",
        "muid_up":"id_up"    ,"muid_down":"id_down",
            "muiso_up":"iso_up"  ,"muiso_down":"iso_down",
                "mupu_up":"pu_up"    ,"mupu_down":"pu_down"}

#single muon control region plots
#                         cuts       mask        loose mus     ak4 jets         leading ak8
def muon_control_region(gghbbcuts,singlemuon_cr, looseMuons, ak4puppijet_pt30, leadingak8jet, plots):
