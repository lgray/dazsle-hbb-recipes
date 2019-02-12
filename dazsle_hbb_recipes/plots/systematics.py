""" enumeration of systematics for easy looping """

jet_syst_list = []
muon_syst_list = []

jet_pt_systs = {"jet_pt":"pt",
                "jet_pt_jes_up":"pt_jes_up","jet_pt_jes_down":"pt_jes_down",
                "jet_pt_jer_up":"pt_jer_up","jet_pt_jer_down":"pt_jer_down"}

jet_weight_systs = {"jetweight":"weight",
                    "jetweight_pu_up":"weight_pu_up",
                    "jetweight_pu_down":"weight_pu_down",
                    "jetweight_trigger_up":"weight_trigger_up",
                    "jetweight_trigger_down":"weight_trigger_down"}

muon_weight_systs = {"muweight":"weight",
                     "muweight_trigger_up":"weight_trigger_up",
                     "muweight_trigger_down":"weight_trigger_down",
                     "muweight_id_up":"weight_id_up"    ,"muweight_id_down":"weight_id_down",
                     "muweight_iso_up":"weight_iso_up"  ,"muweight_iso_down":"weight_iso_down",
                     "muweight_pu_up":"weight_pu_up"    ,"muweight_pu_down":"weight_pu_down"}
