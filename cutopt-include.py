# -*- coding: utf-8 -*-
def run_cutopt(fname, write_report = True):
    model = build_model_from_rootfile(fname)
#    model.scale_predictions(5.0 / 1.1)
    model.fill_histogram_zerobins()
    model.set_signal_processes('zp*')
    for p in model.processes:
        model.add_lognormal_uncertainty('lumi', math.log(1.044), p)
     
    model.add_lognormal_uncertainty('zj_rate', math.log(2.0), 'zjets')
    model.add_lognormal_uncertainty('wj_rate', math.log(1.5), 'wjets')
    model.add_lognormal_uncertainty('qcd_rate', math.log(2.0), 'qcd')
    model.add_lognormal_uncertainty('ttbar_rate', math.log(1.5), 'ttbar')
    model_summary(model)
    
    #dist = model.distribution
    #model.distribution = get_fixed_dist(dist)
    #res = pl_intervals(model, input = 'toys-asimov:0.0', n = 1, write_report = False, nuisance_constraint = dist)
    
    exp, obs = asymptotic_cls_limits(model, use_data=False, signal_process_groups=None, beta_signal_expected=0.0, bootstrap_model=True, input=None, n=1, options=None)
    print "expected limit"
    """    
    for spid in res:
        print spid, res[spid][0.90][0][1]
        
    """
    for x, y in zip(exp.x, exp.y):
        print x, y
       
    print "end"
    if write_report: report.write_html('./htmlout')
