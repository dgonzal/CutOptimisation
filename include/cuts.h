// cuts on the opt_event level; apply reweight instead of cut.

#include "cutopt.cxx"

class OptCut{
public:
    virtual double operator()(const opt_event & evt) const = 0;
};

class PassAll: public OptCut{
public:
    double operator()(const opt_event & evt) const{
        return true;
    }
};

class CountBtags: public OptCut{
public:
    CountBtags(const obs_selection & obs_sel, int nmin_, int nmax_, double threshold): nmin(nmin_), nmax(nmax_), disc_threshold(threshold){
        index_btag1 = obs_sel.indexof("btag1");
        index_btag2 = obs_sel.indexof("btag2");
    }
    virtual double operator()(const opt_event & evt) const {
        int nbtags = 0;
        if(evt.data[index_btag1] > disc_threshold) ++nbtags;
        if(evt.data[index_btag2] > disc_threshold) ++nbtags;
        bool pass = (nbtags >= nmin || nmin < 0) && (nbtags <= nmax || nmax < 0);
        return pass?1.0:0.0;
    }
private:
    double disc_threshold;
    int nmin, nmax;
    size_t index_btag1, index_btag2;
};

class CountToptags: public OptCut{
public:
    CountToptags(const obs_selection & obs_sel, int nmin_, int nmax_, double ptmin_ = 350, double mjmin_ = 140, double mjmax_=250, double mmin_min_=50):
       ptmin(ptmin_), mjmin(mjmin_), mjmax(mjmax_), mmin_min(mmin_min_),
       nmin(nmin_), nmax(nmax_){
         index_ca1m = obs_sel.indexof("ca1m");
         index_ca2m = obs_sel.indexof("ca2m");
         index_ca1mmin = obs_sel.indexof("ca1mmin");
         index_ca2mmin = obs_sel.indexof("ca2mmin");
         index_ca1pt = obs_sel.indexof("ca1pt");
         index_ca2pt = obs_sel.indexof("ca2pt");
    }
    virtual double operator()(const opt_event & evt) const {
        if(evt.event_type == opt_event::signal || evt.event_type==opt_event::ttbar_bkg){
            int ntags = 0;
            if(evt.data[index_ca1pt] > ptmin && evt.data[index_ca1m] > mjmin && evt.data[index_ca1m] < mjmax && evt.data[index_ca1mmin] > mmin_min) ++ntags;
            if(evt.data[index_ca2pt] > ptmin && evt.data[index_ca2m] > mjmin && evt.data[index_ca2m] < mjmax && evt.data[index_ca2mmin] > mmin_min) ++ntags;
            bool pass = (ntags >= nmin || nmin < 0) && (ntags <= nmax || nmax < 0);
            return pass?1.0:0.0;
        }
        else{
            double p1 = toptagger_mistag_rate(evt.data[index_ca1pt]);
            double p2 = toptagger_mistag_rate(evt.data[index_ca2pt]);
            double prob_n0 = (1 - p1) * (1 - p2);
            double prob_n1 = p1 * (1 - p2) + p2 * (1 - p1);
            double prob_n2 = p1 * p2;
            double total_prob = 0;
            if((nmin < 0 || nmin <= 0) && (nmax < 0 || nmax >=0)) total_prob += prob_n0;
            if((nmin < 0 || nmin <= 1) && (nmax < 0 || nmax >=1)) total_prob += prob_n1;
            if((nmin < 0 || nmin <= 2) && (nmax < 0 || nmax >=2)) total_prob += prob_n2;
            return total_prob;
        }
    }
    
private:
    double ptmin, mjmin, mjmax, mmin_min;
    int nmin, nmax;
    size_t index_ca1m, index_ca2m, index_ca1mmin, index_ca2mmin, index_ca1pt, index_ca2pt;
};


class AndCut: public OptCut{
public:
    void add(const boost::shared_ptr<OptCut> & cut){
        cuts.push_back(cut);
    }
    
    virtual double operator()(const opt_event & evt) const{
        double result = 1.0;
        for(size_t i=0; i<cuts.size(); ++i){
            result *= (*cuts[i])(evt);
        }
        return result;
    }
private:
    vector<boost::shared_ptr<OptCut> > cuts;
};

class OrCut: public OptCut{
public:
    void add(const boost::shared_ptr<OptCut> & cut){
        cuts.push_back(cut);
    }
    
    virtual double operator()(const opt_event & evt) const{
        double result = 0.0;
        for(size_t i=0; i<cuts.size(); ++i){
            result = max(result, (*cuts[i])(evt));
        }
        return result;
    }
private:
    vector<boost::shared_ptr<OptCut> > cuts;
};


class CutInfoCut: public OptCut{
public:
    CutInfoCut(const cutinfo & info_, const obs_selection & obs_sel): info(info_) {
       index = obs_sel.indexof(info.obsname);
    }
    
    virtual double operator()(const opt_event & evt) const{
        bool result;
        switch(info.mode){
            case cutinfo::cut_gtr: result = evt.data[index] > *(info.cut_threshold); break;
            case cutinfo::cut_lt: result = evt.data[index] < *(info.cut_threshold); break;
            default:
              throw string("logic error in CutInfoCut::operator()");
        }
        return result?1.0:0.0;
    }
    
private:
    cutinfo info;
    size_t index;
};


