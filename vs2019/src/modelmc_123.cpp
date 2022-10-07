#include "modelmc_123.h"
#include "models/src/state.h"
#include "models/src/convert.h"
#include "version.txt"
#include <algorithm>
#include <limits>
#include <math.h>
#include <stdio.h>

static const double dPa = 101325.0;//标准大气压
static const unsigned long mKurNow = 0x400;  // 0x表示16进制,添加进行D-C模型开发
static const unsigned long mKurPast = 0x800;
static const unsigned long mEtNow = 0x1000;
static const unsigned long mEtPast = 0x2000;
//s->state_ |= tension_now;

// excerpt-export-start
int __stdcall DllMain(void *,unsigned, void *) {
    return 1;
}

extern "C" EXPORT_TAG const char *getName() {
#ifdef MODELDEBUG
    return "cmodelmc_123d";
#else
    return "cmodelmc_123";
#endif
}

extern "C" EXPORT_TAG unsigned getMajorVersion() {
    return MAJOR_VERSION;
}

extern "C" EXPORT_TAG unsigned getMinorVersion() {
    return UPDATE_VERSION;
}

extern "C" EXPORT_TAG void *createInstance() {
    models::ModelMc_123 *m = new models::ModelMc_123();
    return (void *)m;
}
// excerpt-export-end

namespace models {

// excerpt-con-start
    ModelMc_123::ModelMc_123(unsigned short option) : ModelElastic(option) {
    }
// excerpt-con-end

    String ModelMc_123::getName() const { 
#ifdef MODELDEBUG
        return L"mc_123-debug"; 
#else
        return L"mc_123";
#endif
    }

    String ModelMc_123::getFullName() const { 
#ifdef MODELDEBUG
        return L"Mc_123 Debug";
#else
        return L"Mc_123";
#endif
    }

    UInt ModelMc_123::getMinorVersion() const {
        return UPDATE_VERSION;
    }

    String ModelMc_123::getProperties() const {
        return ModelElastic::getProperties() + L",cohesion,friction,dilation,tension,flag-brittle, fricdel, ratiofail, ke, ne, kb, mb, kur";
    }

    String ModelMc_123::getStates() const {
        return L"shear-n,tension-n,shear-p,tension-p,Kur-n, Kur-p, Et-n, Et-p";
    }

    Variant ModelMc_123::getProperty(UInt index) const {
        if (index <= 4)
            return ModelElastic::getProperty(index);
        else {
            switch (index) {     
            case 5: return cohesion_;
            case 6: return friction_; 
            case 7: return dilation_;
            case 8: return tension_;
            case 9: return brittle_;
            case 10: return FricDel_;
            case 11: return FailRatio_;
            case 12: return Ke_;
            case 13: return Ne_;
            case 14: return Kb_;
            case 15: return Mb_;
            case 16: return Kur_;
            }
        }
        return 0.0;
    }

    void ModelMc_123::setProperty(UInt index,const Variant &p,UInt restoreVersion) {
        if (index <= 4)
            ModelElastic::setProperty(index,p,restoreVersion);
        else {
            switch (index) {
            case 5: cohesion_ = p.toDouble(); break;
            case 6: friction_ = p.toDouble(); break;
            case 7: dilation_ = p.toDouble(); break;
            case 8: tension_  = p.toDouble(); break;
            case 9: brittle_  = p.toBool();   break;
            case 10:  FricDel_ = p.toDouble();      break;
            case 11:  FailRatio_ = p.toDouble();    break;
            case 12:  Ke_ = p.toDouble();       break;
            case 13:  Ne_ = p.toDouble();       break;
            case 14:  Kb_ = p.toDouble();       break;
            case 15:  Mb_ = p.toDouble();       break;
            case 16:  Kur_ = p.toDouble();      break;
            }
        }
    }

    bool ModelMc_123::isPropertyAdvanced(UInt i) const {
        if (i <= 4)
            return ModelElastic::isPropertyAdvanced(i);
        else if (i==16) 
            return true;
        return false;
    }

    void ModelMc_123::copy(const ConstitutiveModel *m) {
        const ModelMc_123 *mm = dynamic_cast<const ModelMc_123 *>(m);
        if (!mm) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        //
        ModelElastic::copy(m);
        //
        cohesion_ = mm->cohesion_;
        friction_ = mm->friction_;
        dilation_ = mm->dilation_;
        tension_  = mm->tension_;
        brittle_  = mm->brittle_;
        FricDel_ = mm->FricDel_;
        FailRatio_ = mm->FailRatio_;
        Ke_ = mm->Ke_;
        Ne_ = mm->Ne_;
        Kb_ = mm->Kb_;
        Mb_ = mm->Mb_;
        Kur_ = mm->Kur_;
    }

// excerpt-run-start
    void ModelMc_123::initialize(UByte d,State *s) {
        ConstitutiveModel::initialize(d,s);
        //必须要有E，K和G的初始化，否则会出现刚度为零的报错
        double dEt = 0.25 * Ke_ * dPa * pow(0.02, Ne_);
        bulk_ = dEt / 3;
        shear_ = (9 * bulk_ - dEt <= 0) ? 0 : 3 * bulk_ * dEt / (9 * bulk_ - dEt);

        updateParameters();
    }

    void ModelMc_123::run(UByte d,State *s) {
        ConstitutiveModel::run(d,s);

        if (s->modulus_reduction_factor_ > 0.0)
            moduliReduction(s->modulus_reduction_factor_);
// excerpt-state-start
        if (s->state_ & shear_now) s->state_ |= shear_past;
        s->state_ &= ~shear_now;
        if (s->state_ & tension_now) s->state_ |= tension_past;
        s->state_ &= ~tension_now;
        if (s->state_ & mKurNow) s->state_ |= mKurPast;
        s->state_ &= ~mKurNow;
        if (s->state_ & mEtNow) s->state_ |= mEtPast;
        s->state_ &= ~mEtNow;
// excerpt-state-end
        UInt iPlas = 0;
        ModelElastic::elasticTrial(s);//初始化——更新参数——更新阿尔法——运行——更新应变增量——更新参数
        s->viscous_ = true;
        if (!canFail()) return;

        SymTensorInfo info;
        DVect3 prin = s->stnS_.getEigenInfo(&info);
        //(1)求新的摩擦角，修正摩擦角        
        double dSD = prin.z() - prin.x();//主应力差>0
        
        if (-1.0 * prin.z() > S3Histroy_)       S3Histroy_ = -1.0 * prin.z();
        //考虑摩擦角随围压的影响,内摩擦角修正
        double dlog = (S3Histroy_ != 0.0) ? log10(S3Histroy_ / dPa) : 0.0;//如果delt_f不等于0，那么dlog进行修正
        friction_ = friction_ - FricDel_ * dlog;
        Fcos_ = cos(friction_ * degrad);//能用户自定义的参数加_，否则不加
        Fsin_ = sin(friction_ * degrad);

        //(2)计算应力水平
        double dSF = 2 * cohesion_ * Fcos_ + (-2.0) * prin.z() * Fsin_;//破坏摩尔圆直径
        double dSL = (dSF != 0) ? dSD * (1 - Fsin_) / dSF : 0.0;      //除数(no p)
        if (dSL > SLHistroy_)            SLHistroy_ = dSL;
        if (dSL >= 0.99)                 dSL = 0.99;

        //(3)加卸载判断，计算模量
        double dSPa = S3Histroy_ / dPa;
        double dEi = Ke_ * dPa * pow(dSPa, Ne_);
        double dEt = 0.25 * Ke_ * dPa * pow(0.02, Ne_);

        if (dEt < dEi * (1 - FailRatio_ * dSL) * (1 - FailRatio_ * dSL)) {
            dEt = dEi * (1 - FailRatio_ * dSL) * (1 - FailRatio_ * dSL);//加载
            s->state_ &= ~mEtNow;
        }

        double dSL_ss=dSL* pow(dSPa, 0.25);
        if (dSL_ss > SDHistroy_)                   SDHistroy_ = dSL_ss;
        double dSL_sc = SDHistroy_ / pow(dSPa, 0.25);//峡工程二期围堰低高防渗心墙方案的有限元分析[J]

        if (0.75*dSL_sc > dSL) {
            dEt = Kur_ * dPa * pow(dSPa, Ne_);//卸载
        }

        Double fs  = nph_*prin.z() - prin.x() - csn_;
        Double fsd = fs / rc_;
        Double ftz = prin.z() - tension_;

        if (fsd > 0.0 && fsd >= ftz) {
            shearCorrection(s, &prin, &iPlas, fs);

        }
        else if (ftz > 0.0 && ftz >= fsd) {
            tensionCorrection(s, &prin, &iPlas, ftz, brittle_);

        }
        apexCorrection(friction_, s, &prin, &iPlas, brittle_);
        //根据0<v<0.49确定体积模量的范围0.33~17Et
        bulk_ = Kb_ * dPa * pow(dSPa, Mb_);
        bulk_ = (bulk_ < 0.33 * dEt) ? 0.33333333333333333333333333333333 * dEt : bulk_;
        bulk_ = (bulk_ > 17.0 * dEt) ? 17.0 * dEt : bulk_;
        shear_ = (9 * bulk_ - dEt <= 0) ? 0 : 3 * bulk_ * dEt / (9 * bulk_ - dEt);  // Divide!      

        if (iPlas) {
            s->stnS_ = info.resolve(prin);
            s->viscous_ = false;
        }
    }

    bool ModelMc_123::updateParameters() {
        Double rsin = std::sin(friction_ * degrad);
        nph_ = (1.0 + rsin) / (1.0 - rsin);
        csn_ = 2.0 * cohesion_ * sqrt(nph_);
        rsin = std::sin(dilation_ * degrad);
        nps_ = (1.0 + rsin) / (1.0 - rsin);
        rc_  = std::sqrt(1.0 + nph_*nph_);
        //
        if (friction_) {
            Double apex = cohesion_ / std::tan(friction_ * degrad);
            tension_ = std::min(tension_,apex);
        }
        //
        ModelElastic::updateParameters();
        //
        Double ra = e1_ - nps_ * e2_;
        Double rb = e2_ - nps_ * e1_;
        Double rd = ra - rb * nph_;
        sc1_  = ra / rd;
        sc3_  = rb / rd;
        sc2_  = e2_ * (1.0 - nps_) / rd;
        //
        return true;
    }

    bool ModelMc_123::updateParameters(bool bEUpdated, Double *sf1, Double *sf3) {
        if (cohesion_ < 0.0)
            throw std::runtime_error("Mc_123 model: cohesion is not allowed less than 0.");

        if (friction_ >= 90.0 || friction_ < 0.0)
            throw std::runtime_error("Mc_123 model: friction angle is not in the valid range of 0 to 90.");

        if (dilation_ > friction_)
            throw std::runtime_error("Mc_123 model: dilationn angle is not allowed greater than friction angle.");

        Double rsin = std::sin(friction_ * degrad);
        nph_ = (1.0 + rsin) / (1.0 - rsin);
        csn_ = 2.0 * cohesion_ * sqrt(nph_);
        rsin = std::sin(dilation_ * degrad);
        nps_ = (1.0 + rsin) / (1.0 - rsin);
        rc_  = std::sqrt(1.0 + nph_*nph_);
        //
        if (friction_) {
            Double apex = cohesion_ / std::tan(friction_*degrad);
            tension_ = std::min(tension_,apex);
        }
        //
        if (!bEUpdated) ModelElastic::updateParameters();
        //
        Double ra = e1_ - nps_*e2_;
        Double rb = e2_ - nps_*e1_;
        Double rd = ra  - rb*nph_;
        sc1_ = ra/rd;
        sc3_ = rb/rd;
        sc2_ = e2_*(1.0 - nps_)/rd;
        //
        if (sf1) *sf1 = -1.0/rd;
        if (sf3) *sf3 = nps_/rd;
        //
        return !bEUpdated;
    }
// excerpt-run-end

    Double ModelMc_123::moduliReduction(const Double &factor) {
        Double shear_new = ModelElastic::moduliReduction(factor);
        Double ra = e1_ - nps_*e2_;
        Double rb = e2_ - nps_*e1_;
        Double rd = ra - rb*nph_;
        sc1_  = ra/rd;
        sc3_  = rb/rd;
        sc2_  = e2_*(1.0 - nps_)/rd;
        return shear_new;
    }

    void ModelMc_123::apexCorrection(const Double &fric,State *s,DVect3 *prin,UInt *iPlasticity,bool bBrittle) {
        if (fric > 0.0) {
            Double apex = cohesion_ / tan(fric*degrad);
            if (prin->x()>=apex || prin->y()>=apex || prin->z()>=apex) {
                if(iPlasticity) *iPlasticity = 4;
                s->state_ |= tension_now;
                if (bBrittle) tension_ = 0.0;
                prin->rx() = apex;
                prin->ry() = apex;
                prin->rz() = apex;
            }
        }
    }

    void ModelMc_123::tensionCorrection(State *s,DVect3 *prin,UInt *iPlasticity,const Double &ftz,bool bBrittle) {
        s->state_ |= tension_now;
        if (bBrittle) tension_ = 0.0;
        Double ftx = prin->x() - tension_;
        Double fty = prin->y() - tension_;
        if (ftx > 0.0) {
            if(iPlasticity) *iPlasticity = 4;
            prin->rx()  = tension_;
            prin->ry()  = tension_;
            prin->rz()  = tension_;
        } else if (fty > 0.0) {
            if(iPlasticity) *iPlasticity = 3;
            prin->rx() -= (fty + ftz) * e2_/(e1_+e2_);
            prin->ry()  = tension_;
            prin->rz()  = tension_;
        } else {
            if(iPlasticity) *iPlasticity = 2;
            Double tco = ftz * e2_/e1_;
            prin->rx() -= tco;
            prin->ry() -= tco;
            prin->rz()  = tension_;
        }
    }

    void ModelMc_123::shearCorrection(State *s,DVect3 *prin,UInt *iPlasticity,const Double &fs) {
        if(iPlasticity) *iPlasticity = 1;
// excerpt-state2-start
        s->state_ |= shear_now;
// excerpt-state2-end
        prin->rx() += fs * sc1_;
        prin->ry() += fs * sc2_;
        prin->rz() += fs * sc3_;
    }

    Double ModelMc_123::getStressStrengthRatio(const SymTensor &st) const {
        DVect3 prin = st.getEigenInfo();
        Double rat = 10.0;
        Double tanf = std::tan(friction_*degrad); 
        Double tcut = friction_ ? std::min(tension_,(cohesion_/tanf)) : tension_; 
        if (tcut - prin.z() <= 0.0)
            rat = 0.0;
        else {
            Double sinf = std::sin(friction_*degrad);
            Double denom = 1.0 - sinf;
            Double nph = limits<Double>::max();
            if (denom) nph = (1.0 + sinf) / denom;
            Double sig1f = nph*prin.z() - 2.0*cohesion_*std::sqrt(nph);
            denom = prin.z() - prin.x();
            if (denom) rat = (prin.z() - sig1f) / denom;
        }
        rat = std::min(rat,10.0);
        return rat;
    }

    void ModelMc_123::scaleProperties(const Double &scale,const std::vector<UInt> &props) {
        for (UInt u=0;u<props.size();++u) {
            switch (props[u]) {
            case 5: cohesion_ *= scale;  break;
            case 6: friction_ = std::max(0.0,std::min(89.0,std::atan(std::tan(friction_*degrad)*scale)/degrad));  break;
            case 7: dilation_ = std::max(0.0,std::min(89.0,std::atan(std::tan(dilation_*degrad)*scale)/degrad));  break;
            case 8: tension_  *= scale;  break;
            }
            setValid(0);
        }
    }

} // namespace models
// EOF

// EOF
