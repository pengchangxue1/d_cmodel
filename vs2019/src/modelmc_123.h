#pragma once

#include "models/elastic/src/modelelastic.h"

namespace models {
    class ModelMc_123 : public ModelElastic {
    public:
        ModelMc_123(unsigned short option=0);
        virtual String        getName() const override;
        virtual String        getFullName() const override;
        virtual UInt          getMinorVersion() const override; 
        virtual String        getProperties() const override;
        virtual String        getStates() const override;
        virtual Variant       getProperty(UInt index) const override;
        virtual void          setProperty(UInt index,const Variant &p,UInt restoreVersion=0) override;
        virtual ModelMc_123 *clone() const override { return new ModelMc_123(); }      
        virtual void          copy(const ConstitutiveModel *mod) override;
        virtual void          run(UByte dim,State *s) override; 
        virtual void          initialize(UByte dim,State *s) override;
        // Optional
        virtual bool          isPropertyAdvanced(UInt i) const override;
        virtual bool          supportsStressStrengthRatio() const override { return true; }
        virtual bool          supportsPropertyScaling() const override { return true; }
        virtual Double        getStressStrengthRatio(const SymTensor &st) const override;
        virtual void          scaleProperties(const Double &scale,const std::vector<UInt> &props) override;

    private:
        virtual bool          updateParameters() override;
        virtual bool          updateParameters(bool bEUpdated,Double *sf1=nullptr,Double *sf3=nullptr);
        virtual Double        moduliReduction(const Double &factor) override;
        virtual void          apexCorrection(const Double &fric,State *s,DVect3 *prin,UInt *iPlasticity=nullptr,bool bBrittle=false);
        virtual void          tensionCorrection(State *s,DVect3 *prin,UInt *iPlasticity,const Double &ftz,bool bBrittle=false);
        virtual void          shearCorrection(State *s,DVect3 *prin,UInt *iPlasticity,const Double &fs);
        
        Double cohesion_ = 0.0;
        Double friction_ = 0.0;
        Double dilation_ = 0.0;
        Double tension_ = 0.0;
        bool brittle_ = false;
        Double nph_  = 0.0;
        Double csn_ = 0.0;
        Double nps_ = 0.0;
        Double rc_ = 0.0;
        Double sc1_ = 0.0;
        Double sc2_ = 0.0;
        Double sc3_ = 0.0;


        // 邓肯模型的参数列表
      //分别对应c,fai,Rf破坏比,K初始模量系数,n初始模量指数,
      //Kb,nb为切线泊松比系数，Kur回弹模量
      //G,F初始泊松比系数,D切线泊松比系数,Kur回弹模量
        Double FricDel_ = 0.0;          //delt摩擦角
        Double FailRatio_ = 0.0;        //Rf
        Double Ke_ = 0.0;               //K
        Double Ne_ = 0.0;               //n
        Double Kb_ = 0.0;               //Kb
        Double Mb_ = 0.0;               //m
        Double Kur_ = 0.0;              //Kur
        Double F_ = 0.0;                //修正摩擦角
        Double Fcos_ = 0.0;
        Double Fsin_ = 0.0;
        Double Shear_ = 0.0;
        Double Bulk_ = 0.0;
        Double SDHistroy_ = 0.0;
        Double S3Histroy_ = 0.0;
        Double SLHistroy_ = 0.0;

        //Double dSL = 0;

    };
} // namespace models

// EOF
