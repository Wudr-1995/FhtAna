#ifndef GetFht_h
#define GetFht_h

#include "EvtNavigator/NavBuffer.h"
#include "SniperKernel/AlgFactory.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/RecHeader.h"
#include "Event/CalibHeader.h"
#include "Event/ElecHeader.h"
#include "Event/SimHeader.h"
#include "RootWriter/RootWriter.h"
#include "SniperKernel/AlgBase.h"
#include "EvtNavigator/NavBuffer.h"
#include "RootWriter/RootWriter.h"
#include "Identifier/Identifier.h"
#include "Identifier/CdID.h"
#include "Geometry/RecGeomSvc.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include "PmtProp.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include <algorithm>

using namespace std;

class CdGeom;

class GetFht : public AlgBase {
    public:
		GetFht(const std::string&);
		bool initialize();
		bool execute();
		bool initGeomSvc();
		bool initPmt();
		bool freshPmtData(TH2D*, TH2D*, TH2D*, TH2D*, double&, double&);
		bool finalize();
		bool IfCrossCd(TVector3&, TVector3&, Double_t);
		TVector3 InciOnLS(TVector3&, TVector3&, Double_t);
		TH1D* FixCurve(TH1D*, int);
		int ValleyFinder(TH1D*, int);
		TH1D* GetDiff(TH1D*, int);
		TH1D* GetDiffInte(TH1D*, int);
		TVector3 GetInciPos(TH1D*, TH1D*, int);
		TVector3 GetExitPos(TH1D*, TH1D*, int);
		TVector3 GetChargeCenter();
		TH2D* MapSmooth(TH2D*, int, int, TString);
		TVector3* GetMassPos(TH2D*, int, int, TString);
    private:
		int nSimTrks;
		int m_turn;
		TString m_path;
		JM::NavBuffer* m_buf;
		std::ofstream outFile;
		std::ofstream label;
		std::ofstream trks;
		std::ofstream PhiOut;
		// std::ifstream nEvt;
		// TTree* m_tree;
		// TH2D* m_hist;
		int m_iEvt;
		int m_usedPmtNum;
		Double_t m_LSRadius;
        CdGeom* m_geom;
		PmtTable m_ptab;
		unsigned int totPmtNum;
        Double_t m_3inchRes;
        Double_t m_20inchRes;
		bool m_3inchusedflag;
		bool m_20inchusedflag;
		Double_t m_qcut;
};

#endif
