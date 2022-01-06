#ifndef FhtAna_h
#define FhtAna_h

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
#include "Identifier/WpID.h"
#include "Geometry/RecGeomSvc.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVector2.h"
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
#include "TH1F.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <limits.h>

#define PI TMath::Pi()

using namespace std;

class CdGeom;
class WpGeom;

class FhtAna : public AlgBase {
    public:
		FhtAna(const std::string&);
		bool initialize();
		bool execute();
		bool initGeomSvc();
		bool initPmt();
		bool freshPmtData(TH2D*, TH2D*, TH2D*, TH2D*, TH2D*, double&, double&);
		bool finalize();
		bool IfCrossCd(TVector3&, TVector3&, Double_t);
		TVector3 InciOnLS(TVector3&, TVector3&, Double_t);
		TVector3 PosOnLS(TVector3&, TVector3&, Double_t, int);
		TVector3 GetInciPos(TH1D*, TH1D*, int);
		TVector3 GetExitPos(TH1D*, TH1D*, int);
		TVector3 GetChargeCenter();
		TH2D* MapSmooth(TH2D*, int, int, TString);
		int* GetMassPos(TH2D*, TH2D*, int, int, TH2D*);
		TH2D* PECut(TH2D*, int, int, double);
		void nCorrosion(TH2D*, int, int, int);
		int MarkConnection(TH2D*, int, int, TH2D*, int);
		int AreaCut(TH2D*, TH2D*, int, int, double, bool, bool);
		bool FindTrk(TVector3&, TVector3&, double&, double&, double&, TH2D*, long int*);
		bool FillContent(TH2D*, int, int);
		bool ChooseCut(TH2D*, TH1D*, int, int);
		bool Expansion(TH2D*, int, int);
		bool RMSMap(TH2D*, TH2D*, int, int, int, int, double);
		bool MapExtend(TH2D*, TH2D*, int, int);
		bool Pool(TH2D*, TH2D*, int, int, int);
		bool XOR(TH2D*, TH2D*, int, int);
		TH2D* Combine(TH2D*, TH2D*, int, int);
		bool UnionCut(TH2D*, TH2D*, TH2D*, int, int, double, TH2D*);
		bool AND(TH2D*, TH2D*, int, int);
		long int* GetCenterPos(TH2D*, TH2D*, int, int);
		double FHTPredict(const PmtProp&, TVector3, TVector3, double);
    private:
		char* outPath;
		char* m_name;
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
		WpGeom* m_wpgeom;
		PmtTable m_ptab;
		unsigned int totPmtNum;
        Double_t m_3inchRes;
        Double_t m_20inchRes;
		bool m_3inchusedflag;
		bool m_20inchusedflag;
		Double_t m_qcut;
		void Corrosion(TH2D*, int, int);
		double m_earlist;
};

#endif
