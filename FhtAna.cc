#include "FhtAna.h"
#include "EvtNavigator/NavBuffer.h"
#include "SniperKernel/AlgFactory.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/RecHeader.h"
#include "Event/CalibHeader.h"
#include "Event/ElecHeader.h"
#include "Event/SimHeader.h"
#include "RootWriter/RootWriter.h"
#include "TMath.h"
#include "TArrow.h"

DECLARE_ALGORITHM(FhtAna);

std::ostream& operator << (std::ostream& s, const TVector3& v){
	s << "(" << v.x() <<  "," << v.y() << "," << v.z() << ")";
	return s;
}

FhtAna::FhtAna(const std::string& name)
: AlgBase(name),
m_iEvt(0),
m_buf(0),
m_usedPmtNum(0)
{
	declProp("ChargeCut", m_qcut = 0);
	declProp("Use3inchPmt", m_3inchusedflag = false);
	declProp("Use20inchPmt", m_20inchusedflag = true);
	declProp("Pmt3inchTimeReso", m_3inchRes = 1);
	declProp("Pmt20inchTimeReso", m_20inchRes = 8);
	declProp("LSRadius", m_LSRadius = 17700);
	declProp("FilePath", outPath = "");
	declProp("FileName", m_name = "");
	declProp("FileNumber", m_turn = 0);
}

bool FhtAna::initialize() {
	LogDebug << "Initializing" << std::endl;
	if(not initGeomSvc())
		return false;
	SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
	if (navBuf.invalid()) {
		LogError << "Cannot get the NavBuffer @ /Event" << std::endl;
		return false;
	}
	m_buf = navBuf.data();
	m_path = outPath;
    return true;
}

bool FhtAna::execute() {
	LogDebug << "executing: " << m_iEvt ++ << std::endl;
	if (m_iEvt < 2)
		return true;
	TH2D* FhtDis = new TH2D("FhtDistribution", "FhtDistribution", 100, 0, PI, 500, 0, 500);
	TH2D* FhtPhi = new TH2D("FhtVPhi", "FhtVPhi", 200, -PI, PI, 500, 0, 500);
	TH2D* Fht2D = new TH2D("FhtDistribution2D", "", 100, 0, PI, 200, -PI, PI);
	TH2D* Q2D = new TH2D("ChargeDistribution2D", "", 100, 0, PI, 200, -PI, PI);
	TH2D* nPMT = new TH2D("nPMT", "", 100, 0, PI, 200, - PI, PI);

	TH1F* FhtDiff = new TH1F("FhtDiff", "", 2000, -100, 100);

	JM::SimEvent* simevent = 0;
	JM::EvtNavigator* nav =m_buf->curEvt();
	std::vector<std::string>& paths = nav->getPath();
	JM::SimHeader* simheader = 0;
	for (size_t i = 0; i < paths.size(); ++i) {
		const std::string& path = paths[i];
		if (path == "/Event/SimOrig") {
			simheader = static_cast<JM::SimHeader*>(nav->getHeader("/Event/SimOrig"));
			LogDebug << "SimHeader (/Event/SimOrig): " << simheader << endl;
			if (simheader)
				break;
		}
	}
	if (initPmt())
		LogDebug << "Initializing PMT success" << std::endl;
	else {
		LogError << "Initializing PMT fails" << std::endl;
		return true;
	}
	double InciTheta, InciPhi;
	if (freshPmtData(FhtDis, FhtPhi, Fht2D, Q2D, nPMT, InciTheta, InciPhi))
		LogDebug << "Freshing PMT data success" << std::endl;
	else {
		LogError << "Freshing PMT data fails" << std::endl;
		delete FhtDis;
		delete FhtPhi;
		delete Fht2D;
		delete Q2D;
		return true;
	}

	double tmpN = nPMT->GetMaximum();
	for (int i = 1; i <= 100; i ++)
		for (int j = 1; j <= 200; j ++)
			nPMT->SetBinContent(i, j, nPMT->GetBinContent(i, j) / tmpN);
	
	// FillContent(nPMT, 100, 200);

	TString pdfPath = m_path + "pdf/" + m_name + "_" + m_turn + "_" + m_iEvt + ".pdf";
	TString txtPath = m_path + m_name + "_" + m_turn + "_" + m_iEvt + ".txt";
	ofstream of(txtPath);

	auto c1 = new TCanvas("Fht", "", 800, 800);
	gStyle->SetOptStat(0000);
	gStyle->SetPalette(1);
	c1->SetRightMargin(0.15);
	c1->SetBottomMargin(0.15);
	c1->SetLeftMargin(0.15);
	c1->SetTopMargin(0.15);
	c1->Print(pdfPath + "[");

	TH2D* pmt = (TH2D*)nPMT->Clone("npmt");
	pmt->GetXaxis()->SetTitle("Theta / Radian");
	pmt->GetYaxis()->SetTitle("Phi / Radian");
	pmt->GetXaxis()->SetTitleSize(0.05);
	pmt->GetYaxis()->SetTitleSize(0.05);
	pmt->GetXaxis()->SetLabelSize(0.05);
	pmt->GetYaxis()->SetLabelSize(0.05);
	c1->cd();
	pmt->Draw("colz");
	c1->Print(pdfPath);

	TH2D* orig = (TH2D*)Q2D->Clone("ori");
	orig->GetXaxis()->SetTitle("Theta / Radian");
	orig->GetYaxis()->SetTitle("Phi / Radian");
	orig->GetXaxis()->SetTitleSize(0.05);
	orig->GetYaxis()->SetTitleSize(0.05);
	orig->GetXaxis()->SetLabelSize(0.05);
	orig->GetYaxis()->SetLabelSize(0.05);
	c1->cd();
	orig->Draw("colz");
	c1->Print(pdfPath);

	for (int i = 1; i <= 100; i ++)
		for (int j = 1; j <= 200; j ++)
			if (nPMT->GetBinContent(i, j))
				Q2D->SetBinContent(i, j, Q2D->GetBinContent(i, j) / nPMT->GetBinContent(i, j));

	auto exQ2D = new TH2D("ExQ2D", "", 120, -0.1 * PI, 1.1 * PI, 220, -1.1 * PI, 1.1 * PI);
	MapExtend(exQ2D, Q2D, 100, 200);
	for (int i = 1; i <= 120; i ++) {
		for (int j = 1; j <= 220; j ++) {
			of << exQ2D->GetBinContent(i, j) << "\t";
		}
		of << endl;
	}
	auto exFht2D = new TH2D("ExFht2D", "", 120, -0.1 * PI, 1.1 * PI, 220, -1.1 * PI, 1.1 * PI);
	MapExtend(exFht2D, Fht2D, 100, 200);
	for (int i = 1; i <= 120; i ++) {
		for (int j = 1; j <= 220; j ++) {
			of << exFht2D->GetBinContent(i, j) << "\t";
		}
		of << endl;
	}

	/*
	c1->cd();
	Q2D->GetXaxis()->SetTitle("Theta / Radian");
	Q2D->GetYaxis()->SetTitle("Phi / Radian");
	Q2D->GetXaxis()->SetTitleSize(0.05);
	Q2D->GetYaxis()->SetTitleSize(0.05);
	Q2D->GetXaxis()->SetLabelSize(0.05);
	Q2D->GetYaxis()->SetLabelSize(0.05);
	Q2D->Draw("colz");
	c1->Print(pdfPath);
	*/

	TString na("ChargeSmoothed");
	// FillContent(Q2D, 100, 200);
	for (int i = 0; i < 4; i ++)
		Expansion(Q2D, 100, 200);

	/*
	c1->cd();
	Q2D->GetXaxis()->SetTitle("Theta / Radian");
	Q2D->GetYaxis()->SetTitle("Phi / Radian");
	Q2D->GetXaxis()->SetTitleSize(0.05);
	Q2D->GetYaxis()->SetTitleSize(0.05);
	Q2D->GetXaxis()->SetLabelSize(0.05);
	Q2D->GetYaxis()->SetLabelSize(0.05);
	Q2D->Draw("colz");
	c1->Print(pdfPath);
	*/

	TH2D* Q2Smooth = MapSmooth(Q2D, 100, 200, na);

	/*
	TH2D* step1Q = (TH2D*)Q2Smooth->Clone("Step1Q");
	step1Q->GetXaxis()->SetTitle("Theta / Radian");
	step1Q->GetYaxis()->SetTitle("Phi / Radian");
	step1Q->GetXaxis()->SetTitleSize(0.05);
	step1Q->GetYaxis()->SetTitleSize(0.05);
	step1Q->GetXaxis()->SetLabelSize(0.05);
	step1Q->GetYaxis()->SetLabelSize(0.05);
	c1->cd();
	step1Q->Draw("colz");
	c1->Print(pdfPath);
	*/

	TH2D* Q2Pool = new TH2D("PoolQSmooth", "", 12, - 0.1 * PI, 1.1 * PI, 22, - 1.1 * PI, 1.1 * PI);
	Pool(Q2Smooth, Q2Pool, 120, 220, 10);

	TH2D* RMSPool = new TH2D("RMS", "", 10, 0, PI, 20, -PI, PI);
	RMSMap(Q2Pool, RMSPool, 12, 22, 1, 1, 1.5E6);
	// TH2D* QPCut = PECut(Q2Pool, 12, 22, 0.6);
	delete Q2Pool;
	delete RMSPool;

	TH2D* RMS = new TH2D("RMS", "", 100, 0, PI, 200, - PI, PI);
	RMSMap(Q2Smooth, RMS, 120, 220, 10, 3, 5E4);

	/*
	c1->cd();
	RMS->GetXaxis()->SetTitle("Theta / Radian");
	RMS->GetYaxis()->SetTitle("Phi / Radian");
	RMS->GetXaxis()->SetTitleSize(0.05);
	RMS->GetYaxis()->SetTitleSize(0.05);
	RMS->GetXaxis()->SetLabelSize(0.05);
	RMS->GetYaxis()->SetLabelSize(0.05);
	RMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	TH2D* exRMS = new TH2D("exRMS", "", 120, - PI * 0.1, PI * 1.1, 220, - PI * 1.1, PI * 1.1);
	MapExtend(exRMS, RMS, 100, 200);

	delete RMS;

	/*
	c1->cd();
	exRMS->GetXaxis()->SetTitle("Theta / Radian");
	exRMS->GetYaxis()->SetTitle("Phi / Radian");
	exRMS->GetXaxis()->SetTitleSize(0.05);
	exRMS->GetYaxis()->SetTitleSize(0.05);
	exRMS->GetXaxis()->SetLabelSize(0.05);
	exRMS->GetYaxis()->SetLabelSize(0.05);
	exRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	// TH2D* Q2HCut = PECut(Q2Smooth, 120, 220, 0.8);
	TH2D* R2HCut = PECut(exRMS, 120, 220, 0.8);
	// TH2D* Q2LCut = PECut(Q2Smooth, 120, 220, 0.35);
	TH2D* R2LCut = PECut(exRMS, 120, 220, 0.35);

	/*
	c1->cd();
	R2LCut->GetXaxis()->SetTitle("Theta / Radian");
	R2LCut->GetYaxis()->SetTitle("Phi / Radian");
	R2LCut->GetXaxis()->SetTitleSize(0.05);
	R2LCut->GetYaxis()->SetTitleSize(0.05);
	R2LCut->GetXaxis()->SetLabelSize(0.05);
	R2LCut->GetYaxis()->SetLabelSize(0.05);
	R2LCut->Draw("colz");
	c1->Print(pdfPath);
	*/

	/*
	c1->cd();
	R2HCut->GetXaxis()->SetTitle("Theta / Radian");
	R2HCut->GetYaxis()->SetTitle("Phi / Radian");
	R2HCut->GetXaxis()->SetTitleSize(0.05);
	R2HCut->GetYaxis()->SetTitleSize(0.05);
	R2HCut->GetXaxis()->SetLabelSize(0.05);
	R2HCut->GetYaxis()->SetLabelSize(0.05);
	R2HCut->Draw("colz");
	c1->Print(pdfPath);
	*/

	// TH2D* cHQ = new TH2D("cHQ", "", 120, - 0.1 * PI, 1.1 * PI, 220, - 1.1 * PI, 1.1 * PI);
	// MarkConnection(Q2HCut, 120, 220, cHQ, 20);
	// TH2D* cLQ = new TH2D("cLQ", "", 120, - 0.1 * PI, 1.1 * PI, 220, - 1.1 * PI, 1.1 * PI);
	// MarkConnection(Q2LCut, 120, 220, cLQ, 20);
	TH2D* cHRMS = new TH2D("cHRMS", "", 120, - 0.1 * PI, 1.1 * PI, 220, - 1.1 * PI, 1.1 * PI);
	MarkConnection(R2HCut, 120, 220, cHRMS, 20);
	TH2D* cLRMS = new TH2D("cLRMS", "", 120, - 0.1 * PI, 1.1 * PI, 220, - 1.1 * PI, 1.1 * PI);
	MarkConnection(R2LCut, 120, 220, cLRMS, 20);

	/*
	c1->cd();
	cLRMS->GetXaxis()->SetTitle("Theta / Radian");
	cLRMS->GetYaxis()->SetTitle("Phi / Radian");
	cLRMS->GetXaxis()->SetTitleSize(0.05);
	cLRMS->GetYaxis()->SetTitleSize(0.05);
	cLRMS->GetXaxis()->SetLabelSize(0.05);
	cLRMS->GetYaxis()->SetLabelSize(0.05);
	cLRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	/*
	c1->cd();
	cHRMS->GetXaxis()->SetTitle("Theta / Radian");
	cHRMS->GetYaxis()->SetTitle("Phi / Radian");
	cHRMS->GetXaxis()->SetTitleSize(0.05);
	cHRMS->GetYaxis()->SetTitleSize(0.05);
	cHRMS->GetXaxis()->SetLabelSize(0.05);
	cHRMS->GetYaxis()->SetLabelSize(0.05);
	cHRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	AreaCut(R2HCut, cHRMS, 120, 220, 0.3, false, true);

	/*
	c1->cd();
	cHRMS->GetXaxis()->SetTitle("Theta / Radian");
	cHRMS->GetYaxis()->SetTitle("Phi / Radian");
	cHRMS->GetXaxis()->SetTitleSize(0.05);
	cHRMS->GetYaxis()->SetTitleSize(0.05);
	cHRMS->GetXaxis()->SetLabelSize(0.05);
	cHRMS->GetYaxis()->SetLabelSize(0.05);
	cHRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	TH2D* test1 = new TH2D("test1", "", 120, - 0.1 * PI, 1.1 * PI, 220, - 1.1 * PI, 1.1 * PI);
	if (!UnionCut(cLRMS, cHRMS, R2LCut, 120, 220, 0.75, test1)) {
		LogInfo << "Error in UnionCut()" << endl;
		return true;
	}

	/*
	c1->cd();
	test1->GetXaxis()->SetTitle("Theta / Radian");
	test1->GetYaxis()->SetTitle("Phi / Radian");
	test1->GetXaxis()->SetTitleSize(0.05);
	test1->GetYaxis()->SetTitleSize(0.05);
	test1->GetXaxis()->SetLabelSize(0.05);
	test1->GetYaxis()->SetLabelSize(0.05);
	test1->Draw("colz");
	c1->Print(pdfPath);
	*/

	MarkConnection(R2LCut, 120, 220, cLRMS, 20);

	/*
	c1->cd();
	R2LCut->GetXaxis()->SetTitle("Theta / Radian");
	R2LCut->GetYaxis()->SetTitle("Phi / Radian");
	R2LCut->GetXaxis()->SetTitleSize(0.05);
	R2LCut->GetYaxis()->SetTitleSize(0.05);
	R2LCut->GetXaxis()->SetLabelSize(0.05);
	R2LCut->GetYaxis()->SetLabelSize(0.05);
	R2LCut->Draw("colz");
	c1->Print(pdfPath);
	*/

	/*
	c1->cd();
	cLRMS->GetXaxis()->SetTitle("Theta / Radian");
	cLRMS->GetYaxis()->SetTitle("Phi / Radian");
	cLRMS->GetXaxis()->SetTitleSize(0.05);
	cLRMS->GetYaxis()->SetTitleSize(0.05);
	cLRMS->GetXaxis()->SetLabelSize(0.05);
	cLRMS->GetYaxis()->SetLabelSize(0.05);
	cLRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	// MarkConnection(R2HCut, 120, 220, cHRMS, 20);
	AreaCut(R2HCut, cHRMS, 120, 220, 0.3, true, false);
	AreaCut(R2LCut, cLRMS, 120, 220, 0.3, true, true);

	/*
	c1->cd();
	cLRMS->GetXaxis()->SetTitle("Theta / Radian");
	cLRMS->GetYaxis()->SetTitle("Phi / Radian");
	cLRMS->GetXaxis()->SetTitleSize(0.05);
	cLRMS->GetYaxis()->SetTitleSize(0.05);
	cLRMS->GetXaxis()->SetLabelSize(0.05);
	cLRMS->GetYaxis()->SetLabelSize(0.05);
	cLRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	/*
	c1->cd();
	R2LCut->GetXaxis()->SetTitle("Theta / Radian");
	R2LCut->GetYaxis()->SetTitle("Phi / Radian");
	R2LCut->GetXaxis()->SetTitleSize(0.05);
	R2LCut->GetYaxis()->SetTitleSize(0.05);
	R2LCut->GetXaxis()->SetLabelSize(0.05);
	R2LCut->GetYaxis()->SetLabelSize(0.05);
	R2LCut->Draw("colz");
	c1->Print(pdfPath);
	*/

	/*
	c1->cd();
	cHRMS->GetXaxis()->SetTitle("Theta / Radian");
	cHRMS->GetYaxis()->SetTitle("Phi / Radian");
	cHRMS->GetXaxis()->SetTitleSize(0.05);
	cHRMS->GetYaxis()->SetTitleSize(0.05);
	cHRMS->GetXaxis()->SetLabelSize(0.05);
	cHRMS->GetYaxis()->SetLabelSize(0.05);
	cHRMS->Draw("colz");
	c1->Print(pdfPath);
	*/

	/*
	c1->cd();
	R2HCut->GetXaxis()->SetTitle("Theta / Radian");
	R2HCut->GetYaxis()->SetTitle("Phi / Radian");
	R2HCut->GetXaxis()->SetTitleSize(0.05);
	R2HCut->GetYaxis()->SetTitleSize(0.05);
	R2HCut->GetXaxis()->SetLabelSize(0.05);
	R2HCut->GetYaxis()->SetLabelSize(0.05);
	R2HCut->Draw("colz");
	c1->Print(pdfPath);
	*/

	TH2D* totMark = Combine(cHRMS, cLRMS, 120, 220);

	/*
	c1->cd();
	totMark->GetXaxis()->SetTitle("Theta / Radian");
	totMark->GetYaxis()->SetTitle("Phi / Radian");
	totMark->GetXaxis()->SetTitleSize(0.05);
	totMark->GetYaxis()->SetTitleSize(0.05);
	totMark->GetXaxis()->SetLabelSize(0.05);
	totMark->GetYaxis()->SetLabelSize(0.05);
	totMark->Draw("colz");
	c1->Print(pdfPath);
	*/

	auto Q = new TH1D("Charge", "", 500, 0, 5E6);
	ChooseCut(exRMS, Q, 120, 220);

	delete exRMS;
	delete cLRMS;
	delete cHRMS;
	delete R2HCut;
	delete R2LCut;

	/*
	c1->cd();
	Q->Draw();
	gPad->SetLogy();
	c1->Print(pdfPath);
	gPad->SetLogy(0);
	*/

	if (Q2Smooth == NULL) {
		LogDebug << "Lower than the PE threshold, skip" << endl;
		delete FhtDis;
		delete FhtPhi;
		delete Fht2D;
		delete Q2D;
		delete Q2Smooth;
		return true;
	}

	// PECut(Q2Smooth, 120, 220, 0.35);

	// nCorrosion(Q2Smooth, 120, 220, 2);

	// int nC = MarkConnection(Q2Smooth, 120, 220, conn, 20);

	// nC = AreaCut(Q2Smooth, conn, 120, 220, 0.5);

	long int* mass = new long int[4];
	// TH2D* testM = new TH2D("testMass", "", 1000, 0, 5000, 1000, -400, 400);
	mass = GetCenterPos(Q2Smooth, totMark, 120, 220);
	// mass = GetMassPos(Q2Smooth, totMark, 120, 220, testM);
	
	delete totMark;

	double unit = PI / 100;
	LogInfo << "==================================================" << endl;

	TVector3 rInci, rDir;
	double rDis, rAng, rTi;
	FindTrk(rInci, rDir, rDis, rAng, rTi, Fht2D, mass);
	LogInfo << "PreRec Inci.Theta: " << rInci.Theta() << "\tPhi: " << rInci.Phi() << endl;
	LogInfo << "PreRec Dir.Theta: " << rDir.Theta() << "\tPhi: " << rDir.Phi() << endl;

	LogInfo << "==================================================" << endl;

	TVector3 chargeCenter = GetChargeCenter();

	simevent = dynamic_cast<JM::SimEvent*>(simheader->event());
	if (not simevent) {
		LogInfo << "No sim event" << endl;
		return true;
	}
	LogInfo << "SimEventGot" << std::endl;
	nSimTrks = simevent->getTracksVec().size();
	LogInfo << "Retrieving tracks data" << std::endl;
	short NumCrossCd = 0;
	LogInfo << "Number of Trks: " << nSimTrks << endl;
	double lX = 50000, lY = 50000, lZ = 50000;

	TH2D* exp2D = new TH2D("FhtExp2D", "", 100, 0, PI, 200, -PI, PI);
	TH2D* pos = new TH2D("pos", "", 500, - 25000, 25000, 500, - 25000, 25000);
	TH2D* LiDiff = new TH2D("LiDiff", "", 500, 0, 10000, 200, - 100, 100);
	TH2D* QDiff = new TH2D("QDiff", "", 50, 0, 50, 200, - 100, 100);
	TH2D* TDiff = new TH2D("TDiff", "", 100, 0, 100, 200, - 100, 100);

	for (short i = 1; i <= 1; i ++) {
		JM::SimTrack* strk = simevent->findTrackByTrkID(i);
		LogInfo << "======================================== ID: " << i << endl;
		// if (lX == strk->getInitX() && lY == strk->getInitY() && lZ == strk->getInitZ())
		// 	continue;
		lX = strk->getInitX();
		lY = strk->getInitY();
		lZ = strk->getInitZ();
		TVector3 Inci(strk->getInitX(), strk->getInitY(), strk->getInitZ());
		TVector3 Dir(strk->getInitPx(), strk->getInitPy(), strk->getInitPz());
		TVector3 Exit(strk->getExitX(), strk->getExitY(), strk->getExitZ());
		LogInfo << "Inci Pos: " << Inci << endl;
		if (IfCrossCd(Inci, Dir, m_LSRadius)) {
			NumCrossCd ++;
			TVector3 LSInci = InciOnLS(Inci, Dir, m_LSRadius);
			// IfCrossCd(LSInci, Dir, m_LSRadius);
			TVector3 LSExit = Exit;
			if (Exit.Mag() > m_LSRadius) {
				TVector3 antiDir = - Dir;
				LSExit = InciOnLS(Exit, antiDir, m_LSRadius);
			}
			TVector3 dir = Dir.Unit();
			of << LSInci.Theta() << "\t" << LSInci.Phi() << "\t" << dir.Theta() << "\t" << dir.Phi() << endl;
			LogInfo << "Inci: " << LSInci << endl
					<< "Exit: " << LSExit << endl
					<< "Length: " << (LSInci - LSExit).Mag() << endl
					<< "Edep: " << strk->getEdep() << endl
					<< "QEdep: " << strk->getQEdep() << endl;
			LogInfo << "Inci.Theta: " << LSInci.Theta() << "\tPhi: " << LSInci.Phi() << "\tMag: " << LSInci.Mag() << endl;
			LogInfo << "Exit.Theta: " << LSExit.Theta() << "\tPhi: " << Exit.Phi() << endl;
			LogInfo << "Dir.Theta: " << dir.Theta() << "\tPhi: " << dir.Phi() << endl;

			int nPMTs = m_ptab.size();
			Dir = Dir.Unit();
			for (int i = 0; i < nPMTs; i ++) {
				if (m_ptab[i].q < 1 || m_ptab[i].fht > 90)
					continue;

				double nLS = 1.485;
				double cLight = 299.;
				double vMuon = 299.;
				double nW = 1.34;
				double ti = 0;

				double tan = TMath::Sqrt(nW * nW - 1);

				TVector3 perp = Inci + Dir * (m_ptab[i].pos - Inci) * Dir;

				TVector3 liSource = perp - (m_ptab[i].pos - perp).Mag() / tan * Dir;

				TVector3 liDir = m_ptab[i].pos - liSource;

				double dis = TMath::Sqrt(TMath::Power(m_ptab[i].pos.Mag(), 2) - TMath::Power(m_ptab[i].pos * liDir.Unit(), 2));

				double expFht = ti + (liSource - Inci) * Dir / vMuon + (m_ptab[i].pos - liSource).Mag() * nW / cLight;

				int binx = m_ptab[i].pos.Theta() / PI * 100 + 1;
				int biny = m_ptab[i].pos.Phi() / PI * 100 + 101;
				exp2D->SetBinContent(binx, biny, TMath::Abs(expFht - m_ptab[i].fht));
				FhtDiff->Fill(expFht - m_ptab[i].fht);

				if (ti + (liSource - Inci) * Dir / vMuon < 0)
					continue;
				binx = TMath::Sqrt(TMath::Power(liSource.X(), 2) + TMath::Power(liSource.Y(), 2)) / 100 + 251;
				binx = liSource.X() / 100 + 251;
				biny = liSource.Z() / 100 + 251;
				pos->SetBinContent(binx, biny, ti + (liSource - Inci) * Dir / vMuon);

				binx = TMath::Sqrt(TMath::Power(m_ptab[i].pos.X(), 2) + TMath::Power(m_ptab[i].pos.Y(), 2)) / 100 + 251;
				binx = m_ptab[i].pos.X() / 100 + 251;
				biny = m_ptab[i].pos.Z() / 100 + 251;
				pos->SetBinContent(binx, biny, TMath::Abs(expFht - m_ptab[i].fht));

				LiDiff->Fill((m_ptab[i].pos - liSource).Mag(), expFht - m_ptab[i].fht);
				QDiff->Fill(m_ptab[i].q, expFht - m_ptab[i].fht);
				TDiff->Fill(m_ptab[i].fht, expFht - m_ptab[i].fht);
			}
		}
	}
	c1->cd();
	FhtDiff->SetTitle("");
	FhtDiff->GetXaxis()->SetTitleSize(0.05);
	FhtDiff->GetYaxis()->SetTitleSize(0.05);
	FhtDiff->GetXaxis()->SetLabelSize(0.05);
	FhtDiff->GetYaxis()->SetLabelSize(0.05);
	FhtDiff->GetXaxis()->SetTitle("(exp - truth) / ns");
	FhtDiff->GetYaxis()->SetTitle("Count");
	FhtDiff->Draw();
	c1->Print(pdfPath);

	c1->cd();
	exp2D->SetTitle("");
	exp2D->GetXaxis()->SetTitleSize(0.05);
	exp2D->GetYaxis()->SetTitleSize(0.05);
	exp2D->GetXaxis()->SetLabelSize(0.05);
	exp2D->GetYaxis()->SetLabelSize(0.05);
	exp2D->GetXaxis()->SetTitle("Theta / Radian");
	exp2D->GetYaxis()->SetTitle("Phi / Radian");
	exp2D->Draw("colz");
	c1->Print(pdfPath);

	c1->cd();
	pos->SetTitle("");
	pos->GetXaxis()->SetTitleSize(0.05);
	pos->GetYaxis()->SetTitleSize(0.05);
	pos->GetXaxis()->SetLabelSize(0.02);
	pos->GetYaxis()->SetLabelSize(0.05);
	pos->GetXaxis()->SetTitle("sqrt(x^{2} + y^{2}) / mm");
	pos->GetYaxis()->SetTitle("z / mm");
	pos->Draw("colz");
	c1->Print(pdfPath);

	c1->cd();
	LiDiff->SetTitle("");
	LiDiff->GetXaxis()->SetTitleSize(0.05);
	LiDiff->GetYaxis()->SetTitleSize(0.05);
	LiDiff->GetXaxis()->SetLabelSize(0.05);
	LiDiff->GetYaxis()->SetLabelSize(0.05);
	LiDiff->GetXaxis()->SetTitle("Light route / mm");
	LiDiff->GetYaxis()->SetTitle("(exp - truth) / ns");
	LiDiff->Draw("colz");
	c1->Print(pdfPath);

	c1->cd();
	TDiff->SetTitle("");
	TDiff->GetXaxis()->SetTitleSize(0.05);
	TDiff->GetYaxis()->SetTitleSize(0.05);
	TDiff->GetXaxis()->SetLabelSize(0.05);
	TDiff->GetYaxis()->SetLabelSize(0.05);
	TDiff->GetXaxis()->SetTitle("Light route / mm");
	TDiff->GetYaxis()->SetTitle("(exp - truth) / ns");
	TDiff->Draw("colz");
	c1->Print(pdfPath);

	// c1->cd();
	// QDiff->SetTitle("");
	// QDiff->GetXaxis()->SetTitleSize(0.05);
	// QDiff->GetYaxis()->SetTitleSize(0.05);
	// QDiff->GetXaxis()->SetLabelSize(0.05);
	// QDiff->GetYaxis()->SetLabelSize(0.05);
	// QDiff->GetXaxis()->SetTitle("nPE / p.e.");
	// QDiff->GetYaxis()->SetTitle("(exp - truth) / ns");
	// QDiff->Draw("colz");
	// c1->Print(pdfPath);

	c1->cd();
	Fht2D->SetTitle("");
	Fht2D->GetXaxis()->SetTitleSize(0.05);
	Fht2D->GetYaxis()->SetTitleSize(0.05);
	Fht2D->GetXaxis()->SetLabelSize(0.05);
	Fht2D->GetYaxis()->SetLabelSize(0.05);
	Fht2D->GetXaxis()->SetTitle("Theta / Radian");
	Fht2D->GetYaxis()->SetTitle("Phi / Radian");
	Fht2D->Draw("colz");
	c1->Print(pdfPath);

	// c1->cd();
	c1->Print(pdfPath + "]");

	of.close();

	delete FhtDis;
	delete FhtPhi;
	delete Fht2D;
	delete Q2D;
	delete Q2Smooth;
	delete FhtDiff;
	delete pos;
	delete LiDiff;
	delete QDiff;
	delete c1;
	// delete testM;
	// outFile.close();
	// PhiOut.close();
	// label.close();
	// trks.close();
	LogDebug << "Executed" << endl;
	return true;
}

bool FhtAna::initGeomSvc() {
	SniperPtr<RecGeomSvc> rgSvc(getParent(), "RecGeomSvc");
	if (rgSvc.invalid()) {
		LogError << "Failed to get RecGeomSvc instance" << std::endl;
		return false;
	}
	m_geom = rgSvc->getCdGeom();
	m_wpgeom = rgSvc->getWpGeom();
	return true;
}

bool FhtAna::initPmt() {
	LogDebug << "Initializing PMTs" << std::endl;
	totPmtNum = 0;
	totPmtNum = m_wpgeom->getPmtNum();
	if (!totPmtNum) {
		LogError << "Wrong PMT Number" << std::endl;
		return false;
	}
	LogDebug << "PMT Number Got" << std::endl;
	m_ptab.reserve(totPmtNum);
	m_ptab.resize(totPmtNum);
	for (unsigned int pid = 0; pid < totPmtNum; pid ++) {
		Identifier Id = Identifier(WpID::id(pid, 0));
		PmtGeom* pmt = m_wpgeom->getPmt(Id);
		if (!pmt) {
			LogError << "Wrong PMT ID" << std::endl;
			return false;
		}
		TVector3 pmtCenter = pmt->getCenter();
		m_ptab[pid].pos = pmtCenter;
		// if (WpID::is3inch(Id)) {
		// 	m_ptab[pid].res = m_3inchRes;
		// 	m_ptab[pid].type = _PMTINCH3;
		// }
		if (WpID::is20inch(Id)) {
			m_ptab[pid].res = m_20inchRes;
			m_ptab[pid].type = _PMTINCH20;
		}
		else {
			LogError << "Pmt[" << pid << "] is neither 3-inch or 20-inch" << std::endl;
			return false;
		}
		m_ptab[pid].q = -1;
		m_ptab[pid].fht = 99999;
		m_ptab[pid].used = false;
	}
	return true;
}

bool FhtAna::freshPmtData(TH2D* ht, TH2D* pht, TH2D *h2d, TH2D *q2d, TH2D* nPMT, double &theta, double &phi) {
	JM::EvtNavigator* nav = m_buf->curEvt();
	if (not nav) {
		LogError << "Cannot retrieve current navigator" << std::endl;
		return false;
	}
	JM::CalibHeader* calibheader = (JM::CalibHeader*)nav->getHeader("/Event/Calib");
	if (not calibheader) {
		LogError << "Cannot retrieve '/Event/Calib'" << std::endl;
		return false;
	}
	const std::list<JM::CalibPMTChannel*>& chhlist = calibheader->event()->calibPMTCol();
	std::list<JM::CalibPMTChannel*>::const_iterator chit = chhlist.begin();
	if (chit == chhlist.end()) {
		LogDebug << "chit == chhlist.end()" << std::endl;
		return false;
	}
	double earliest = 1000;
	while (chit != chhlist.end()) {
		JM::CalibPMTChannel* calib = *chit ++;
		Identifier id = Identifier(calib->pmtId());
		Identifier::value_type value = id.getValue();
		unsigned int pid;
		if (not ((value & 0xFF000000) >> 24 == 0x20)) {
			continue;
		}
		if (pid > totPmtNum) {
			LogError << "Data/Geometry Mis-Match : PmtId(" << pid << ") >= the number of PMTs." << std::endl;
			return false;
		}
		pid = WpID::module(id);
		m_ptab[pid].q = calib->nPE();
		m_ptab[pid].fht = calib->firstHitTime();
		if ((WpID::is20inch(id) && m_20inchusedflag)) {
			if (earliest > m_ptab[pid].fht) {
				earliest = m_ptab[pid].fht;
				theta = m_ptab[pid].pos.Theta();
				phi = m_ptab[pid].pos.Phi();
			}
			m_ptab[pid].loc = 1;
			m_ptab[pid].used = true;
			m_usedPmtNum ++;
			ht->Fill(m_ptab[pid].pos.Theta(), m_ptab[pid].fht);
			pht->Fill(m_ptab[pid].pos.Phi(), m_ptab[pid].fht);
			int binx = m_ptab[pid].pos.Theta() / (PI / 100);
			int biny = (m_ptab[pid].pos.Phi() + TMath::Pi()) / (PI / 100);
			if (m_ptab[pid].fht < 100)
				h2d->SetBinContent(binx, biny,
								   m_ptab[pid].fht < h2d->GetBinContent(binx, biny) || h2d->GetBinContent(binx, biny) == 0 ?
								   m_ptab[pid].fht : h2d->GetBinContent(binx, biny));
			// q2d->SetBinContent(binx, biny, q2d->GetBinContent(binx, biny) + m_ptab[pid].q);
			q2d->Fill(m_ptab[pid].pos.Theta(), m_ptab[pid].pos.Phi(), m_ptab[pid].q);
			// nPMT->SetBinContent(binx, biny, nPMT->GetBinContent(binx, biny) + 1);
			nPMT->Fill(m_ptab[pid].pos.Theta(), m_ptab[pid].pos.Phi(), 1);
			// LogDebug << m_ptab[pid].fht << std::endl;
			// LogDebug << h2d->GetBinContent(binx, biny) << std::endl;
		}
	}
	LogDebug << "Loading calibration data done" << std::endl;
	return true;
}

bool FhtAna::IfCrossCd(TVector3& Inci, TVector3& Dir, Double_t R) {
	TVector3 dir = Dir.Unit();
	Double_t Dis = TMath::Sqrt(Inci.Mag() * Inci.Mag() - fabs(Inci * dir) * fabs(Inci * dir));
	// LogDebug << "Distance to center: " << Dis << endl;
	if (R < Dis)
		return false;
	return true;
}

TVector3 FhtAna::InciOnLS(TVector3& Inci, TVector3& Dir, Double_t R) {
	TVector3 dir = Dir.Unit();
	Double_t Dis2 = Inci.Mag() * Inci.Mag() - fabs(Inci * dir) * fabs(Inci * dir);
	TVector3 LSInci = Inci + dir * (fabs(Inci * dir) - TMath::Sqrt(R * R - Dis2));
	return LSInci;
}

TVector3 FhtAna::PosOnLS(TVector3& pos, TVector3& Dir, Double_t R, int co) {
	TVector3 dir = Dir.Unit();
	Double_t Dis2 = pos.Mag2() - fabs(pos * dir) * fabs(pos * dir);
	TVector3 LSpos = pos + co * dir * (TMath::Sqrt(R * R - Dis2) + co * ((pos * dir > 0) ? -1 : 1) * fabs(pos * dir));
	return LSpos;
}

bool FhtAna::finalize() {
	LogDebug << "Finalizing" << std::endl;
	return true;
}

TVector3 FhtAna::GetInciPos(TH1D *ThevFht, TH1D *ThevPhi, int n) {
	int thetaBin;
	int fht = 1000;
	for (int i = 1; i < n - 3; i ++) {
		if (ThevFht->GetBinContent(i) < fht) {
			thetaBin = i;
			fht = ThevFht->GetBinContent(i);
		}
	}
	LogInfo << "Minimum Fht: " << ThevFht->GetBinContent(thetaBin) << endl;
	double phi = ThevPhi->GetBinContent(thetaBin);
	phi = phi * TMath::Pi() / 100 - TMath::Pi();
	double theta = (double)thetaBin * TMath::Pi() / n;
	TVector3 InciPos;
	InciPos.SetMagThetaPhi(17700, theta, phi);
	LogInfo << "Theta: " << theta << "\t"
			<< "Phi: " << phi << endl;
	return InciPos;
}

TVector3 FhtAna::GetExitPos(TH1D *ThevFht, TH1D *ThevPhi, int n) {
	int fht = 1000;
	int thetaBin;
	for (int i = 1; i < n - 3; i ++) {
		if (ThevFht->GetBinContent(i) < fht) {
			thetaBin = i;
			fht = ThevFht->GetBinContent(i);
		}
	}
	double theta = (double)thetaBin * TMath::Pi() / n;
	for (int i = 5; i <= n - 4; i ++) {
		double Vi = ThevFht->GetBinContent(i);
		bool breakFlag = false;
		for (int j = i - 4; j <= i + 4; j ++) {
			if (j == i)
				continue;
			if (ThevFht->GetBinContent(j) < Vi) {
				breakFlag = true;
				break;
			}
		}
		if (!breakFlag)
			if(TMath::Abs(i - thetaBin) < 10)
				continue;
			else {
				thetaBin = i;
				break;
			}
	}
	if (thetaBin == ThevFht->GetMinimumBin()) {
		return TVector3(0, 0, 0);
	}

	theta = (double)thetaBin * TMath::Pi() / n;
		
	double phi = ThevPhi->GetBinContent(thetaBin);
	phi = phi * TMath::Pi() / 100 - TMath::Pi();
	TVector3 ExitPos;
	ExitPos.SetMagThetaPhi(19400, theta, phi);
	LogInfo << "Theta: " << theta << "\t"
			<< "Phi: " << phi << endl;
	return ExitPos;
}

TVector3 FhtAna::GetChargeCenter() {
	int n = m_ptab.size();
	double totCharge = 0;
	TVector3 totChaPos(0, 0, 0);
	for (int i = 0; i < n; i ++) {
		if (!m_ptab[i].used)
			continue;
		totCharge += m_ptab[i].q;
		totChaPos += m_ptab[i].q * m_ptab[i].pos;
	}
	totChaPos *= 1 / totCharge;
	return totChaPos;
}

TH2D* FhtAna::MapSmooth(TH2D* ori, int nx, int ny, TString name) {
	TH2D* h = (TH2D*)ori->Clone("tmph");

	// Extend edge of the map
	double unit = PI / nx;
	TH2D* ret = new TH2D(name, "", nx + 20, 0 - 10 * unit, PI + 10 * unit, ny + 20, -PI - 10 * unit, PI + 10 * unit);
	MapExtend(ret, h, nx, ny);
	delete h;

	// Smooth process
	TH2D* bc = (TH2D*)ret->Clone("BackH");
	for (int i = 3; i <= nx + 18; i ++) {
		for (int j = 3; j <= ny + 18; j ++) {
			double tmp = 0;
			for (int k = i - 2; k <= i + 2; k ++)
				for (int l = j - 2; l <= j + 2; l ++)
					tmp += bc->GetBinContent(k, l);
			tmp /= 25;
			ret->SetBinContent(i, j, tmp);
		}
	}

	delete bc;
	return ret;
}

bool FhtAna::FillContent(TH2D* h, int nx, int ny) {
	for (int i = 0; i < 2; i ++);
		Expansion(h, nx, ny);
	return true;
}

TH2D* FhtAna::PECut(TH2D* ori, int nx, int ny, double thr) {
	// LogDebug << "ori address: " << ori << endl;
	TString name("Q2D");
	name += thr;
	name += ori->GetName();
	TH2D* ret = new TH2D(name, "", nx, - 0.1 * PI, 1.1 * PI, ny, - 1.1 * PI, 1.1 * PI);
	double peak = ori->GetBinContent(ori->GetMaximumBin());
	// if (peak > 12000)
	// 	thr *= peak;
	// else
	// 	thr = (thr + 0.1) * peak;
	thr = thr * peak;
	LogDebug << "Threshold: " << thr << endl;
	// if (thr < 1000) {
	// 	LogDebug << "No track in CD" << endl;
	// 	return NULL;
	// }
	for (int i = 1; i <= nx; i ++)
		for (int j = 1; j <= ny; j ++) {
			double tmp = ori->GetBinContent(i, j);
			// ori->SetBinContent(i, j, tmp > thr ? tmp : 0);
			ret->SetBinContent(i, j, tmp > thr ? tmp : 0);
		}
	return ret;
}

void FhtAna::nCorrosion(TH2D* ori, int nx, int ny, int nturn) {
	for (int i = 0; i < nturn; i ++)
		Corrosion(ori, nx, ny);
}

void FhtAna::Corrosion(TH2D* ori, int nx, int ny) {
	TH2D* bc = (TH2D*)ori->Clone("oriBack");
	for (int i = 2; i <= nx - 8; i ++) {
		for (int j = 2; j <= ny - 8; j ++) {
			int matchFlag = 0;
			for (int k = i - 1; k <= i + 1; k ++) {
				for (int l = j - 1; l <= j + 1; l ++) {
					if (bc->GetBinContent(k, l))
						matchFlag ++;
				}
			}
			if (matchFlag < 6)
				ori->SetBinContent(i, j, 0);
		}
	}
	delete bc;
}

int FhtAna::AreaCut(TH2D* ori, TH2D* mark, int nx, int ny, double thr, bool cutOut, bool cutIn) {
	TH2D* bc = (TH2D*)ori->Clone("oriBack");
	struct Area {
		int area = 0;
		double aIn = 0;
		double aOut = 0;
		double max = 0;
		double min = 1E9;
		TVector2 st;
		TVector2 ed;
	};
	map<int, struct Area> mArea;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmp = mark->GetBinContent(i, j);
			if (tmp) {
				double val = ori->GetBinContent(i, j);
				mArea[(int)tmp].area ++;
				if (i >= 11 && i <= nx - 10 && j >= 11 && j <= ny - 10)
					mArea[(int)tmp].aIn ++;
				else
					mArea[(int)tmp].aOut ++;
				if (mArea[(int)tmp].area == 1) {
					mArea[(int)tmp].st.Set((double)i, (double)j);
					mArea[(int)tmp].ed.Set((double)i, (double)j);
					mArea[(int)tmp].max = val;
				}
				else {
					mArea[(int)tmp].st.Set(i < mArea[(int)tmp].st.X() ? i : mArea[(int)tmp].st.X(), j < mArea[(int)tmp].st.Y() ? j : mArea[(int)tmp].st.Y());
					mArea[(int)tmp].ed.Set(i > mArea[(int)tmp].ed.X() ? i : mArea[(int)tmp].ed.X(), j > mArea[(int)tmp].ed.Y() ? j : mArea[(int)tmp].ed.Y());
					mArea[(int)tmp].max = val > mArea[(int)tmp].max ? val : mArea[(int)tmp].max;
					mArea[(int)tmp].min = val < mArea[(int)tmp].min ? val : mArea[(int)tmp].min;
				}
			}
		}
	}

	map<int, Area>::iterator it = mArea.begin();
	bool overArea = false;
	LogInfo << mArea.size() << endl;
	while (it != mArea.end() && cutIn) {
		if ((it->second).area > 200) {
			overArea = true;
			double th = thr * ((it->second).max - (it->second).min) + (it->second).min;
			// if ((it->second).max < 1E6)
			// 	th = 1.1E6;
			LogInfo << "Threshold: " << th << endl;
			for (int i = ((it->second).st.X()); i <= ((it->second).ed.X()); i ++) {
				for (int j = (it->second).st.Y(); j <= (it->second).ed.Y(); j ++) {
					double tmp = ori->GetBinContent(i, j);
					int ma = (int)mark->GetBinContent(i, j);
					if (tmp < th && ma == it->first) {
						ori->SetBinContent(i, j, 0);
						mark->SetBinContent(i, j, 0);
					}
				}
			}
		}
		it ++;
	}

	if (overArea) {
		for (int i = 1; i <= nx; i ++)
			for (int j = 1; j <= ny; j ++)
				mark->SetBinContent(i, j, 0);
		MarkConnection(ori, nx, ny, mark, 10);
	}

	if (!cutOut)
		return mArea.size();

	map<int, Area> areas;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmp = mark->GetBinContent(i, j);
			if (tmp) {
				double val = ori->GetBinContent(i, j);
				areas[(int)tmp].area ++;
				if (i >= 11 && i <= nx - 10 && j >= 11 && j <= ny - 10)
					areas[(int)tmp].aIn ++;
				else
					areas[(int)tmp].aOut ++;
				if (areas[(int)tmp].area == 1) {
					areas[(int)tmp].st.Set((double)i, (double)j);
					areas[(int)tmp].ed.Set((double)i, (double)j);
					areas[(int)tmp].max = val;
				}
				else {
					areas[(int)tmp].st.Set(i < areas[(int)tmp].st.X() ? i : areas[(int)tmp].st.X(), j < areas[(int)tmp].st.Y() ? j : areas[(int)tmp].st.Y());
					areas[(int)tmp].ed.Set(i > areas[(int)tmp].ed.X() ? i : areas[(int)tmp].ed.X(), j > areas[(int)tmp].ed.Y() ? j : areas[(int)tmp].ed.Y());
					areas[(int)tmp].max = val > areas[(int)tmp].max ? val : areas[(int)tmp].max;
					areas[(int)tmp].min = val < areas[(int)tmp].min ? val : areas[(int)tmp].min;
				}
			}
		}
	}

	overArea = false;
	it = areas.begin();
	while (it != areas.end()) {
		LogInfo << "AreaID: " << it->first << endl;
		LogInfo << "AreaIn: " << (it->second).aIn << endl;
		LogInfo << "AreaOut: " << (it->second).aOut << endl;
		if ((it->second).aIn < (it->second).aOut) {
			overArea = true;
			for (int i = ((it->second).st.X()); i <= ((it->second).ed.X()); i ++) {
				for (int j = (it->second).st.Y(); j <= (it->second).ed.Y(); j ++) {
					int ma = (int)mark->GetBinContent(i, j);
					if (ma == it->first) {
						ori->SetBinContent(i, j, 0);
						mark->SetBinContent(i, j, 0);
					}
				}
			}
		}
		it ++;
	}
	if (overArea) {
		for (int i = 1; i <= nx; i ++)
			for (int j = 1; j <= ny; j ++)
				mark->SetBinContent(i, j, 0);
		return MarkConnection(ori, nx, ny, mark, 10);
	}
	return mArea.size();
}

int* FhtAna::GetMassPos(TH2D* ori, TH2D* mark, int nx, int ny, TH2D* test) {
	struct PosPE
	{
		int pos;
		double PE;
	};

	vector<PosPE> ct;
	
	// for (int i = 11; i <= nx - 10; i ++) {
	// 	for (int j = 11; j <= ny - 10; j ++) {
	// 		double tmp = 0;
	// 		for (int k = i - 5; k <= i + 5; k ++)
	// 			for (int l = j - 5; l <= j + 5; l ++)
	// 				tmp += ori->GetBinContent(k, l);
	// 		tmp /= 121;
	// 		if (tmp > 7000) {
	// 			PosPE pp;
	// 			pp.pos = (i - 10) * 1000 + (j - 10);
	// 			pp.PE = tmp;
	// 			ct.push_back(pp);
	// 		}
	// 	}
	// }

	// double x[225];
	map<int, TVector3> qp;
	map<int, double> q;
	map<int, int> area;
	double unit = TMath::Pi() / 100;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmp = mark->GetBinContent(i, j);
			if (tmp > 0) {
				// LogDebug << "tmp: " << tmp << endl;
				TVector3 p;
				double phi;
				double the;
				if (j < 11) {
					phi = (190 + j) * unit - TMath::Pi();
					if (i < 11) {
						the = (11 - i) * unit;
						phi = phi > 0 ? phi - TMath::Pi() : phi == 0 ? 0 : phi + TMath::Pi();
					}
					else if (i > 110) {
						the = (210 - i) * unit;
						phi = phi > 0 ? phi - TMath::Pi() : phi == 0 ? 0 : phi + TMath::Pi();
					}
					else {
						the = (i - 10) * unit;
					}
				}
				else if (j > 210) {
					phi = (j - 210) * unit - TMath::Pi();
					if (i < 11) {
						the = (11 - i) * unit;
						phi = phi > 0 ? phi - TMath::Pi() : phi == 0 ? 0 : phi + TMath::Pi();
					}
					else if (i > 110) {
						the = (210 - i) * unit;
						phi = phi > 0 ? phi - TMath::Pi() : phi == 0 ? 0 : phi + TMath::Pi();
					}
					else {
						the = (i - 10) * unit;
					}
				}
				else {
					the = (i - 10) * unit;
					phi = (j - 10) * unit - TMath::Pi();
				}
				p.SetMagThetaPhi(m_LSRadius, the, phi);
				// LogDebug << the << ", " << phi << ", " << p << endl;
				qp[(int)tmp] += (p * ori->GetBinContent(i, j));
				// LogDebug << qp[(int)tmp] << endl;
				q[(int)tmp] += ori->GetBinContent(i, j);
				area[(int)tmp] ++;
			}
		}
	}

	static int mass[4] = {0};
	map<int, TVector3>::iterator qpIt = qp.begin();
	map<int, double>::iterator qIt = q.begin();
	map<int, int>::iterator aIt = area.begin();
	map<int, TVector3> rec;
	int m = 0;
	while (qpIt != qp.end()) {
		LogDebug << "nIterator: " << qpIt->first << endl;
		TVector3 p = 1 / qIt->second * qpIt->second;
		LogDebug << p << endl;
		if (m < 4) {
			map<int, TVector3>::iterator recIt = rec.begin();
			int i = 0;
			bool breakFlag = false;
			while (recIt != rec.end()) {
				if ((p - recIt->second).Mag() < 3000 &&
					(p.Theta() < 0.314 && (recIt->second).Theta() < 0.314 ||
					 p.Theta() > 2.826 && (recIt->second).Theta() > 2.826 ||
					 p.Phi() < -2.826 && (recIt->second).Phi() < -2.826 ||
					 p.Phi() > 2.826 && (recIt->second).Phi() > 2.826)) {
					if (area[recIt->first] < area[qpIt->first]) {
						mass[i] = (int)(p.Theta() / unit) * 1000 + (int)((p.Phi() + TMath::Pi()) / unit);
						rec.erase(recIt);
					}
					breakFlag = true;
					break;
				}
				i ++;
				recIt ++;
			}
			if (!breakFlag) {
				int pp = (int)(p.Theta() / unit) * 1000 + (int)((p.Phi() + TMath::Pi()) / unit);
				mass[m] = pp;
				m ++;
			}
		}
		rec[qpIt->first] = p;
		qpIt ++;
		qIt ++;
	}

	for (int i = 0; i < 4; i ++)
		LogDebug << "mass[" << i << "]: " << mass[i] << endl;

	qp.clear();
	q.clear();
	area.clear();
	rec.clear();

	return mass;
}

int FhtAna::MarkConnection(TH2D* ori, int nx, int ny, TH2D* mark, int thr) {
	vector<int> st;
	vector<int> ed;
	double ID = 0;
	map<int, int> parent;
	for (int i = 1; i <= nx; i ++) {
		double last = 0;
		vector<int> start;
		vector<int> end;
		for (int j = 1; j <= ny; j ++) {
			if (ori->GetBinContent(i, j)) {
				if (last == 0) {
					ID ++;
					start.push_back(j);
					parent[ID] = ID;
				}
				if (j == ny) {
					end.push_back(j);
					// LogDebug << "st.size: " << start.size() << " ed.size: " << end.size() << endl;
					int k = 0;
					for (int& u : ed) {
						if (start.back() <= u && end.back() >= st[k]) {
							double tmpID = mark->GetBinContent(i - 1, u);
							while (tmpID != parent[tmpID])
								tmpID = parent[tmpID];
							if (parent[ID] < parent[tmpID])
								parent[tmpID] = parent[ID];
							else
								parent[ID] = parent[tmpID];
							// LogDebug << "InIf: " << tmpID << endl;
						}
						k ++;
					}
				}
				mark->SetBinContent(i, j, ID);
			}
			else {
				if (last) {
					end.push_back(j - 1);
					int k = 0;
					for (int& u : ed) {
						if (start.back() <= u && end.back() >= st[k]) {
							double tmpID = mark->GetBinContent(i - 1, u);
							while (tmpID != parent[tmpID])
								tmpID = parent[tmpID];
							double bcID = ID;
							while (bcID != parent[bcID])
								bcID = parent[bcID];
							if (parent[bcID] < parent[tmpID])
								parent[tmpID] = parent[bcID];
							else
								parent[bcID] = parent[tmpID];
							// LogDebug << "InIf: " << tmpID << endl;
						}
						k ++;
					}
					// LogDebug << "PosID: " << i << ", " << j << ", " << parent[ID] << ", " << k << endl;
				}
			}
			last = ori->GetBinContent(i, j);
		}
		st.swap(start);
		ed.swap(end);
	}

	int ret = 0;
	map<int, int> count;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmp = mark->GetBinContent(i, j);
			if (tmp) {
				while (tmp != parent[tmp])
					tmp = parent[tmp];
				mark->SetBinContent(i, j, parent[tmp]);
				count[parent[tmp]] ++;
			}
		}
	}

	map<int, int>::iterator it = count.begin();
	ID = 1;
	while (it != count.end()) {
		if (it->second < thr) {
			it->second = 0;
		}
		else {
			it->second = (int)ID;
			ID += 1.0;
		}
		it ++;
	}

	ret = (int)ID;

	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmp = mark->GetBinContent(i, j);
			if (tmp) {
				mark->SetBinContent(i, j, count[parent[tmp]]);
				if (!count[parent[tmp]])
					ori->SetBinContent(i, j, 0);
			}
		}
	}

	// ret = AreaCut(ori, mark, nx, ny, 0.5);

	return ret;
}

bool FhtAna::FindTrk(TVector3& inci, TVector3& dir, double& dis, double& ang, double& ti, TH2D* tMap, long int* mass) {
	struct posFht {
		double theta;
		double phi;
		double fht;
		double mag;
		double z;
	};
	vector<posFht> points;
	double unit = PI / 100;
	int nMass = 0;
	for (int i  = 0; i < 4; i ++) {
		if (mass[i])
			nMass ++;
		double mag = (int)(mass[i] / 1000000);
		mass[i] = (int)(mass[i] % 1000000);
		double the = (int)(mass[i] / 1000);
		double phi = (int)(mass[i] % 1000);
		double fht = tMap->GetBinContent(the, phi);
		the *= unit;
		phi = phi * unit - PI;
		LogInfo << "mag: " << mag << "\ttheta: " << the << "\tphi: " << phi << "\tfht: " << fht << endl;
		posFht pf;
		pf.phi = phi;
		pf.theta = the;
		pf.fht = mass[i] == 0 ? 10000 : fht;
		pf.mag = m_LSRadius;
		TVector3 tmp;
		tmp.SetMagThetaPhi(mag, the, phi);
		pf.z = mass[i] == 0 ? INT_MIN : tmp.Z();
		points.push_back(pf);
		mass[i] = 0;
	}
	sort(points.begin(), points.end(), [](posFht p1, posFht p2) {
		return p1.z > p2.z;
	});
	vector<double> angs;
	TVector3 p1, p2, p3, p4;
	int ID = 0;
	p1.SetMagThetaPhi(points[0].mag, points[0].theta, points[0].phi);
	p2.SetMagThetaPhi(points[1].mag, points[1].theta, points[1].phi);
	p3.SetMagThetaPhi(points[2].mag, points[2].theta, points[2].phi);
	p4.SetMagThetaPhi(points[3].mag, points[3].theta, points[3].phi);

	if (nMass == 0) {
		LogInfo << "No Track" << endl;
		return true;
	}	
	else if (nMass == 1) {
		TVector3 tmp = GetChargeCenter();
		dir = (tmp.Z() < p1.Z()) ? (tmp - p1).Unit() : (p1 - tmp).Unit();
		// inci = (tmp.Z() < p1.Z()) ? PosOnLS(p1, dir, m_LSRadius, -1) : PosOnLS(tmp, dir, m_LSRadius, -1);
		inci = PosOnLS(p1, dir, m_LSRadius, -1);
		dis = 0;
		ang = 0;
		ti = tMap->GetBinContent(inci.Theta() / unit, (inci.Phi() + PI) / unit);
		return true;
	}
	else if (nMass == 2) {
		dir = (p1.Z() < p2.Z()) ? (p1 - p2).Unit() : (p2 - p1).Unit();
		// inci = p1.Z() < p2.Z() ? PosOnLS(p2, dir, m_LSRadius, -1) : PosOnLS(p1, dir, m_LSRadius, -1);
		inci = PosOnLS(p1, dir, m_LSRadius, -1);
		TVector3 tmp = GetChargeCenter();
		tmp = tmp - (inci + dir * (tmp - inci) * inci);
		dis = tmp.Mag() * 2;
		TVector3 ori(0, dir.Z(), -dir.Y());
		ang = tmp.Angle(ori);
		ori.Rotate(ang, dir);
		if (tmp.Angle(ori) > 0.2)
			ang = 2 * PI - ang;
		ti = tMap->GetBinContent(inci.Theta() / unit, (inci.Phi() + PI) / unit);
		return true;
	}
	else if (nMass == 3) {
		TVector3 d1, d2, d3;
		d1 = p1 - p2;
		d2 = p2 - p3;
		d3 = p3 - p1;
		double a1, a2, a3;
		a1 = fabs(d1.Angle(TVector3(0, 0, -1)));
		a1 = (a1 > PI / 2) ? PI - a1 : a1;
		a2 = fabs(d2.Angle(TVector3(0, 0, -1)));
		a2 = (a2 > PI / 2) ? PI - a2 : a2;
		a3 = fabs(d3.Angle(TVector3(0, 0, -1)));
		a3 = (a3 > PI / 2) ? PI - a3 : a3;
		int i;
		TVector3 tmp;
		if (a1 < a2 && a1 < a3) {
			dir = (d1 * TVector3(0, 0, -1) > 0) ? d1 : - d1;
			dir = dir.Unit();
			inci = PosOnLS(p1, dir, m_LSRadius, -1);
			tmp = p3 - (inci + dir * (p3 - inci) * dir);
		}
		else if (a2 < a1 && a2 < a3) {
			dir = (d2 * TVector3(0, 0, -1) > 0) ? d2 : - d2;
			dir = dir.Unit();
			inci = PosOnLS(p3, dir, m_LSRadius, -1);
			tmp = p1 - (inci + dir * (p1 - inci) * dir);
		}
		else if (a3 < a1 && a3 < a2) {
			dir = (d3 * TVector3(0, 0, -1) > 0) ? d3 : - d3;
			dir = dir.Unit();
			inci = PosOnLS(p1, dir, m_LSRadius, -1);
			tmp = p2 - (inci + dir * (p2 - inci) * dir);
		}
		dis = tmp.Mag();
		ti = tMap->GetBinContent(inci.Theta() / unit, (inci.Phi() + PI) / unit);
		TVector3 ori(0, dir.Z(), - dir.Y());
		ang = tmp.Angle(ori);
		ori.Rotate(ang, dir);
		if (tmp.Angle(ori) > 0.2)
			ang = 2 * PI - ang;
		return true;
	}
	else if (nMass == 4) {
		double a1, a2, a3;
		TVector3 tmp;
		a1 = fabs((p1 - p2).Angle(p3 - p4));
		a1 = (a1 > PI / 2) ? PI - a1 : a1;
		a2 = fabs((p1 - p3).Angle(p2 - p4));
		a2 = (a2 > PI / 2) ? PI - a2 : a2;
		a3 = fabs((p1 - p4).Angle(p2 - p3));
		a3 = (a3 > PI / 2) ? PI - a3 : a3;
		if (a1 < a2 && a1 < a3) {
			dir = p1 - p2;
			dir = (dir * TVector3(0, 0, -1) > 0) ? dir : - dir;
			dir = dir.Unit();
			inci = PosOnLS(p1, dir, m_LSRadius, -1);
			tmp = p3 - (inci + dir * (p3 - inci) * dir);
		}
		else if (a2 < a1 && a2 < a3) {
			dir = p1 - p3;
			dir = (dir * TVector3(0, 0, -1) > 0) ? dir : - dir;
			dir = dir.Unit();
			inci = PosOnLS(p1, dir, m_LSRadius, -1);
			tmp = p2 - (inci + dir * (p3 - inci) * dir);
		}
		else if (a3 < a1 && a3 < a2) {
			dir = p1 - p4;
			dir = (dir * TVector3(0, 0, -1) > 0) ? dir : - dir;
			dir = dir.Unit();
			inci = PosOnLS(p1, dir, m_LSRadius, -1);
			tmp = p2 - (inci + dir * (p2 - inci) * dir);
		}
		dis = tmp.Mag();
		ti = tMap->GetBinContent(inci.Theta() / unit, (inci.Phi() + PI) / unit);
		TVector3 ori(0, dir.Z(), - dir.Y());
		ang = tmp.Angle(ori);
		ori.Rotate(ang, dir);
		if (tmp.Angle(ori) > 0.2)
			ang = 2 * PI - ang;
		return true;
	}
	else {
		LogInfo << "Find trk fail" << endl;
		return false;
	}
}

bool FhtAna::ChooseCut(TH2D* ori, TH1D* q, int nx, int ny) {
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			q->Fill(ori->GetBinContent(i, j));
		}
	}
	return true;
}

bool FhtAna::Expansion(TH2D* ori, int nx, int ny) {
	if (ori == NULL) {
		LogInfo << "The map is NULL" << endl;
		return false;
	}
	// LogInfo << "------------------------------ Expanding ------------------------------" << endl;
	TH2D* cp = (TH2D*)ori->Clone("tmpCp");
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			if (cp->GetBinContent(i, j))
				continue;
			double sum = 0;
			double max = 0;
			int n = 0;
			if (i > 1 && cp->GetBinContent(i - 1, j)) {
				sum += cp->GetBinContent(i - 1, j);
				max = cp->GetBinContent(i - 1, j) > max ? cp->GetBinContent(i - 1, j) : max;
				n ++;
			}
			if (i < nx && cp->GetBinContent(i + 1, j)) {
				sum += cp->GetBinContent(i + 1, j);
				max = cp->GetBinContent(i + 1, j) > max ? cp->GetBinContent(i + 1, j) : max;
				n ++;
			}
			if (j > 1 && cp->GetBinContent(i, j - 1)) {
				sum += cp->GetBinContent(i, j - 1);
				max = cp->GetBinContent(i, j - 1) > max ? cp->GetBinContent(i, j - 1) : max;
				n ++;
			}
			if (j < ny && cp->GetBinContent(i, j + 1)) {
				sum += cp->GetBinContent(i, j + 1);
				max = cp->GetBinContent(i, j + 1) > max ? cp->GetBinContent(i, j + 1) : max;
				n ++;
			}
			if (i > 1 && j > 1 && cp->GetBinContent(i - 1, j - 1)) {
				sum += cp->GetBinContent(i - 1, j - 1);
				max = cp->GetBinContent(i - 1, j - 1) > max ? cp->GetBinContent(i - 1, j - 1) : max;
				n ++;
			}
			if (i > 1 && j < ny && cp->GetBinContent(i - 1, j + 1)) {
				sum += cp->GetBinContent(i - 1, j + 1);
				max = cp->GetBinContent(i - 1, j + 1) > max ? cp->GetBinContent(i - 1, j + 1) : max;
				n ++;
			}
			if (i < nx && j > 1 && cp->GetBinContent(i + 1, j - 1)) {
				sum += cp->GetBinContent(i + 1, j - 1);
				max = cp->GetBinContent(i + 1, j - 1) > max ? cp->GetBinContent(i + 1, j - 1) : max;
				n ++;
			}
			if (i < nx && j < ny && cp->GetBinContent(i + 1, j + 1)) {
				sum += cp->GetBinContent(i + 1, j + 1);
				max = cp->GetBinContent(i + 1, j + 1) > max ? cp->GetBinContent(i + 1, j + 1) : max;
				n ++;
			}
			if (n) {
				// ori->SetBinContent(i, j, max);
				ori->SetBinContent(i, j, sum / n);
			}
		}
	}
	return true;
}
bool FhtAna::RMSMap(TH2D* ori, TH2D* rms, int nx, int ny, int u, int len, double thr) {
	if (ori == NULL) {
		LogInfo << "The map is NULL" << endl;
		return false;
	}
	/*
	for (int i = 11; i <= nx - 10; i ++) {
		for (int j = 11; j <= ny - 10; j ++) {
			double sum = ori->Integral(i - 5, i + 5, j - 5, j + 5);
			double ave = sum / 121;
			sum = 0;
			for (int k = i - 5; k <= i + 5; k ++)
				for (int l = j - 5; l <= j + 5; l ++)
					sum += TMath::Power(ave - ori->GetBinContent(k, l), 2);
			sum /= 121;
			rms->SetBinContent(i - 10, j - 10, TMath::Sqrt(sum));
		}
	}
	*/
	for (int i = u + 1; i <= nx - u; i ++) {
		for (int j = u + 1; j <= ny - u; j ++) {
			double sum = 0;
			TVector2 ave;
			double p = 0;
			double sig = 0;
			double tmp = ori->GetBinContent(i, j);
			// sum += ori->GetBinContent(i, j) - ori->GetBinContent(i - 5, j);
			// sum += ori->GetBinContent(i, j) - ori->GetBinContent(i + 5, j);
			// sum += ori->GetBinContent(i, j) - ori->GetBinContent(i, j - 5);
			// sum += ori->GetBinContent(i, j) - ori->GetBinContent(i, j + 5);
/*
			for (int k = i - len; k <= i + len; k ++) {
				for (int l = j - len; l <= j + len; l ++) {
					double coe = len - TMath::Abs(i - k) + 1;
					if (coe > (len - TMath::Abs(j - l) + 1))
						coe = len - TMath::Abs(j - l) + 1;
					sum += ori->GetBinContent(k, l) * coe;
					coes += coe;
				}
			}
			sum /= coes;
*/
			for (int k = i - len; k <= i + len; k ++) {
				for (int l = j - len; l <= j + len; l ++) {
					ave += ori->GetBinContent(k, l) * TVector2(k, l);
					sum += ori->GetBinContent(k, l);
				}
			}
			ave /= sum;
			for (int k = i - len; k <= i + len; k ++) {
				for (int l = j - len; l <= i + len; l ++) {
					TVector2 tmpv(k, l);
					p += ori->GetBinContent(k, l) * TMath::Power((tmpv - ave).Mod(), 4);
					sig += ori->GetBinContent(k, l) * TMath::Power((tmpv - ave).Mod(), 2);
				}
			}
			p /= sum;
			sig = TMath::Power(sig / sum, 2);

			rms->SetBinContent(i - u, j - u, sum);
			// rms->SetBinContent(i - u, j - u, sum > thr ? sum : 0);
			// if (ori->GetBinContent(i, j) == 200000)
			// 	rms->SetBinContent(i - 10, j - 10, 0);
		}
	}
	return true;
}

bool FhtAna::MapExtend(TH2D* ret, TH2D* h, int nx, int ny) {
	if (h == NULL || ret == NULL) {
		LogInfo << "The map is NULL" << endl;
		return false;
	}
	double unit = PI / nx;
	for (int i = 1; i <= nx; i ++)
		for (int j = 1; j <= ny; j ++)
			ret->SetBinContent(i + 10, j + 10, h->GetBinContent(i, j));
	for (int i = 1; i <= ny; i ++) {
		int pos = (i == 100 ? 1 : (i < 100 ? i + 100 : i - 100));
		for (int j = 1; j <= 10; j ++)
			ret->SetBinContent(j, i + 10, h->GetBinContent(11 - j, pos));
		for (int j = nx + 11; j <= nx + 20; j ++)
			ret->SetBinContent(j, i + 10, h->GetBinContent(nx - (j - nx - 11), pos));
	}
	for (int i = 1; i <= nx + 20; i ++) {
		for (int j = 1; j <= 10; j ++) {
			ret->SetBinContent(i, j, ret->GetBinContent(i, j + ny));
			ret->SetBinContent(i, ny + 10 + j, ret->GetBinContent(i, j + 10));
		}
	}
	return true;
}

bool FhtAna::Pool(TH2D* ori, TH2D* pool, int nx, int ny, int s) {
	if (ori == NULL) {
		LogInfo << "The map is NULL" << endl;
		return false;
	}
	nx /= s;
	ny /= s;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double sum = 0;
			for (int k = (i - 1) * s + 1; k <= i * s; k ++) {
				for (int l = (j - 1) * s + 1; l <= j * s; l ++) {
					sum += ori->GetBinContent(k, l);
				}
			}
			pool->SetBinContent(i, j, sum);
		}
	}
	return true;
}

bool FhtAna::XOR(TH2D* a, TH2D* b, int nx, int ny) {
	// a is the map of low threshold, b is the map of high threshold
	if (!a || !b) {
		LogInfo << "Input map is NULL" << endl;
		return false;
	}
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmpa = a->GetBinContent(i, j);
			double tmpb = b->GetBinContent(i, j);
			a->SetBinContent(i, j, ((!tmpa && tmpb) || (tmpa && !tmpb)) ? (tmpa ? tmpa : tmpb) : 0);
		}
	}
	return a;
}

TH2D* FhtAna::Combine(TH2D* a, TH2D* b, int nx, int ny) {
	// Combine the map a & b, if the connection area partially overlap, perform AND, or perform OR
	if (!a || !b) {
		LogInfo << "Input map is NULL" << endl;
		return NULL;
	}

	// map<int, int> marks;
	// for (int i = 1; i <= nx; i ++) {
	// 	for (int j = 1; j <= ny; j ++) {
	// 		double tmpa = a->GetBinContent(i, j);
	// 		double tmpb = b->GetBinContent(i, j);
	// 		if (tmpa && tmpb) {
	// 			marks[tmpa] = 1;
	// 		}
	// 	}
	// }

	// TString name("Combine");
	// name += a->GetName();
	// name += b->GetName();
	// TH2D* ret = new TH2D(name, "", nx, 0, 1.1 * PI, ny, -1.1 * PI, 1.1 * PI);
	// for (int i = 1; i <= nx; i ++) {
	// 	for (int j = 1; j <= ny; j ++) {
	// 		double tmpa = a->GetBinContent(i, j);
	// 		double tmpb = b->GetBinContent(i, j);
	// 		if (marks[tmpa])
	// 			ret->SetBinContent(i, j, tmpa && tmpb ? tmpa : 0);
	// 		else
	// 			ret->SetBinContent(i, j, tmpa || tmpb ? (tmpb > 0 ? tmpb : tmpa) : 0);
	// 	}
	// }
	// return ret;

	TString name("Combine");
	name += a->GetName();
	name += b->GetName();
	TH2D* ret = new TH2D(name, "", nx, 0, 1.1 * PI, ny, -1.1 * PI, 1.1 * PI);
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmpa = a->GetBinContent(i, j);
			double tmpb = b->GetBinContent(i, j);
			a->SetBinContent(i, j, tmpa + tmpb);
		}
	}
	MarkConnection(a, nx, ny, ret, 5);
	return ret;
}

bool FhtAna::AND(TH2D* ori, TH2D* co, int nx, int ny) {
	if (!ori || !co) {
		LogInfo << "Input map is NULL" << endl;
		return false;
	}
	TString name("Mask");
	name += ori->GetName();
	TH2D* ret = new TH2D(name, "", nx, 0, 1.1 * PI, ny, -1.1 * PI, 1.1 * PI);
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmpa = ori->GetBinContent(i, j);
			double tmpb = co->GetBinContent(i, j);
			ori->SetBinContent(i, j, (tmpa && tmpb) ? tmpa : 0);
		}
	}
	return true;
}

bool FhtAna::UnionCut(TH2D* l, TH2D* h, TH2D* ori, int nx, int ny, double thr, TH2D* test1) {
	if (!l || !h) {
		LogInfo << "Input map is NULL" << endl;
		return false;
	}
	TH2D* H = (TH2D*)h->Clone("CloneH");
	for (int i = 0; i < 13; i ++)
		Expansion(H, 120, 220);
	if (!Expansion(H, 120, 220)) {
		LogInfo << "Error in Expansion()" << endl;
		return false;
	}

	for (int i = 1; i <= nx; i ++)
		for (int j = 1; j <= ny; j ++)
			test1->SetBinContent(i, j, H->GetBinContent(i, j));

	struct Area {
		int area = 0;
		double aIn = 0;
		double aOut = 0;
		double max = 0;
		TVector2 st;
		TVector2 ed;
		double nOL = 0;
		double lastMark = 0;
		double hMax = 0;
	};
	LogInfo << "Checking..." << endl;
	map<int, struct Area> areas;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmpl = l->GetBinContent(i, j);
			if (tmpl) {
				double tmph = H->GetBinContent(i, j);
				double val = ori->GetBinContent(i, j);
				areas[(int)tmpl].area ++;
				if (i >= 11 && i <= nx - 10 && j >= 11 && j <= ny - 10)
					areas[(int)tmpl].aIn ++;
				else
					areas[(int)tmpl].aOut ++;
				if (areas[(int)tmpl].area == 1) {
					areas[(int)tmpl].st.Set((double)i, (double)j);
					areas[(int)tmpl].ed.Set((double)i, (double)j);
					areas[(int)tmpl].max = val;
				}
				else {
					TVector2 ST = areas[(int)tmpl].st;
					TVector2 ED = areas[(int)tmpl].ed;
					areas[(int)tmpl].st.Set(i < ST.X() ? i : ST.X(), j < ST.Y() ? j : ST.Y());
					areas[(int)tmpl].ed.Set(i > ED.X() ? i : ED.X(), j > ED.Y() ? j : ED.Y());
					if (!tmph)
						areas[(int)tmpl].max = val > areas[(int)tmpl].max ? val : areas[(int)tmpl].max;
					areas[(int)tmpl].hMax = val > areas[(int)tmpl].hMax ? val : areas[(int)tmpl].hMax;
				}
				if (tmph && tmph != areas[(int)tmpl].lastMark)
					areas[(int)tmpl].nOL ++;
				if (tmph)
					areas[(int)tmpl].lastMark = tmph;
			}
		}
	}

	LogInfo << "Processing..." << endl;
	map<int, Area>::iterator it = areas.begin();
	bool overArea = false;
	while (it != areas.end()) {
		// XOR & AreaCut
		LogInfo << "n overlap: " << (it->second).nOL << endl;
		if ((it->second).area > 200 && (it->second).nOL >= 2) {
			overArea = true;
			double th = thr * (it->second).max;
			if ((it->second).max < 0.7 * (it->second).hMax)
				th = (it->second).hMax * 0.7;
			LogInfo << "Threshold: " << th << endl;
			for (int i = (it->second).st.X(); i <= (it->second).ed.X(); i ++) {
				for (int j = (it->second).st.Y(); j <= (it->second).ed.Y(); j ++) {
					double tmpl = l->GetBinContent(i, j);
					double tmph = H->GetBinContent(i, j);
					if (tmpl && tmph) {
						l->SetBinContent(i, j, 0);
						ori->SetBinContent(i, j, 0);
					}
					double tmp = ori->GetBinContent(i, j);
					if (tmp && tmp < th && tmpl == it->first) {
						ori->SetBinContent(i, j, 0);
						l->SetBinContent(i, j, 0);
					}
				}
			}
		}
		if ((it->second).area > 300 && (it->second).nOL == 1) {
			overArea = true;
			double th = thr * (it->second).max;
			if ((it->second).max < 0.65 * (it->second).hMax)
				th = (it->second).hMax * 0.5;
			LogInfo << "Threshold: " << th << endl;
			for (int i = (it->second).st.X(); i <= (it->second).ed.X(); i ++) {
				for (int j = (it->second).st.Y(); j <= (it->second).ed.Y(); j ++) {
					double tmpl = l->GetBinContent(i, j);
					double tmph = H->GetBinContent(i, j);
					if (tmpl && tmph) {
						l->SetBinContent(i, j, 0);
						ori->SetBinContent(i, j, 0);
					}
					double tmp = ori->GetBinContent(i, j);
					if (tmp && tmp < th && tmpl == it->first) {
						ori->SetBinContent(i, j, 0);
						l->SetBinContent(i, j, 0);
					}
				}
			}
		}
		it ++;
	}

	if (overArea) {
		for (int i = 1; i <= nx; i ++)
			for (int j = 1; j <= ny; j ++)
				l->SetBinContent(i, j, 0);
		MarkConnection(ori, nx, ny, l, 10);
	}
	else
		return true;

	map<int, struct Area> Areas;
	for (int i = 1; i <= nx; i ++) {
		for (int j = 1; j <= ny; j ++) {
			double tmpl = l->GetBinContent(i, j);
			if (tmpl) {
				double tmph = H->GetBinContent(i, j);
				double val = ori->GetBinContent(i, j);
				Areas[(int)tmpl].area ++;
				if (i >= 11 && i <= nx - 10 && j >= 11 && j <= ny - 10)
					Areas[(int)tmpl].aIn ++;
				else
					Areas[(int)tmpl].aOut ++;
				if (Areas[(int)tmpl].area == 1) {
					Areas[(int)tmpl].st.Set((double)i, (double)j);
					Areas[(int)tmpl].ed.Set((double)i, (double)j);
					Areas[(int)tmpl].max = val;
				}
				else {
					TVector2 ST = Areas[(int)tmpl].st;
					TVector2 ED = Areas[(int)tmpl].ed;
					Areas[(int)tmpl].st.Set(i < ST.X() ? i : ST.X(), j < ST.Y() ? j : ST.Y());
					Areas[(int)tmpl].ed.Set(i > ED.X() ? i : ED.X(), j > ED.Y() ? j : ED.Y());
					if (!tmph)
						Areas[(int)tmpl].max = val > Areas[(int)tmpl].max ? val : Areas[(int)tmpl].max;
				}
				if (tmph != Areas[(int)tmpl].lastMark)
					Areas[(int)tmpl].nOL ++;
				Areas[(int)tmpl].lastMark = tmph;
			}
		}
	}

	overArea = false;
	it = Areas.begin();
	while (it != Areas.end()) {
		// Delete the area near the edge
		if ((it->second).aIn < (it->second).aOut) {
			overArea = true;
			LogInfo << "Cut outside..." << endl;
			for (int i = (it->second).st.X(); i <= (it->second).ed.X(); i ++) {
				for (int j = (it->second).st.Y(); j <= (it->second).ed.Y(); j ++) {
					double tmpl = l->GetBinContent(i, j);
					if (tmpl == it->first) {
						ori->SetBinContent(i, j, 0);
						l->SetBinContent(i, j, 0);
					}
				}
			}
		}
		it ++;
	}
	LogInfo << "Marking..." << endl;
	if (overArea) {
		for (int i = 1; i <= nx; i ++)
			for (int j = 1; j <= ny; j ++)
				l->SetBinContent(i, j, 0);
	}
	return true;
}

long int* FhtAna::GetCenterPos(TH2D* ori, TH2D* mark, int nx, int ny) {
	static long int mass[4] = {0};
	if (!ori || !mark) {
		LogInfo << "Input map is NULL" << endl;
		return mass;
	}
	map<int, TVector3> qp;
	map<int, double> q;
	map<int, int> area;
	double unit = PI / 100;
	for (int i = 0; i < totPmtNum; i ++) {
		if (!m_ptab[i].used)
			continue;
		double the = m_ptab[i].pos.Theta();
		double phi = m_ptab[i].pos.Phi();
		int x = the / unit + 11;
		int y = (phi + PI) / unit + 11;
		if (mark->GetBinContent(x, y)) {
			double tmp = mark->GetBinContent(x, y);
			qp[(int)tmp] += ori->GetBinContent(x, y) * m_ptab[i].pos;
			q[(int)tmp] += ori->GetBinContent(x, y);
			area[(int)tmp] ++;
		}
		if (x <= 20) {
			if ((y > 20 && y <= 100) || (y > 120 && y <= 200)) {
				y = y < 110 ? y + 100 : y - 100;
				x = 21 - x;
				double tmp = mark->GetBinContent(x, y);
				if (tmp) {
					qp[(int)tmp] += ori->GetBinContent(x, y) * m_ptab[i].pos;
					q[(int)tmp] += ori->GetBinContent(x, y);
					area[(int)tmp] ++;
				}
			}
			else {
				double y1 = y <= 110 ? y + 100 : y - 100;
				double x1 = 21 - x;
				double tmp = mark->GetBinContent(x1, y1);
				if (tmp) {
					qp[(int)tmp] += ori->GetBinContent(x1, y1) * m_ptab[i].pos;
					q[(int)tmp] += ori->GetBinContent(x1, y1);
					area[(int)tmp] ++;
				}
				if (y <= 20) {
					y1 = y + 200;
					x1 = x;
				}
				else if (y <= 110 && y > 100) {
					y1 -= 200;
				}
				else if (y <= 120 && y > 110) {
					y1 = y1 + 200;
				}
				else {
					y1 = y - 200;
					x1 = x;
				}
				tmp = mark->GetBinContent(x1, y1);
				if (tmp) {
					qp[(int)tmp] += ori->GetBinContent(x1, y1) * m_ptab[i].pos;
					q[(int)tmp] += ori->GetBinContent(x1, y1);
					area[(int)tmp] ++;
				}
			}
		}
		else if (x > 100) {
			if ((y > 20 && y <= 90) || (y > 110 && y <= 200)) {
				y = y < 110 ? y + 100 : y - 100;
				x = 221 - x;
				double tmp = mark->GetBinContent(x, y);
				if (tmp) {
					qp[(int)tmp] += ori->GetBinContent(x, y) * m_ptab[i].pos;
					q[(int)tmp] += ori->GetBinContent(x, y);
					area[(int)tmp] ++;
				}
			}
			else {
				double y1 = y <= 110 ? y + 100 : y - 100;
				double x1 = 221 - x;
				double tmp = mark->GetBinContent(x1, y1);
				if (tmp) {
					qp[(int)tmp] += ori->GetBinContent(x1, y1) * m_ptab[i].pos;
					q[(int)tmp] += ori->GetBinContent(x1, y1);
					area[(int)tmp] ++;
				}
				if (y <= 20) {
					y1 = y + 200;
					x1 = x;
				}
				else if (y <= 110 && y > 100) {
					y1 -= 200;
				}
				else if (y <= 120 && y > 110) {
					y1 = y1 + 200;
				}
				else {
					y1 = y - 200;
					x1 = x;
				}
				tmp = mark->GetBinContent(x1, y1);
				if (tmp) {
					qp[(int)tmp] += ori->GetBinContent(x1, y1) * m_ptab[i].pos;
					q[(int)tmp] += ori->GetBinContent(x1, y1);
					area[(int)tmp] ++;
				}
			}
		}
		else if (y <= 20) {
			double y1 = y + 200;
			double x1 = x;
			double tmp = mark->GetBinContent(x1, y1);
			if (tmp) {
				qp[(int)tmp] += ori->GetBinContent(x1, y1) * m_ptab[i].pos;
				q[(int)tmp] += ori->GetBinContent(x1, y1);
				area[(int)tmp] ++;
			}
		}
		else if (y > 200) {
			double y1 = y - 200;
			double x1 = x;
			double tmp = mark->GetBinContent(x1, y1);
			if (tmp) {
				qp[(int)tmp] += ori->GetBinContent(x1, y1) * m_ptab[i].pos;
				q[(int)tmp] += ori->GetBinContent(x1, y1);
				area[(int)tmp] ++;
			}
		}
	}

	map<int, TVector3>::iterator qpIt = qp.begin();
	map<int, double>::iterator qIt = q.begin();
	map<int, int>::iterator aIt = area.begin();
	map<int, TVector3> rec;
	int m = 0;
	while (qpIt != qp.end()) {
		LogInfo << "The " << qpIt->first << "th mass." << endl;
		TVector3 p = 1 / qIt->second * qpIt->second;
		if (m < 4) {
			map<int, TVector3>::iterator recIt = rec.begin();
			int i = 0;
			bool breakFlag = false;
			while (recIt != rec.end()) {
				if ((p - recIt->second).Mag() < 3000 &&
					(p.Theta() < 0.314 && (recIt->second).Theta() < 0.314 ||
					 p.Theta() > 2.826 && (recIt->second).Theta() > 2.826)) {
					if (area[recIt->first] < area[qpIt->first]) {
						mass[i] = (long int)p.Mag() * 1E6 + (long int)(p.Theta() / unit) * 1000 + (long int)((p.Phi() + TMath::Pi()) / unit);
						rec.erase(recIt);
					}
					breakFlag = true;
					break;
				}
				i ++;
				recIt ++;
			}
			if (!breakFlag) {
				long int pp = (long int)p.Mag() * 1E6 + (long int)(p.Theta() / unit) * 1000 + (long int)((p.Phi() + TMath::Pi()) / unit);
				mass[m] = pp;
				m ++;
			}
		}
		rec[qpIt->first] = p;
		qpIt ++;
		qIt ++;
	}

	for (int i = 0; i < 4; i ++)
		LogDebug << "mass[" << i << "]: " << mass[i] << endl;

	qp.clear();
	q.clear();
	area.clear();
	rec.clear();

	return mass;
}

double FhtAna::FHTPredict(const PmtProp& pmt, TVector3 inci, TVector3 dir, double ti) {
	double nLS = 1.485;
	double cLight = 299.;
	double vMuon = 299.;
	double nW = 1.34;

	double tan = TMath::Sqrt(nW * nW - 1);

	TVector3 perp = inci + dir * (pmt.pos - inci) * dir;

	TVector3 liSource = perp - (pmt.pos - perp).Mag() / tan * dir;

	TVector3 liDir = pmt.pos - liSource;

	double dis = TMath::Sqrt(TMath::Power(pmt.pos.Mag(), 2) - TMath::Power(pmt.pos * liDir.Unit(), 2));

	return ti + (liSource - inci) * dir / vMuon + (pmt.pos - liSource).Mag() * nW / cLight;
}
