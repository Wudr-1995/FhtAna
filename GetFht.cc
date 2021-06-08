#include "GetFht.h"
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

DECLARE_ALGORITHM(GetFht);

GetFht::GetFht(const std::string& name)
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
	declProp("FileNumber", m_turn = 0);
}

bool GetFht::initialize() {
	LogDebug << "Initializing" << std::endl;
	if(not initGeomSvc())
		return false;
	SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
	if (navBuf.invalid()) {
		LogError << "Cannot get the NavBuffer @ /Event" << std::endl;
		return false;
	}
	m_buf = navBuf.data();
	// SniperPtr<RootWriter> m_rw(*getParent(), "RootWriter");
	// if (!m_rw.valid()) {
	// 	LogError << "Cannot get the RootWriter instance" << std::endl;
	// 	return false;
	// }
	// m_tree = m_rw->bookTree("Trks", "Trks");
	// m_hist = new TH2D("FhtVTheta", "FhtVTheta", )
	m_path = "/junofs/users/wudr/Muon_Sim/Real_Muon/Rec_job/";
    return true;
}

bool GetFht::execute() {
	LogDebug << "executing: " << m_iEvt ++ << std::endl;
	TH2D* FhtDis = new TH2D("FhtDistribution", "FhtDistribution", 100, 0, 3.14, 500, 0, 500);
	TH2D* FhtPhi = new TH2D("FhtVPhi", "FhtVPhi", 200, -3.14, 3.14, 500, 0, 500);
	TH2D *Fht2D = new TH2D("FhtDistribution2D", "", 100, 0, 3.14, 200, -3.14, 3.14);
	TH2D *Q2D = new TH2D("ChargeDistribution2D", "", 100, 0, 3.14, 200, -3.14, 3.14);
	// nEvt.open(m_path + "nEvt.txt");
	// int globalCount;
	// if (!nEvt.good()) {
	// 	outFile.open(m_path + "nEvt.txt");
	// 	outFile << 1 << endl;
	// 	outFile.close();
	// 	nEvt.close();
	// 	globalCount = 1;
	// }
	// else {
	// 	nEvt >> globalCount;
	// 	globalCount ++;
	// 	nEvt.close();
	// 	outFile.open(m_path + "nEvt.txt");
	// 	outFile << globalCount << endl;
	// 	outFile.close();
	// }

	JM::SimEvent* simevent = 0;
	JM::EvtNavigator* nav =m_buf->curEvt();
	std::vector<std::string>& paths = nav->getPath();
	JM::SimHeader* simheader = 0;
	for (size_t i = 0; i < paths.size(); ++i) {
		const std::string& path = paths[i];
		if (path == "/Event/SimOrig") {
			simheader = static_cast<JM::SimHeader*>(nav->getHeader("/Event/SimOrig"));
			std::cout << "SimHeader (/Event/SimOrig): " << simheader << std::endl;
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
	if (freshPmtData(FhtDis, FhtPhi, Fht2D, Q2D, InciTheta, InciPhi))
		LogDebug << "Freshing PMT data success" << std::endl;
	else {
		LogError << "Freshing PMT data fails" << std::endl;
		delete FhtDis;
		delete FhtPhi;
		delete Fht2D;
		delete Q2D;
		return true;
	}

	// TString path = m_path + "Sample_" + m_turn + "_" + m_iEvt + ".txt";
	// outFile.open(path);
	// path = m_path + "SamplePhi_" + m_turn + "_" + m_iEvt + ".txt";
	// LogDebug << path << std::endl;
	// PhiOut.open(path);
	// path = m_path + "Label_" + m_turn + "_" + m_iEvt + ".txt";
	// label.open(path);
	// path = m_path + "Trks_" + m_turn + "_" + m_iEvt + ".txt";
	// trks.open(path);

	simevent = (JM::SimEvent*)simheader->event();
	LogDebug << "SimEventGot" << std::endl;
	nSimTrks = simevent->getTracksVec().size();
	LogDebug << "Retrieving tracks data" << std::endl;
	short NumCrossCd = 0;
	for (short i = 0; i < nSimTrks; i ++) {
		JM::SimTrack* strk = simevent->findTrackByTrkID(i);
		TVector3 Inci(strk->getInitX(), strk->getInitY(), strk->getInitZ());
		TVector3 Dir(strk->getInitPx(), strk->getInitPy(), strk->getInitPz());
		TVector3 Exit(strk->getExitX(), strk->getExitY(), strk->getExitZ());
		if (IfCrossCd(Inci, Dir, m_LSRadius)) {
			NumCrossCd ++;
			TVector3 LSInci = InciOnLS(Inci, Dir, m_LSRadius);
			TVector3 LSExit = Exit;
			if (Exit.Mag() > m_LSRadius) {
				TVector3 antiDir = - Dir;
				LSExit = InciOnLS(Exit, antiDir, m_LSRadius);
			}
			TVector3 dir = Dir.Unit();
			// trks << LSInci.X() << "\t" << LSInci.Y() << "\t" << LSInci.Z() << "\t"
			// 	<< dir.X() << "\t" << dir.Y() << "\t" << dir.Z() << endl;
			trks << LSInci.Theta() << "\t" << LSInci.Phi() << "\t"
				<< LSExit.Theta() << "\t" << LSExit.Phi() << "\t" << LSExit.Mag() << "\t"
				<< dir.Theta() << "\t" << dir.Phi() << endl;
		}
	}
	// label << NumCrossCd << endl;
	// for (short i = 1; i <= 100; i ++) {
		// for (short j = 1; j <= 500; j ++) {
			// outFile << FhtDis->GetBinContent(i, j) << "\t";
		// }
		// outFile << endl;
	// }
	// for (short i = 1; i <= 200; i ++) {
		// for (short j = 1; j <= 500; j ++) {
			// PhiOut << FhtPhi->GetBinContent(i, j) << "\t";
		// }
		// PhiOut << endl;
	// }
	TString pdfPath = m_path + "Fht2D_" + m_turn + "_" + m_iEvt + ".pdf";
	auto c1 = new TCanvas("Fht", "", 800, 1200);
	c1->Print(pdfPath + "[");
	TArrow *a1;
	TArrow *a2;
	for (short i = 0; i < nSimTrks; i ++) {
		JM::SimTrack* strk = simevent->findTrackByTrkID(i);
		TVector3 Inci(strk->getInitX(), strk->getInitY(), strk->getInitZ());
		TVector3 Dir(strk->getInitPx(), strk->getInitPy(), strk->getInitPz());
		TVector3 Exit(strk->getExitX(), strk->getExitY(), strk->getExitZ());
		if (IfCrossCd(Inci, Dir, m_LSRadius)) {
			TVector3 LSInci = InciOnLS(Inci, Dir, m_LSRadius);
			TVector3 LSExit = InciOnLS(Exit, Dir, m_LSRadius);
			if (i == 0) {
				a1 = new TArrow(LSInci.Theta(), LSInci.Phi(), LSExit.Theta(), LSExit.Phi(), 0.5, "|->");
			}
			else {
				a2 = new TArrow(LSInci.Theta(), LSInci.Phi(), LSExit.Theta(), LSExit.Phi(), 0.5, "|->");
			}
			
		}
	}
	// delete p1;
	// delete p2;
	c1->cd();
	TH1D *proj;
	TH1D *Fht1D = new TH1D("Fht1D", "", 100, 0, 100);
	TH1D *FhtvThe = new TH1D("FhtVTheta", "", 100, 0, 100);
	TH1D *Q1D = new TH1D("Q1D", "", 100, 0, 100);
	TH1D *QvThe = new TH1D("QVTheta", "", 100, 0, 100);
	for (int i = 1; i <= 100; i ++) {
		proj = Fht2D->ProjectionY("tmp", i, i, "");
		double tmpFht = 1000, FhtPhi;
		for (int j = 1; j <= 200; j ++) {
			double t = proj->GetBinContent(j);
			if (t && t < tmpFht) {
				FhtPhi = j;
				tmpFht = t;
			}
		}
		Fht1D->SetBinContent(i, FhtPhi);
		FhtvThe->SetBinContent(i, proj->GetBinContent(FhtPhi));
		proj->Delete();
		proj = NULL;
		proj = Q2D->ProjectionY("tmp", i, i, "");
		Q1D->SetBinContent(i, proj->GetMaximumBin());
		QvThe->SetBinContent(i, proj->GetBinContent(proj->GetMaximumBin()));
	}
	TH1D *FhtSmoth = FixCurve(FhtvThe, 100);
	TVector3 ExitPos = GetExitPos(FhtSmoth, Fht1D, 100);
	TVector3 InciPos;
	TArrow *recA;
	if (ExitPos[0] != 0 && ExitPos[1] != 0) {
		InciPos = GetInciPos(FhtSmoth, Fht1D, 100);
		// recA = new TArrow(InciPos.Theta(), InciPos.Phi(), ExitPos.Theta(), ExitPos.Phi(), 0.5, "|->");
		recA = new TArrow(InciTheta, InciPhi, ExitPos.Theta(), ExitPos.Phi(), 0.5, "|->");
		LogInfo << "Inci theta: " << InciPos.Theta() << ". Inci phi: " << InciPos.Phi() << endl;
		LogInfo << "Exit theta: " << ExitPos.Theta() << ". Exit phi: " << ExitPos.Phi() << endl;
	}
	TPad *p1 = new TPad("FTheta", "", 0, 0.5, 1, 1);
	// p1->SetTopMargin(0.1);
	// p1->SetBottomMargin(0.11);
	// p1->SetRightMargin(0.13);
	// p1->SetLeftMargin(0.15);
	p1->SetTopMargin(0.2);
	p1->SetBottomMargin(0);
	p1->SetRightMargin(0.2);
	p1->SetLeftMargin(0.2);
	p1->Draw();
	TPad *p2 = new TPad("FTheta", "", 0, 0, 1, 0.5);
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0.2);
	p2->SetRightMargin(0.2);
	p2->SetLeftMargin(0.2);
	p2->Draw();
	p1->cd();
	gStyle->SetOptStat(0000);
	gStyle->SetPalette(1);
	Fht2D->SetTitle("");
	Fht2D->GetXaxis()->SetTitleSize(0.05);
	Fht2D->GetYaxis()->SetTitleSize(0.05);
	Fht2D->GetXaxis()->SetLabelSize(0.05);
	Fht2D->GetYaxis()->SetLabelSize(0.05);
	Fht2D->Draw("colz");
	a1->SetLineColor(kRed);
	a1->SetFillColor(kRed);
	a1->SetLineWidth(2);
	a1->Draw("same");
	a2->SetLineColor(kRed);
	a2->SetFillColor(kRed);
	a2->SetLineWidth(2);
	a2->Draw("same");
	if (recA != NULL) {
		recA->SetLineColor(kBlue);
		recA->SetFillColor(kBlue);
		recA->SetLineWidth(2);
		recA->Draw("same");
	}
	// c1->SaveAs(pdfPath);

	// pdfPath = m_path + "Q2D_" + m_turn + "_" + m_iEvt + ".pdf";
	// auto c2 = new TCanvas("Fht", "", 800, 600);
	// c2->cd();
	p2->cd();
	gStyle->SetOptStat(0000);
	gStyle->SetPalette(1);
	Q2D->SetTitle("");
	Q2D->GetXaxis()->SetTitleSize(0.05);
	Q2D->GetYaxis()->SetTitleSize(0.05);
	Q2D->GetXaxis()->SetLabelSize(0.05);
	Q2D->GetYaxis()->SetLabelSize(0.05);
	Q2D->Draw("colz");
	// c1->SaveAs(pdfPath);
	c1->Print(pdfPath);
	// FhtvThe->Smooth();
	// TH1D *FhtSmoth = GetDiffInte(FhtvThe, 100);
	int nVs = ValleyFinder(FhtSmoth, 100);
	LogInfo << "The number of valleys: " << nVs << endl;
	TH1D *FhtDif = GetDiff(FhtvThe, 100);
	c1->cd();
	p1->cd();
	// Fht1D->Draw();
	FhtvThe->Draw();
	p2->cd();
	FhtSmoth->Draw();
	// FhtvThe->Draw();
	// FhtSmoth->Draw("same");
	c1->Print(pdfPath);
	p1->cd();
	Q1D->Draw();
	p2->cd();
	QvThe->Draw();
	c1->Print(pdfPath);
	c1->Print(pdfPath + "]");

	delete FhtDis;
	delete FhtPhi;
	delete Fht2D;
	delete Q2D;
	delete Fht1D;
	delete Q1D;
	delete FhtvThe;
	delete QvThe;
	delete FhtSmoth;
	delete FhtDif;
	// outFile.close();
	// PhiOut.close();
	// label.close();
	// trks.close();
	LogDebug << "Executed" << endl;
	return true;
}

bool GetFht::initGeomSvc() {
	SniperPtr<RecGeomSvc> rgSvc(getParent(), "RecGeomSvc");
	if (rgSvc.invalid()) {
		LogError << "Failed to get RecGeomSvc instance" << std::endl;
		return false;
	}
	m_geom = rgSvc->getCdGeom();
	return true;
}

bool GetFht::initPmt() {
	LogDebug << "Initializing PMTs" << std::endl;
	totPmtNum = 0;
	totPmtNum = m_geom->getPmtNum();
	if (!totPmtNum) {
		LogError << "Wrong PMT Number" << std::endl;
		return false;
	}
	LogDebug << "PMT Number Got" << std::endl;
	m_ptab.reserve(totPmtNum);
	m_ptab.resize(totPmtNum);
	for (unsigned int pid = 0; pid < totPmtNum; pid ++) {
		Identifier Id = Identifier(CdID::id(pid, 0));
		PmtGeom* pmt = m_geom->getPmt(Id);
		if (!pmt) {
			LogError << "Wrong PMT ID" << std::endl;
			return false;
		}
		TVector3 pmtCenter = pmt->getCenter();
		m_ptab[pid].pos = pmtCenter;
		if (CdID::is3inch(Id)) {
			m_ptab[pid].res = m_3inchRes;
			m_ptab[pid].type = _PMTINCH3;
		}
		else if (CdID::is20inch(Id)) {
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

bool GetFht::freshPmtData(TH2D* ht, TH2D* pht, TH2D *h2d, TH2D *q2d, double &theta, double &phi) {
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
		if (not (value & 0xFF000000) >> 24 == 0x10)
			continue;
		unsigned int pid = CdID::module(id);
		if (pid > totPmtNum) {
			LogError << "Data/Geometry Mis-Match : PmtId(" << pid << ") >= the number of PMTs." << std::endl;
			return false;
		}
		m_ptab[pid].q = calib->nPE();
		m_ptab[pid].fht = calib->firstHitTime();
		if ((CdID::is3inch(id) && m_3inchusedflag) || (CdID::is20inch(id) && m_20inchusedflag)) {
			if (earliest > m_ptab[pid].fht) {
				earliest = m_ptab[pid].fht;
				theta = m_ptab[pid].pos.Theta();
				phi = m_ptab[pid].pos.Phi();
			}
			m_ptab[pid].used = true;
			m_usedPmtNum ++;
			ht->Fill(m_ptab[pid].pos.Theta(), m_ptab[pid].fht);
			pht->Fill(m_ptab[pid].pos.Phi(), m_ptab[pid].fht);
			int binx = m_ptab[pid].pos.Theta() / 0.0314;
			int biny = (m_ptab[pid].pos.Phi() + TMath::Pi()) / 0.0314;
			h2d->SetBinContent(binx, biny,
							   m_ptab[pid].fht < h2d->GetBinContent(binx, biny) || h2d->GetBinContent(binx, biny) == 0 ?
							   m_ptab[pid].fht : h2d->GetBinContent(binx, biny));
			q2d->SetBinContent(binx, biny, q2d->GetBinContent(binx, biny) + m_ptab[pid].q);
			// LogDebug << m_ptab[pid].fht << std::endl;
			// LogDebug << h2d->GetBinContent(binx, biny) << std::endl;
		}
	}
	LogDebug << "Loading calibration data done" << std::endl;
	return true;
}

bool GetFht::IfCrossCd(TVector3& Inci, TVector3& Dir, Double_t R) {
	TVector3 dir = Dir.Unit();
	Double_t Dis = TMath::Sqrt(Inci.Mag() * Inci.Mag() - fabs(Inci * dir) * fabs(Inci * dir));
	LogDebug << "Distance to center: " << Dis << endl;
	if (R < Dis)
		return false;
	return true;
}

TVector3 GetFht::InciOnLS(TVector3& Inci, TVector3& Dir, Double_t R) {
	TVector3 dir = Dir.Unit();
	Double_t Dis2 = Inci.Mag() * Inci.Mag() - fabs(Inci * dir) * fabs(Inci * dir);
	TVector3 LSInci = Inci + (Inci * dir < 0 ? 1 : -1) * dir * (TMath::Sqrt(Inci.Mag() * Inci.Mag() - Dis2) - TMath::Sqrt(R * R - Dis2));
	return LSInci;
}

bool GetFht::finalize() {
	LogDebug << "Finalizing" << std::endl;
	return true;
}

TH1D* GetFht::FixCurve(TH1D *h, int n) {
	LogInfo << "In FixCurve" << endl;
	if (h == NULL)
		return NULL;
	double p1, p2, p3, p4;
	TH1D *tmph = new TH1D("tmpHist", "", n, 0, n);
	for (int i = 1; i <= n; i ++) {
		tmph->SetBinContent(i, h->GetBinContent(i));
	}
	TH1D *diff = GetDiff(tmph, 100);
	// tmph->Smooth();
	for (int i = 1; i <= n / 2;) {
		if (diff->GetBinContent(i) < -6) {
			LogInfo << "Start Content: " << diff->GetBinContent(i) << endl;
			int endEdge = 0;
			for (int j = i + 1; j <= n / 2; j ++) {
				if (diff->GetBinContent(j) > 5)
					endEdge = j;
				if ((j - endEdge > 10 || (diff->GetBinContent(j) < -6 && j - endEdge > 2)) && endEdge != 0)
					break;
			}
			if (endEdge == 0 || endEdge == n / 2) {
				i ++;
				continue;
			}
			LogInfo << "i + 1: " << i + 1 << endl;
			LogInfo << "endEdge: " << endEdge << endl;
			p1 = tmph->GetBinContent(i);
			p2 = tmph->GetBinContent(endEdge + 1);
			for (int j = i + 1; j <= endEdge; j ++) {
				LogInfo << "Content set: " << p1 + (p2 - p1) * (j - i) / (endEdge + 1 - i) << endl;
				LogInfo << "Content before & after: " << tmph->GetBinContent(j - 1) << ", " << tmph->GetBinContent(j + 1) << endl;
				tmph->SetBinContent(j, p1 + (p2 - p1) * (j - i) / (endEdge + 1 - i));
			}
			// i = endEdge + 1;
			diff->Delete();
			diff = GetDiff(tmph, 100);
		}
		else
			i ++;
		
	}
	// tmph->Smooth();
	/*
	for (int i = 1; i <= n - 2; i ++) {
		p1 = tmph->GetBinContent(i);
		p2 = tmph->GetBinContent(i + 1);
		p3 = tmph->GetBinContent(i + 2);
		if (TMath::Abs(p3 - p2) > 5) {
			double p3cp;
			int j;
			for (j = i + 3; j <= n - 4; j ++) {
				p4 = tmph->GetBinContent(j);
				p3cp = tmph->GetBinContent(j - 1);
				if (TMath::Abs(p4 - p3cp) > 5 &&
					diff->GetBinContent(j) < 5 &&
					diff->GetBinContent(j + 1) < 5 &&
					diff->GetBinContent(j + 2) < 5 &&
					diff->GetBinContent(j + 3) < 5)
					break;
			}
			for (int k = i + 2; k < j; k ++) {
				tmph->SetBinContent(k,
									p2 + (p4 - p2) / (j - i - 1) * (k - i - 1));
			}
			i = j;
		}
	}
	*/
	return tmph;
}

int GetFht::ValleyFinder(TH1D *h, int n) {
	int nVs = 0;
	for (int i = 5; i <= n - 4; i ++) {
		double Vi = h->GetBinContent(i);
		bool breakFlag = false;
		for (int j = i - 4; j <= i + 4; j ++) {
			if (j == i)
			continue;
			if (h->GetBinContent(j) < Vi) {
				breakFlag = true;
				break;
			}
		}
		if (!breakFlag)
			nVs ++;
	}
	return nVs;
}

TH1D* GetFht::GetDiff(TH1D *h, int n) {
	if (h == NULL) {
		LogInfo << "The input histgram is NULL." << endl;
		return NULL;
	}
	TH1D *tmph = new TH1D("DiffH", "", n, 0, n);
	for (int i = 1; i <= n; i ++)
		tmph->SetBinContent(i, h->GetBinContent(i));
	for (int i = 1; i < n; i ++) {
		tmph->SetBinContent(i, tmph->GetBinContent(i + 1) - tmph->GetBinContent(i));
	}
	tmph->SetBinContent(n, 0);
	return tmph;
}

TH1D* GetFht::GetDiffInte(TH1D *h, int n) {
	if (h == NULL) {
		LogInfo << "The input histgram is NULL." << endl;
		return NULL;
	}
	TH1D *tmph = GetDiff(h, n);
	tmph->SetBinContent(1, TMath::Abs(tmph->GetBinContent(1)));
	for (int i = 2; i <= n; i ++) {
		tmph->SetBinContent(i, tmph->GetBinContent(i - 1) + TMath::Abs(tmph->GetBinContent(i)));
	}
	return tmph;
}

TVector3 GetFht::GetInciPos(TH1D *ThevFht, TH1D *ThevPhi, int n) {
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

TVector3 GetFht::GetExitPos(TH1D *ThevFht, TH1D *ThevPhi, int n) {
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
	ExitPos.SetMagThetaPhi(17700, theta, phi);
	LogInfo << "Theta: " << theta << "\t"
			<< "Phi: " << phi << endl;
	return ExitPos;
}

TVector3 GetFht::GetChargeCenter() {
	
}
